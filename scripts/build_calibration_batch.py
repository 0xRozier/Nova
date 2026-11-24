"""The batch that asks whether our measurement reads binding at all.

Every GPU batch so far compared our candidates to the field's and asked which
ranks higher. None of them ever checked the premise: that BoltzGen's iiptm,
computed with the epoch's binding site, is larger for a sequence that really
binds that site than for one that does not. We had no sequence of known truth
to check it with.

FINDINGS 73 produced five. The PDB holds five single-domain binders of this
week's target (IL-6, P05231), and their epitopes overlap the config's eight
residues to very different degrees:

    61H7       (4O9H)  8/8   the structure the config itself names
    VHH6       (5FUC)  7/8   a true VHH, but junctional: needs IL-6R-alpha
    68F2       (4ZS7)  3/8
    I9H        (8YWQ)  2/8
    Olokizumab (4CNI)  0/8   binds IL-6 elsewhere entirely

That is a calibration ladder: same target, same scoring, known and graded
truth. It costs 15 of the 87 sequences here and it is the most valuable part
of the batch, because the answer changes how every earlier measurement reads.

The second question is FINDINGS 73 against FINDINGS 68. 68 measured
"donor CDR3 onto our framework" and got 0.338, the worst arm we have run. 73
proposes the opposite -- keep the binder's framework, redesign only the CDR3 --
and 61H7 is the one binder for which that is not absurd: its paratope is only
25% CDR3 (CDR1 31%, CDR2 38%, framework 6%), so three quarters of its contact
surface survives a CDR3 swap. The other four are 40-62% CDR3.

So the batch crosses framework against CDR3 and keeps the controls that make
the cross readable:

    ref_intact    the six crystallographic binders           6
    camelised     61H7 made submittable, and the steps       3
    scramble      camelised framework, CDR3 interior shuffled 5
    keep_fw       camelised framework, archive CDR3s         32
    move_cdr3     61H7's CDR3 on other frameworks             9
    near_paratope archive sequences nearest 61H7's CDR1/2    15
    archive_far   ordinary archive candidates, the mode      15

There are no replicates, and not by choice: scripts.boltzgen_score names each
input by the sequence's SHA-256, so the same sequence submitted twice is one
file and one measurement. The sigma of 0.017-0.028 that FINDINGS 49 measured
cannot be resampled inside one run, which is why precision here comes from the
number of DISTINCT sequences in an arm rather than from repeats.

That shapes how the ladder is read. Five anchors alone are too few for a
Spearman estimate to mean much, so the ladder's floor is widened with genuine
non-binders of IL-6: `archive_far` holds 15 sequences designed against the
previous target (PD-L1), and 8AOK's VHH binds PD-L1 crystallographically. If
61H7 does not clear that floor, no ordering among the five matters.

READINGS FIXED BEFORE THE RUN. These are the point of the batch and must not be
renegotiated after seeing it:

  - If `ref_intact` does not order by epitope overlap -- Spearman rho under
    0.3 across the five, or 61H7 at or below Olokizumab -- then iiptm computed
    at this binding site does not read binding, and every GPU hour since
    FINDINGS 49 measured something else. That outcome is more important than
    any ranking this batch could produce, and it is the one I expect least.
  - If `camelised` falls below `ref_intact`'s 61H7 by more than the replicate
    spread, the three mutations that make it submittable also cost the binding,
    and the submittable sequence is not the binder. That would close this route
    as cleanly as FINDINGS 68 closed the graft, and it is the risk this arm is
    here to price: IMGT 49 and 50 are where a VHH fills the space a light chain
    used to occupy, so camelising a light-chain-paired VH is not free in the
    literature either.
  - If `scramble` matches the camelised scaffold within noise, the metric
    ignores the CDR3 and `keep_fw` versus `move_cdr3` cannot be read from it.
  - FINDINGS 73 is supported only if `keep_fw` beats `move_cdr3` by more than
    the replicate spread AND `keep_fw` reaches `ref_intact`'s 61H7. Beating
    `archive_far` alone is not enough: that is the comparison FINDINGS 59
    already showed is confounded by scaffold choice.

Usage:
    python -m scripts.build_calibration_batch --out inputs/kaggle_calibration.csv
"""
import argparse
import random
from pathlib import Path

import pandas as pd

from nova_lab import DATA_DIR, REPO_ROOT
from nova_lab.config import MIN_HUMAN_J_IDENTITY
from nova_lab.nanobodies import check_sequence, germline_identity, regions

# The VH domains of FINDINGS 73, read out of the mmCIF files, truncated at the
# end of FR4 (the W-G-x-G-T...VTVSS motif). Provenance in the docstring.
BINDERS = {
    "61H7": ("EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYRMYWVRQPPGKGLEWVSAISAGGGSTYYGDSVKGRFTIS"
             "RDNAKNTVYLQMNSLKPEDTAVYYCANRAGWGMGDYWGQGTQVTVSS"),
    "VHH6": ("DVQFVESGGGSVHAGGSLRLNCATSGYIYSTYCMGWFRQAPGKEREGVAHIYTNSGRTYYADSVKGRFTIS"
             "QDNAKNTVYLQMNSLKPEDTAIYYCAARPSIRCASFSATEYKDWGQGTQVTVSS"),
    "68F2": ("EVQLQESGPGLVKPSQTLSLTCTVSGGSITTRYYAWSWIRQPPGKGLEWMGVIDYDGDTYYSPSLKSRTSI"
             "SWDTSKNQFSLQLSSVTPEDTAVYYCARDPDVVTGFHYDYWGQGTQVTVSS"),
    "I9H": ("QITLKESGPALVKPTQTLTLTCTFSGFSLSTSGMGVGWIRQPPGKALEWLAHIWWDDDKYYNPSLKSRLTI"
            "SKDTSKNQVVLTMTNMDPVDTATYYCARRANYGTSYDYGMDYWGQGTLVTVSS"),
    # 8AOK. A true VHH that binds PD-L1, the PREVIOUS target, with its epitope
    # verified 21/21 in FINDINGS 66. Against IL-6 it is a known non-binder, and
    # it passes the hallmark rule that 61H7 fails -- so it anchors the floor
    # with a sequence that is camelid, submittable and simply aimed elsewhere.
    "8AOK-VHH": ("GEVQLVESGGGLVQAGGSLRLSCAASGRTFSNYHMAWFRQAPGKEREFVAGISWTGRGTYYTDSVKGRF"
                 "TISRDNAKNTVYLQMNSLKPEDTAVYYCAAEGTLYGSGGRTHQSAYDYWGQGTQVTVSS"),
    # 4CNI. The bottom rung: Olokizumab binds IL-6 at a site that shares none
    # of the config's eight residues, so a metric that reads the binding site
    # must score it below the others. Inlined rather than parsed from the
    # mmCIF, which upstream does not carry and which is 1.2 MB.
    "Olokizumab": ("EVQLVESGGGLVQPGGSLRLSCAASGFNFNDYFMNWVRQAPGKGLEWVAQMRNKNYQYGTYYAESLEGRF"
                   "TISRDDSKNSLYLQMNSLKTEDTAVYYCARESYYGFTSYWGQGTLVTVSS"),
}
# 61H7 made submittable. The live upstream config enforces three rules it
# fails: enforce_vhh_hallmarks (E/Q at IMGT 49, R at 50), min_human_v_identity
# 82 and min_human_j_identity 86. Raw 61H7 has H49=G, H50=L and J identity
# 85.7, so it is refused twice over. Three point mutations fix both -- G49E,
# L50R, and Q->L in FR4 towards human IGHJ4 -- and NONE of the three is a
# contact: the paratope sits in CDR1 (T28 S30 S31 Y32 R33), CDR2 (S52 A53 G54
# G55 G56 S57), CDR3 (R99 A100 G101 W102) and Y59, and all four regions come
# out byte-identical. Measured after the change: VHH hallmarks present, V 86.7
# on IGHV3-23, J 92.9, nativeness 0.657, humanness 0.937, Layer 1 clean.
# Whether those three mutations cost binding is what the `camelised` arm asks.
CAMELISED = ("EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYRMYWVRQPPGKEREWVSAISAGGGSTYYGDSVKGRFTIS"
             "RDNAKNTVYLQMNSLKPEDTAVYYCANRAGWGMGDYWGQGTLVTVSS")

# Epitope overlap with the config's eight residues, measured in FINDINGS 73.
# This is the ladder `ref_intact` is read against.
OVERLAP = {"61H7": 8, "VHH6": 7, "68F2": 3, "I9H": 2, "Olokizumab": 0,
           "8AOK-VHH": 0}

ARCHIVE_REGIONS = DATA_DIR / "archive" / "Q9NZQ7_regions.csv"

# 61H7's CDR1/CDR2 paratope positions, 0-based in its own CDR strings.
PARATOPE_CDR1, PARATOPE_CDR2 = [2, 4, 5, 6, 7], [1, 2, 3, 4, 5, 6]
EPOCH_CDR3_LENGTHS = [11, 13, 14, 19]   # 61H7's own, then lengths epochs asked for
SEED = 73



def swap_cdr3(sequence: str, cdr3: str) -> str | None:
    """`sequence` with its CDR3 replaced. None when the CDR3 is not placeable."""
    r = regions(sequence)
    if not r or sequence.count(r["CDR3"]) != 1:
        return None
    return sequence.replace(r["CDR3"], cdr3)


def j_compatible(sequence: str) -> bool:
    """Does this built sequence still clear the validator's J floor?

    Swapping a CDR3 moves the J-region alignment, because the CDR3's last
    residues are J-encoded: a foreign CDR3 whose tail does not match human
    IGHJ4 drops J identity from 92.9 to 78.6 and the submission is refused
    (the live config sets min_human_j_identity 86; FINDINGS 7ter measured 58%
    of candidates lost this way). Unfiltered, two thirds of `keep_fw` would be
    sequences we could never submit, which is not the arm we want to price.
    """
    g = germline_identity(sequence)
    return g is not None and g["J"]["identity_pct"] >= MIN_HUMAN_J_IDENTITY


def paratope_score(cdr1: str, cdr2: str, ref1: str, ref2: str) -> int:
    """How many of 61H7's 11 CDR1/CDR2 contact positions a sequence reproduces."""
    return (sum(1 for i in PARATOPE_CDR1 if i < len(cdr1) and cdr1[i] == ref1[i])
            + sum(1 for i in PARATOPE_CDR2 if i < len(cdr2) and cdr2[i] == ref2[i]))


def build() -> pd.DataFrame:
    rng = random.Random(SEED)
    binders = dict(BINDERS)
    ref = regions(binders["61H7"])
    if not ref:
        raise SystemExit("61H7 n'est pas numerotable : rien d'autre ne tient")
    rows: list[dict] = []

    def add(arm: str, label: str, seq: str, note: str = ""):
        rows.append({"sequence": seq, "arm": arm, "label": label, "note": note})

    # 1. the ladder
    for name, seq in binders.items():
        add("ref_intact", name, seq, f"epitope {OVERLAP[name]}/8")

    # 2. what the three submission-enabling mutations cost
    add("camelised", "61H7 camelise+FR4", CAMELISED,
        "la version soumettable : G49E, L50R, FR4 Q->L")
    add("camelised", "61H7 camelise", CAMELISED.replace("WGQGTLVTVSS", "WGQGTQVTVSS"),
        "hallmarks seuls ; reste refuse sur J 85,7 < 86")
    add("camelised", "61H7 FR4 humain", binders["61H7"].replace("WGQGTQVTVSS", "WGQGTLVTVSS"),
        "FR4 seul ; reste refuse sur les hallmarks")

    # 3. negative control at constant framework
    # Only the middle is shuffled. The first three and last three residues of a
    # CDR3 are V- and J-encoded, so shuffling them would change J identity and
    # make the control fail on a rule rather than on binding -- the comparison
    # with `keep_fw` has to hold everything but the loop's interior constant.
    head, middle, tail = ref["CDR3"][:3], list(ref["CDR3"][3:-3]), ref["CDR3"][-3:]
    seen = {ref["CDR3"]}
    for _ in range(5):
        for _try in range(200):
            rng.shuffle(middle)
            shuffled = head + "".join(middle) + tail
            if shuffled not in seen:
                break
        seen.add(shuffled)
        built = swap_cdr3(CAMELISED, shuffled)
        if built:
            add("scramble", f"camelise+{shuffled}", built,
                "CDR3 melange, composition conservee")

    # 4 and 5: the cross, from real archive material
    arch = pd.read_csv(ARCHIVE_REGIONS).dropna(subset=["CDR1", "CDR2", "CDR3"])
    arch = arch.drop_duplicates("sequence")
    arch["p"] = [paratope_score(a, b, ref["CDR1"], ref["CDR2"])
                 for a, b in zip(arch.CDR1, arch.CDR2)]

    for length in EPOCH_CDR3_LENGTHS:
        pool = sorted({c for c in arch.CDR3 if len(c) == length})
        if not pool:
            continue
        taken = 0
        for cdr3 in rng.sample(pool, len(pool)):
            if taken == 8:
                break
            built = swap_cdr3(CAMELISED, cdr3)
            if built and j_compatible(built):
                add("keep_fw", f"camelise+{cdr3}", built, f"CDR3 d'archive L={length}")
                taken += 1
        if taken < 8:
            print(f"  note : seulement {taken} CDR3 d'archive de longueur {length} "
                  f"restent au-dessus du plancher J")

    hosts = [(n, binders[n]) for n in ("VHH6", "68F2", "I9H")]
    hosts += [(f"arch{i}", s) for i, s in
              enumerate(rng.sample(sorted(arch.nlargest(200, "p").sequence), 6), 1)]
    for name, host in hosts:
        built = swap_cdr3(host, ref["CDR3"])
        if built:
            add("move_cdr3", f"{name}+61H7cdr3", built,
                "le bras du §68 : CDR3 valide sur une charpente etrangere")

    # 6 and 7: archive sequences near and far from 61H7's paratope
    for seq, p in arch.nlargest(15, "p")[["sequence", "p"]].itertuples(index=False):
        add("near_paratope", f"p{p}", seq, f"paratope CDR1+CDR2 {p}/11")
    mode = arch[arch.p.between(4, 5)]
    for seq, p in mode.sample(15, random_state=SEED)[["sequence", "p"]].itertuples(index=False):
        add("archive_far", f"p{p}", seq, f"paratope CDR1+CDR2 {p}/11")

    frame = pd.DataFrame(rows)
    frame["group"] = "candidate"
    frame["epoch"] = pd.NA
    frame["true_rank_sum"] = pd.NA
    # Recorded, never used to drop a row: a calibration anchor is wanted even
    # when it would fail the validator -- VHH6 is refused on humanness and is
    # still the 7/8 rung of the ladder.
    frame["layer1"] = [", ".join(check_sequence(s)) or "ok" for s in frame.sequence]
    # Depth-major so an interrupted run still covers every arm.
    frame["depth"] = frame.groupby("arm").cumcount()
    return frame.sort_values(["depth", "arm"]).drop(columns="depth").reset_index(drop=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path,
                        default=REPO_ROOT / "inputs" / "kaggle_calibration.csv")
    args = parser.parse_args()

    frame = build()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.out, index=False)

    print(f"{len(frame)} sequences -> {args.out}")
    counts = frame.arm.value_counts()
    for arm in ("ref_intact", "camelised", "scramble", "keep_fw", "move_cdr3",
                "near_paratope", "archive_far"):
        if arm in counts:
            ok = (frame[frame.arm == arm].layer1 == "ok").sum()
            print(f"  {arm:<14} {counts[arm]:3d}   {ok} passent la couche 1 hors epoque")
    print(f"  {len(frame) * 1.4 / 60:.1f} h sur 2x T4 a 1,4 min/sequence")
    print("  les lectures sont fixees dans le docstring du module, avant le run")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
