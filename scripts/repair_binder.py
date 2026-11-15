"""Build candidates by repairing an archive binder in place, not by regrafting it.

## Why this operation and not humanisation

`scripts.humanise` moved a good binder's CDRs onto a fresh human germline
framework. Measured on GPU (FINDINGS 49, 50): 28 % median retention of the
donor's binding, 20 % when the CDR3 was kept residue for residue, and no
correlation between retention and the number of framework substitutions
(Spearman -0.13, p = 0.49, median 12). The conclusion recorded there is that
the binding of these sequences does not travel with the motif -- it depends on
the framework context humanisation replaces.

So this script replaces nothing. It takes the donor sequence as it is and
applies the *smallest* set of substitutions that the validator's two walls
force, with two constraints:

  1. every substitution moves a residue **towards a human germline residue**
     (an IGHV allele of the epoch's gene, or an IGHJ allele), so compliance can
     only improve;
  2. every substitution sits in the **CDR3 stalk (IMGT 105-106), the CDR3 tail
     (115-117) or FR4 (118-128)** -- the C-terminal strand that carries the J
     region. CDR1, CDR2 and the CDR3 apex (107-114) are never touched, and
     neither are FR1, FR2 and FR3.

The arithmetic that forces the size of the set (docs/REPAIR_PISTE.md):

  * novelty: the donor is itself archived at identity 1.0, so `whole_identity`
    has to fall under 0.95, which for a 120-mer costs 7 substitutions (8 to be
    safe); and because that leaves whole identity above 0.85, the CDR3 branch
    of `is_duplicate` turns into a rejection unless CDR3 similarity falls under
    0.85, which costs floor(0.15 L) + 1 substitutions inside the CDR3;
  * compliance: J identity is counted over 14 positions (IMGT 115-128), so the
    floor of 86 means at most one mismatch there.

The two walls are not in conflict for a donor that fails J: repairing J to a
germline allele *is* where the novelty substitutions come from. That is the
whole point -- for donors that already pass J there is nothing germline-directed
left to move, and the script refuses them rather than inventing apex mutations.

**What is not proven: that the binding survives.** Nothing here is measured on
a structure. The donor's design_iiptm belongs to the donor; the candidate has
to be scored on GPU. Every emitted sequence is checked against
`nova_lab.nanobodies.check_sequence` with the full archive, so compliance and
novelty are measured, and only those.

Usage:
    python -m scripts.repair_binder --cdr3-length 13 --allowed-v IGHV3-74 --n 20 \
        --out inputs/repaired.csv
    python -m scripts.repair_binder --cdr3-length 13 --allowed-v IGHV3-74 \
        --cutoff 26299        # honest replay: only donors that existed then
"""
import argparse
import json
from pathlib import Path

from anarci.germlines import all_germlines

from nova_lab import DATA_DIR
from nova_lab.config import NANOBODY_TARGET
from nova_lab.nanobodies import NanobodyArchive, check_sequence, germline_identity

VREF = all_germlines["V"]["H"]["human"]
JREF = all_germlines["J"]["H"]["human"]
J_WINDOW = tuple(range(115, 129))
TAIL = (115, 116, 117)          # the CDR3 positions the J reference covers
STALK = (105, 106)              # the CDR3 positions the V reference covers
APEX = tuple(range(107, 115))   # never touched
FR1 = tuple(range(1, 27))       # only the `germline` strategy reaches these
FR3 = tuple(range(66, 105))
PROTECT = {23, 49, 50, 104}     # canonical cysteines and the VHH FR2 hallmarks


def cdr3_mutations_needed(length: int, threshold: float = 0.85) -> int:
    """Substitutions inside a CDR3 of `length` that put cdr_similarity under
    `threshold`. `cdr_similarity` is 1 - mismatches / length for equal lengths."""
    return next(m for m in range(length + 1) if 1 - m / length < threshold)


def imgt_map(seq: str) -> dict | None:
    """{IMGT position: index in seq}, from ANARCI's own state vector -- the same
    numbering `germline_identity` scores against."""
    from anarci.anarci import check_for_j, run_hmmer
    records = [("query", seq)]
    try:
        alignments = run_hmmer(records, ncpu=1, hmmer_species=None)
        check_for_j(records, alignments, "imgt")
    except Exception:  # noqa: BLE001 - an unnumberable donor is skipped
        return None
    _, state_vectors, details = alignments[0]
    if len(state_vectors) != 1 or len(details) != 1 or details[0]["chain_type"] != "H":
        return None
    states = dict(state_vectors[0])
    return {p: i for p in range(1, 129) if (i := states.get((p, "m"))) is not None}


def _positions(reference) -> dict:
    return {p: aa for p, aa in enumerate(reference, 1) if aa != "-"}


def repair(seq: str, allowed_v: str, cdr3_need: int, budget: int = 8,
           pos_map: dict | None = None, strategy: str = "confined"
           ) -> tuple[str | None, dict | str]:
    """The repaired sequence, or (None, why not).

    `cdr3_need` substitutions have to land inside the CDR3 and `budget` in
    total; both come from the duplicate rule, not from taste.

    `strategy` says where the substitutions the novelty budget still needs may
    go once the J window is exhausted:

    * `confined` -- nowhere else. Only a donor whose J window is far from every
      human allele can pay for its own novelty, which is exactly the donor the
      J floor refuses. FR1, FR2 and FR3 are never touched.
    * `germline` -- FR1 then FR3, at positions where the donor already differs
      from the V allele of the epoch's gene, so V identity rises too. This is
      the only route for a donor that already passes J (its J window is already
      germline, so there is nothing left to move there), and it is the route
      FINDINGS 50 warns about: humanisation moved a median of 12 framework
      residues and kept 28 % of the binding. It moves far fewer, and FR2 (the
      former VL interface, which carries the VHH hallmarks) is left alone.
    """
    pos_map = pos_map if pos_map is not None else imgt_map(seq)
    if pos_map is None:
        return None, "not numberable"
    if not all(p in pos_map for p in J_WINDOW):
        return None, "IMGT 115-128 not fully numbered"

    def mismatches(ref: dict) -> list[int]:
        return [p for p in sorted(ref)
                if p in pos_map and p not in PROTECT and seq[pos_map[p]] != ref[p]]

    v_options = [(len(mismatches(r)), a, r)
                 for a, r in ((a, _positions(r)) for a, r in VREF.items())
                 if a.split("*")[0] == allowed_v]
    if not v_options:
        return None, f"no human V allele for {allowed_v}"
    _, v_allele, vref = min(v_options, key=lambda o: o[0])
    stalk = [p for p in STALK if p in mismatches(vref)]

    j_options = []
    for allele, reference in JREF.items():
        jref = _positions(reference)
        mism = mismatches(jref)
        j_options.append((allele, jref, mism, [p for p in mism if p in TAIL]))
    feasible = [o for o in j_options if len(o[3]) + len(stalk) >= cdr3_need]
    if not feasible:
        best = max(len(o[3]) for o in j_options) + len(stalk)
        return None, (f"CDR3 novelty needs {cdr3_need} germline-directed substitutions "
                      f"in the stalk and tail, at most {best} are available")
    # the allele that offers the most substitutions inside the window, so none has
    # to be spent in FR1/FR2/FR3
    feasible.sort(key=lambda o: (-len(o[2]), -len(o[3])))
    j_allele, jref, j_mism, j_tail = feasible[0]

    muts: dict[int, str] = {p: jref[p] for p in j_tail}
    for p in stalk:
        if sum(105 <= q <= 117 for q in muts) >= cdr3_need:
            break
        muts[p] = vref[p]
    if sum(105 <= q <= 117 for q in muts) < cdr3_need:
        return None, "CDR3 quota unreachable without an apex substitution"
    for p in j_mism:                       # the rest of the J window: J goes to 100 %
        muts.setdefault(p, jref[p])
    for p in stalk:                        # last resort inside the window
        if len(muts) >= budget:
            break
        muts.setdefault(p, vref[p])
    if len(muts) < budget and strategy == "germline":
        for p in FR1 + FR3:
            if len(muts) >= budget:
                break
            if p not in muts and seq[pos_map[p]] != vref.get(p, seq[pos_map[p]]):
                muts[p] = vref[p]
    if len(muts) < budget:
        where = ("the CDR3 stalk and FR4" if strategy == "confined"
                 else "the CDR3 stalk, FR4, FR1 and FR3")
        return None, (f"only {len(muts)} germline-directed substitutions inside "
                      f"{where}, the novelty budget is {budget}")

    out = list(seq)
    for p, aa in muts.items():
        out[pos_map[p]] = aa
    return "".join(out), {
        "v_allele": v_allele, "j_allele": j_allele,
        "mutations": {p: f"{seq[pos_map[p]]}{p}{aa}" for p, aa in sorted(muts.items())},
        "n_mut": len(muts),
        "n_cdr3": sum(105 <= p <= 117 for p in muts),
        "n_apex": sum(p in APEX for p in muts),
        "n_fr4": sum(118 <= p <= 128 for p in muts),
        "n_framework": sum(p < 105 for p in muts),
    }


def donors(target: str, cdr3_length: int, floor: float, allowed_v: str | None = None,
           cutoff: int | None = None):
    """The archive's good binders of the right CDR3 length, best first.

    `cutoff` is an honest replay: a donor is kept only if it carries an epoch
    before the cutoff, or if `archive_ledger` says we had already seen it. The
    archive has no submission dates, so an undated sequence is only admissible
    when the ledger stamps it early enough.
    """
    import pandas as pd
    df = pd.read_parquet(DATA_DIR / "dataset" / f"{target}_nanobodies.parquet")
    df = df[(df["design_iiptm"] >= floor) & (df["cdr3"].str.len() == cdr3_length)]
    if allowed_v:
        df = df[df["v_genes"].fillna("").str.split(",").apply(lambda g: allowed_v in g)]
    if cutoff is not None:
        # the same snapshot `scripts.propose_nanobody` replays from: the ledger
        # when it exists, the frozen prebuilt file otherwise. An undated row
        # outside it is a sequence the archive gained since and cannot be placed
        # in time, so it is not admissible -- this is the hole that
        # tests/test_humanise_cutoff.py pins open in `scripts.humanise`.
        from scripts.propose_nanobody import dated_archive_snapshot
        seen = dated_archive_snapshot(target, cutoff)
        dated = df["epoch"].notna() & (df["epoch"] < cutoff)
        df = df[dated | (df["epoch"].isna() & df["sequence"].isin(seen))]
    return df.sort_values("design_iiptm", ascending=False)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--target", default=NANOBODY_TARGET[0])
    parser.add_argument("--cdr3-length", type=int, required=True,
                        help="the epoch's CDR3 length; substitutions do not change it")
    parser.add_argument("--allowed-v", required=True, help="the epoch's V gene")
    parser.add_argument("--floor", type=float, default=0.72,
                        help="donor design_iiptm floor; 0.72 is the field's level")
    parser.add_argument("--n", type=int, default=20, help="candidates wanted")
    parser.add_argument("--budget", type=int, default=8,
                        help="substitutions for the whole-sequence novelty wall")
    parser.add_argument("--strategy", choices=("confined", "germline"),
                        default="confined",
                        help="confined: every substitution in the CDR3 stalk or FR4; "
                             "germline: FR1/FR3 may also take germline-directed ones")
    parser.add_argument("--cutoff", type=int,
                        help="replay: keep only donors that existed before this epoch")
    parser.add_argument("--allow-warnings", action="store_true",
                        help="accept a candidate the archive rule only warns about")
    parser.add_argument("--out", type=Path, help="CSV for scripts.boltzgen_score")
    args = parser.parse_args()

    pool = donors(args.target, args.cdr3_length, args.floor, args.allowed_v, args.cutoff)
    print(f"{len(pool)} donors bind at {args.floor}+ with a CDR3 of "
          f"{args.cdr3_length} and call {args.allowed_v}"
          + (f", seen before epoch {args.cutoff}" if args.cutoff else ""))
    if pool.empty:
        return 1

    need = cdr3_mutations_needed(args.cdr3_length,
                                 0.85 if args.allow_warnings else 0.80)
    print(f"CDR3 similarity has to fall under {0.85 if args.allow_warnings else 0.80}: "
          f"{need} substitutions inside a CDR3 of {args.cdr3_length}")
    archive = NanobodyArchive(args.target)
    rows, refused = [], {}
    for row in pool.itertuples():
        if len(rows) >= args.n:
            break
        seq, info = repair(row.sequence, args.allowed_v, need, args.budget,
                           strategy=args.strategy)
        if seq is None:
            refused[info] = refused.get(info, 0) + 1
            continue
        warnings: list[str] = []
        problems = check_sequence(seq, allowed_v=args.allowed_v,
                                  allowed_cdr3_length=args.cdr3_length,
                                  archive=archive, warnings=warnings)
        if problems or (warnings and not args.allow_warnings):
            key = (problems + warnings)[0].split(":")[0]
            refused[key] = refused.get(key, 0) + 1
            continue
        g = germline_identity(seq)
        rows.append({"sequence": seq, "donor": row.sequence,
                     "donor_design_iiptm": row.design_iiptm,
                     "v_identity": round(g["V"]["identity_pct"], 1),
                     "j_identity": round(g["J"]["identity_pct"], 1),
                     "n_mut": info["n_mut"], "n_cdr3": info["n_cdr3"],
                     "n_apex": info["n_apex"], "n_fr4": info["n_fr4"],
                     "n_framework": info["n_framework"],
                     "warnings": "; ".join(warnings),
                     "mutations": json.dumps(info["mutations"])})
    print(f"{len(rows)} candidates pass check_sequence against the whole archive")
    for why, n in sorted(refused.items(), key=lambda kv: -kv[1])[:6]:
        print(f"  {n:4d} refused: {why}")
    if not rows:
        return 1
    import pandas as pd
    out = pd.DataFrame(rows)
    print(f"\ndonor design_iiptm: median {out.donor_design_iiptm.median():.3f}, "
          f"max {out.donor_design_iiptm.max():.3f} -- that is the DONOR's score. "
          f"What these read is UNTESTED and has to be measured on GPU.")
    print(f"substitutions: median {out.n_mut.median():.0f}, "
          f"of which {out.n_cdr3.median():.0f} in the CDR3 stalk and tail, "
          f"{out.n_apex.sum()} in the CDR3 apex, {out.n_framework.sum()} in FR1/FR2/FR3")
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(args.out, index=False)
        print(f"\nsaved to {args.out} -- screen with scripts.preflight, "
              f"score with scripts.boltzgen_score")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
