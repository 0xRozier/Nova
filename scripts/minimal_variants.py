"""Perturb ONE good parent as little as the rules allow, instead of mixing two.

The diagnosis this comes from (FINDINGS 15). Our generator takes a framework
from one archived sequence and a CDR3 from another, and the result caps at
design_iiptm 0.674 while both donors score 0.73-0.79. The parts do not
transfer: a nanobody's binding depends on the specific pairing of a framework
with the loop it presents. Four other explanations were measured and refuted --
the GPU (our distribution is bounded at 0.674, so throughput cannot help), the
donors (they are excellent and from the current regime), the distance to the
archive (rho = -0.07), and too little design surface (three CDRs scored worse
than one).

So: keep the pairing. Take one parent whole and change only what the rules
force, which is far less than a graft changes:

    our generator     16 of 16 CDR3 residues  (a whole foreign loop, then mutated)
    minimum legal      3 of 16                (+ 4 elsewhere)

The rules, from the validator's own is_duplicate (nova-upstream
utils/nanobodies.py:164, a hard rejection there):

    whole-sequence identity >= 0.95        -> duplicate
    CDR3 identity >= 0.85                  -> duplicate (same clonotype)
    CDR3 >= 0.80 with CDR1 and CDR2 >= 0.90 -> duplicate

The third tier binds here, because our variants keep CDR1 and CDR2 whole and
so sit at 100% on both. The CDR3 must therefore come in under 0.80, not 0.85:
for 16 residues that is 4 substitutions (12/16 = 75%), not 3. Learned the hard
way -- the first run built its variants with 3, read 81.25%, and had all 24
refused as duplicates by tier 3. Plus at least 7 substitutions in all
(116/123 = 94% < 95%), so 3 more go in the framework.

Where the other 4 go matters. Framework positions mostly do not contact the
target, and moving them TOWARD the archive consensus -- which is close to the
human germline the validator scores against -- creates the required difference
while pushing V identity up rather than down. The CDR3's first two and last
three residues are left alone: they are V- and J-encoded (FINDINGS 7ter).

Feasibility, checked before writing this: of 45 (V gene, CDR3 length)
combinations, 41 have parents, typically dozens, all with a measured
design_iiptm, the best per combination between 0.68 and 0.79.

This is the fifth hypothesis of the investigation and the four before it were
refuted, so the point of the script is to be measurable, not to be believed:
it emits a batch whose parents' scores are known, so the variants can be read
as a delta against them.

Usage:
    python -m scripts.minimal_variants --allowed-v IGHV3-74 --cdr3 16 --parents 8
"""
from __future__ import annotations

import argparse
import random
from collections import Counter
from pathlib import Path

import pandas as pd

from nova_lab import DATA_DIR
from nova_lab.cdr3_generators import J_TAIL
from nova_lab.config import (MIN_HUMAN_J_IDENTITY, MIN_HUMAN_V_IDENTITY,
                             NANOBODY_TARGET)
from nova_lab.nanobodies import (NanobodyArchive, check_sequence,
                                 germline_identity, regions)

# From is_duplicate: a submission must differ from every archived sequence by
# more than 5% of the whole sequence and more than 15% of the CDR3.
WHOLE_IDENTITY_CAP = 0.95
CDR3_IDENTITY_CAP = 0.85
# Tier 3: a CDR3 at or above this with CDR1 and CDR2 conserved is a duplicate.
# Binding for us, since we keep CDR1 and CDR2 whole.
CDR3_CONSERVED_CDR_CAP = 0.80
V_HEAD = 2   # V-encoded, left alone
J_FIXED = 3  # J-encoded, left alone
AAS = "ACDEFGHIKLMNPQRSTVWY"
# PSH is the hydrophobic patch score of the CDR SURFACE, and substituting from
# a uniform draw over all twenty puts a hydrophobic residue in the binding loop
# two times in five. Measured: variants built that way read PSH 100.6 median
# against 89.4 for their accepted parents, +11.2 points, and 72% breach the
# measured break of 96 against 25% of parents. So substitute from the
# hydrophilic residues only. C is excluded anyway (the canonical cysteines are
# structural and filter_cysteine rejects extras), and G and P are kept out of
# this set because they change the loop's backbone rather than its surface.
HYDROPHILIC = "DEKNQRSTH"
REGION_ORDER = ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")


def required_changes(sequence_length: int, cdr3_length: int) -> tuple[int, int]:
    """(CDR3 substitutions, total substitutions) the rules force, at minimum.

    is_duplicate has THREE tiers, and the third is the one that matters when
    CDR1 and CDR2 are left alone:

        identity >= 0.95                              -> duplicate
        CDR3 >= 0.85                                  -> duplicate
        CDR3 >= 0.80 AND CDR1 >= 0.90 AND CDR2 >= 0.90 -> duplicate

    Our variants keep CDR1 and CDR2 whole, so they sit at 100% on both and the
    third tier binds: the CDR3 must come in under 0.80, not 0.85. Measured the
    hard way -- with 3 substitutions of 16 the CDR3 reads 81.25%, clears tier 2
    and is caught by tier 3, and all 24 variants of the first run were refused
    as duplicates. Four substitutions give 75%.
    """
    in_cdr3 = int(cdr3_length * (1 - CDR3_CONSERVED_CDR_CAP)) + 1
    total = int(sequence_length * (1 - WHOLE_IDENTITY_CAP)) + 1
    return in_cdr3, max(total, in_cdr3)


def changes_for(sequence_length: int, cdr3_length: int,
                in_cdr3: int | None = None) -> tuple[int, int]:
    """Same, but allowing MORE CDR3 substitutions than the rules force.

    The minimum turned out to be the wrong place to stand (FINDINGS 19: 20 of
    20 minimal variants scored below their parent, median -0.172). The epochs'
    actual winners sit at a CDR3 similarity of 0.667 to their nearest archived
    ancestor -- roughly a third of the loop changed, where the minimum changes
    a quarter -- and gain +0.074 over it (FINDINGS 20). So the substitution
    count has to be a dial, not a constant, or that band cannot be tested.
    """
    minimum, total = required_changes(sequence_length, cdr3_length)
    if in_cdr3 is None:
        return minimum, total
    if in_cdr3 < minimum:
        raise ValueError(
            f"{in_cdr3} CDR3 substitutions would leave the loop at "
            f"{1 - in_cdr3 / cdr3_length:.0%} identity, and tier 3 of "
            f"is_duplicate refuses anything at or above 80% while CDR1 and "
            f"CDR2 are untouched: the minimum here is {minimum}"
        )
    return in_cdr3, max(total, in_cdr3)


def consensus(sequences: list[str], region: str) -> dict[int, str]:
    """Most common residue per position of a region, over the archive.

    Used to move framework positions toward what the archive agrees on, which
    is close to the human germline the validator scores identity against. Only
    sequences whose region has the same length contribute, since positions are
    only comparable then.
    """
    by_length: dict[int, list[str]] = {}
    for sequence in sequences:
        r = regions(sequence)
        if not r:
            continue
        part = r[region]
        by_length.setdefault(len(part), []).append(part)
    if not by_length:
        return {}
    length = max(by_length, key=lambda k: len(by_length[k]))
    parts = by_length[length]
    return {i: Counter(p[i] for p in parts).most_common(1)[0][0]
            for i in range(length)}


def parents(target: str, allowed_v: str, cdr3_length: int,
            n: int) -> pd.DataFrame:
    """Archived sequences that already have this epoch's V gene and CDR3 length.

    Ordered by their own measured design_iiptm: the point of the strategy is to
    start from something known to bind.
    """
    df = pd.read_parquet(DATA_DIR / "dataset" / f"{target}_nanobodies.parquet")
    ok = df[(df["v_identity"] >= MIN_HUMAN_V_IDENTITY)
            & (df["j_identity"] >= MIN_HUMAN_J_IDENTITY)
            & (df["cdr3"].str.len() == cdr3_length)
            & df["v_genes"].str.split(",").apply(lambda g: allowed_v in g)
            & df["design_iiptm"].notna()]
    return ok.nlargest(n, "design_iiptm")[
        ["sequence", "cdr3", "design_iiptm", "v_identity", "j_identity", "epoch"]]


def variant(parent: str, cdr3_length: int, fr_consensus: dict,
            rng: random.Random, cdr3_subs: int | None = None) -> tuple[str, dict] | None:
    """One perturbed child of this parent, or None if it cannot be built.

    cdr3_subs=None keeps the legal minimum; a larger number walks out to the
    band the winners occupy.
    """
    r = regions(parent)
    if not r:
        return None
    in_cdr3, total = changes_for(len(parent), cdr3_length, cdr3_subs)
    loop = list(r["CDR3"])
    free = list(range(V_HEAD, len(loop) - J_FIXED))
    if len(free) < in_cdr3:
        return None
    for position in rng.sample(free, in_cdr3):
        loop[position] = rng.choice([a for a in HYDROPHILIC if a != loop[position]])

    # The rest go in the framework, toward the archive consensus: those
    # positions do not contact the target and the move raises V identity
    # instead of lowering it.
    parts = {key: r[key] for key in REGION_ORDER}
    parts["CDR3"] = "".join(loop)
    candidates = []
    for region in ("FR1", "FR2", "FR3"):
        for i, residue in enumerate(parts[region]):
            wanted = fr_consensus.get(region, {}).get(i)
            if wanted and wanted != residue:
                candidates.append((region, i, wanted))
    rng.shuffle(candidates)
    changed = []
    for region, i, wanted in candidates[:max(total - in_cdr3, 0)]:
        text = list(parts[region])
        text[i] = wanted
        parts[region] = "".join(text)
        changed.append(f"{region}{i + 1}{wanted}")

    child = "".join(parts[key] for key in REGION_ORDER)
    return child, {"cdr3_substitutions": in_cdr3,
                   "framework_moves": len(changed),
                   "framework_detail": ",".join(changed)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default=NANOBODY_TARGET[0])
    parser.add_argument("--allowed-v", required=True)
    parser.add_argument("--cdr3", type=int, required=True)
    parser.add_argument("--parents", type=int, default=8)
    parser.add_argument("--per-parent", type=int, default=3)
    parser.add_argument("--cdr3-subs", type=int,
                        help="CDR3 substitutions per variant (default: the legal "
                             "minimum). Higher walks toward the winners' band.")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", type=Path,
                        default=Path("inputs/kaggle_minimal_variants.csv"))
    args = parser.parse_args()

    pool = parents(args.target, args.allowed_v, args.cdr3, args.parents)
    if pool.empty:
        raise SystemExit(f"no scored parent with {args.allowed_v} and a "
                         f"{args.cdr3}-residue CDR3")
    print(f"{len(pool)} parents, their own design_iiptm "
          f"{pool.design_iiptm.min():.3f}-{pool.design_iiptm.max():.3f}")
    in_cdr3, total = changes_for(len(pool.sequence.iat[0]), args.cdr3, args.cdr3_subs)
    minimum, _ = required_changes(len(pool.sequence.iat[0]), args.cdr3)
    print(f"{in_cdr3} CDR3 substitutions and {total} in all "
          f"(the rules force {minimum}); a graft changes all {args.cdr3}")
    print(f"  CDR3 identity to the parent: {1 - in_cdr3 / args.cdr3:.0%}  "
          f"(winners sit near 67%, the minimum here is {1 - minimum / args.cdr3:.0%})")

    archive = NanobodyArchive(args.target)
    fr_consensus = {region: consensus(list(archive.sequences), region)
                    for region in ("FR1", "FR2", "FR3")}
    rng = random.Random(args.seed)
    rows, refused = [], Counter()
    for parent in pool.itertuples():
        made = 0
        for _ in range(args.per_parent * 40):
            if made >= args.per_parent:
                break
            built = variant(parent.sequence, args.cdr3, fr_consensus, rng,
                            args.cdr3_subs)
            if not built:
                continue
            child, how = built
            # The cheap gates first, and the archive rule only on what survives
            # them. The rule costs 5.6 s per candidate -- a difflib alignment
            # against each of 14424 archived sequences -- so running it inside
            # the retry loop made this script take 90 minutes and get killed.
            # FINDINGS 7bis had already fixed exactly this in the miner; this is
            # the same defect written again in a new file.
            problems = check_sequence(child, args.allowed_v, args.cdr3)
            if problems:
                refused[problems[0][:44]] += 1
                continue
            identity = germline_identity(child)
            r = regions(child)
            rows.append({
                "sequence": child, "predicted_rank_sum": "", "stratum": 0,
                "cdr3": r["CDR3"] if r else "", "group": "candidate",
                "true_rank_sum": "", "repo": "",
                "v_identity": round(identity["V"]["identity_pct"], 2) if identity else "",
                "parent_design_iiptm": round(parent.design_iiptm, 4),
                "parent_cdr3": parent.cdr3,
                # The parent's own sequence, so the pair can be recovered
                # exactly. Matching on a rounded design_iiptm failed: the
                # batch stores 4 decimals and the dataset keeps full
                # precision, so every join came back empty.
                "parent_sequence": parent.sequence, **how})
            made += 1
        print(f"  parent at {parent.design_iiptm:.3f}: {made} variants")
    if refused:
        print("\nrefused by our own gates:")
        for cause, n in refused.most_common(5):
            print(f"  {n:>4}  {cause}")
    if not rows:
        raise SystemExit("no variant cleared the gates")

    # Now the archive rule, once per accepted variant. By construction each
    # differs from its parent by at least the 7 residues the rule requires, so
    # this should pass -- but "should" is not "does", and a duplicate is a
    # certain rejection.
    print(f"\narchive rule on the {len(rows)} accepted variants "
          f"({len(rows) * 5.6 / 60:.0f} min)...")
    kept = []
    for row in rows:
        problems = check_sequence(row["sequence"], args.allowed_v, args.cdr3, archive)
        if problems:
            refused[f"[archive] {problems[0][:34]}"] += 1
        else:
            kept.append(row)
    print(f"  {len(kept)} of {len(rows)} are new")
    rows = kept
    if not rows:
        raise SystemExit("every variant duplicates the archive")
    out = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"\n{len(out)} variants -> {args.out}")
    print(f"  parents' design_iiptm: {out.parent_design_iiptm.min():.3f}-"
          f"{out.parent_design_iiptm.max():.3f}")
    print(f"  V identity of variants: {out.v_identity.min():.1f}-{out.v_identity.max():.1f}"
          f"  (floor {MIN_HUMAN_V_IDENTITY})")
    print("\nRead the result as a delta: each variant against its own parent's")
    print("score. That is what tells whether the minimum legal perturbation")
    print("preserves binding or destroys it.")


if __name__ == "__main__":
    main()
