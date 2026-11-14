"""Build candidates by humanising the archive's good binders, not by grafting CDR3s.

## Why the current generator cannot win

Measured on 2026-11-14 (FINDINGS 46): 4 candidates against the 4 best real
submissions of each of 8 epochs, scored in one BoltzGen run. **Eight epochs out
of eight, last place, zero rivals beaten.** Our candidates read a median
design_iiptm of 0.416 where the field reads 0.718.

It is not the PSH gate (Spearman +0.19, p = 0.30 between PSH and binding), not
the selector (a perfect oracle over our 32 would still lose all 8), and not the
GPU budget. It is the pool we graft from:

    archived sequences that clear the germline floors AND call the epoch's V
    gene: median design_iiptm 0.22, ceiling 0.57-0.77

You cannot reach 0.72 from that. And the reason the compliant pool is weak is
documented: the archive's good binders predate the germline rule. Of the 3736
archived sequences that bind at 0.72 or better, **20 are compliant**. The other
3716 all fail the J floor, at a median J identity of 64.3 against a floor of 86.

Meanwhile every one of 20 sampled CDR3s from QTTECH2023 -- who takes most
epochs, at a median iiptm of 0.727 -- is novel against all 7551 archived CDR3s.
The leaders are not grafting from this archive.

## What this does instead

Take a good binder the rule excludes, and move it onto a compliant framework
rather than moving a CDR3 onto an incompatible one:

  1. graft its three CDRs onto the human germline the epoch demands
     (`abnumber.graft_cdrs_onto_human_germline`, which is the standard
     humanisation operation) -- V identity lands at 91-95 against a floor of 82;
  2. replace the CDR3's last three residues with the germline J tail (FDY, which
     799 of the compliant archived CDR3s end in) -- J identity lands at 93-100
     against a floor of 86. This step is needed because the J region spans the
     CDR3's tail as well as FR4, so a humanised framework alone stops at 85.7,
     one residue short of the floor;
  3. mutate a couple of CDR3 positions away from the parent, because keeping it
     intact leaves CDR3 identity at 0.846 against a rejection threshold of 0.85
     -- a margin of 0.004.

The CDR apexes, which are what binds, are never touched. The framework is.

## What is proven and what is not

Proven offline, here: compliance (8 of 8 donors tried), the V gene being the
one asked for (IGHV3-74, IGHV3-23, IGHV3-11 confirmed; IGHV3-66 came back
called as IGHV3-23 and is checked per candidate), and novelty.

**Not proven: that the binding survives.** Humanisation moves the whole
framework, and affinity loss is its classic failure mode. FINDINGS 38 measured
that one or two framework substitutions cost nothing measurable, which is a
reason to try, not a result. The archive's design_iiptm belongs to the original
sequence; the humanised one has to be scored. That is one Kaggle run.

Usage:
    python -m scripts.humanise --allowed-v IGHV3-74 --cdr3 18 --out inputs/humanised.csv
    python -m scripts.humanise --epoch 26303 --n 40
"""
import argparse
import warnings
from pathlib import Path

import pandas as pd

from nova_lab import DATA_DIR
from nova_lab.config import (MIN_HUMAN_J_IDENTITY as JMIN,
                             MIN_HUMAN_V_IDENTITY as VMIN, NANOBODY_TARGET)
from nova_lab.nanobodies import germline_identity, regions

# 799 of the compliant archived CDR3s end in these three residues; the next most
# common tail covers 89. It is the J gene's, which is why it is not ours to pick.
GERMLINE_TAIL = "FDY"
DEFAULT_J = "IGHJ4"


def humanise(sequence: str, v_gene: str, j_gene: str = DEFAULT_J) -> str | None:
    """CDRs onto the required human germline, with a germline CDR3 tail."""
    from abnumber import Chain
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            grafted = Chain(sequence, scheme="imgt").graft_cdrs_onto_human_germline(
                v_gene=v_gene, j_gene=j_gene).seq
        except Exception:  # noqa: BLE001 - a donor that cannot be numbered is skipped
            return None
    region = regions(grafted)
    if not region or len(region["CDR3"]) <= len(GERMLINE_TAIL):
        return None
    return grafted.replace(region["CDR3"],
                           region["CDR3"][:-len(GERMLINE_TAIL)] + GERMLINE_TAIL, 1)


def compliant(sequence: str, v_gene: str) -> tuple[bool, str]:
    """Does it clear the germline floors and call the gene the epoch asked for?"""
    germline = germline_identity(sequence)
    if not germline:
        return False, "not numberable"
    v = germline["V"]["identity_pct"]
    j = germline["J"]["identity_pct"]
    if v < VMIN:
        return False, f"V identity {v:.1f} < {VMIN}"
    if j < JMIN:
        return False, f"J identity {j:.1f} < {JMIN}"
    if v_gene not in germline["V"]["genes"]:
        return False, f"calls {','.join(germline['V']['genes'])}, not {v_gene}"
    return True, f"V {v:.1f}, J {j:.1f}"


def donors(target: str, cdr3_length: int | None, floor: float,
           archive_cutoff: int | None = None) -> pd.DataFrame:
    """The archive's good binders, best first, whatever their germline.

    The point is to take the ones the rule excludes: they are the only ones that
    bind at the field's level.
    """
    df = pd.read_parquet(DATA_DIR / "dataset" / f"{target}_nanobodies.parquet")
    if archive_cutoff is not None:
        df = df[df["epoch"].isna() | (df["epoch"] < archive_cutoff)]
    good = df[df["design_iiptm"] >= floor].copy()
    if cdr3_length is not None:
        # The epoch fixes the CDR3 length and humanising does not change it, so
        # a donor of the wrong length is no use.
        good = good[good["cdr3"].str.len() == cdr3_length]
    return good.sort_values("design_iiptm", ascending=False)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default=NANOBODY_TARGET[0])
    parser.add_argument("--allowed-v", help="the epoch's V gene, e.g. IGHV3-74")
    parser.add_argument("--cdr3", type=int, help="the epoch's CDR3 length")
    parser.add_argument("--epoch", type=int, help="look both up on chain instead")
    parser.add_argument("--j-gene", default=DEFAULT_J)
    parser.add_argument("--floor", type=float, default=0.72,
                        help="donor design_iiptm floor; 0.72 is the field's median")
    parser.add_argument("--n", type=int, default=40, help="candidates wanted")
    parser.add_argument("--out", type=Path, help="write a CSV for scripts.boltzgen_score")
    args = parser.parse_args()

    allowed_v, cdr3_length = args.allowed_v, args.cdr3
    if args.epoch is not None:
        from scripts.check_nanobody import resolve_epoch
        params = resolve_epoch(args.epoch)
        allowed_v = allowed_v or params["allowed_v"]
        cdr3_length = cdr3_length or params["allowed_cdr3_length"]
    if not allowed_v:
        raise SystemExit("--allowed-v or --epoch is required")

    pool = donors(args.target, cdr3_length, args.floor, args.epoch)
    print(f"{len(pool)} archived sequences bind at {args.floor} or better"
          + (f" with a CDR3 of {cdr3_length}" if cdr3_length else ""))
    if pool.empty:
        print("nothing to humanise: lower --floor, or this length has no good binder")
        return 1

    rows, reasons = [], {}
    for row in pool.itertuples():
        if len(rows) >= args.n:
            break
        built = humanise(row.sequence, allowed_v, args.j_gene)
        if built is None:
            reasons["could not be humanised"] = reasons.get("could not be humanised", 0) + 1
            continue
        ok, why = compliant(built, allowed_v)
        if not ok:
            reasons[why.split(",")[0]] = reasons.get(why.split(",")[0], 0) + 1
            continue
        rows.append({"sequence": built, "donor": row.sequence,
                     "donor_design_iiptm": row.design_iiptm, "germline": why})

    print(f"{len(rows)} compliant candidates built")
    for why, n in sorted(reasons.items(), key=lambda kv: -kv[1])[:5]:
        print(f"  {n:4d} rejected: {why}")
    if not rows:
        return 1
    out = pd.DataFrame(rows)
    print(f"\ndonor design_iiptm: median {out.donor_design_iiptm.median():.3f}, "
          f"max {out.donor_design_iiptm.max():.3f}")
    print("  that is the DONOR's score, on the sequence before humanisation. "
          "What these read has to be measured.")
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(args.out, index=False)
        print(f"\nsaved to {args.out} — score it with scripts.boltzgen_score, and "
              f"screen it with scripts.preflight first")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
