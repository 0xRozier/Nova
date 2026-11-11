"""Parent -> variant deltas on everything measurable without a GPU.

The minimal variants of `minimal_variants.py` exist to answer one question:
does the smallest perturbation the uniqueness rule allows preserve binding?
That answer costs a GPU session. What costs nothing is the rest of the
validator's verdict -- PSH, Aggrescan3D, V and J identity -- and a variant
whose developability collapses against its own parent is not submittable
whatever its design_iiptm turns out to be. So we price that first.

PSH comes from the preflight logs (3 folds per sequence, already paid for).
Aggrescan3D is recovered from the preflight residue cache rather than
re-folded: `data/preflight/run/a3d{1,2,3}/residues/<sha256>/residue-scores.csv`
holds the per-residue scores, and the average and the peak are the mean and
the maximum of that column -- verified against `summary/aggrescan3d.csv` to
seven decimals. V and J identity are recomputed with ANARCI, which is cheap.

The parent of each variant is carried in the batch CSV (`parent_sequence`),
because matching on a rounded design_iiptm gave empty joins.
"""

from __future__ import annotations

import csv
import hashlib
import re
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from nova_lab.nanobodies import germline_identity  # noqa: E402

A3D_AVERAGE_CAP = -0.680
A3D_PEAK_CAP = 1.580
PSH_BREAK = 96.0  # measured in FINDINGS 3ter, not the advertised 89

PREFLIGHT_RUN = ROOT / "data" / "preflight" / "run"
MEASUREMENTS = ROOT / "data" / "measurements"
BATCHES = [
    ("IGHV3-74", 16),
    ("IGHV3-74", 18),
    ("IGHV3-66", 12),
    ("IGHV3-64", 18),
]

PSH_LINE = re.compile(r"PSH ([\d.]+) over (\d+) folds")
HEADER = re.compile(r"^### (?:V2|PARENTS) (IGHV3-\d+) CDR3=(\d+)")


def parse_preflight(path: Path) -> dict[tuple[str, int], list[tuple[float | None, str]]]:
    """Group -> per-sequence verdicts, in the order the pipeline screened them.

    The logs print only the first 34 residues, so the join is positional: the
    screening walks the batch CSV top to bottom. Each entry is (psh, verdict),
    with psh None for a sequence the deterministic gates or A3D refused.
    """
    out: dict[tuple[str, int], list[tuple[float | None, str]]] = {}
    key = None
    for line in path.read_text().splitlines():
        head = HEADER.match(line)
        if head:
            key = (head.group(1), int(head.group(2)))
            out[key] = []
            continue
        if key is None or "..." not in line:
            continue
        if "REFUSED" in line:
            out[key].append((None, line.split("REFUSED — ", 1)[-1]))
            continue
        hit = PSH_LINE.search(line)
        out[key].append((float(hit.group(1)), "gates pass") if hit else (None, line))
    return out


def a3d_from_cache(sequence: str) -> tuple[float | None, float | None, int]:
    """Average and peak over whichever of the three fold replicates are cached."""
    digest = hashlib.sha256(sequence.encode()).hexdigest()
    averages, peaks = [], []
    for replicate in ("a3d1", "a3d2", "a3d3"):
        scores = PREFLIGHT_RUN / replicate / "residues" / digest / "residue-scores.csv"
        if not scores.exists():
            continue
        values = [float(r["score_full_precision"]) for r in csv.DictReader(scores.open())]
        if not values:
            continue
        averages.append(statistics.fmean(values))
        peaks.append(max(values))
    if not averages:
        return None, None, 0
    return statistics.fmean(averages), statistics.fmean(peaks), len(averages)


def identities(sequence: str) -> tuple[float | None, float | None]:
    germline = germline_identity(sequence)
    if germline is None:
        return None, None
    return germline["V"]["identity_pct"], germline["J"]["identity_pct"]


def fmt(value, spec="+.1f"):
    return "    -" if value is None else format(value, spec)


def main() -> int:
    variant_psh = parse_preflight(MEASUREMENTS / "minimal2_preflight.log")
    parent_psh = parse_preflight(MEASUREMENTS / "parents_preflight.log")

    parent_measured: dict[str, float] = {}
    for gene, length in BATCHES:
        sequences = (MEASUREMENTS / f"parents_{gene}_{length}.txt").read_text().split()
        verdicts = parent_psh[(gene, length)]
        if len(sequences) != len(verdicts):
            raise SystemExit(
                f"parents {gene}/{length}: {len(sequences)} sequences, {len(verdicts)} verdicts"
            )
        for sequence, (psh, _) in zip(sequences, verdicts):
            if psh is not None:
                parent_measured[sequence] = psh

    rows = []
    for gene, length in BATCHES:
        batch = list(csv.DictReader((ROOT / "inputs" / f"minimal2_{gene}_{length}.csv").open()))
        verdicts = variant_psh[(gene, length)]
        if len(batch) != len(verdicts):
            raise SystemExit(
                f"{gene}/{length}: {len(batch)} variants, {len(verdicts)} verdicts -- "
                "the positional join is not valid, do not read the deltas"
            )
        for row, (psh, verdict) in zip(batch, verdicts):
            variant, parent = row["sequence"], row["parent_sequence"]
            v_avg, v_peak, folds = a3d_from_cache(variant)
            p_avg, p_peak, _ = a3d_from_cache(parent)
            v_v, v_j = identities(variant)
            p_v, p_j = identities(parent)
            rows.append(
                {
                    "group": f"{gene}/{length}",
                    "sequence": variant,
                    "parent_sequence": parent,
                    "parent_design_iiptm": float(row["parent_design_iiptm"]),
                    "verdict": verdict,
                    "psh": psh,
                    "parent_psh": parent_measured.get(parent),
                    "a3d_average": v_avg,
                    "a3d_peak": v_peak,
                    "a3d_folds": folds,
                    "parent_a3d_average": p_avg,
                    "parent_a3d_peak": p_peak,
                    "v_identity": v_v,
                    "j_identity": v_j,
                    "parent_v_identity": p_v,
                    "parent_j_identity": p_j,
                    "cdr3_substitutions": int(row["cdr3_substitutions"]),
                    "framework_moves": int(row["framework_moves"]),
                }
            )

    for row in rows:
        row["d_psh"] = None if None in (row["psh"], row["parent_psh"]) else row["psh"] - row["parent_psh"]
        row["d_a3d_average"] = (
            None
            if None in (row["a3d_average"], row["parent_a3d_average"])
            else row["a3d_average"] - row["parent_a3d_average"]
        )
        row["d_a3d_peak"] = (
            None
            if None in (row["a3d_peak"], row["parent_a3d_peak"])
            else row["a3d_peak"] - row["parent_a3d_peak"]
        )
        row["d_v"] = None if None in (row["v_identity"], row["parent_v_identity"]) else row["v_identity"] - row["parent_v_identity"]
        row["d_j"] = None if None in (row["j_identity"], row["parent_j_identity"]) else row["j_identity"] - row["parent_j_identity"]
        row["submittable"] = bool(
            row["psh"] is not None
            and row["psh"] < PSH_BREAK
            and row["a3d_average"] is not None
            and row["a3d_average"] <= A3D_AVERAGE_CAP
            and row["a3d_peak"] is not None
            and row["a3d_peak"] <= A3D_PEAK_CAP
        )

    out = ROOT / "data" / "measurements" / "minimal2_deltas.csv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"{len(rows)} variants -> {out.relative_to(ROOT)}\n")
    print("group        n   soumissibles   dPSH median   dA3D avg   dA3D peak   dV    dJ")
    for gene, length in BATCHES:
        group = [r for r in rows if r["group"] == f"{gene}/{length}"]
        kept = [r for r in group if r["submittable"]]
        def median(field):
            values = [r[field] for r in group if r[field] is not None]
            return statistics.median(values) if values else None
        print(
            f"{gene}/{length:<4} {len(group):>3}   {len(kept):>12}   "
            f"{fmt(median('d_psh')):>11}   {fmt(median('d_a3d_average'), '+.3f'):>8}   "
            f"{fmt(median('d_a3d_peak'), '+.3f'):>9}   {fmt(median('d_v')):>5} {fmt(median('d_j')):>5}"
        )

    kept = [r for r in rows if r["submittable"]]
    paired = [r for r in rows if r["d_psh"] is not None]
    print(
        f"\nensemble : {len(kept)}/{len(rows)} soumissibles  "
        f"(PSH < {PSH_BREAK:g} et A3D sous les deux plafonds)"
    )
    if paired:
        worse = sum(r["d_psh"] > 0 for r in paired)
        print(
            f"dPSH apparie (n={len(paired)}) : median {statistics.median(r['d_psh'] for r in paired):+.1f}, "
            f"{worse}/{len(paired)} degrades, pire {max(r['d_psh'] for r in paired):+.1f}, "
            f"meilleur {min(r['d_psh'] for r in paired):+.1f}"
        )
    refused = [r for r in rows if r["verdict"] != "gates pass"]
    if refused:
        print(f"\nrefusees aux portes deterministes / A3D : {len(refused)}")
        for row in refused:
            print(f"  {row['group']:<13} {row['verdict']}")

    print("\nretenues pour la mesure GPU, par parent :")
    for row in sorted(kept, key=lambda r: (r["group"], r["psh"])):
        print(
            f"  {row['group']:<13} PSH {row['psh']:6.1f} ({fmt(row['d_psh'])})  "
            f"A3D {row['a3d_average']:+.3f}/{row['a3d_peak']:.3f}  "
            f"V {row['v_identity']:.1f} J {row['j_identity']:.1f}  "
            f"parent iiptm {row['parent_design_iiptm']:.3f}"
        )
    print(f"\nparents distincts a embarquer comme controles : {len({r['parent_sequence'] for r in kept})}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
