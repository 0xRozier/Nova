"""Assemble the GPU batch: the submittable minimal variants and their parents.

The parents are the point. Each one has an archived design_iiptm measured in an
earlier Kaggle session; if they do not re-score near that value here, the run
measured the session rather than the sequences and the variant numbers mean
nothing. So every parent of a kept variant rides in the same batch, and each
row carries the archived value to compare against.

Rows are ordered parent-first, then that parent's variants, so `TOP = 2` in the
notebook scores one pair and nothing else -- a useful smoke test.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DELTAS = ROOT / "data" / "measurements" / "minimal2_deltas.csv"
OUT = ROOT / "inputs" / "minimal2_gpu.csv"

FIELDS = [
    "sequence",
    "role",
    "group",
    "parent_sequence",
    "parent_design_iiptm",
    "psh",
    "d_psh",
    "a3d_average",
    "a3d_peak",
    "v_identity",
    "j_identity",
]


def main() -> int:
    rows = [r for r in csv.DictReader(DELTAS.open()) if r["submittable"] == "True"]
    if not rows:
        raise SystemExit(f"{DELTAS}: no submittable variant, nothing to measure")

    by_parent: dict[str, list[dict]] = {}
    for row in rows:
        by_parent.setdefault(row["parent_sequence"], []).append(row)

    out = []
    for parent, variants in sorted(by_parent.items(), key=lambda kv: -float(kv[1][0]["parent_design_iiptm"])):
        first = variants[0]
        out.append(
            {
                "sequence": parent,
                "role": "parent",
                "group": first["group"],
                "parent_sequence": "",
                "parent_design_iiptm": first["parent_design_iiptm"],
                "psh": first["parent_psh"],
                "d_psh": "",
                "a3d_average": first["parent_a3d_average"],
                "a3d_peak": first["parent_a3d_peak"],
                "v_identity": first["parent_v_identity"],
                "j_identity": first["parent_j_identity"],
            }
        )
        for variant in sorted(variants, key=lambda r: float(r["psh"])):
            out.append({"role": "variant", **{k: variant.get(k, "") for k in FIELDS if k != "role"}})

    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(out)

    sequences = [r["sequence"] for r in out]
    if len(set(sequences)) != len(sequences):
        raise SystemExit("duplicate sequence in the batch: BoltzGen would score it once")

    parents = [r for r in out if r["role"] == "parent"]
    variants = [r for r in out if r["role"] == "variant"]
    print(f"{len(out)} sequences -> {OUT.relative_to(ROOT)}")
    print(f"  {len(variants)} variantes minimales, {len(parents)} parents comme controles internes")
    print(f"  ~{len(out) * 1.52 / 60:.1f} h de GPU a 1,5 min/sequence (mesure sur 84)")
    print("\n  parent (iiptm archive)   variantes   PSH des variantes")
    for parent in parents:
        kin = [v for v in variants if v["parent_sequence"] == parent["sequence"]]
        psh = ", ".join(f"{float(v['psh']):.1f}" for v in kin)
        print(f"  {parent['group']:<13} {float(parent['parent_design_iiptm']):.3f}   {len(kin):>9}   {psh}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
