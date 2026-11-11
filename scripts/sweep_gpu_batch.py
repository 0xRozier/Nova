"""Assemble the distance sweep: one variable, the CDR3 substitution count.

FINDINGS 19 measured the minimum legal perturbation (4 substitutions of 18,
CDR3 identity 0.78 to the parent) and found 20 of 20 variants below their
parent. FINDINGS 20 then measured what the epochs' winners actually do: they
sit at 0.667 to their nearest archived ancestor and gain +0.074 over it. The
two numbers only mean something together if nothing else changes, so this
batch holds the same parents and the same pipeline, and moves one dial:

    4 substitutions -> CDR3 identity 0.78   already measured, the control arm
    6                               0.667   the winners' band
    8                               0.556
   10                               0.444

Two rules the batch must satisfy, both asserted rather than assumed:

* every variant's parent has a design_iiptm measured in THIS pipeline -- from
  the FINDINGS 19 run, or measured here for a parent that run never scored.
  Without it a variant is an absolute number with nothing to subtract.
* three parents the FINDINGS 19 run already scored ride along again. That run
  validated the chain to within 0.022, but this is a different session and a
  control is not something one assumes twice.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from nova_lab.nanobodies import cdr_similarity, regions  # noqa: E402

ARMS = (6, 8, 10)
RECONTROLS = 3
PREVIOUS = ROOT / "data" / "measurements" / "minimal2_gpu_scores.csv"
OUT = ROOT / "inputs" / "sweep18_gpu.csv"
COLUMNS = ["sequence", "role", "arm", "cdr3_substitutions", "cdr3_sim_parent",
           "parent_sequence", "parent_design_iiptm", "v_identity"]


def main() -> int:
    previous = pd.read_csv(PREVIOUS)
    scored = previous[previous.role == "parent"].drop_duplicates("sequence")
    measured = dict(zip(scored.sequence, scored.design_iiptm))

    frames = []
    for n in ARMS:
        path = ROOT / "inputs" / f"sweep18_s{n}.csv"
        df = pd.read_csv(path)
        if "cdr3_sim_parent" not in df.columns:
            df["cdr3_sim_parent"] = [
                cdr_similarity(regions(r.sequence), regions(r.parent_sequence))["CDR3"]
                for r in df.itertuples()
            ]
            df.to_csv(path, index=False)
        df["arm"] = f"s{n}"
        df["cdr3_substitutions"] = n
        df["role"] = "variant"
        frames.append(df[COLUMNS])
    sweep = pd.concat(frames, ignore_index=True)

    # Parents: the best already-scored ones as a re-control, plus every parent
    # the previous run never scored -- those are not optional, their variants
    # are unreadable without them.
    needed = set(sweep.parent_sequence)
    unscored = sorted(needed - set(measured))
    recontrol = [s for s in scored.nlargest(RECONTROLS, "design_iiptm").sequence if s in needed]
    rows = [
        {"sequence": s, "role": "parent", "arm": "recontrol", "cdr3_substitutions": 0,
         "cdr3_sim_parent": 1.0, "parent_sequence": "",
         "parent_design_iiptm": measured[s], "v_identity": ""}
        for s in recontrol
    ] + [
        {"sequence": s, "role": "parent", "arm": "unscored", "cdr3_substitutions": 0,
         "cdr3_sim_parent": 1.0, "parent_sequence": "",
         "parent_design_iiptm": "", "v_identity": ""}
        for s in unscored
    ]
    batch = pd.concat([pd.DataFrame(rows, columns=COLUMNS), sweep], ignore_index=True)

    if not batch.sequence.is_unique:
        raise SystemExit("sequence en double : BoltzGen ne la noterait qu une fois")
    readable = set(measured) | set(batch[batch.role == "parent"].sequence)
    orphans = sorted(needed - readable)
    if orphans:
        raise SystemExit(
            f"{len(orphans)} parents sans score mesure et absents du lot : "
            "leurs variantes seraient des nombres absolus sans rien a soustraire"
        )

    batch.to_csv(OUT, index=False)
    print(f"{len(batch)} sequences -> {OUT.relative_to(ROOT)}  "
          f"(~{len(batch) * 1.52 / 60:.1f} h de GPU a 1,5 min/sequence)")
    print(batch.groupby("arm", sort=False)
          .agg(n=("sequence", "size"), identite_CDR3_au_parent=("cdr3_sim_parent", "median"))
          .to_string())
    print(f"\n  {len(recontrol)} parents en re-controle (deja mesures : "
          f"{', '.join(f'{measured[s]:.3f}' for s in recontrol)})")
    print(f"  {len(unscored)} parent(s) jamais mesure(s), embarque(s) pour que leurs "
          "variantes soient lisibles")
    print("\n  bras temoin a 4 substitutions (identite 0,78) : deja mesure au §19,")
    print("  il ne coute pas une minute de GPU ici.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
