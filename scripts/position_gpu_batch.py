"""Le test de position, a nombre de substitutions constant.

Le §21 a mesure des ecarts par position du CDR3 allant jusqu a 0,43 et a refuse
de les reporter : les positions etaient tirees au hasard, chaque variante en
touchait 4 a 10, donc les appartenances etaient correlees et confondues avec
"quelle variante s est effondree".

Ici une seule chose change entre les deux bras, et c est l ENSEMBLE des
positions, a 4 substitutions dans les deux cas :

    peripherie  positions 2, 3, 13, 14  -- les bords de la fenetre libre,
                                           contre la tete V et la queue J figees
    milieu      positions 7, 8, 9, 10   -- le centre de la boucle

Lecture decidee d avance, et c est une comparaison de deux echantillons, pas la
recherche d un effet par position : medianes des deux bras, taux
d effondrement des deux bras, et Mann-Whitney sur les scores. Avec les
effectifs disponibles le test voit un effet large (d ~ 1,3) et manquera un
effet modere -- c est dit avant, pas apres.

Si la peripherie preserve et que le milieu detruit, le generateur doit cesser
de tirer ses positions au hasard, et tout le plan du §29 gagne le meme facteur.
Si les deux bras se valent, la position n explique rien de plus que la distance
(§21) et il ne reste que le volume de selection.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

MEASURED = ["minimal2_gpu_scores.csv", "sweep18_gpu_scores.csv",
            "tail18_s4_scores.csv", "tail18b_s4_scores.csv"]
ARMS = {"periph": ("inputs/pos18_periph.csv", "inputs/pos18_periph2.csv"),
        "milieu": ("inputs/pos18_milieu.csv", "inputs/pos18_milieu2.csv")}
OUT = ROOT / "inputs" / "pos18_gpu.csv"
COLUMNS = ["sequence", "role", "arm", "cdr3_positions", "cdr3_substitutions",
           "parent_sequence", "parent_design_iiptm", "previous_design_iiptm", "v_identity"]


def main() -> int:
    measured = pd.concat([pd.read_csv(ROOT / "data" / "measurements" / n) for n in MEASURED],
                         ignore_index=True)
    scored = measured[measured.role == "parent"].drop_duplicates("sequence")
    parent_score = dict(zip(scored.sequence, scored.design_iiptm))

    rows, counts = [], {}
    for arm, paths in ARMS.items():
        frames = [pd.read_csv(ROOT / p) for p in paths if (ROOT / p).exists()]
        if not frames:
            raise SystemExit(f"bras {arm} : aucun fichier genere")
        df = pd.concat(frames, ignore_index=True).drop_duplicates("sequence")
        counts[arm] = len(df)
        for r in df.itertuples():
            rows.append({"sequence": r.sequence, "role": "variant", "arm": arm,
                         "cdr3_positions": r.cdr3_positions,
                         "cdr3_substitutions": r.cdr3_substitutions,
                         "parent_sequence": r.parent_sequence,
                         "parent_design_iiptm": r.parent_design_iiptm,
                         "previous_design_iiptm": "", "v_identity": r.v_identity})

    variants = pd.DataFrame(rows, columns=COLUMNS).drop_duplicates("sequence")
    # un bras dont les positions ne sont pas toutes identiques n est pas un bras
    for arm in ARMS:
        sets = variants[variants.arm == arm].cdr3_positions.unique()
        if len(sets) != 1:
            raise SystemExit(f"bras {arm} : plusieurs ensembles de positions {sets}")

    # Le temoin : un parent deja note, le meilleur disponible. La chaine a ete
    # validee trois fois (0,022 / 0,018 / 0,003), un seul suffit.
    needed = set(variants.parent_sequence)
    unscored = sorted(needed - set(parent_score))
    control = [{"sequence": s, "role": "parent", "arm": "unscored", "cdr3_positions": "",
                "cdr3_substitutions": 0, "parent_sequence": "", "parent_design_iiptm": "",
                "previous_design_iiptm": "", "v_identity": ""} for s in unscored]
    best = max((s for s in needed if s in parent_score), key=lambda s: parent_score[s],
               default=None)
    if best:
        control.append({"sequence": best, "role": "parent", "arm": "recontrol",
                        "cdr3_positions": "", "cdr3_substitutions": 0, "parent_sequence": "",
                        "parent_design_iiptm": parent_score[best],
                        "previous_design_iiptm": parent_score[best], "v_identity": ""})

    batch = pd.concat([pd.DataFrame(control, columns=COLUMNS), variants],
                      ignore_index=True).drop_duplicates("sequence", keep="first")
    readable = set(parent_score) | set(batch[batch.role == "parent"].sequence)
    orphans = sorted(needed - readable)
    if orphans:
        raise SystemExit(f"{len(orphans)} parents sans score mesure et absents du lot")

    batch.to_csv(OUT, index=False)
    print(f"{len(batch)} sequences -> {OUT.relative_to(ROOT)}  "
          f"(~{len(batch) * 1.5 / 60:.1f} h de GPU a 1,5 min/sequence)")
    print(batch.groupby("arm", sort=False).size().to_string())
    for arm, paths in ARMS.items():
        pos = variants[variants.arm == arm].cdr3_positions.iat[0]
        print(f"\n  {arm:<8} positions {pos}, n={counts[arm]}, "
              f"{variants[variants.arm == arm].parent_sequence.nunique()} parents")
    n1, n2 = counts["periph"], counts["milieu"]
    print(f"\n  Puissance : {n1} contre {n2}. Mann-Whitney voit un effet large")
    print("  (d ~ 1,3) et manquera un effet modere. Dit avant, pas apres.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
