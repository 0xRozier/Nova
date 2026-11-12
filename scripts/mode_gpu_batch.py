"""Positions imposees contre positions tirees, DANS LA MEME SESSION.

Le §31 avait observe 0 effondrement sur 22 variantes a positions imposees,
contre 22 % historiquement a positions tirees. Le §32 a montre que ce 0/22 est
confondu avec la session : le taux d effondrement varie de 0 % a 31 % d une
session a l autre (chi2 p = 0,0029). Les deux explications -- mode de generation
et session -- sont donc inseparables dans les donnees existantes, et le seul
dispositif qui les separe met les deux modes dans le meme lot.

Une seule variable : les positions sont tirees, ou imposees a 2,3,13,14. Memes
parents, meme nombre de substitutions (4), meme chaine, meme session.

Et trois replicats, parce que le §28 a montre ce que vaut un maximum non
replique : le 0,770 du §26 a re-note 0,756 et les quatre meilleures ont toutes
baisse. Les trois meilleures du §31 (0,783 / 0,757 / 0,735) repassent donc ici.
Le nom de run doit etre neuf, sans quoi BoltzGen rendrait leurs anciens
resultats au lieu de recalculer -- ce qui annulerait precisement ce controle.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from nova_lab.nanobodies import regions  # noqa: E402

MEASURED = ["minimal2_gpu_scores.csv", "sweep18_gpu_scores.csv", "tail18_s4_scores.csv",
            "tail18b_s4_scores.csv", "pos18_scores.csv"]
IMPOSED_POSITIONS = (2, 3, 13, 14)
REPLICATES = 3
OUT = ROOT / "inputs" / "mode18_gpu.csv"
COLUMNS = ["sequence", "role", "arm", "cdr3_positions", "cdr3_substitutions",
           "parent_sequence", "parent_design_iiptm", "previous_design_iiptm", "v_identity"]


def moved_positions(sequence: str, parent: str) -> tuple[int, ...] | None:
    a, b = regions(sequence), regions(parent)
    if not a or not b or len(a["CDR3"]) != len(b["CDR3"]):
        return None
    return tuple(i for i, (x, y) in enumerate(zip(a["CDR3"], b["CDR3"])) if x != y)


def main() -> int:
    measured = pd.concat([pd.read_csv(ROOT / "data" / "measurements" / n) for n in MEASURED],
                         ignore_index=True)
    scored_parents = measured[measured.role == "parent"].drop_duplicates("sequence")
    parent_score = dict(zip(scored_parents.sequence, scored_parents.design_iiptm))
    already = set(measured.sequence)

    rows = []
    for arm, path in (("aleatoire", "inputs/rand18_c.csv"),
                      ("impose", "inputs/pos18_periph2.csv")):
        df = pd.read_csv(ROOT / path).drop_duplicates("sequence")
        df = df[~df.sequence.isin(already)]          # rien de deja note
        sets = set()
        for r in df.itertuples():
            pos = moved_positions(r.sequence, r.parent_sequence)
            if pos is None or len(pos) != 4:
                continue
            sets.add(pos)
            rows.append({"sequence": r.sequence, "role": "variant", "arm": arm,
                         "cdr3_positions": ",".join(map(str, pos)),
                         "cdr3_substitutions": 4, "parent_sequence": r.parent_sequence,
                         "parent_design_iiptm": r.parent_design_iiptm,
                         "previous_design_iiptm": "", "v_identity": r.v_identity})
        # le bras impose doit avoir UN seul ensemble, l aleatoire plusieurs :
        # sinon ce ne sont pas deux modes differents et le lot ne teste rien.
        if arm == "impose" and sets != {IMPOSED_POSITIONS}:
            raise SystemExit(f"bras impose : ensembles {sorted(sets)}, attendu "
                             f"{[IMPOSED_POSITIONS]}")
        if arm == "aleatoire" and len(sets) < 5:
            raise SystemExit(f"bras aleatoire : seulement {len(sets)} ensembles distincts, "
                             "ce n est pas un tirage")

    variants = pd.DataFrame(rows, columns=COLUMNS)
    counts = variants.arm.value_counts().to_dict()

    # les replicats : les trois meilleures variantes du §31
    pos18 = pd.read_csv(ROOT / "data" / "measurements" / "pos18_scores.csv")
    top = pos18[pos18.role == "variant"].nlargest(REPLICATES, "design_iiptm")
    fixed = [{"sequence": r.sequence, "role": "replicate", "arm": "replicate_top",
              "cdr3_positions": r.cdr3_positions, "cdr3_substitutions": 4,
              "parent_sequence": r.parent_sequence,
              "parent_design_iiptm": r.parent_design_iiptm,
              "previous_design_iiptm": r.design_iiptm, "v_identity": ""}
             for r in top.itertuples()]
    best = max(parent_score, key=parent_score.get)
    fixed.append({"sequence": best, "role": "parent", "arm": "recontrol",
                  "cdr3_positions": "", "cdr3_substitutions": 0, "parent_sequence": "",
                  "parent_design_iiptm": parent_score[best],
                  "previous_design_iiptm": parent_score[best], "v_identity": ""})

    batch = pd.concat([pd.DataFrame(fixed, columns=COLUMNS), variants],
                      ignore_index=True).drop_duplicates("sequence", keep="first")
    orphans = sorted(set(variants.parent_sequence) - set(parent_score)
                     - set(batch[batch.role == "parent"].sequence))
    if orphans:
        raise SystemExit(f"{len(orphans)} parents sans score mesure et absents du lot")

    batch.to_csv(OUT, index=False)
    print(f"{len(batch)} sequences -> {OUT.relative_to(ROOT)}  "
          f"(~{len(batch) * 1.6 / 60:.1f} h de GPU)")
    print(batch.groupby("arm", sort=False).size().to_string())
    print(f"\n  bras aleatoire {counts.get('aleatoire', 0)} contre "
          f"bras impose {counts.get('impose', 0)}, memes parents, 4 substitutions")
    print("  replicats, valeur du §31 a retrouver ou non :")
    for r in batch[batch.role == "replicate"].itertuples():
        print(f"    {float(r.previous_design_iiptm):.5f}")
    print("\n  Predictions enregistrees au §32 :")
    print("   1. les positions imposees suppriment l effondrement -> impose ~0%, aleatoire ~20%")
    print("   2. c etait la session -> les deux bras au meme niveau, et 22% du budget restent perdus")
    print("   3. le 0,783 baisse de ~0,015 et atterrit vers 0,768")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
