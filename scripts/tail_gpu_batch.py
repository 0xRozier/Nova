"""Le lot de queue haute, plus un test de reproductibilite de l effondrement.

Le §22 isole la derniere asymetrie non refutee : combien de candidats sont notes
avant de choisir. Le bras a 4 substitutions est le seul dont les survivantes
atteignent 0,707, mais il n en comptait que 4, donc son maximum de 0,718 est le
maximum d un echantillon de 4 et non le plafond de la bande. D ou 37 variantes
a la meme distance, pour lire la queue haute la ou elle compte.

Et une question qui vaut six sequences. L effondrement (design_iiptm ~ 0,13,
interaction_pae > 22) frappe une variante sur deux sans dependre de la distance
(42-56 %, sans tendance) ni du nombre de substitutions (Mann-Whitney p = 0,81).
Or l etape de design de BoltzGen est stochastique -- le carnet porte un
REPLICATE pour cela. Donc : l effondrement est-il une propriete de la sequence
ou du tirage ? Trois survivantes et trois effondrees deja mesurees repassent
ici. Si une effondree re-note haut, le « 20 sur 20 en dessous du parent » du
§19 est en partie du bruit et la selection n a qu a re-noter ; si elle
s effondre encore, c est la sequence.

BoltzGen REPREND un run existant, donc ces six sequences ne redonneraient leurs
anciens resultats que sous l ancien nom de run. Le nom est neuf (tail18_s4) et
le carnet refuse de demarrer si le repertoire du run existe deja.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

MEASURED = [ROOT / "data" / "measurements" / name for name in
            ("minimal2_gpu_scores.csv", "sweep18_gpu_scores.csv")]
VARIANTS = ROOT / "inputs" / "tail18_s4.csv"
OUT = ROOT / "inputs" / "tail18_gpu.csv"
COLLAPSE = 0.30
REPLICATES_PER_SIDE = 3
MAX_SEQUENCES = 70  # ~1h50 de GPU ; on coupe dans les variantes, jamais ailleurs
COLUMNS = ["sequence", "role", "arm", "cdr3_substitutions", "parent_sequence",
           "parent_design_iiptm", "previous_design_iiptm", "v_identity"]


def row(sequence, role, arm, subs=4, parent="", parent_score="", previous="", v=""):
    return {"sequence": sequence, "role": role, "arm": arm,
            "cdr3_substitutions": subs, "parent_sequence": parent,
            "parent_design_iiptm": parent_score, "previous_design_iiptm": previous,
            "v_identity": v}


def main() -> int:
    measured = pd.concat([pd.read_csv(p) for p in MEASURED], ignore_index=True)
    parent_score = dict(zip(measured[measured.role == "parent"].sequence,
                            measured[measured.role == "parent"].design_iiptm))
    scored_variants = measured[measured.role == "variant"].drop_duplicates("sequence")

    variants = pd.read_csv(VARIANTS)
    rows = []

    # 1. les parents : ceux qui n ont pas de score mesure sont obligatoires,
    #    sinon leurs variantes n ont rien a quoi se comparer.
    needed = set(variants.parent_sequence)
    unscored = sorted(needed - set(parent_score))
    for s in unscored:
        rows.append(row(s, "parent", "unscored", 0))
    # un seul re-controle suffit ici : la chaine a deja ete validee deux fois
    # (§19 intra-session a 0,022, §21 inter-session a 0,018).
    recontrol = max((s for s in needed if s in parent_score),
                    key=lambda s: parent_score[s], default=None)
    if recontrol:
        rows.append(row(recontrol, "parent", "recontrol", 0,
                        parent_score=parent_score[recontrol],
                        previous=parent_score[recontrol]))

    # 2. les replicats : trois survivantes et trois effondrees, prises aux
    #    extremes pour que le verdict ne tienne pas a la zone grise.
    collapsed = scored_variants.nsmallest(REPLICATES_PER_SIDE, "design_iiptm")
    survivors = scored_variants.nlargest(REPLICATES_PER_SIDE, "design_iiptm")
    if collapsed.design_iiptm.max() >= COLLAPSE or survivors.design_iiptm.min() < COLLAPSE:
        raise SystemExit("les replicats ne sont pas des deux cotes de la coupure")
    for df, arm in ((survivors, "replicate_survivor"), (collapsed, "replicate_collapsed")):
        for r in df.itertuples():
            rows.append(row(r.sequence, "replicate", arm, r.cdr3_substitutions,
                            parent=r.parent_sequence, previous=r.design_iiptm))

    fixed = pd.DataFrame(rows, columns=COLUMNS)
    budget = MAX_SEQUENCES - len(fixed)
    if len(variants) > budget:
        print(f"  {len(variants)} variantes ramenees a {budget} pour tenir le budget")
        variants = variants.head(budget)
    sweep = pd.DataFrame([
        row(r.sequence, "variant", "s4", 4, r.parent_sequence,
            r.parent_design_iiptm, "", r.v_identity)
        for r in variants.itertuples()
    ], columns=COLUMNS)

    batch = pd.concat([fixed, sweep], ignore_index=True)
    batch = batch.drop_duplicates("sequence", keep="first")
    readable = set(parent_score) | set(batch[batch.role == "parent"].sequence)
    orphans = sorted(set(sweep.parent_sequence) - readable)
    if orphans:
        raise SystemExit(f"{len(orphans)} parents sans score mesure et absents du lot")
    if len(batch) > MAX_SEQUENCES:
        raise SystemExit(f"{len(batch)} sequences, au-dela du budget de {MAX_SEQUENCES}")

    batch.to_csv(OUT, index=False)
    print(f"{len(batch)} sequences -> {OUT.relative_to(ROOT)}  "
          f"(~{len(batch) * 1.6 / 60:.1f} h de GPU a 1,6 min/sequence)")
    print(batch.groupby("arm", sort=False).size().to_string())
    print("\n  replicats, valeur deja mesuree a retrouver ou non :")
    for r in batch[batch.role == "replicate"].itertuples():
        print(f"    {r.arm:<20} {float(r.previous_design_iiptm):.3f}")
    print("\n  Lecture decidee d avance, non parametrique : combien de variantes")
    print("  depassent 0,72 / 0,74 / 0,75, et le maximum observe. Pas d ajustement")
    print("  de loi (le §12bis a produit un 1,042 sur une grandeur bornee par 1).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
