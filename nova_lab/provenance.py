"""Quelle cible a produit chaque fichier de scores ? Et le refus de melanger.

Les CSV de `results/` n'enregistrent pas la cible contre laquelle ils ont ete
notes. Tant qu'il n'y avait qu'une cible, l'omission etait invisible. La bascule
du 20 novembre 2026 (Q9NZQ7 / PD-L1 / 8AOK -> P05231 / IL-6 / 4O9H) l'a rendue
couteuse : le lot `scramble` (FINDINGS 90) a tire ses dix donneurs « hauts »
d'une liste de fichiers nommee `LIVE_TARGET_RUNS`, ecrite avant la bascule et
devenue fausse apres. Les sequences choisies pour leur iiptm de 0,715-0,795
notaient 0,145-0,224 sur la cible vivante, et 40 des 56 sequences du lot ont ete
payees sur une premisse vide.

**La date de commit ne dit pas la cible** : `scores_nanobody_xtal_26301.csv` est
commite le 23 novembre, trois jours apres la bascule, mais note sur l'ancienne
cible. C'est pourquoi la table ci-dessous est explicite, et appuyee sur deux
signatures qui concordent sur les quinze fichiers :

1. la colonne `epoch`. Tous les runs de l'ancienne cible en portent une (25996 a
   26347) ; aucun run de la cible vivante n'en a ;
2. l'echelle de l'`interaction_pae`. Ancienne cible : mediane 10,9 a 21,9 et
   groupe `real` a 8,8-15,0. Cible vivante : mediane 23,2 a 24,2, et les
   sequences du champ qui notaient 0,70+ sur PD-L1 tombent dans la masse.

Les runs a venir n'ont plus besoin de cette table : `scripts.boltzgen_score`
ecrit desormais une colonne `target`, qui est lue en priorite.

Usage :
    from nova_lab.provenance import assert_current_target
    assert_current_target(paths)   # leve avant de depenser du GPU
"""
from pathlib import Path

import pandas as pd

from nova_lab import REPO_ROOT
from nova_lab.config import NANOBODY_TARGET, PREVIOUS_NANOBODY

CURRENT = NANOBODY_TARGET[0]
PREVIOUS = PREVIOUS_NANOBODY["target"]

# Etabli par les deux signatures du docstring, fichier par fichier. Un fichier
# absent d'ici est inconnu, pas suppose a jour : `target_of` rend None et
# `assert_current_target` refuse.
RUN_TARGET: dict[str, str] = {
    # Ancienne cible (Q9NZQ7, PD-L1, 8AOK) : colonne epoch presente, pae basse
    "scores_nanobody_arms_26299_26306.csv": PREVIOUS,
    "scores_nanobody_deep_26304_26345_26347.csv": PREVIOUS,
    "scores_nanobody_design3arms_26345.csv": PREVIOUS,
    "scores_nanobody_k1sample_26345.csv": PREVIOUS,
    "scores_nanobody_wide4arms_26345.csv": PREVIOUS,
    "scores_nanobody_repaired_26299_26306.csv": PREVIOUS,
    "scores_nanobody_observation_combined.csv": PREVIOUS,
    "scores_nanobody_observation_cutoff.csv": PREVIOUS,
    "scores_nanobody_observation_depth4.csv": PREVIOUS,
    # Commite le 23 novembre, mais note sur l'ancienne cible : epoch 26301, et
    # ses quatre `real` a 0,700 / pae 9,19. La date de commit ne dit rien.
    "scores_nanobody_xtal_26301.csv": PREVIOUS,
    # Sans colonne epoch, mais commites le 16 novembre : la cible vivante
    # n'existait pas encore dans la config, donc ils ne peuvent pas la porter.
    # C'est le seul sens ou la date de commit prouve quelque chose.
    "scores_nanobody_nanobody_batch3.csv": PREVIOUS,
    "scores_nanobody_nanobody_batch4.csv": PREVIOUS,
    # Cible vivante (P05231, IL-6, 4O9H) : pas de colonne epoch, pae 23-24
    "scores_nanobody_calibration.csv": CURRENT,
    "scores_nanobody_learning.csv": CURRENT,
    "scores_nanobody_transfer.csv": CURRENT,
    "scores_nanobody_ptm.csv": CURRENT,
    "scores_nanobody_scramble.csv": CURRENT,
}


def target_of(path: str | Path) -> str | None:
    """La cible d'un fichier de scores, ou None si elle n'est pas etablie.

    La colonne `target` du fichier gagne sur la table : les runs posterieurs a
    FINDINGS 90 l'ecrivent eux-memes.
    """
    path = Path(path)
    try:
        cols = pd.read_csv(path, nrows=1).columns
        if "target" in cols:
            values = set(pd.read_csv(path, usecols=["target"])["target"].dropna())
            if len(values) == 1:
                return str(values.pop())
            if len(values) > 1:
                raise ValueError(f"{path.name} melange plusieurs cibles : {sorted(values)}")
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return None
    return RUN_TARGET.get(path.name)


def assert_current_target(paths, what: str = "ces fichiers de scores") -> None:
    """Lever si l'un des fichiers n'a pas ete note sur la cible vivante.

    A appeler dans tout constructeur de lot qui puise des sequences, des scores
    ou des ancres dans l'historique. Le transfert d'une cible a l'autre vaut
    rho = +0,003 pour les metriques d'interface (FINDINGS 80) : un score de
    l'ancienne cible ne dit rien du score sur celle-ci, et ne vaut donc pas
    d'etre interroge sur GPU.
    """
    wrong, unknown = [], []
    for p in paths:
        p = Path(p)
        if not p.is_absolute():
            p = REPO_ROOT / p
        t = target_of(p)
        if t is None:
            unknown.append(p.name)
        elif t != CURRENT:
            wrong.append(f"{p.name} ({t})")
    if wrong or unknown:
        msg = [f"{what} ne portent pas tous la cible vivante {CURRENT} :"]
        if wrong:
            msg.append(f"  notes sur une autre cible : {', '.join(sorted(wrong))}")
        if unknown:
            msg.append(f"  cible inconnue : {', '.join(sorted(unknown))} — "
                       f"l'etablir et l'inscrire dans nova_lab/provenance.RUN_TARGET")
        raise RuntimeError("\n".join(msg))


def runs_on_current_target() -> list[str]:
    """Les fichiers de results/ notes sur la cible vivante, verifies sur disque."""
    out = []
    for name, target in sorted(RUN_TARGET.items()):
        if target == CURRENT and (REPO_ROOT / "results" / name).exists():
            out.append(name)
    return out
