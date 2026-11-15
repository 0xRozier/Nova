"""La coupure de `scripts.humanise` porte encore le trou corrige dans `propose`.

`scripts/humanise.py:120` filtre sur la seule colonne `epoch` :

    df = df[df["epoch"].isna() | (df["epoch"] < archive_cutoff)]

C'est mot pour mot la ligne que le commit 7a90b69 a remplacee dans
`scripts/propose_nanobody.py`, et elle est ici sur le chemin de la piste 2
(bras `humanised` des §47 et §49). Le vivier d'`humanise` est precisement la
population que ce commit a mesuree comme entierement future : les donneurs non
dates qui lient a 0,72 ou mieux.

Mesure du 15 novembre, sur les huit epoques rejouees (26299-26306) et a la
longueur de CDR3 que chacune impose : 117 donneurs avec le filtre en place,
**0** avec la coupure corrigee. Tous les 117 sont hors de l'instantane date.

Le premier test epingle le trou (il passe aujourd'hui, et c'est le probleme).
Le second dit ce que la fonction devrait faire ; il est `xfail` jusqu'a ce que
`humanise.donors` passe par `propose_nanobody.dated_archive_snapshot` comme
`propose` le fait maintenant.

    python -m pytest tests/test_humanise_cutoff.py -q
"""
import pandas as pd
import pytest

import scripts.humanise as mod
import scripts.propose_nanobody as propose_mod


@pytest.fixture
def pool(tmp_path, monkeypatch):
    """Quatre donneurs, et un instantane qui n'en date qu'un seul.

    `NEW_ARCHIVE` est une ligne que l'archive publique a gagnee depuis
    l'instantane : elle ne porte aucune date, donc elle ne peut pas etre situee
    avant l'epoque rejouee.
    """
    monkeypatch.setattr(mod, "DATA_DIR", tmp_path)
    monkeypatch.setattr(propose_mod, "DATA_DIR", tmp_path)
    (tmp_path / "dataset").mkdir(parents=True)
    (tmp_path / "prebuilt").mkdir(parents=True)
    pd.DataFrame({"sequence": ["OLD_ARCHIVE"]}).to_csv(
        tmp_path / "prebuilt" / "T_regions.csv", index=False)
    pd.DataFrame({
        "sequence": ["OLD_ARCHIVE", "NEW_ARCHIVE", "SUBMITTED_BEFORE", "SUBMITTED_AFTER"],
        "epoch": [None, None, 26200.0, 26400.0],
        "design_iiptm": [0.73, 0.78, 0.74, 0.80],
        "cdr3": ["AAAAAAAAAAAAA"] * 4,
    }).to_parquet(tmp_path / "dataset" / "T_nanobodies.parquet")
    return tmp_path


def donor_sequences(cutoff):
    return set(mod.donors("T", cdr3_length=None, floor=0.72,
                          archive_cutoff=cutoff)["sequence"])


def test_the_hole_is_still_open(pool):
    """Ce test passe, et c'est exactement le defaut.

    `NEW_ARCHIVE` n'a pas d'epoque, donc la coupure la garde, alors qu'elle a
    pu etre soumise apres l'epoque rejouee -- et sur les vraies donnees, les
    197 lignes de cette forme sont toutes posterieures a l'instantane.
    """
    assert "NEW_ARCHIVE" in donor_sequences(26300), (
        "si cette assertion echoue, humanise.py a ete corrige : "
        "supprimer ce test et retirer le xfail du suivant")


def test_the_cutoff_keeps_what_existed_and_drops_the_rest(pool):
    """Ce que la coupure devrait faire : les deux memes regles que `propose`."""
    assert "SUBMITTED_AFTER" not in donor_sequences(26300), (
        "une soumission datee de l'epoque rejouee ou apres doit tomber")
    assert "SUBMITTED_BEFORE" in donor_sequences(26300)
    assert "OLD_ARCHIVE" in donor_sequences(26300), (
        "une ligne non datee presente dans l'instantane est admissible")


@pytest.mark.xfail(reason="scripts/humanise.py:120 filtre sur la seule colonne "
                          "epoch : une ligne non datee passe n'importe quelle "
                          "coupure. Voir docs/REPLAY_AUDIT.md section 6a.",
                   strict=True)
def test_an_undated_row_the_archive_gained_since_should_be_refused(pool):
    """Le trou lui-meme, exprime comme l'exigence qu'il viole."""
    assert "NEW_ARCHIVE" not in donor_sequences(26300)


def test_without_a_cutoff_nothing_is_dropped(pool):
    """Le mineur en direct n'a pas de fuite par construction : tout ce qui est
    dans le vivier d'aujourd'hui est bien dans le passe."""
    assert donor_sequences(None) == {
        "OLD_ARCHIVE", "NEW_ARCHIVE", "SUBMITTED_BEFORE", "SUBMITTED_AFTER"}
