"""Tests du lot a budget egal.

Deux choses a garantir. Que le lot, s'il existe, a bien la forme appariee qu'il
annonce. Et surtout que `analyse`/`verdict` ne peuvent pas etre lus apres coup :
chaque verdict possible est verifie sur un aller-retour synthetique ou la
reponse est plantee a l'avance (§7).
"""
import numpy as np
import pandas as pd
import pytest

from nova_lab import REPO_ROOT
from scripts.build_equalbudget_batch import (
    CALIB_BANDS,
    CALIB_OFFSET_DEPRESSED,
    CALIB_OFFSET_OK,
    CALIB_RHO,
    INCUMBENT,
    JACCARD_CEILING,
    K,
    MIN_CDR3,
    N_SLOTS,
    RIVAL_GAP,
    SIGMA_2,
    WINS_FOR_FRESH,
    WINS_FOR_PERMUT,
    analyse,
    verdict,
)
from scripts.build_scramble_batch import jaccard_5mer

BATCH_PATH = REPO_ROOT / "inputs" / "kaggle_equalbudget.csv"


# ---------------------------------------------------------------------------
# les seuils sont pre-enregistres : ils ne doivent pas bouger en silence
# ---------------------------------------------------------------------------
def test_thresholds_are_the_ones_written_in_the_docstring():
    assert (RIVAL_GAP, SIGMA_2, INCUMBENT) == (0.014, 0.027, 0.3869)
    assert (N_SLOTS, K, MIN_CDR3) == (8, 3, 12)
    assert (WINS_FOR_PERMUT, WINS_FOR_FRESH, JACCARD_CEILING) == (5, 3, 0.90)
    assert (CALIB_RHO, CALIB_OFFSET_OK, CALIB_OFFSET_DEPRESSED) == (0.80, 0.05, -0.10)
    assert len(CALIB_BANDS) == 7 and CALIB_BANDS[0][0] == 0.15 and CALIB_BANDS[-1][1] == 0.90


# ---------------------------------------------------------------------------
# lecture 0 : l'instrument. Les trois etats doivent etre retrouves.
# ---------------------------------------------------------------------------
def _calib(tmp_path, offset: float, scramble: bool = False):
    """14 sequences d'archive, notees chez nous a `offset` de leur valeur.

    `scramble` casse la correlation : notre valeur devient independante de la
    leur, ce qui est le cas « l'instrument ne reproduit pas ».
    """
    rng = np.random.default_rng(1)
    rows, scores = [], []
    arch_values = np.linspace(0.16, 0.86, 14)
    for i, a in enumerate(arch_values):
        s = f"CALIB{i:02d}" + "A" * 100
        rows.append(dict(sequence=s, arm="calib", slot=-1, label=f"calib#{i}",
                         group="calibration", previous_iiptm=a, previous_run="archive",
                         archive_iiptm=a, archive_pae=12.0, band_low=0.15, band_high=0.90,
                         v_gene=None, cdr3_length=np.nan, cdr3_original="", cdr3_now="",
                         v_identity=85.0, j_identity=92.0, layer1="ok",
                         jaccard_max_batch=0.5))
        ours = (0.2 + 0.1 * rng.random()) if scramble else a + offset
        scores.append(dict(sequence=s, design_iiptm=ours, design_ptm=0.80))
    bp, sp = tmp_path / "cbatch.csv", tmp_path / "cscores.csv"
    pd.DataFrame(rows).to_csv(bp, index=False)
    pd.DataFrame(scores).to_csv(sp, index=False)
    return str(sp), str(bp)


def test_reading0_calls_a_faithful_instrument_faithful(tmp_path):
    sp, bp = _calib(tmp_path, offset=0.001)
    res = analyse(sp, bp, log=lambda *a, **kw: None)
    assert res["reading0_state"] == "fidele"
    assert res["reading0_rho"] > 0.99
    assert any("INSTRUMENT FIDELE" in line for line in verdict(res))


def test_reading0_catches_a_depressed_instrument(tmp_path):
    """Le cas qui nous pend au nez : monotone, mais 0,45 trop bas."""
    sp, bp = _calib(tmp_path, offset=-0.45)
    res = analyse(sp, bp, log=lambda *a, **kw: None)
    assert res["reading0_state"] == "monotone mais deprime"
    assert res["reading0_rho"] > 0.99
    assert any("DEPRIME" in line for line in verdict(res))


def test_reading0_catches_an_instrument_that_does_not_reproduce(tmp_path):
    sp, bp = _calib(tmp_path, offset=0.0, scramble=True)
    res = analyse(sp, bp, log=lambda *a, **kw: None)
    assert res["reading0_state"] == "ne reproduit pas"
    assert any("NE REPRODUIT PAS" in line for line in verdict(res))


def test_reading0_does_not_hide_an_intermediate_offset(tmp_path):
    """Entre -0,10 et -0,05 : ni fidele ni deprime, et il faut le dire."""
    sp, bp = _calib(tmp_path, offset=-0.08)
    res = analyse(sp, bp, log=lambda *a, **kw: None)
    assert res["reading0_state"] == "monotone, decalage intermediaire"
    assert any("decalage intermediaire" in line for line in verdict(res))


# ---------------------------------------------------------------------------
# aller-retour synthetique : on plante la reponse, analyse doit la retrouver
# ---------------------------------------------------------------------------
def _synthetic(tmp_path, permut_offset: float, n_slots: int = N_SLOTS,
               donor_drift: float = 0.0, break_layer1: bool = False,
               jaccard: float = 0.5, wins: int | None = None):
    """Un lot et ses scores, ou la permutation bat les neufs de `permut_offset`.

    `wins` limite le nombre de creneaux ou l'avantage est applique ; les autres
    recoivent l'avantage oppose, pour planter une majorite sans planter une
    unanimite.
    """
    rows, scores = [], []
    rng = np.random.default_rng(0)
    for slot in range(n_slots):
        base = 0.20 + 0.01 * slot
        advantage = permut_offset if (wins is None or slot < wins) else -permut_offset
        seq = f"DONOR{slot:02d}" + "A" * 100
        rows.append(dict(sequence=seq, arm="donneur", slot=slot, label=f"donneur#{slot:02d}",
                         group="candidate", previous_iiptm=base, previous_run="r.csv",
                         v_gene="IGHV3-74", cdr3_length=14, cdr3_original="AR" + "Y" * 12,
                         cdr3_now="AR" + "Y" * 12, v_identity=85.0, j_identity=92.0,
                         layer1="ok", jaccard_max_batch=jaccard))
        scores.append(dict(sequence=seq, design_iiptm=base + donor_drift, design_ptm=0.78))
        for arm, bump in (("permut", advantage), ("neuf", 0.0)):
            for i in range(K):
                s = f"{arm.upper()}{slot:02d}{i}" + "A" * 100
                bad = break_layer1 and arm == "permut" and slot == 0 and i == 0
                rows.append(dict(sequence=s, arm=arm, slot=slot, label=f"{arm}#{slot:02d}.{i}",
                                 group="candidate", previous_iiptm=np.nan, previous_run="",
                                 v_gene="IGHV3-74", cdr3_length=14,
                                 cdr3_original="AR" + "Y" * 12, cdr3_now="AR" + "X" * 12,
                                 v_identity=79.0 if bad else 85.0, j_identity=92.0,
                                 layer1="human V identity 79.0 < 82" if bad else "ok",
                                 jaccard_max_batch=jaccard))
                # bruit minuscule pour que les maxima ne soient pas ex aequo
                scores.append(dict(sequence=s, design_iiptm=base + bump + 1e-4 * rng.random(),
                                   design_ptm=0.78))
    bp, sp = tmp_path / "batch.csv", tmp_path / "scores.csv"
    pd.DataFrame(rows).to_csv(bp, index=False)
    pd.DataFrame(scores).to_csv(sp, index=False)
    return str(sp), str(bp)


def test_analyse_finds_a_planted_permutation_win(tmp_path):
    sp, bp = _synthetic(tmp_path, permut_offset=0.05, wins=7)
    res = analyse(sp, bp, log=lambda *a, **kw: None)
    assert res["reading1_wins"] == 7
    assert res["reading1_mean"] > RIVAL_GAP
    assert any("LA PERMUTATION GAGNE" in line for line in verdict(res))


def test_analyse_finds_a_planted_fresh_win(tmp_path):
    sp, bp = _synthetic(tmp_path, permut_offset=-0.05, wins=7)
    res = analyse(sp, bp, log=lambda *a, **kw: None)
    assert res["reading1_wins"] <= WINS_FOR_FRESH
    assert res["reading1_mean"] < -RIVAL_GAP
    assert any("LES CANDIDATS NEUFS GAGNENT" in line for line in verdict(res))


def test_analyse_calls_a_small_difference_no_difference(tmp_path):
    """Un avantage sous l'ecart entre rivaux ne doit pas etre appele."""
    sp, bp = _synthetic(tmp_path, permut_offset=0.004, wins=8)
    res = analyse(sp, bp, log=lambda *a, **kw: None)
    assert res["reading1_wins"] == 8          # unanime...
    assert abs(res["reading1_mean"]) < RIVAL_GAP
    lines = verdict(res)
    assert any("AUCUNE DIFFERENCE" in line for line in lines)   # ...mais trop petit
    assert not any("GAGNE" in line for line in lines)


def test_the_control_catches_a_drifting_run(tmp_path):
    sp, bp = _synthetic(tmp_path, permut_offset=0.05, wins=7, donor_drift=0.10)
    res = analyse(sp, bp, log=lambda *a, **kw: None)
    assert res["reading4_pass"] is False
    assert any("TEMOIN EN ECHEC" in line for line in verdict(res))


def test_a_clean_run_passes_the_control(tmp_path):
    sp, bp = _synthetic(tmp_path, permut_offset=0.05, wins=7, donor_drift=0.005)
    res = analyse(sp, bp, log=lambda *a, **kw: None)
    assert res["reading4_pass"] is True


def test_a_permutation_that_breaks_layer_one_is_reported(tmp_path):
    sp, bp = _synthetic(tmp_path, permut_offset=0.05, wins=7, break_layer1=True)
    res = analyse(sp, bp, log=lambda *a, **kw: None)
    assert res["reading5"]["permut"]["layer1_ok"] == K * N_SLOTS - 1
    assert any("COUCHE 1 CASSEE" in line for line in verdict(res))


def test_a_jaccard_above_the_ceiling_is_reported(tmp_path):
    sp, bp = _synthetic(tmp_path, permut_offset=0.05, wins=7, jaccard=0.95)
    res = analyse(sp, bp, log=lambda *a, **kw: None)
    assert res["reading5_diversity_pass"] is False
    assert any("DIVERSITE EN ECHEC" in line for line in verdict(res))


def test_the_high_tail_reading_counts_against_the_incumbent(tmp_path):
    """Une seule sequence au-dessus du record soumettable doit etre comptee."""
    sp, bp = _synthetic(tmp_path, permut_offset=0.05, wins=8)
    scores = pd.read_csv(sp)
    scores.loc[scores["sequence"].str.startswith("PERMUT000"), "design_iiptm"] = INCUMBENT + 0.01
    scores.to_csv(sp, index=False)
    res = analyse(sp, bp, log=lambda *a, **kw: None)
    assert res["reading3_permut_above"] == 1
    assert res["reading3_neuf_above"] == 0
    assert any("nouveau record soumettable" in line for line in verdict(res))


# ---------------------------------------------------------------------------
# le lot lui-meme, s'il a ete construit
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def batch():
    if not BATCH_PATH.exists():
        pytest.skip(f"{BATCH_PATH} pas encore construit")
    return pd.read_csv(BATCH_PATH)


def test_the_budget_is_equal_between_the_two_arms(batch):
    counts = batch["arm"].value_counts()
    assert counts["permut"] == counts["neuf"], "budget inegal : la comparaison serait truquee"
    assert counts["donneur"] == batch["slot"].nunique()
    for slot, sub in batch.groupby("slot"):
        assert (sub["arm"] == "permut").sum() == (sub["arm"] == "neuf").sum() == K, \
            f"creneau {slot} desequilibre"


def test_no_sequence_is_queried_twice(batch):
    """BoltzGen nomme ses entrees par l'empreinte de la sequence (§6)."""
    assert batch["sequence"].is_unique


def test_each_slot_is_matched_on_v_gene_and_cdr3_length(batch):
    """Le bras neuf doit tomber dans le meme creneau, sinon V ou L explique tout."""
    for slot, sub in batch.groupby("slot"):
        assert sub["v_gene"].nunique() == 1, f"creneau {slot} : plusieurs genes V"
        assert sub["cdr3_length"].nunique() == 1, f"creneau {slot} : plusieurs longueurs"
        assert (sub["cdr3_now"].str.len() == sub["cdr3_length"]).all(), \
            f"creneau {slot} : une sequence n'a pas la longueur de CDR3 du creneau"


def test_every_sequence_is_submittable(batch):
    """Un bras dont rien ne peut etre soumis ne mesure rien d'operationnel."""
    bad = batch[batch["layer1"] != "ok"]
    assert bad.empty, f"{len(bad)} sequences ne passent pas la couche 1 : {list(bad['label'])}"


def test_the_donors_cdr3_are_long_enough_to_permute(batch):
    donors = batch[batch["arm"] == "donneur"]
    assert (donors["cdr3_length"] >= MIN_CDR3).all()


def test_diversity_stays_under_the_validator_ceiling(batch):
    assert batch["jaccard_max_batch"].max() < JACCARD_CEILING


def test_permutations_keep_the_donor_cdr3_composition(batch):
    """Une permutation est un reordonnancement : meme multiset de residus."""
    for slot, sub in batch.groupby("slot"):
        donor = sub[sub["arm"] == "donneur"].iloc[0]
        for p in sub[sub["arm"] == "permut"].itertuples():
            assert sorted(p.cdr3_now) == sorted(donor["cdr3_original"]), \
                f"{p.label} n'est pas une permutation du CDR3 du donneur"


def test_fresh_candidates_are_not_near_twins_of_the_donor(batch):
    """Sinon le bras neuf serait une permutation deguisee."""
    for slot, sub in batch.groupby("slot"):
        donor = sub[sub["arm"] == "donneur"].iloc[0]["sequence"]
        for n in sub[sub["arm"] == "neuf"].itertuples():
            j = jaccard_5mer(n.sequence, donor)
            assert j < JACCARD_CEILING, f"{n.label} est a Jaccard {j:.3f} de son donneur"


def test_the_arms_are_interleaved(batch):
    """Un run interrompu a mi-chemin doit garder les deux bras equilibres."""
    half = batch.iloc[: len(batch) // 2]["arm"].value_counts()
    assert abs(half.get("permut", 0) - half.get("neuf", 0)) <= 2


def test_the_donors_were_scored_on_the_live_target(batch):
    from nova_lab.provenance import runs_on_current_target
    live = set(runs_on_current_target())
    donors = batch[batch["arm"] == "donneur"]
    assert donors["previous_iiptm"].notna().all(), "un donneur sans valeur precedente"
    assert set(donors["previous_run"]) <= live, \
        f"donneurs venant d'un run hors cible vivante : {set(donors['previous_run']) - live}"
