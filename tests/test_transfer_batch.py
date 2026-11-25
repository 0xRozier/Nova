"""The transfer batch, and the reading written before its data exists.

This batch answers one question -- do scores against the previous target
predict scores against the live one -- and the answer decides whether 1 490
paid GPU queries are an asset or history. Two things can quietly ruin it:

  range restriction, which would attenuate a transfer that is real and talk us
  out of reusing the queries;

  a verdict chosen after seeing the number, which is why `verdict()` is tested
  against its own thresholds rather than described in prose.
"""
import numpy as np
import pandas as pd
import pytest

from nova_lab.config import PREVIOUS_NANOBODY
from scripts.build_transfer_batch import (BANDS, PER_BAND, analyse, build,
                                          repeat_rows, transfer_rows, verdict)
from scripts.train_score_surrogate import load


@pytest.fixture(scope="module")
def batch():
    return build()


def test_no_sequence_is_queried_twice(batch):
    """BoltzGen keys by sequence hash: a repeat is paid for and not measured."""
    assert batch.sequence.duplicated().sum() == 0


def test_both_arms_are_present(batch):
    counts = batch.arm.value_counts()
    assert counts["transfer"] >= 2 * BANDS
    assert counts["repeat"] >= 3, "sans repetitions, un rho faible est ambigu"


def test_every_transfer_sequence_has_a_previous_score_over_several_runs(batch):
    arm = batch[batch.arm == "transfer"]
    previous = arm.previous_iiptm.astype(float)
    assert previous.notna().all()
    assert (arm.previous_runs.astype(float) > 1).all(), (
        "un point dont l'abscisse vient d'une seule requete porte le bruit "
        "complet sur l'axe qu'on croit precis")


def test_the_sample_is_not_range_restricted(batch):
    """The whole design rests on this: a narrow sample hides a real transfer."""
    pool = load(PREVIOUS_NANOBODY["structure"])
    pool = pool[pool.runs > 1]
    arm = batch[batch.arm == "transfer"].previous_iiptm.astype(float)
    span_pool = pool.y.max() - pool.y.min()
    assert (arm.max() - arm.min()) / span_pool > 0.85
    # and spread, not clumped at the ends
    edges = np.linspace(pool.y.min(), pool.y.max(), BANDS + 1)
    occupied = np.unique(np.clip(np.digitize(arm, edges[1:-1]), 0, BANDS - 1))
    assert len(occupied) == BANDS, f"bandes vides : {set(range(BANDS)) - set(occupied)}"


def test_each_band_is_sampled_evenly():
    arm = transfer_rows()
    sizes = arm.groupby("band").size()
    assert sizes.min() >= min(PER_BAND, sizes.max())


def test_the_repeats_span_the_calibration_range():
    again = repeat_rows()
    values = again.calibration_iiptm.astype(float)
    assert values.max() - values.min() > 0.05, (
        "des repetitions prises dans un coin ne mesurent pas le bruit de la plage")


def test_the_batch_is_deterministic():
    assert list(build().sequence) == list(build().sequence)


def test_the_verdict_follows_the_pre_registered_thresholds():
    assert "des maintenant" in verdict({"rho": 0.75})
    assert "partiel" in verdict({"rho": 0.45})
    assert "propre a la cible" in verdict({"rho": 0.10})
    assert "pas assez" in verdict({})


def test_analyse_recovers_a_planted_correlation(tmp_path, batch):
    """A synthetic return trip, so the reader is exercised before the GPU runs."""
    batch_path = tmp_path / "batch.csv"
    batch.to_csv(batch_path, index=False)
    rng = np.random.default_rng(0)
    previous = pd.to_numeric(batch.previous_iiptm, errors="coerce").fillna(0.3).to_numpy()
    scores = pd.DataFrame({
        "sequence": batch.sequence,
        # a real transfer plus noise, and a fold confidence unrelated to it
        "design_iiptm": 0.1 + 0.4 * previous + rng.normal(0, 0.01, len(batch)),
        "design_ptm": rng.normal(0.8, 0.02, len(batch)),
    })
    scores_path = tmp_path / "scores.csv"
    scores.to_csv(scores_path, index=False)

    result = analyse(scores_path, batch_path)
    assert result["n_transfer"] == int((batch.arm == "transfer").sum())
    assert result["rho"] > 0.8
    assert result["rho_sans_ptm"] > 0.8, "un ptm non lie ne doit rien retirer"
    assert "des maintenant" in verdict(result)


def test_analyse_reports_no_transfer_when_there_is_none(tmp_path, batch):
    batch_path = tmp_path / "batch.csv"
    batch.to_csv(batch_path, index=False)
    rng = np.random.default_rng(1)
    scores = pd.DataFrame({
        "sequence": batch.sequence,
        "design_iiptm": rng.normal(0.2, 0.03, len(batch)),
        "design_ptm": rng.normal(0.8, 0.02, len(batch)),
    })
    scores_path = tmp_path / "scores.csv"
    scores.to_csv(scores_path, index=False)
    result = analyse(scores_path, batch_path)
    assert abs(result["rho"]) < 0.35
    assert "propre a la cible" in verdict(result) or "partiel" in verdict(result)
