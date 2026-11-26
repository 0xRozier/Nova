"""The ptm batch's design guarantees, and its reading written before the data.

Two ways this batch could lie. The bands could differ on CDR3 length, and a
rise in iiptm would then be the loop's length rather than the fold confidence
-- so the match is tested, not assumed. And the verdict could be chosen after
seeing the number, so `verdict` is tested against its own thresholds.
"""
import numpy as np
import pandas as pd
import pytest

from scripts.build_ptm_batch import (BANDS, NOISE, PER_BAND, analyse, build,
                                     training_set, verdict)

SMALL = dict(pool_size=600, per_band=5)


@pytest.fixture(scope="module")
def batch():
    return build(**SMALL)


def test_no_sequence_is_queried_twice(batch):
    """BoltzGen keys by sequence hash: a repeat is paid for and not measured."""
    assert batch.sequence.duplicated().sum() == 0


def test_none_of_them_has_been_queried_on_the_live_target(batch):
    from nova_lab.config import NANOBODY_STRUCTURE
    from scripts.train_score_surrogate import load
    assert not set(batch.sequence) & set(load(NANOBODY_STRUCTURE).sequence)


def test_the_bands_are_ordered_and_span_the_threshold(batch):
    medians = [batch[batch.band == label].predicted_ptm.median()
               for _, _, label in BANDS if (batch.band == label).any()]
    assert medians == sorted(medians), "les bandes doivent monter"
    assert min(medians) < 0.749 < max(medians), (
        "le lot doit traverser le seuil du §87, sinon il ne teste rien")


def test_cdr3_length_is_matched_across_bands():
    """The confound: unmatched, a rise in iiptm could be the loop's length.

    This one needs a real pool, and that is part of the guarantee rather than a
    concession. The matching fills each band with the CDR3 lengths in the
    proportions the whole pool shows; with the 5-per-band configuration the
    other tests use, the low band has too few sequences of each length to fill
    the quota and falls back to the band at large, which lands its median at 18
    against 13 elsewhere. Measured: the gap between band medians is 6 at a pool
    of 3 000, 2 at 6 000 with 10 per band, and 0 at the shipped 6 000 with 15.
    So this asserts the SHIPPED configuration -- the only one whose matching is
    a claim about the batch that actually ran -- and pays the minute it costs.
    """
    matched = build()
    medians = matched.groupby("band").cdr3_length.median()
    assert medians.max() - medians.min() <= 2, medians.to_dict()


def test_the_training_set_pools_both_targets():
    """FINDINGS 87: this metric transfers, which is why pooling is legitimate."""
    d = training_set()
    assert len(d) > 1000
    assert d.y.between(0.6, 0.9).all()


def test_the_batch_is_deterministic():
    assert list(build(**SMALL).sequence) == list(build(**SMALL).sequence)


def test_layer_one_is_recorded_not_enforced(batch):
    assert "layer1" in batch.columns
    assert (batch.layer1 != "ok").any(), (
        "si tout passe, le lot a ete filtre et ne couvre plus les bandes voulues")


def test_the_verdict_follows_its_own_thresholds():
    closed = {"rho_predictor": 0.1, "rho_dose": 0.9, "ecart": 1.0, "seuil_ecart": 2 * NOISE}
    assert "ne generalise pas" in verdict(closed), "le critere 1 doit dominer"
    good = {"rho_predictor": 0.6, "rho_dose": 0.5, "ecart": 0.05, "seuil_ecart": 2 * NOISE}
    assert "utilisable" in verdict(good)
    half = {"rho_predictor": 0.6, "rho_dose": 0.5, "ecart": 0.001, "seuil_ecart": 2 * NOISE}
    assert "un seul des deux" in verdict(half)
    dead = {"rho_predictor": 0.6, "rho_dose": 0.0, "ecart": 0.001, "seuil_ecart": 2 * NOISE}
    assert "inutile comme procuration" in verdict(dead)
    assert "pas assez" in verdict({})


def test_analyse_recovers_a_planted_threshold(tmp_path, batch):
    """A synthetic return trip, so the reader runs before the GPU does."""
    bp = tmp_path / "batch.csv"; batch.to_csv(bp, index=False)
    rng = np.random.default_rng(0)
    ptm = batch.predicted_ptm.to_numpy() + rng.normal(0, 0.004, len(batch))
    scores = pd.DataFrame({"sequence": batch.sequence, "design_ptm": ptm,
                           "design_iiptm": 0.15 + 2.0 * np.maximum(ptm - 0.749, 0)
                           + rng.normal(0, 0.005, len(batch))})
    sp = tmp_path / "scores.csv"; scores.to_csv(sp, index=False)
    r = analyse(sp, bp)
    assert r["rho_predictor"] > 0.8
    assert r["rho_dose"] > 0.5
    assert r["ecart"] > r["seuil_ecart"]
    assert "utilisable" in verdict(r)


def test_analyse_reports_the_clean_negative(tmp_path, batch):
    """Predictor works, iiptm does not follow: the FINDINGS 87 lever closes."""
    bp = tmp_path / "batch.csv"; batch.to_csv(bp, index=False)
    rng = np.random.default_rng(1)
    ptm = batch.predicted_ptm.to_numpy() + rng.normal(0, 0.004, len(batch))
    scores = pd.DataFrame({"sequence": batch.sequence, "design_ptm": ptm,
                           "design_iiptm": rng.normal(0.19, 0.01, len(batch))})
    sp = tmp_path / "scores.csv"; scores.to_csv(sp, index=False)
    r = analyse(sp, bp)
    assert r["rho_predictor"] > 0.8
    assert abs(r["rho_dose"]) < 0.3
    assert "inutile comme procuration" in verdict(r)
