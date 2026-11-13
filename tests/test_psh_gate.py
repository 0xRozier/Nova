"""The PSH gate, and the calibration it now rests on.

PSH > 89 is the one validator rule that refused every production submission,
so these tests guard two different things:

  - the arithmetic of `nova_lab.psh_calibration`, which turns our fold-averaged
    reading into the probability that the validator's own draw exceeds 89;
  - the evidence underneath it, `results/psh_calibration_epoch26174.csv`. That
    file is the 36 submissions the validators accepted plus the 1 they refused,
    folded three times each on the install they run. If a change to the
    toolchain ever puts an accepted submission above 89 again, the scale has
    drifted off theirs and the gate is measuring something else -- which is
    exactly the failure that cost twelve epochs, and it was invisible for a
    month because nothing checked it.

    python -m pytest tests/test_psh_gate.py -q
"""
import ast
import math
from pathlib import Path

import pytest

from nova_lab.config import PSH_STRICT_REJECTION_THRESHOLD
from nova_lab.psh_calibration import (SIGMA, THRESHOLD, expected_value, rank,
                                      refusal_probability, verdict)

CALIBRATION = Path(__file__).resolve().parent.parent / "results" / "psh_calibration_epoch26174.csv"


def test_threshold_is_the_validators_own_number():
    """89, from upstream config.yaml, not a margin we chose."""
    assert THRESHOLD == float(PSH_STRICT_REJECTION_THRESHOLD) == 89.0


def test_refusal_probability_is_a_half_at_the_threshold():
    assert refusal_probability(89.0, folds=3) == pytest.approx(0.5)


def test_refusal_probability_rises_with_the_reading():
    probabilities = [refusal_probability(p, folds=3) for p in (80, 85, 88, 89, 92, 99)]
    assert probabilities == sorted(probabilities)
    assert probabilities[0] < 0.001 and probabilities[-1] > 0.999


def test_more_folds_sharpen_the_verdict():
    """Averaging folds cuts our own noise; the validator's single draw remains.

    So the probability moves towards 0 or 1 with k, but never past the
    validator's own spread -- a point 2 below the threshold cannot be made
    certain by folding.
    """
    assert refusal_probability(87.0, folds=1) > refusal_probability(87.0, folds=3)
    assert refusal_probability(91.0, folds=1) < refusal_probability(91.0, folds=3)
    assert refusal_probability(87.0, folds=99) > 0.1  # 2 / 1.74 is only 1.15 sigma


def test_an_unmeasured_candidate_has_no_probability():
    assert math.isnan(refusal_probability(float("nan")))
    assert math.isnan(refusal_probability(float("nan"), folds=3))


def test_expected_value_discounts_quality_by_the_risk():
    safe = expected_value(0.776, 84.4, folds=3)
    risky = expected_value(0.784, 88.6, folds=3)
    # The decision this changed: the better binder is worth less once the 42%
    # chance of being thrown out before scoring is priced in.
    assert safe > risky
    assert expected_value(0.5, float("nan")) == 0.5  # unmeasured keeps its quality


def test_rank_puts_a_measured_candidate_ahead_of_an_unmeasured_one():
    """An unscreened sequence is what the production run submitted eleven times."""
    quality = {"SAFE": 0.70, "RISKY": 0.80, "UNKNOWN": 0.99}
    psh = {"SAFE": 84.0, "RISKY": 91.0}
    assert rank(quality, psh, folds=3) == ["SAFE", "RISKY", "UNKNOWN"]


def test_verdict_says_when_nothing_was_measured():
    assert "not measured" in verdict(float("nan"))
    assert "89" not in verdict(float("nan"))


@pytest.fixture(scope="module")
def calibration():
    if not CALIBRATION.exists():
        pytest.skip(f"{CALIBRATION} is missing; it is the evidence under the gate")
    import pandas as pd
    return pd.read_csv(CALIBRATION)


def test_every_accepted_submission_reads_below_the_threshold(calibration):
    """The scale is the validator's, or the gate is meaningless.

    On the install whose DSSP step was failing, 44-62% of these same accepted
    submissions read above 89 and the project concluded the threshold could not
    be used. With mkdssp 3.1.4 and the filter's patch in place, 36 of 36 read
    below it.
    """
    accepted = calibration[calibration["truth"] == "accepted"]["psh_mean"].dropna()
    assert len(accepted) >= 36, "the calibration set lost sequences — check for failed folds"
    over = accepted[accepted > THRESHOLD]
    assert over.empty, (
        f"{len(over)} of {len(accepted)} submissions the validators ACCEPTED read "
        f"above {THRESHOLD} here (max {accepted.max():.1f}). The toolchain has "
        f"drifted off theirs: run python -m scripts.doctor.")


def test_the_refused_submission_reads_near_the_threshold(calibration):
    """Their log says 89.57 for it. We read ~88, i.e. about 1.5 low."""
    refused = calibration[calibration["truth"] == "refused"]["psh_mean"].dropna()
    assert len(refused) == 1
    assert 85.0 < refused.iloc[0] < 92.0, (
        f"we read {refused.iloc[0]:.1f} where the validator logged 89.57")


def test_the_fold_noise_matches_the_sigma_the_gate_uses(calibration):
    """SIGMA is pooled from this file; a drift means the gate is miscalibrated."""
    import numpy as np
    folds = [ast.literal_eval(s) for s in calibration["psh_folds"]]
    triples = [[float(x) for x in f] for f in folds if len(f) >= 3]
    assert len(triples) >= 30, "too few three-fold sequences to pool a deviation"
    pooled = math.sqrt(np.mean([np.var(t, ddof=1) for t in triples]))
    assert pooled == pytest.approx(SIGMA, abs=0.5), (
        f"one measurement now deviates by {pooled:.2f}, not {SIGMA}; "
        f"refusal_probability is calibrated on the old figure")
