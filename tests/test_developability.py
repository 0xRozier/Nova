"""The validator's developability bounds, and the two false refusals that found them.

Of the 36 submissions the validators accepted in epoch 26174, our gate refused
two. Both were our fault, in different ways, and the fix for each is pinned
here with its real numbers:

  - one read PPC 1.2858 on the first fold and was refused, while its three
    folds mean 1.1687 and pass. We were averaging PSH and Aggrescan3D over
    folds and taking every other flag from fold 1 -- on a metric that deviates
    by 0.037, and by 0.23 at its worst.
  - the other read PSH 70.3 on all three folds, under the validator's lower
    bound of 73.4, so it was not noise: it is our scale sitting 2 to 3 points
    below theirs. Shifted onto their scale it lands on the bound rather than
    under it.

The second is also why the bounds are checked on a shifted reading. A scale
that reads low is wrong at both ends and in opposite directions -- too
permissive at the top, too strict at the bottom -- and only the top end costs
us an epoch, which is why the shift is the conservative choice.

    python -m pytest tests/test_developability.py -q
"""
import pytest

from nova_lab.developability import (BOUNDS, PSH_BIAS, SIGMA, band, borderline,
                                     breaches, on_their_scale)

# The two accepted submissions our gate used to refuse (epoch 26174).
PPC_FOLDS = [1.2858, 1.1078, 1.1125]      # mean 1.1687
PSH_LOW_FOLDS = [69.4961, 71.5934, 69.7235]  # mean 70.2710


def mean(values):
    return sum(values) / len(values)


def test_the_bounds_are_the_filters_not_tnps():
    """metanano/config.py, rounded to two decimals, which is what decides.

    TNP's own red thresholds differ in the fourth decimal and in both
    directions (PPC 1.1853 against 1.18, PSH max 155.4699 against 155.47), and
    the validator uses its own, so these are the numbers that matter.
    """
    assert BOUNDS["PSH"] == (73.4, 155.47)
    assert BOUNDS["PPC"] == (None, 1.18)
    assert BOUNDS["PNC"] == (None, 1.88)
    assert BOUNDS["C"] == (0.56, 1.61)
    assert BOUNDS["L"] == (20, 39) and BOUNDS["L3"] == (5, 23)


def test_the_ppc_false_refusal_is_gone_once_folds_are_averaged():
    assert breaches({"PPC": PPC_FOLDS[0]}, folds=1), "fold 1 alone did breach"
    assert not breaches({"PPC": mean(PPC_FOLDS)}, folds=3)


def test_the_low_psh_false_refusal_is_gone_once_the_scale_is_shifted():
    """Three folds all read under the bound, so folding cannot fix this one."""
    assert all(v < BOUNDS["PSH"][0] for v in PSH_LOW_FOLDS)
    assert not breaches({"PSH": mean(PSH_LOW_FOLDS)}, folds=3)
    assert borderline({"PSH": mean(PSH_LOW_FOLDS)}, folds=3), (
        "it should still be reported as sitting on the bound")


def test_the_shift_applies_to_psh_and_nothing_else():
    assert on_their_scale("PSH", 80.0) == 80.0 + PSH_BIAS
    for metric in ("C", "PPC", "PNC", "L", "L3"):
        assert on_their_scale(metric, 1.0) == 1.0


def test_it_still_refuses_what_the_validator_refuses():
    assert breaches({"PSH": 60.0}, folds=3)      # far under the floor
    assert breaches({"PSH": 170.0}, folds=3)     # far over the ceiling
    assert breaches({"PPC": 1.5}, folds=3)
    assert breaches({"PNC": 2.5}, folds=3)
    assert breaches({"C": 1.8}, folds=3)
    assert breaches({"L3": 24}, folds=3) and breaches({"L3": 4}, folds=3)
    assert breaches({"L": 45}, folds=3)


def test_a_metric_we_did_not_measure_is_not_a_metric_that_passed():
    """It is simply not judged, and the caller has to notice its absence.

    Treating a missing measurement as a pass is how eleven submissions went out
    unscreened; treating it as a failure would throw away every candidate whose
    model TNP could not build.
    """
    assert breaches({}, folds=3) == []
    assert breaches({"PSH": float("nan")}, folds=3) == []
    assert breaches({"PSH": 60.0, "PPC": float("nan")}, folds=3) == [
        b for b in breaches({"PSH": 60.0}, folds=3)]


def test_the_band_shrinks_with_folds_and_is_zero_for_the_integers():
    assert band("PSH", 1) > band("PSH", 3) > band("PSH", 9)
    assert band("PSH", 4) == pytest.approx(band("PSH", 1) / 2)
    # L and L3 come off the IMGT numbering, not the structure, so they cannot
    # move between folds and get no tolerance at all.
    assert band("L", 3) == 0.0 and band("L3", 3) == 0.0
    assert SIGMA["L"] == 0.0 and SIGMA["L3"] == 0.0


def test_a_value_exactly_on_a_bound_is_not_a_breach():
    """The validator refuses `psh > max`, strictly, and `psh < min`, strictly."""
    low, high = BOUNDS["PSH"]
    assert not breaches({"PSH": high - PSH_BIAS}, folds=3)
    assert not breaches({"PSH": low - PSH_BIAS}, folds=3)
    assert not breaches({"PPC": BOUNDS["PPC"][1]}, folds=3)
