"""What our PSH reading actually says about the validator's verdict.

`config.PSH_STRICT_REJECTION_THRESHOLD` is 89 and that is the validator's own
number, applied to *its* measurement. Ours is a different measurement of the
same quantity, and treating 89 as our threshold was wrong in a way that cost us
a third of the calendar.

The evidence, from 300 submissions the validators ruled on, profiled by our own
TNP (`results/psh_by_length.csv`, `results/psh_accepted_vs_refused.csv`):

  - Sequences the validators ACCEPTED read above 89 in our pipeline 28 % of the
    time at CDR3 5-8, 33 % at 9-12, 60 % at 13-16 and **76 % at 17+**, median
    92.0 there. Reading 92 at a long CDR3 is what an accepted sequence reads.
  - Splitting accepted against refused at the same lengths: at 17+ a reading of
    96 or less is accepted 68 % of the time [57-77] and above 96, 28 %
    [18-42] (Fisher p = 2e-05). At 5-8 the two groups are indistinguishable --
    accepted median 84.1, refused 83.7, AUC 0.46 -- so at short lengths our
    number carries no information about the verdict at all.

Two consequences for the miner:

  1. Refusing an epoch because our reading exceeds 89 discards candidates that
     would have been accepted two times in three. The old rule was stricter
     still (89 minus a 4.6 margin, so 84.4), which at CDR3 17+ refused every
     candidate in all nine epochs replayed offline -- the same nine we ranked
     first in six of.
  2. A refused submission costs the epoch, which is lost anyway when we submit
     nothing. So the bar for submitting is not "certain to pass", it is "better
     than not submitting": any positive acceptance probability clears it. The
     guard that remains is there for the pathological case, not the ordinary
     one.

The curve is empirical on purpose. A four-parameter fit of the measurement
model (our bias, the validator's noise) was tried and withdrawn: a quarter of
the refusals at 17+ read BELOW the accepted median, so they were refused for
something our PSH does not see, the model had no term for it, and its
parameters ran to whatever bound the search allowed. See FINDINGS 3ter.
"""
from __future__ import annotations

import bisect

# P(accepted | our PSH in this band), reweighted to the field's real base rate.
# Each entry is (upper bound of the band, probability, n). None as a bound means
# "everything above". Bands are as measured; do not interpolate inside one.
LONG_CURVE = [(88, 0.66, 25), (92, 0.74, 26), (96, 0.62, 23),
              (100, 0.25, 25), (None, 0.31, 27)]

# At CDR3 5-8 the accepted and refused distributions overlap (AUC 0.46), so
# there is no curve to give: the reading does not predict the verdict. What was
# measured is kept here as a record -- 77 % of 38 at or below 88, then 64 % of
# 8, 69 % of 7, and one single sequence above 96 that happened to be accepted.
# Those last three bands are too thin to read anything into, and the last would
# otherwise hand the miner a 100 % from n = 1, so the short band answers with
# its base rate instead (see acceptance_probability).
SHORT_CURVE = [(88, 0.77, 38), (92, 0.64, 8), (96, 0.69, 7), (None, 1.00, 1)]

# Share of submissions the validators accepted, per band, over the whole paired
# field: 651 of 869 at CDR3 5-8, 741 of 1433 at 17+.
BASE_RATE_SHORT = 0.75
BASE_RATE_LONG = 0.52

# Below this many sequences a band's own rate is not reported as a probability.
MIN_BAND_N = 10

LONG_CDR3 = 17  # where the long-band measurement applies

# Above this reading, acceptance at a long CDR3 drops from ~68% to ~28%
# (Fisher p = 2e-05). Used for reporting and for ordering candidates, not as a
# refusal: 28% beats the 0% of submitting nothing.
LONG_BREAK = 96

# Refuse the epoch only below this estimated acceptance probability. The
# measured floor is ~25%, so this does not fire on ordinary epochs by design;
# it exists for a candidate far outside anything we have measured.
MIN_ACCEPTANCE = 0.10


def curve_for(cdr3_length: int) -> list[tuple[int | None, float, int]]:
    return LONG_CURVE if cdr3_length >= LONG_CDR3 else SHORT_CURVE


def acceptance_probability(psh: float, cdr3_length: int) -> float:
    """P(the validators accept it | our PSH reads `psh`), as measured.

    A step function over the measured bands, because that is what the data
    supports: at CDR3 17+ the three bands below 96 (66 %, 74 %, 62 %) have
    overlapping intervals, so there is no gradient to read inside them, and
    smoothing one in would invent precision we do not have.
    """
    if not discriminates(cdr3_length):
        # The reading says nothing at this length, so the honest answer is the
        # band's base rate, not a per-bin figure drawn from a handful of
        # sequences -- one of those bins holds a single sequence.
        return BASE_RATE_SHORT
    _, probability, n = _band(psh, cdr3_length)
    return probability if n >= MIN_BAND_N else BASE_RATE_LONG


def _band(psh: float, cdr3_length: int) -> tuple[int | None, float, int]:
    curve = curve_for(cdr3_length)
    bounds = [b for b, _, _ in curve[:-1]]
    return curve[bisect.bisect_left(bounds, psh)]


def band_label(psh: float, cdr3_length: int) -> str:
    """The measured band a reading falls in, for logs and the journal."""
    curve = curve_for(cdr3_length)
    low = None
    for upper, probability, n in curve:
        if upper is None or psh <= upper:
            span = (f"> {low}" if upper is None
                    else f"<= {upper}" if low is None else f"{low}-{upper}")
            if not discriminates(cdr3_length):
                return (f"{span}: at this length our PSH does not predict the "
                        f"verdict (AUC 0.46), base rate {BASE_RATE_SHORT:.0%}")
            if n < MIN_BAND_N:
                return f"{span}: only {n} measured, falling back to the base rate"
            return f"{span} (accepted {probability:.0%} of {n} measured)"
        low = upper
    return "unmeasured"


def discriminates(cdr3_length: int) -> bool:
    """Whether our PSH predicts the verdict at this length at all.

    At CDR3 5-8 it does not (AUC 0.46 on 54 sequences), so a PSH-based refusal
    there rests on nothing. Selecting the lowest reading is still harmless.
    """
    return cdr3_length >= LONG_CDR3
