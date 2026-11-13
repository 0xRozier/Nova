"""The validator's developability bounds, applied to metrics that are noisy.

## Whose thresholds

TNP prints its own verdict per metric -- "(RED flag)" -- and this project used
those. The validator does not: `NOVA-nanobody-filter`'s DevelopabilityService
takes TNP's *values* and compares them against its own configuration
(`metanano/filters/developability.py`, `metanano/config.py`), and
`neurons/validator/nanobody_validity.py` then drops any sequence carrying a
flag. Those numbers are TNP's rounded to two decimals, and the rounding goes
both ways:

    metric   TNP's own red       the filter's config     who is stricter
    C        0.5634 - 1.6132     0.56 - 1.61             filter, at the top
    PSH      73.4039 - 155.4699  73.4 - 155.47           TNP, at the top
    PPC      > 1.1853            > 1.18                  filter
    PNC      > 1.8838            > 1.88                  filter

Hundredths, but the gate is a cliff, so the values below are the filter's.

## Why a band and not a comparison

Every one of these is computed on a structure NanoBodyBuilder2 predicted, and
that prediction is not deterministic. Measured over 37 sequences folded three
times each (`results/psh_calibration_epoch26174.csv`), the deviation of ONE
measurement is:

    L, L3   0           (integers off the numbering, not the structure)
    C       0.005
    PSH     2.7         (five estimates from 1.74 to 3.29; see
                         nova_lab.psh_calibration.SIGMA)
    PPC     0.037       but 0.23 at its worst
    PNC     0.023       but 0.14 at its worst

So a single reading next to a bound says nothing. That was not academic: of the
36 submissions the validators accepted in epoch 26174, our gate refused one on
PPC whose three folds read 1.2858, 1.1078 and 1.1125 -- the first fold alone
breached 1.18, the mean is 1.1687 and passes. We were averaging PSH and
Aggrescan3D over folds while taking the other flags from fold 1.

A breach therefore counts only beyond two deviations of the fold mean, which is
the rule Aggrescan3D already used here.

## The bias we appear to carry

Three independent readings say our PSH sits **2 to 3 points below theirs**:

  - the one sequence whose log we have: we read 88.0 over three folds where
    their log says 89.57 (-1.6);
  - an accepted submission we read at 70.3 over three folds, which their own
    lower bound of 73.4 means they read at least 73.4 (<= -3.1);
  - the threshold that best separates their 296 verdicts on our scale
    (epochs 26299-26306, Youden) is 86.5 against their 89 (-2.5).

Three weak estimates agreeing in sign and magnitude, on an axis where being
wrong is expensive in one direction only: under-reading makes us think a
candidate is safe when the validator refuses it. PSH_BIAS is therefore applied
when estimating the risk of refusal, and it explains the second false refusal
too -- 70.3 plus the bias sits on the lower bound rather than below it.

It is an estimate, not a calibration. Narrowing it needs sequences whose
validator-side reading we can see, which means the W&B logs.

A scale that reads low is wrong at **both** of PSH's bounds, and in opposite
directions: at the top it makes us too permissive, at the bottom too strict. So
the bias is added to the PSH reading wherever it meets one of the validator's
numbers -- the upper bound here, the strict 89 in nova_lab.psh_calibration --
and the two never apply to the same comparison. The other metrics get no such
correction: we have no evidence of a bias on them, and inventing one would be
worse than leaving them raw.
"""
from __future__ import annotations

# metanano/config.py: (minimum, maximum); None means unbounded on that side.
BOUNDS: dict[str, tuple[float | None, float | None]] = {
    "L": (20, 39),
    "L3": (5, 23),
    "C": (0.56, 1.61),
    "PSH": (73.4, 155.47),
    "PPC": (None, 1.18),
    "PNC": (None, 1.88),
}

# Deviation of ONE measurement, pooled over 37 sequences folded three times.
SIGMA: dict[str, float] = {
    "L": 0.0, "L3": 0.0, "C": 0.005, "PSH": 2.7, "PPC": 0.037, "PNC": 0.023,
}

# How many deviations of the fold mean a breach must clear to count.
BANDS = 2.0

# Our PSH appears to read this far below the validator's (see the module note).
PSH_BIAS = 2.0


def band(metric: str, folds: int = 1) -> float:
    """How far past a bound a reading must sit before it means anything."""
    return BANDS * SIGMA.get(metric, 0.0) / max(folds, 1) ** 0.5


def on_their_scale(metric: str, value: float) -> float:
    """Our reading, shifted to what the validator would read for it."""
    return value + PSH_BIAS if metric == "PSH" else value


def breaches(values: dict[str, float], folds: int = 1) -> list[str]:
    """Which of the validator's bounds these fold-averaged metrics breach.

    `values` maps TNP's metric names to the mean over `folds` folds. A metric
    that is absent is not judged: a measurement we do not have is not a
    measurement that passed, and the caller has to say so itself.
    """
    out = []
    for metric, (low, high) in BOUNDS.items():
        raw = values.get(metric)
        if raw is None or raw != raw:
            continue
        value = on_their_scale(metric, raw)
        margin = band(metric, folds)
        if low is not None and value < low - margin:
            out.append(f"{metric} {value:.4g} below {low} (band {margin:.3g})")
        if high is not None and value > high + margin:
            out.append(f"{metric} {value:.4g} above {high} (band {margin:.3g})")
    return out


def borderline(values: dict[str, float], folds: int = 1) -> list[str]:
    """Bounds a reading crosses but not past the noise -- worth logging, not refusing."""
    out = []
    for metric, (low, high) in BOUNDS.items():
        raw = values.get(metric)
        if raw is None or raw != raw:
            continue
        value = on_their_scale(metric, raw)
        margin = band(metric, folds)
        if low is not None and low - margin <= value < low:
            out.append(f"{metric} {value:.4g} just under {low}")
        if high is not None and high < value <= high + margin:
            out.append(f"{metric} {value:.4g} just over {high}")
    return out
