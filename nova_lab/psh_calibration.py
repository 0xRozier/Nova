"""What our PSH reading says about the validator's verdict, measured twice.

The validator's rule is one line of its own code
(`neurons/validator/nanobody_validity.py`, upstream db9bf19):

    psh_strict_rejection = [... if result["surface_hydrophobic_patches"]
                            > config["psh_strict_rejection_threshold"]]

with the threshold at 89. A sequence above it is dropped before anything is
scored, and the epoch is lost. There is no band, no curve and no partial
credit: the only question is whether the validator's own reading exceeds 89.

## Why this file was rewritten

The version it replaces fitted an empirical acceptance curve -- bands at 88,
92, 96, a "break" at 96, and the advice that a reading of 92 at a long CDR3 is
"what an accepted sequence reads". Every number in it came from a local TNP
whose DSSP step was failing. TNP prints one line when that happens and carries
on with a PSH computed another way, so the readings looked plausible and were
measurements of something else. On that install, 44-62 % of the submissions the
validators had ACCEPTED read above 89, which is what made the threshold look
unusable and the curve look necessary.

With the install the validator actually runs -- TNP at 29dcac7, the filter's
process_pdb patch, mkdssp 3.1.4, FreeSASA 2.0.1, ImmuneBuilder's refine.py typo
fixed; `nova_lab.toolchain` checks all of it -- the same 36 accepted
submissions of epoch 26174, folded three times each, read:

    minimum 70.3, median 83.5, maximum 88.9, and NONE above 89.

36 of 36 below the threshold the validator applies. Our reading is the
validator's reading. The one submission it refused in that epoch reads 88.0
here against the 89.57 in its log, so if anything we read about 1.5 low on that
single point of comparison -- a reason to keep margin, not to distrust the
scale.

## The noise, re-measured

NanoBodyBuilder2's refinement is not deterministic, so one measurement is a
draw. Over the 37 calibration sequences folded three times each, the pooled
deviation of ONE measurement is **1.74** (median per-sequence 0.95, worst
3.88). The 4.01 this project carried came from the broken path and was more
than twice too large.

Folding k times gives a mean with deviation 1.74/sqrt(k); the validator still
takes a single draw of its own. So from our k-fold mean m:

    P(the validator reads above 89)
        = Phi((m + bias - 89) / (2.7 * sqrt(1 + 1/k)))

where `bias` is how far below theirs our scale appears to read, estimated at
2 to 3 points by three independent routes (nova_lab.developability). For k = 3
and the default bias of 2.0 that gives 12 % at m = 83.4, 35 % at 86.3 and 63 %
at 88.6 -- against 1 %, 9 % and 42 % under the first repaired model, which
ignored the bias and took sigma at its low end. Both errors ran the same way:
they made a candidate look safer than the validator treats it.

## What it is for

Not a verdict -- a cost. A refused submission loses the epoch, and so does
submitting nothing, so the candidate to send is the one with the best expected
score, `predicted quality * P(not refused)`, not the one under some ceiling.
That ordering is what `expected_value` gives, and on our own 20 measured
candidates it changes the answer: the variant this project had settled on reads
88.6 (42 % chance of refusal) while FR211Q->E reads 84.4 (1 %) for the same
binding score within noise.
"""
from __future__ import annotations

from math import erf, sqrt

# The validator's own threshold, applied to its own measurement.
THRESHOLD = 89.0

# Deviation of ONE PSH measurement on the faithful install. Five estimates:
#
#   the 36 accepted submissions, run 1        1.74
#   the same 36, run 2                        2.49
#   our 20 v8/v9 candidates, run 1            3.29
#   the same 20, run 2                        2.70
#   between those two candidate runs          2.85   (paired, so independent
#                                                     of the within-run folds)
#
# The between-run figure matching the within-run ones is what settles it: the
# folds are genuinely independent, and 2.15 was the low end of a noisy
# estimate. Our own candidates are noisier than the accepted field, which fits
# -- they sit higher on PSH, on less stable structures. 2.7 is the conservative
# reading of the five.
#
# The consequence is worth stating plainly: with three folds the mean deviates
# by 1.6, so this gate separates "safe" from "doomed" and cannot separate 84
# from 86. Two independent three-fold runs of the same 20 candidates agree at
# Spearman 0.64, with one sequence moving 4.8 points. Do not read an ordering
# into gaps of a point or two.
SIGMA = 2.7

# What the accepted field reads, for context in logs and reports: the 36
# submissions the validators accepted in epoch 26174, three folds each.
ACCEPTED_MEDIAN = 83.5
ACCEPTED_MAX = 88.9


def refusal_probability(psh_mean: float, folds: int = 1,
                        bias: float | None = None) -> float:
    """P(the validator's own reading exceeds 89), given our k-fold mean.

    Both readings are draws around the sequence's true PSH, so the spread that
    matters combines our mean's deviation with the validator's single draw.

    `bias` is how far below theirs our scale appears to read, and it is added
    to our mean before the comparison. Three weak estimates put it at 2 to 3
    points and agree in sign; the default is the conservative end of that, and
    nova_lab.developability documents the evidence. Pass 0.0 for the raw
    comparison -- which is what this function did until the evidence turned up,
    and which makes a candidate look safer than it is.
    """
    if psh_mean != psh_mean:  # NaN: never measured, so nothing is known
        return float("nan")
    if bias is None:
        from nova_lab.developability import PSH_BIAS
        bias = PSH_BIAS
    spread = SIGMA * sqrt(1.0 + 1.0 / max(folds, 1))
    return 0.5 * (1.0 + erf((psh_mean + bias - THRESHOLD) / (spread * sqrt(2.0))))


def expected_value(quality: float, psh_mean: float, folds: int = 1) -> float:
    """`quality` discounted by the chance the submission is refused on PSH.

    `quality` is whatever we believe about the candidate on the merits -- the
    surrogate's predicted design_iiptm in the miner, a measured one when we
    have it. A candidate that cannot be measured keeps its quality untouched
    and is ranked below any measured candidate by `rank`, since an unscreened
    sequence is exactly what the production run submitted eleven times.
    """
    p = refusal_probability(psh_mean, folds)
    return quality if p != p else quality * (1.0 - p)


def rank(candidates: dict[str, float], psh: dict[str, float],
         folds: int = 1) -> list[str]:
    """Sequences best first, by expected value; measured ones ahead of not.

    `candidates` maps sequence -> quality, `psh` maps sequence -> our k-fold
    mean for those that were measured.
    """
    def key(sequence: str):
        measured = sequence in psh
        value = expected_value(candidates[sequence],
                               psh.get(sequence, float("nan")), folds)
        return (0 if measured else 1, -value)

    return sorted(candidates, key=key)


def verdict(psh_mean: float, folds: int = 1) -> str:
    """One line for a log: the reading, the risk, and what the field reads."""
    if psh_mean != psh_mean:
        return "PSH not measured — the gate that refused every production submission"
    p = refusal_probability(psh_mean, folds)
    where = ("below" if psh_mean <= ACCEPTED_MEDIAN else
             "above" if psh_mean > ACCEPTED_MAX else "within")
    return (f"PSH {psh_mean:.1f} over {folds} fold(s) -> {p:.0%} chance the "
            f"validator refuses it ({where} the accepted field, median "
            f"{ACCEPTED_MEDIAN}, max {ACCEPTED_MAX})")
