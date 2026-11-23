"""The validator's germline gate, as it is rather than as this repo guessed it.

What was wrong. `config.MIN_HUMAN_V_IDENTITY = 82` and
`MIN_HUMAN_J_IDENTITY = 86` are this lab's reconstruction. The filter
(`metanano/utils/igblast_nativeness.py`) hard-rejects on V identity under 0.70,
V coverage under 0.80 and a CDR3 outside 6-35, and has **no J floor at all**
(FINDINGS 67). Measured against sequences whose verdict we know (FINDINGS 69):

  - the hard rejects never fire -- 0 of 60 accepted submissions, 0 of 60 archive
    binders;
  - the nativeness threshold of 0.80 is **not** enforced: 0 of 60 ACCEPTED
    submissions reach it (median 0.507). Applying it would refuse everything the
    validators took;
  - the humanness threshold of 0.75 **is** the one that discriminates: 40 of 40
    accepted submissions pass, against 21 % of the P05231 archive's binders.

So the real germline gate is humanness >= 0.75, and nothing else.

Why a proxy. Humanness needs IgBLAST, about a second per sequence: twelve hours
for the 42 328 sequences of the P05231 archive. But abnumber's V identity, which
this lab already computes in milliseconds, tracks it at **Spearman +0.93** (J
identity does not: -0.15, which is why our J floor was both invented and
useless). Calibrated on 190 sequences spanning V 76.5 to 89.8:

| seuil V | precision | rappel | exactitude |
|---|---|---|---|
| 78.0 | 0.53 | 1.00 | 0.62 |
| **79.0** | **0.88** | **0.96** | **0.93** |
| 80.0 | 0.96 | 0.77 | 0.88 |
| 82.0 (le notre) | 0.94 | **0.57** | 0.79 |

V >= 79 is the best trade, and our 82 was discarding **43 %** of the sequences
the real gate accepts. The proxy is a screen, not a verdict: `verify` runs the
filter itself on what survives, which is the two-stage shape this lab already
uses for PSH.
"""
from __future__ import annotations

# Calibrated on 190 sequences against the filter's own humanness score.
HUMANNESS_THRESHOLD = 0.75
V_IDENTITY_PROXY = 79.0
# What the filter hard-rejects on. Reported for completeness: none of these has
# ever fired on a sequence this lab measured.
HARD_MIN_V_IDENTITY = 70.0
HARD_MIN_V_COVERAGE = 80.0
HARD_CDR3_RANGE = (6, 35)


def screens(v_identity: float | None) -> bool:
    """The fast screen: does this sequence plausibly clear humanness >= 0.75?

    96 % recall, 88 % precision. False positives are expected and are what
    `verify` is for; a false negative is a donor lost for nothing, which is why
    recall is what this threshold optimises.
    """
    return v_identity is not None and v_identity >= V_IDENTITY_PROXY


def hard_rejected(v_identity: float | None, v_coverage: float | None,
                  cdr3_length: int | None) -> str | None:
    """The filter's own hard gates, or None. Kept separate from the threshold:
    a hard reject is a malformed sequence, not a borderline one."""
    if v_identity is None or v_identity < HARD_MIN_V_IDENTITY:
        return "low_v_identity"
    if v_coverage is not None and v_coverage < HARD_MIN_V_COVERAGE:
        return "low_v_coverage"
    if cdr3_length is None:
        return "missing_cdr3"
    lo, hi = HARD_CDR3_RANGE
    if not lo <= cdr3_length <= hi:
        return "cdr3_length_out_of_range"
    return None


def verify(sequences: list[str]) -> dict[str, bool]:
    """Run the filter itself on a shortlist. {sequence: passes humanness}.

    Raises if the filter is not installed, rather than returning False for
    everything: `NativenessFilter.analyze` swallows its own exceptions and
    returns None, and reading that None as a refusal is what first scored four
    known-accepted submissions as refused (FINDINGS 69).
    """
    from scripts.nativeness import available, score
    usable, why = available()
    if not usable:
        raise RuntimeError(f"le filtre n'est pas executable : {why}")
    frame = score(sequences)
    return {row.sequence: (row.humanness is not None
                           and row.humanness >= HUMANNESS_THRESHOLD)
            for row in frame.itertuples()}
