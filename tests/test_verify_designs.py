"""The guard that was missing for two design runs.

210 designs over two runs were read as evidence that structure-conditioned
design does not work. None of them had a single cysteine: the canonical
Cys23-Cys104 disulfide had been stripped, so they could not fold like a
nanobody and the validator would have refused them on min_cysteines. The low
scores measured a broken pipeline and were published as a refutation.

    python -m pytest tests/test_verify_designs.py -q
"""
import pandas as pd

from scripts.verify_designs import check, cysteines, frozen_intact

GOOD = ("EVQLVESGGGLVQPGGSLRLSCAASGRTFSYNPMGWFRQAPGKGRELVAAISRTGGSTYYPDSVEGRFT"
        "ISRDNAKRMVYLQMNSLRAEDTAVYYCAAAGVRAEDGRVRTLPSEYTFWGQGTQVTVSS")
# what the vendored BoltzGen actually returned: Cys23 -> Glu, no cysteine left
BROKEN = GOOD.replace("C", "E")


def test_it_catches_the_run_that_was_published_as_a_refutation():
    report = check(pd.DataFrame({"s": [BROKEN] * 3}), "s")
    assert report["ok"] is False
    assert report["no_cysteine"] == 3
    assert report["below_validator_floor"] == 3


def test_a_real_nanobody_passes():
    report = check(pd.DataFrame({"s": [GOOD] * 3}), "s")
    assert report["ok"] is True
    assert report["median_cysteines"] == cysteines(GOOD) == 2


def test_one_cysteine_clears_the_validator_but_not_the_canonical_pair():
    """min_cysteines is typically 1, but Cys23-Cys104 needs two. A design with
    one would be accepted and still fold wrong, so the guard refuses it."""
    single = GOOD.replace("C", "E", 1)
    report = check(pd.DataFrame({"s": [single]}), "s")
    assert report["below_validator_floor"] == 0
    assert report["below_canonical_pair"] == 1
    assert report["ok"] is False


def test_it_notices_a_framework_the_spec_had_frozen_drifting():
    drifted = "QEQ" + GOOD[3:]
    report = check(pd.DataFrame({"s": [GOOD, drifted]}), "s", scaffolds=[GOOD])
    assert report["frozen_broken"] == 1 and report["ok"] is False


def test_several_fr1_variants_are_reported_even_without_a_scaffold():
    """The published run had four FR1 variants where the spec imposed one, and
    nobody looked. Reported whether or not a scaffold is passed."""
    report = check(pd.DataFrame({"s": [GOOD, "QEQ" + GOOD[3:], "EVK" + GOOD[3:]]}), "s")
    assert report["distinct_fr1"] == 3


def test_frozen_intact_compares_only_what_the_spec_froze():
    """The CDR3 is designed, so a different CDR3 must not count as drift."""
    other_cdr3 = GOOD[:100] + "WWWWWWWW" + GOOD[108:]
    assert frozen_intact(other_cdr3, GOOD) is True
