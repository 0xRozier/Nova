"""The one place a silent bug would flatter the ESM comparison.

`pools` builds two representations, one over the whole sequence and one over
the CDR3 alone. If a sequence whose CDR3 cannot be located fell back to the
whole-sequence mean, the `cdr3` column would quietly become a copy of the
`whole` column for those rows, and the two representations would agree because
they were the same numbers -- not because the model saw the same thing.
"""
import numpy as np
import pytest

from scripts.esm_surrogate import pools

GOOD = ("EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYRMYWVRQPPGKEREWVSAISAGGGSTYYGDSVKGRFTIS"
        "RDNAKNTVYLQMNSLKPEDTAVYYCANRAGWGMGDYWGQGTLVTVSS")
JUNK = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"


def fake_reps(sequences):
    """Row i of a sequence's representation is i, so a mean is checkable."""
    return {s: np.arange(len(s), dtype=np.float32).repeat(4).reshape(len(s), 4)
            for s in sequences}


def test_the_cdr3_pool_is_the_cdr3_and_not_the_whole_sequence():
    from nova_lab.nanobodies import regions
    reps = fake_reps([GOOD])
    out = pools([GOOD], reps)
    cdr3 = regions(GOOD)["CDR3"]
    at = GOOD.find(cdr3)
    assert out["cdr3"][0][0] == pytest.approx(np.arange(at, at + len(cdr3)).mean())
    assert out["whole"][0][0] == pytest.approx(np.arange(len(GOOD)).mean())
    assert out["cdr3"][0][0] != pytest.approx(out["whole"][0][0])


def test_an_unlocatable_cdr3_gives_zeros_not_a_silent_fallback():
    reps = fake_reps([JUNK])
    out = pools([JUNK], reps)
    assert np.all(out["cdr3"][0] == 0.0)
    assert out["whole"][0][0] > 0.0


def test_both_is_the_two_pools_side_by_side():
    reps = fake_reps([GOOD])
    out = pools([GOOD], reps)
    width = out["whole"].shape[1]
    assert out["both"].shape[1] == 2 * width
    assert np.allclose(out["both"][0, :width], out["whole"][0])
    assert np.allclose(out["both"][0, width:], out["cdr3"][0])
