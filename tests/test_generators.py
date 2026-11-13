"""The CDR3 generators' contract, and the one new generator's claim.

Every generator must return exactly `cdr3_length` residues and leave the last
`j_tail` of them alone: the end of the CDR3 is coded by the J gene, and
mutating it failed 285 candidates out of 800 before it was protected.

HydrophilicGenerator was built to attack one measured problem -- at CDR3 18 a
fresh chimeric candidate reads a median PSH of 100.6 where the submissions the
validators accept read 83.5. What these tests check is that it does what it
claims on the sequence. Whether that moves the measured PSH was a separate
experiment and the answer was no: 100 candidates per arm give 14/100 against
13/100 under the threshold (results/psh_generator_ab_cdr3_18.csv). The tests
below therefore guard a mechanism that works and a result that does not, which
is on purpose -- if someone changes the bias and the hydropathy claim breaks,
that is a bug; if the PSH claim is ever met, that is a discovery and it belongs
in the A/B, not here.

    python -m pytest tests/test_generators.py -q
"""
import random
import statistics

import pytest

from nova_lab.cdr3_generators import (AAS, GENERATORS, HYDROPATHY,
                                      ChimericGenerator, HydrophilicGenerator,
                                      make_generator)

DONOR = "ARVVLLIIFFWWYYSSGGFDY"  # deliberately greasy, so a bias has room to show
J_TAIL = 4


def hydropathy(cdr3: str) -> float:
    return sum(HYDROPATHY[a] for a in cdr3)


@pytest.mark.parametrize("name", sorted(GENERATORS - {"iglm", "pssm", "pssm_denovo"}))
@pytest.mark.parametrize("length", [5, 8, 18, 23])
def test_length_is_exact_and_the_j_tail_survives(name, length):
    generator = make_generator(name)
    rng = random.Random(0)
    for _ in range(20):
        out = generator.generate_cdr3(DONOR, length, J_TAIL, rng)
        assert len(out) == length
        expected_tail = DONOR[:length].ljust(length, "G")[length - J_TAIL:]
        assert out[length - J_TAIL:] == expected_tail
        assert set(out) <= set(AAS)


def test_beta_zero_draws_the_same_distribution_as_the_baseline():
    """The knob is honest: at 0 the bias is gone.

    Not the same *stream* -- weighted draws consume the generator differently,
    so the sequences differ -- but the same distribution, which is what matters:
    it means a measured PSH difference between the two comes from the bias and
    not from a different number of mutations or a different position rule.
    """
    plain, unbiased = ChimericGenerator(), HydrophilicGenerator(beta=0.0)
    rng_a, rng_b = random.Random(7), random.Random(7)
    base = [hydropathy(plain.generate_cdr3(DONOR, 18, J_TAIL, rng_a)) for _ in range(600)]
    same = [hydropathy(unbiased.generate_cdr3(DONOR, 18, J_TAIL, rng_b)) for _ in range(600)]
    # Both means sit on the same value; 1.5 is well inside the sampling error of
    # 600 draws (the bias at beta=0.4 moves the median by more than 2).
    assert abs(statistics.mean(base) - statistics.mean(same)) < 1.5

    # And the positions it touches are uniform again, not pulled to the greasy end.
    donor = "VVVVVLLLLL" + "SSSSSS" + "GFDY"
    rng = random.Random(5)
    greasy = total = 0
    for _ in range(300):
        out = unbiased.generate_cdr3(donor, 20, J_TAIL, rng)
        for position, (before, after) in enumerate(zip(donor, out)):
            if before != after:
                total += 1
                greasy += position < 10
    assert 0.4 < greasy / total < 0.75  # 10 of 16 mutable positions are greasy


def test_it_lowers_cdr3_hydropathy_against_the_baseline():
    """The claim, on the sequence: less of what PSH counts."""
    plain = ChimericGenerator()
    biased = HydrophilicGenerator()
    rng_a, rng_b = random.Random(3), random.Random(3)
    base = [hydropathy(plain.generate_cdr3(DONOR, 18, J_TAIL, rng_a)) for _ in range(300)]
    new = [hydropathy(biased.generate_cdr3(DONOR, 18, J_TAIL, rng_b)) for _ in range(300)]
    assert statistics.median(new) < statistics.median(base) - 2.0, (
        f"median hydropathy {statistics.median(new):.1f} against "
        f"{statistics.median(base):.1f}: the bias is not biting")


def test_it_spends_its_mutations_on_the_greasy_positions():
    """Where the mutations land, not just what replaces them.

    A donor whose first half is hydrophobic and second half polar: the biased
    generator should touch the first half more often. With beta at 0 this is a
    coin toss, which is what makes it a test of the position rule.
    """
    donor = "VVVVVLLLLL" + "SSSSSS" + "GFDY"  # 10 greasy, 6 polar, then the J tail
    biased = HydrophilicGenerator()
    rng = random.Random(5)
    greasy_hits = 0
    total = 0
    for _ in range(300):
        out = biased.generate_cdr3(donor, 20, J_TAIL, rng)
        for position, (before, after) in enumerate(zip(donor, out)):
            if before != after:
                total += 1
                greasy_hits += position < 10
    assert total > 0
    assert greasy_hits / total > 0.65, (
        f"only {greasy_hits}/{total} mutations landed on the hydrophobic half; "
        f"the position bias is not working")


def test_an_all_polar_donor_is_left_roughly_alone_in_hydropathy():
    """The bias cannot make a polar CDR3 much worse, which is the safety case."""
    donor = "SSTTNNQQDDEEGGSSFDY"
    biased = HydrophilicGenerator()
    rng = random.Random(11)
    values = [hydropathy(biased.generate_cdr3(donor, 18, J_TAIL, rng)) for _ in range(200)]
    assert statistics.median(values) < hydropathy(donor[:18]) + 4.0


def test_unknown_generator_names_are_refused():
    with pytest.raises(ValueError, match="unknown generator"):
        make_generator("hydrophobic-by-mistake")
