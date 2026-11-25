"""The learning batch's design guarantees.

Two ways to waste 2.6 h of GPU here: query a sequence whose answer we already
have, or let the "coverage" arm be a random sample in disguise. A random sample
would re-draw the archive's own shape instead of correcting the lopsidedness of
our live-target set, where 30% of 132 sequences share one framework.
"""
import numpy as np
import pandas as pd
import pytest

from nova_lab.config import NANOBODY_STRUCTURE
from scripts.build_learning_batch import (SEED, build, candidate_pool,
                                          farthest_points, project)
from scripts.train_score_surrogate import load

# The shipped batch uses a 4 000-sequence pool; the tests exercise the same
# code paths on a tenth of it, because the farthest-point pass is quadratic and
# a six-minute test file is a tax on every future run.
POOL, N_COVER, N_EXPLOIT = 400, 12, 6
SMALL = dict(pool_size=POOL, n_cover=N_COVER, n_exploit=N_EXPLOIT)


@pytest.fixture(scope="module")
def known():
    return load(NANOBODY_STRUCTURE)


@pytest.fixture(scope="module")
def batch():
    return build(**SMALL)


def test_no_sequence_is_queried_twice(batch, known):
    """A repeat is paid for and not measured: BoltzGen keys by sequence hash."""
    assert batch.sequence.duplicated().sum() == 0
    assert not set(batch.sequence) & set(known.sequence)


def test_both_arms_are_present_and_sized(batch):
    counts = batch.arm.value_counts()
    assert counts["cover"] == N_COVER
    assert counts["exploit"] == N_EXPLOIT


def test_the_exploit_arm_is_predicted_clearly_higher(batch):
    """Without separation there is no prospective test of the surrogate."""
    cover = batch[batch.arm == "cover"].predicted
    exploit = batch[batch.arm == "exploit"].predicted
    assert exploit.min() > cover.median()


def test_coverage_picks_are_further_out_than_a_random_sample(known):
    """The whole point of farthest-point sampling, measured rather than assumed."""
    pool = candidate_pool(set(known.sequence), pool=POOL)
    space = project(list(known.sequence) + list(pool.sequence))
    seeds = len(known)
    chosen = farthest_points(space, seeds, N_COVER)

    to_known = lambda idx: np.linalg.norm(
        space[seeds + np.asarray(idx), None, :] - space[None, :seeds, :],
        axis=2).min(axis=1)
    picked = to_known(chosen).mean()
    rng = np.random.default_rng(SEED)
    random = np.mean([to_known(rng.choice(len(pool), N_COVER, replace=False)).mean()
                      for _ in range(5)])
    assert picked > 1.5 * random, f"couverture {picked:.2f} contre hasard {random:.2f}"


def test_farthest_points_never_returns_a_seed():
    """Seeds only push the picks away; returning one would re-query a known answer."""
    rng = np.random.default_rng(0)
    space = rng.normal(size=(60, 4))
    chosen = farthest_points(space, n_seeds=20, take=10)
    assert len(chosen) == len(set(chosen)) == 10
    assert all(0 <= i < 40 for i in chosen)


def test_the_batch_is_deterministic():
    assert list(build(**SMALL).sequence) == list(build(**SMALL).sequence)


def test_layer_one_is_recorded_but_not_used_to_drop_rows(batch):
    """A refused sequence is still a valid measurement of the oracle."""
    assert "layer1" in batch.columns
    assert (batch.layer1 != "ok").any(), (
        "si tout passe, le lot a ete filtre et ne couvre plus l'espace voulu")
