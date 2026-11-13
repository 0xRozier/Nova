"""The two numbers every epoch rests on, and the change announced for 26331.

`challenge_params` mirrors upstream's `get_challenge_params_from_blockhash`. If
it drifts, every candidate is built for the wrong challenge and every
submission is refused -- which looks exactly like a lost epoch, so nothing
would say so.

The subnet's Discord announced on 2026-11-13 that PR 102 merges at epoch 26331
and warned nanobody miners that it "will also affect how valid germlines and
CDR3 lengths are chosen each epoch due to RNG state consumption". The code on
the rxn-group-filtering branch does not: both the old and the new reaction
selection are arithmetic on the seed and never touch the `rng` object, which is
consumed afterwards by the V-gene and CDR3-length draws, in the same order.

These tests pin that down, so that if the merged version differs the suite says
which half moved. What they cannot settle is which source is right -- for that
there is `scripts.verify_challenge`, which asks the accepted submissions
themselves and needs no trust in either.

    python -m pytest tests/test_challenge_params.py -q
"""
import random
import secrets

from nova_lab.config import (CDR3_LENGTH_RANGE, RXN_GROUPS,
                             RXN_GROUP_FILTERING_FROM_EPOCH,
                             V_NEIGHBORHOOD_CATEGORIES)
from scripts.epoch_info import challenge_params

N_REACTIONS = 44
BEFORE = RXN_GROUP_FILTERING_FROM_EPOCH - 1
AFTER = RXN_GROUP_FILTERING_FROM_EPOCH


def hashes(n: int, seed: int = 0) -> list[str]:
    rng = random.Random(seed)
    return ["0x" + "".join(rng.choice("0123456789abcdef") for _ in range(64))
            for _ in range(n)]


def test_pr102_leaves_the_nanobody_parameters_alone():
    """The claim in the announcement, checked on 20 000 block hashes.

    If this ever fails, the merged PR consumes the RNG differently from the
    branch we read and every nanobody parameter we derive is wrong from that
    epoch on. It is the single most expensive thing in this file.
    """
    for block_hash in hashes(20_000):
        old = challenge_params(block_hash, N_REACTIONS, epoch=BEFORE)
        new = challenge_params(block_hash, N_REACTIONS, epoch=AFTER)
        assert old["allowed_v"] == new["allowed_v"]
        assert old["allowed_cdr3_length"] == new["allowed_cdr3_length"]


def test_the_reaction_rule_does_change():
    """And the half that does move, moves: one reaction becomes several."""
    block_hash = hashes(1)[0]
    before = challenge_params(block_hash, N_REACTIONS, epoch=BEFORE)
    after = challenge_params(block_hash, N_REACTIONS, epoch=AFTER)
    assert isinstance(before["allowed_reaction"], str)
    assert isinstance(after["allowed_reaction"], list)
    assert len(after["allowed_reaction"]) == len(RXN_GROUPS)


def test_without_an_epoch_the_old_rule_applies():
    """So every replay of a past epoch stays reproducible."""
    block_hash = hashes(1)[0]
    assert (challenge_params(block_hash, N_REACTIONS)["allowed_reaction"]
            == challenge_params(block_hash, N_REACTIONS, epoch=BEFORE)["allowed_reaction"])


def test_each_group_contributes_exactly_one_reaction():
    for block_hash in hashes(500, seed=1):
        picked = challenge_params(block_hash, N_REACTIONS, epoch=AFTER)["allowed_reaction"]
        ids = [int(name.removeprefix("rxn:")) for name in picked]
        assert len(ids) == len(RXN_GROUPS)
        for reaction, group in zip(ids, RXN_GROUPS):
            assert reaction in group


def test_the_groups_cover_every_reaction_exactly_once():
    """The branch's comment says "together these lists cover reactions 1-44"."""
    flat = [reaction for group in RXN_GROUPS for reaction in group]
    assert sorted(flat) == list(range(1, 45))


def test_the_parameters_stay_inside_what_the_config_allows():
    genes = {gene for group in V_NEIGHBORHOOD_CATEGORIES.values() for gene in group}
    for block_hash in hashes(2_000, seed=2):
        for epoch in (BEFORE, AFTER):
            params = challenge_params(block_hash, N_REACTIONS, epoch=epoch)
            assert params["allowed_v"] in genes
            assert params["allowed_cdr3_length"] in CDR3_LENGTH_RANGE


def test_the_same_hash_always_gives_the_same_answer():
    """It is derived, not drawn: two calls must agree, or replays are nonsense."""
    for block_hash in hashes(200, seed=3):
        assert (challenge_params(block_hash, N_REACTIONS, epoch=AFTER)
                == challenge_params(block_hash, N_REACTIONS, epoch=AFTER))
