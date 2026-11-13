"""Show the challenge parameters of the current (or a given) NOVA epoch.

Mirrors utils/challenge.get_challenge_params_from_blockhash() upstream: the hash
of the epoch's first block seeds everything, so the parameters are known as soon
as the epoch starts.

Usage:
    python -m scripts.epoch_info              # current epoch
    python -m scripts.epoch_info --epoch 26005
    python -m scripts.epoch_info --last 24    # table of the last 24 epochs
"""
import argparse
import random

import bittensor as bt

from nova_lab.chain import current_block
from nova_lab.config import (CDR3_LENGTH_RANGE, EPOCH_LENGTH, NO_SUBMISSION_BLOCKS,
                             RXN_GROUPS, RXN_GROUP_FILTERING_FROM_EPOCH,
                             V_NEIGHBORHOOD_CATEGORIES)
from nova_lab.reactions import reaction_ids

NETUID = 68


def challenge_params(block_hash: str, n_reactions: int,
                     epoch: int | None = None) -> dict:
    """The epoch's imposed reaction(s), V gene and CDR3 length, from one hash.

    The order matters and is upstream's: the reaction is derived by arithmetic
    on the seed, then `rng` is consumed by the V-gene draw and the CDR3-length
    draw. Nothing before those two draws touches `rng`, which is why PR 102
    changes the reaction and leaves the nanobody parameters alone -- see the
    note on RXN_GROUP_FILTERING_FROM_EPOCH in nova_lab.config.

    `epoch` selects which reaction rule applies. Without it the pre-26331 rule
    is used, which keeps every replay of a past epoch reproducible; pass the
    epoch to get the rule in force then.
    """
    seed = int(block_hash[2:], 16)
    rng = random.Random(seed)

    if epoch is not None and epoch >= RXN_GROUP_FILTERING_FROM_EPOCH:
        # One reaction from each group, so an epoch allows several.
        allowed_reaction = [f"rxn:{group[seed % len(group)]}" for group in RXN_GROUPS]
    else:
        allowed_reaction = f"rxn:{seed % n_reactions + 1}"

    genes = V_NEIGHBORHOOD_CATEGORIES[rng.choice(list(V_NEIGHBORHOOD_CATEGORIES))]
    allowed_v = rng.choice(genes) if len(genes) > 1 else genes[0]
    return {
        "allowed_reaction": allowed_reaction,
        "allowed_v": allowed_v,
        "allowed_cdr3_length": rng.choice(list(CDR3_LENGTH_RANGE)),
    }


def epoch_params(subtensor, epoch: int, n_reactions: int) -> dict:
    start = epoch * EPOCH_LENGTH
    params = challenge_params(subtensor.determine_block_hash(start), n_reactions,
                              epoch=epoch)
    return {"epoch": epoch, "start_block": start, **params}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--epoch", type=int)
    parser.add_argument("--last", type=int)
    parser.add_argument("--network", default="finney")
    args = parser.parse_args()

    st = bt.Subtensor(args.network)
    block = current_block(st)
    current = block // EPOCH_LENGTH
    n_reactions = len(reaction_ids())

    if args.last:
        for epoch in range(current - args.last + 1, current + 1):
            p = epoch_params(st, epoch, n_reactions)
            print(f"{p['epoch']}  {p['allowed_reaction']:>7}  {p['allowed_v']:<9}  CDR3={p['allowed_cdr3_length']}")
        return

    epoch = args.epoch if args.epoch is not None else current
    p = epoch_params(st, epoch, n_reactions)
    end = (epoch + 1) * EPOCH_LENGTH
    print(f"current block       {block}")
    print(f"epoch               {epoch} (blocks {p['start_block']} -> {end})")
    if epoch == current:
        left = end - NO_SUBMISSION_BLOCKS - block
        print(f"submission window   {left} blocks left (~{left * 12 // 60} min)")
    print(f"allowed reaction    {p['allowed_reaction']}")
    print(f"allowed V gene      {p['allowed_v']}")
    print(f"allowed CDR3 length {p['allowed_cdr3_length']}")


if __name__ == "__main__":
    main()
