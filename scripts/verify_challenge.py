"""Is our (V gene, CDR3 length) derivation still the validator's? Ask the field.

Everything this miner does in an epoch rests on two numbers derived from the
epoch's first block hash. If that derivation drifts from upstream's, every
candidate is built for the wrong challenge and every submission is refused --
silently, since a refused submission looks exactly like a lost epoch.

We normally guard that by reading upstream's code (`scripts.check_config`). That
is not enough twice over: the values can live on a branch that `main` does not
show yet, and the announcement can be wrong about its own effects. On
2026-11-13 the subnet's Discord warned nanobody miners that PR 102 would change
how germlines and CDR3 lengths are drawn "due to RNG state consumption"; the
branch's code does no such thing, and 200 000 simulated hashes agree it does
not. One of those two sources is mistaken and we cannot tell which from the
outside.

This script needs neither. **A nanobody submission is only accepted if its V
gene and CDR3 length match the epoch's**, so every accepted submission is a
witness to the real parameters. If our derivation says IGHV3-66/13 and thirty
accepted sequences say IGHV3-74/18, the derivation is wrong and that is the end
of the argument.

What it cannot do is see an epoch before it has been scored, so it is a trailing
check: run it after `scripts.fetch_submissions` and `scripts.analyze_nanobodies`,
and it tells you within one epoch of a change. Run it before engaging TAO after
any upstream change, and at 26331 in particular.

Usage:
    python -m scripts.verify_challenge               # every epoch we have
    python -m scripts.verify_challenge --last 10
"""
import argparse
import json
import time
from collections import Counter

import pandas as pd

from nova_lab import DATA_DIR
from nova_lab.config import EPOCH_LENGTH, NANOBODY_TARGET
from nova_lab.nanobodies import germline_identity, regions, seq_hash


def accepted_by_epoch(target: str) -> dict[int, list[str]]:
    """The sequences the validator scored, per epoch: our witnesses."""
    results = DATA_DIR / "submissions" / f"{target}_nanobody_results.csv"
    decrypted = DATA_DIR / "submissions" / "submissions.jsonl"
    if not results.exists() or not decrypted.exists():
        raise SystemExit(
            f"missing {results.name} or {decrypted.name}: run "
            f"scripts.fetch_submissions then scripts.analyze_nanobodies first")
    scored = set(pd.read_csv(results)["hash"])
    out: dict[int, list[str]] = {}
    for line in decrypted.open():
        row = json.loads(line)
        if row["sequences"] and seq_hash(row["sequences"][0]) in scored:
            out.setdefault(int(row["epoch"]), []).append(row["sequences"][0])
    return out


def observed(sequences: list[str]) -> tuple[Counter, Counter]:
    """What the accepted sequences themselves say the parameters were."""
    genes, lengths = Counter(), Counter()
    for sequence in sequences:
        region = regions(sequence)
        germline = germline_identity(sequence)
        if region:
            lengths[len(region.get("CDR3", ""))] += 1
        if germline and germline["V"]["genes"]:
            # A sequence can match several germlines; any of them may be the
            # epoch's, so all are counted and the check asks for an overlap.
            for gene in germline["V"]["genes"]:
                genes[gene] += 1
    return genes, lengths


def derived(epoch: int) -> dict:
    """Our own derivation, with retries: the archive node answers 429 under load."""
    import bittensor as bt

    from nova_lab.reactions import reaction_ids
    from scripts.epoch_info import challenge_params
    for attempt in range(5):
        try:
            st = bt.Subtensor("wss://archive.chain.opentensor.ai:443")
            block_hash = st.determine_block_hash(epoch * EPOCH_LENGTH)
            return challenge_params(block_hash, len(reaction_ids()), epoch=epoch)
        except Exception as error:  # noqa: BLE001
            if attempt == 4:
                raise
            print(f"  epoch {epoch}: {type(error).__name__}, retrying", flush=True)
            time.sleep(5 * 2 ** attempt)
    raise AssertionError("unreachable")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default=NANOBODY_TARGET[0])
    parser.add_argument("--last", type=int, help="only the N most recent epochs")
    args = parser.parse_args()

    witnesses = accepted_by_epoch(args.target)
    epochs = sorted(witnesses)
    if args.last:
        epochs = epochs[-args.last:]
    if not epochs:
        raise SystemExit("no scored submissions on disk to check against")

    print(f"checking our derivation against {sum(len(witnesses[e]) for e in epochs)} "
          f"accepted submissions over {len(epochs)} epochs\n")
    print(f"{'epoch':>6}  {'ours':<22} {'the field':<26} verdict")
    bad = []
    for epoch in epochs:
        params = derived(epoch)
        genes, lengths = observed(witnesses[epoch])
        if not genes or not lengths:
            print(f"{epoch:>6}  {'(not numberable)':<22}")
            continue
        gene_ok = params["allowed_v"] in genes
        length_ok = params["allowed_cdr3_length"] in lengths
        theirs = (f"{max(genes, key=genes.get)}/"
                  f"{max(lengths, key=lengths.get)} on {len(witnesses[epoch])}")
        mine = f"{params['allowed_v']}/{params['allowed_cdr3_length']}"
        verdict = "ok" if gene_ok and length_ok else "MISMATCH"
        if not (gene_ok and length_ok):
            bad.append((epoch, mine, dict(genes), dict(lengths)))
        print(f"{epoch:>6}  {mine:<22} {theirs:<26} {verdict}")

    print()
    if bad:
        print(f"{len(bad)} of {len(epochs)} epochs DISAGREE. Our derivation is not "
              f"the validator's any more, so every candidate we build is for the "
              f"wrong challenge. Do not submit. Re-read "
              f"utils/challenge.get_challenge_params_from_blockhash upstream and "
              f"nova_lab.config.RXN_GROUP_FILTERING_FROM_EPOCH.")
        for epoch, mine, genes, lengths in bad[:5]:
            print(f"  epoch {epoch}: ours {mine}, accepted sequences carry "
                  f"{genes} and CDR3 lengths {lengths}")
        return 1
    print(f"all {len(epochs)} epochs agree: the derivation is still theirs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
