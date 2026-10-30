"""Check the offline nanobody checker against the validators' own verdicts.

A submitted sequence has published BoltzGen metrics for its epoch only if the
validator accepted it, so the decrypted submissions give us labelled data.

Our checker runs without structure prediction, so it cannot see the
developability rejections: it may accept a sequence the validator refused
(false accept), but it must never reject one the validator accepted.

    python -m pytest tests/test_nanobody_parity.py -q

Needs data/submissions/Q9NZQ7_nanobody_results.csv (scripts.analyze_nanobodies)
and the decrypted submissions (scripts.fetch_submissions).
"""
import json
import random

import pytest

from nova_lab import DATA_DIR
from nova_lab.nanobodies import NanobodyArchive, check_sequence, seq_hash

SAMPLE = 100


def epoch_params(epoch: int) -> dict:
    """The epoch's allowed V gene and CDR3 length, from its first block hash."""
    from scripts.epoch_info import challenge_params
    import bittensor as bt
    from nova_lab.config import EPOCH_LENGTH
    from nova_lab.reactions import reaction_ids
    st = bt.Subtensor("wss://archive.chain.opentensor.ai:443")
    return challenge_params(st.determine_block_hash(epoch * EPOCH_LENGTH), len(reaction_ids()))


@pytest.fixture(scope="module")
def labelled():
    """[(sequence, epoch, accepted_by_validator)], one row per submission."""
    import pandas as pd
    path = DATA_DIR / "submissions" / "Q9NZQ7_nanobody_results.csv"
    if not path.exists():
        pytest.skip("run scripts.analyze_nanobodies first")
    valid = set(pd.read_csv(path)["hash"])
    rows = []
    for line in (DATA_DIR / "submissions" / "submissions.jsonl").open():
        s = json.loads(line)
        if s["sequences"]:
            h = seq_hash(s["sequences"][0])
            rows.append((s["sequences"][0], s["epoch"], h in valid))
    random.Random(0).shuffle(rows)
    return rows[:SAMPLE]


@pytest.fixture(scope="module")
def params(labelled):
    return {e: epoch_params(e) for e in {e for _, e, _ in labelled}}


@pytest.fixture(scope="module")
def archive():
    return NanobodyArchive()


@pytest.fixture(scope="module")
def exclude_by_epoch():
    """Sequences submitted at or after an epoch: already archived today, but not
    when that epoch ran."""
    by_epoch = {}
    for line in (DATA_DIR / "submissions" / "submissions.jsonl").open():
        s = json.loads(line)
        if s["sequences"]:
            by_epoch.setdefault(s["epoch"], set()).add(s["sequences"][0])
    return lambda epoch: {s for e, seqs in by_epoch.items() if e >= epoch for s in seqs}


def test_never_rejects_a_sequence_the_validator_accepted(labelled, params, archive, exclude_by_epoch):
    wrong = []
    for seq, epoch, accepted in labelled:
        if not accepted:
            continue
        p = params[epoch]
        problems = check_sequence(seq, p["allowed_v"], p["allowed_cdr3_length"],
                                  archive, exclude_by_epoch(epoch), warnings=[])
        if problems:
            wrong.append((epoch, seq[:30], problems))
    assert not wrong, f"{len(wrong)} false rejections, e.g. {wrong[:3]}"


def test_catches_the_rejections_it_can_see(labelled, params, archive, exclude_by_epoch):
    """Structure-based developability rejections are invisible to us; the
    sequence-level ones (germline, V gene, CDR3 length) must be caught."""
    rejected = [(s, e) for s, e, accepted in labelled if not accepted]
    caught = sum(bool(check_sequence(s, params[e]["allowed_v"], params[e]["allowed_cdr3_length"],
                                     archive, exclude_by_epoch(e), warnings=[]))
                 for s, e in rejected)
    print(f"\ncaught {caught}/{len(rejected)} of the validator's rejections "
          f"(the rest need structure-based developability)")
    assert caught >= 0.15 * len(rejected)
