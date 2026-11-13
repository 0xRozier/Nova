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
from collections import Counter

import pytest

from nova_lab import DATA_DIR
from nova_lab.config import NANOBODY_TARGET
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

    # submissions.jsonl only covers the epochs we decrypted, so it cannot
    # exclude everything the archive gained since. The archive carries no date
    # column, but data/prebuilt/<target>_regions.csv is keyed by sequence and
    # was built when the archive held 9466 entries, covering up to epoch 26025
    # — so it is a dated snapshot. Anything outside it arrived later and must
    # not count as a duplicate when replaying an earlier epoch. Without this,
    # refreshing the archive to 12842 made the checker reject 4 of 85 sequences
    # the validators had accepted.
    snapshot, snapshot_epoch = set(), 26025
    snap = DATA_DIR / "prebuilt" / f"{NANOBODY_TARGET[0]}_regions.csv"
    if snap.exists():
        import pandas as pd
        frame = pd.read_csv(snap, keep_default_na=False)
        if "sequence" in frame.columns:
            snapshot = set(frame["sequence"])

    current = set(NanobodyArchive(with_regions=False).sequences)
    arrived_later = (current - snapshot) if snapshot else set()

    def exclude(epoch: int) -> set:
        later = {s for e, seqs in by_epoch.items() if e >= epoch for s in seqs}
        if epoch <= snapshot_epoch:
            later |= arrived_later
        return later

    return exclude


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


def test_the_germline_floor_is_still_enforced(archive):
    """The one rejection kind this checker reproduces, probed directly.

    This replaces an assertion that at least 8 % of the validator's rejections
    in a sample of 100 were caught here. That number was a proxy for "the
    germline check still works", and it was measured on epochs 25956-26025,
    where germline failures were common. It is not a property of our code: on
    epochs 26299-26306 the field has cleaned up — 66 rejections in 296
    submissions, and we catch 1 — so the proxy failed while the check it stood
    for was fine. Today's rejections are PSH and developability, which
    check_sequence deliberately does not model because they need a structure;
    scripts.preflight is where those are screened.

    So the mechanism is tested on a sequence built to fail it, which cannot
    drift with the field.
    """
    from nova_lab.nanobodies import germline_identity
    parent = next(s for s in archive.sequences if germline_identity(s))
    mutated = list(parent)
    for i in range(30, min(len(mutated), 80), 3):
        mutated[i] = "W"  # wreck the framework, leave the ends alone
    mutated = "".join(mutated)

    allowed_v = germline_identity(parent)["V"]["genes"][0]
    problems = check_sequence(mutated, allowed_v, None, archive, set(), warnings=[])
    assert any("identity" in p for p in problems), (
        f"a sequence mutated far off the human germline was not refused: "
        f"{problems or 'no problem reported at all'}. The germline check is the "
        f"one hard gate this checker reproduces; if it is silent, ANARCI is "
        f"returning None (run python -m scripts.doctor: without hmmscan it does "
        f"that without raising).")


def test_it_flags_nothing_the_validator_was_happy_with(labelled, params, archive,
                                                       exclude_by_epoch):
    """Precision, which is the direction that costs us epochs.

    A false accept wastes one submission; a false refusal makes the miner skip
    a candidate, or an epoch, for nothing. The recall figure is printed for
    information: what we catch now is whatever the field happens to be getting
    wrong that week, and it is no longer a property of this code.
    """
    rejected = [(s, e) for s, e, accepted in labelled if not accepted]
    caught = []
    for seq, epoch in rejected:
        problems = check_sequence(seq, params[epoch]["allowed_v"],
                                  params[epoch]["allowed_cdr3_length"],
                                  archive, exclude_by_epoch(epoch), warnings=[])
        if problems:
            caught.append(problems[0])
    print(f"\ncaught {len(caught)}/{len(rejected)} of the validator's rejections:")
    for kind, n in Counter(r.split("(")[0].strip()[:40] for r in caught).most_common():
        print(f"  {n:4d}  {kind}")
    print("  (the rest are refused on grounds that need a structure — chiefly "
          "PSH > 89, which scripts.preflight screens)")

    accepted = [(s, e) for s, e, ok in labelled if ok]
    wrong = []
    for seq, epoch in accepted:
        problems = check_sequence(seq, params[epoch]["allowed_v"],
                                  params[epoch]["allowed_cdr3_length"],
                                  archive, exclude_by_epoch(epoch), warnings=[])
        if problems:
            wrong.append((epoch, problems[0]))
    assert not wrong, f"{len(wrong)}/{len(accepted)} false refusals, e.g. {wrong[:3]}"
