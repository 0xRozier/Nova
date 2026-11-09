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


def test_catches_the_rejections_it_can_see(labelled, params, archive, exclude_by_epoch):
    """What we catch is germline identity, and that is most of what we can.

    Measured over 147 of the validator's rejections: we catch 18, and all 18 are
    human V identity below the floor. The rest are refused for reasons we do not
    model as hard gates — above all PSH, where our own measurement deviates by
    4.01 against a threshold of 89, so treating it as a verdict would refuse
    sequences the validator accepts (39% of the accepted set reads above 89
    here).

    The floor was 15% when this test was written, but that number was inflated
    by the anachronistic duplicates fixed in exclude_by_epoch: the checker was
    "catching" sequences that were novel at their epoch and only look archived
    today. 12% is the honest rate, so the floor sits at 8% — low enough not to
    be brittle, high enough that a broken germline check (which would drop it to
    roughly zero) still fails the test.
    """
    rejected = [(s, e) for s, e, accepted in labelled if not accepted]
    reasons = []
    for seq, epoch in rejected:
        problems = check_sequence(seq, params[epoch]["allowed_v"],
                                  params[epoch]["allowed_cdr3_length"],
                                  archive, exclude_by_epoch(epoch), warnings=[])
        if problems:
            reasons.append(problems[0])
    caught = len(reasons)
    kinds = Counter(r.split("(")[0].strip()[:40] for r in reasons)
    print(f"\ncaught {caught}/{len(rejected)} of the validator's rejections:")
    for kind, n in kinds.most_common():
        print(f"  {n:4d}  {kind}")
    print("  (the rest are refused on grounds we deliberately do not gate on, "
          "chiefly PSH)")

    assert caught >= 0.08 * len(rejected), (
        f"only {caught}/{len(rejected)} caught; the germline check may be broken")
    assert any("identity" in r for r in reasons), (
        "no germline-identity rejection caught at all — that check is the one "
        "kind we reproduce, so its absence means the reproduction is broken")
