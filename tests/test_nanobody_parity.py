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

# data/prebuilt/<target>_regions.csv was frozen with the archive at 9466
# entries, which is every sequence up to this epoch.
SNAPSHOT_EPOCH = 26025


PARAMS_CACHE = DATA_DIR / "submissions" / "epoch_params_cache.json"


def epoch_params(epoch: int) -> dict:
    """The epoch's allowed V gene and CDR3 length, from its first block hash.

    Cached in `data/submissions/epoch_params_cache.json`, because a block hash
    is immutable: once an epoch is read it never needs the chain again. Without
    the cache this file alone took 23 minutes and failed at random on the public
    archive node's HTTP 429 -- the same rate limit that killed the production
    miner for five hours. A test that fails on someone else's quota teaches
    nothing.

    A cache miss still goes to the chain, with backoff, and is written back.
    """
    import time

    cache = {}
    if PARAMS_CACHE.exists():
        try:
            cache = json.loads(PARAMS_CACHE.read_text())
        except json.JSONDecodeError:  # a truncated write is not worth failing on
            cache = {}
    if str(epoch) in cache:
        return cache[str(epoch)]

    from scripts.epoch_info import challenge_params
    import bittensor as bt
    from nova_lab.config import EPOCH_LENGTH
    from nova_lab.reactions import reaction_ids

    last = None
    for attempt in range(5):
        try:
            st = bt.Subtensor("wss://archive.chain.opentensor.ai:443")
            block_hash = st.determine_block_hash(epoch * EPOCH_LENGTH)
            params = challenge_params(block_hash, len(reaction_ids()))
            cache[str(epoch)] = {"allowed_v": params["allowed_v"],
                                 "allowed_cdr3_length": params["allowed_cdr3_length"]}
            PARAMS_CACHE.parent.mkdir(parents=True, exist_ok=True)
            PARAMS_CACHE.write_text(json.dumps(cache, indent=1, sort_keys=True))
            return cache[str(epoch)]
        except Exception as error:  # noqa: BLE001 - any transport failure retries
            last = error
            time.sleep(5 * 2 ** attempt)
    pytest.skip(f"the archive node would not answer for epoch {epoch}, and it is "
                f"not in {PARAMS_CACHE.name}: {type(last).__name__}: {last}")


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
    """Sequences already in today's archive that did not exist at an epoch.

    The validator's novelty gate compares a submission against the archive **as
    it stood then**, so replaying an epoch means subtracting everything that
    arrived later. Three sources date the archive, and together they do not
    cover it:

    - `submissions.jsonl`, for the epochs we decrypted;
    - `data/prebuilt/<target>_regions.csv`, keyed by sequence and frozen when
      the archive held 9466 entries, i.e. everything up to epoch 26025. What is
      outside it arrived later;
    - `<target>_nanobody_results.csv`, which carries an epoch per scored
      submission and dates the 230 that FINDINGS 48 merged into the pool.

    What remains undated is about 7000 sequences that the archive fetch brought
    in from epochs after 26025, with no date anywhere in the data. So the
    archive **cannot be reconstructed for an epoch after 26025**, and the
    novelty gate cannot be replayed there -- see `novelty_is_replayable`.
    """
    by_epoch = {}
    for line in (DATA_DIR / "submissions" / "submissions.jsonl").open():
        s = json.loads(line)
        if s["sequences"]:
            by_epoch.setdefault(s["epoch"], set()).add(s["sequences"][0])

    import pandas as pd

    snapshot = set()
    snap = DATA_DIR / "prebuilt" / f"{NANOBODY_TARGET[0]}_regions.csv"
    if snap.exists():
        frame = pd.read_csv(snap, keep_default_na=False)
        if "sequence" in frame.columns:
            snapshot = set(frame["sequence"])

    current = set(NanobodyArchive(with_regions=False).sequences)
    arrived_later = (current - snapshot) if snapshot else set()

    dated = {}
    results = DATA_DIR / "submissions" / f"{NANOBODY_TARGET[0]}_nanobody_results.csv"
    if results.exists():
        frame = pd.read_csv(results)
        if {"sequence", "epoch"} <= set(frame.columns):
            dated = frame.groupby("sequence")["epoch"].min().to_dict()

    def exclude(epoch: int) -> set:
        later = {s for e, seqs in by_epoch.items() if e >= epoch for s in seqs}
        if epoch <= SNAPSHOT_EPOCH:
            later |= arrived_later
        later |= {s for s, e in dated.items() if e >= epoch}
        return later

    return exclude


def novelty_is_replayable(epoch: int) -> bool:
    """Whether today's archive can be rewound to this epoch.

    Only up to the dated snapshot. Past it, roughly 7000 archive entries carry
    no date in any of our data, so a sequence that was novel at epoch 26300 may
    sit 95 % away from something that arrived at 26400 -- and the gate would
    refuse it for a reason that did not exist then. Asserting on novelty there
    would test the age of our archive, not this code.
    """
    return epoch <= SNAPSHOT_EPOCH


def without_novelty(problems: list[str]) -> list[str]:
    return [p for p in problems if "archived sequence" not in p]


def test_never_rejects_a_sequence_the_validator_accepted(labelled, params, archive,
                                                        exclude_by_epoch):
    """No false refusal on a submission the validators took.

    Novelty is asserted only where the archive can be rewound
    (`novelty_is_replayable`); elsewhere a duplicate verdict says our archive
    is younger than the epoch, not that the gate is wrong.
    """
    wrong, unverifiable = [], 0
    for seq, epoch, accepted in labelled:
        if not accepted:
            continue
        p = params[epoch]
        problems = check_sequence(seq, p["allowed_v"], p["allowed_cdr3_length"],
                                  archive, exclude_by_epoch(epoch), warnings=[])
        if not novelty_is_replayable(epoch):
            unverifiable += len(problems) - len(without_novelty(problems))
            problems = without_novelty(problems)
        if problems:
            wrong.append((epoch, seq[:30], problems))
    if unverifiable:
        print(f"\n{unverifiable} novelty verdicts left unasserted: their epochs "
              f"are past {SNAPSHOT_EPOCH} and the archive cannot be rewound there")
    assert not wrong, f"{len(wrong)} false rejections, e.g. {wrong[:3]}"


def test_the_novelty_gate_is_exact_where_the_archive_can_be_rewound(
        labelled, params, archive, exclude_by_epoch):
    """The gate itself, on the epochs the snapshot dates.

    This is the half of the guarantee the test above gives up. If no sampled
    epoch is old enough there is nothing to assert, and saying so is better
    than passing silently.
    """
    old = [(s, e) for s, e, ok in labelled if ok and novelty_is_replayable(e)]
    if not old:
        pytest.skip(f"no accepted submission in the sample is at or before "
                    f"epoch {SNAPSHOT_EPOCH}")
    wrong = []
    for seq, epoch in old:
        p = params[epoch]
        problems = check_sequence(seq, p["allowed_v"], p["allowed_cdr3_length"],
                                  archive, exclude_by_epoch(epoch), warnings=[])
        duplicate = [x for x in problems if "archived sequence" in x]
        if duplicate:
            wrong.append((epoch, seq[:30], duplicate))
    assert not wrong, (f"{len(wrong)}/{len(old)} accepted submissions called "
                       f"duplicates on a datable epoch: {wrong[:3]}")


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
        if not novelty_is_replayable(epoch):
            problems = without_novelty(problems)
        if problems:
            wrong.append((epoch, problems[0]))
    assert not wrong, f"{len(wrong)}/{len(accepted)} false refusals, e.g. {wrong[:3]}"
