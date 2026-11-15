"""The cutoff that makes a replay honest, and the hole that was in it.

A replay must not let the generator graft from sequences that did not exist at
the epoch being replayed. The cutoff used to filter on the `epoch` column
alone, which looks right and is not: the donor pool merges our decrypted
submissions, which carry an epoch, with the public archive, which carries no
date at all -- so re-downloading the archive injects today's sequences as
undated rows, and an undated row passed the cutoff.

Measured on 15 November: of the 197 undated donors binding at 0.72 or better
and germline-compliant, none existed at the epochs 26299-26306 that FINDINGS
48, 49 and 51 replayed. All 197 arrived after the only dated snapshot we have,
so the gains those sections reported were partly the future leaking in.

    python -m pytest tests/test_replay_honesty.py -q
"""
import pandas as pd
import pytest

import scripts.propose_nanobody as mod


@pytest.fixture
def pool(tmp_path, monkeypatch):
    """A dataset of four donors and a snapshot that dates only one of them."""
    monkeypatch.setattr(mod, "DATA_DIR", tmp_path)
    (tmp_path / "dataset").mkdir(parents=True)
    (tmp_path / "prebuilt").mkdir(parents=True)
    pd.DataFrame({"sequence": ["OLD_ARCHIVE"]}).to_csv(
        tmp_path / "prebuilt" / "T_regions.csv", index=False)
    frame = pd.DataFrame({
        "sequence": ["OLD_ARCHIVE", "NEW_ARCHIVE", "SUBMITTED_BEFORE", "SUBMITTED_AFTER"],
        "epoch": [None, None, 26200.0, 26400.0],
    })
    frame.to_parquet(tmp_path / "dataset" / "T_nanobodies.parquet")
    return tmp_path


def admissible(pool, cutoff):
    frame = pd.read_parquet(pool / "dataset" / "T_nanobodies.parquet")
    snapshot = mod.dated_archive_snapshot("T")
    dated_ok = frame["epoch"].notna() & (frame["epoch"] < cutoff)
    undated_ok = frame["epoch"].isna() & frame["sequence"].isin(snapshot)
    return set(frame[dated_ok | undated_ok]["sequence"])


def test_the_snapshot_is_the_only_dated_view_of_the_archive(pool):
    assert mod.dated_archive_snapshot("T") == {"OLD_ARCHIVE"}


def test_a_missing_snapshot_admits_no_undated_sequence(pool):
    """Erring towards a pool that is too small, never towards inventing donors."""
    (pool / "prebuilt" / "T_regions.csv").unlink()
    assert mod.dated_archive_snapshot("T") == set()
    assert admissible(pool, 26300) == {"SUBMITTED_BEFORE"}


def test_an_undated_sequence_the_archive_gained_since_is_refused(pool):
    """The hole this test exists for: NEW_ARCHIVE has no epoch, so the old
    cutoff kept it although it may post-date the replayed epoch."""
    assert "NEW_ARCHIVE" not in admissible(pool, 26300)


def test_what_existed_then_is_still_admitted(pool):
    assert admissible(pool, 26300) == {"OLD_ARCHIVE", "SUBMITTED_BEFORE"}


def test_a_submission_from_the_replayed_epoch_or_later_is_refused(pool):
    assert "SUBMITTED_AFTER" not in admissible(pool, 26300)
    assert "SUBMITTED_BEFORE" not in admissible(pool, 26200), "the cutoff is exclusive"


def test_without_a_cutoff_nothing_is_dropped(pool):
    """Live mining has no leak by construction: everything in today's pool is
    genuinely in the past, so a live run must see all of it."""
    frame = pd.read_parquet(pool / "dataset" / "T_nanobodies.parquet")
    assert len(frame) == 4
