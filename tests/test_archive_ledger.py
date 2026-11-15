"""Stamping the archive, so a replay stops grafting from the future.

The public archive is a flat list with no submission date. Re-downloading it
injects today's sequences, and before this a replay could not tell them from
six-week-old ones: of the 197 undated donors binding at 0.72 or better and
clearing the germline floors, none existed at the epochs FINDINGS 48, 49 and 51
replayed.

The stamp is an upper bound on the real submission epoch -- a sequence was
submitted at or before the moment we first saw it -- which is the safe
direction: a replay may miss a donor that existed, never invent one.

    python -m pytest tests/test_archive_ledger.py -q
"""
import pytest

from nova_lab import archive_ledger as led


@pytest.fixture
def data(tmp_path, monkeypatch):
    monkeypatch.setattr(led, "DATA_DIR", tmp_path)
    (tmp_path / "prebuilt").mkdir(parents=True)
    (tmp_path / "prebuilt" / "T_regions.csv").write_text("sequence\nOLD\n")
    return tmp_path


def test_without_a_ledger_it_falls_back_to_the_frozen_snapshot(data):
    assert led.bootstrap("T") == {"OLD": led.SNAPSHOT_EPOCH}
    assert led.seen_before("T", 26300) == {"OLD"}
    assert led.seen_before("T", led.SNAPSHOT_EPOCH) == set(), "the bound is exclusive"


def test_new_sequences_are_stamped_at_the_epoch_we_saw_them(data):
    assert led.record("T", ["OLD", "FRESH"], 26340) == 1, "only FRESH is new"
    assert led.load("T") == {"OLD": led.SNAPSHOT_EPOCH, "FRESH": 26340}


def test_a_stamp_is_never_moved_by_seeing_the_sequence_again(data):
    """Seeing it again says nothing new, and moving the stamp forward would
    hide a donor that a later replay is entitled to."""
    led.record("T", ["FRESH"], 26340)
    assert led.record("T", ["FRESH"], 26400) == 0
    assert led.load("T")["FRESH"] == 26340


def test_a_replay_sees_only_what_was_stamped_before_it(data):
    led.record("T", ["A", "B"], 26340)
    led.record("T", ["C"], 26350)
    assert led.seen_before("T", 26345) == {"OLD", "A", "B"}
    assert "C" not in led.seen_before("T", 26345)
    assert led.seen_before("T", 26360) == {"OLD", "A", "B", "C"}


def test_recording_nothing_new_leaves_the_file_alone(data):
    led.record("T", ["A"], 26340)
    before = led.path("T").read_text()
    assert led.record("T", ["A"], 26350) == 0
    assert led.path("T").read_text() == before


def test_the_write_is_atomic(data):
    """A killed rebuild must leave the previous ledger, not a truncated one."""
    led.record("T", ["A"], 26340)
    led.record("T", ["B"], 26350)
    assert set(led.load("T")) == {"OLD", "A", "B"}
    assert not list(data.glob("archive/*.tmp")), "no temporary file left behind"


def test_empty_sequences_are_not_stamped(data):
    assert led.record("T", ["", None], 26340) == 0
