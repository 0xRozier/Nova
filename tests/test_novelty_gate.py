"""The archive-novelty gate, on a synthetic archive.

This gate used to be covered only through `test_nanobody_parity.py`, by
replaying real epochs. That coverage is gone and cannot come back: the
validator compares a submission against the archive **as it stood then**, and
about 7000 of today's 16700 entries carry no date in any of our data, so an
epoch after 26025 cannot be rewound (see `novelty_is_replayable` there).

Replaying is also the wrong tool for the rule itself. A hand-built archive
pins the thresholds of `is_duplicate` exactly -- 95 % whole sequence, 85 % CDR3,
and the 80/90/90 combination -- and says which verdicts are rejections rather
than warnings, which is what decides whether the miner drops a candidate.

    python -m pytest tests/test_novelty_gate.py -q
"""
import pandas as pd
import pytest

from nova_lab.nanobodies import NanobodyArchive, is_duplicate, seq_hash

# caplacizumab, whose numbering is known to work
PROBE = ("EVQLVESGGGLVQPGGSLRLSCAASGRTFSYNPMGWFRQAPGKGRELVAAISRTGGSTYYPDSVEGRFT"
         "ISRDNAKRMVYLQMNSLRAEDTAVYYCAAAGVRAEDGRVRTLPSEYTFWGQGTQVTVSS")


def archive_of(sequences, tmp_path, monkeypatch):
    """A NanobodyArchive over exactly these sequences, numbered for real."""
    directory = tmp_path / "archive"
    directory.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"sequence": sequences,
                  "sequence_hash": [seq_hash(s) for s in sequences]}
                 ).to_csv(directory / "T_nanobodies.csv", index=False)
    import nova_lab.nanobodies as mod
    monkeypatch.setattr(mod, "DATA_DIR", tmp_path)
    return NanobodyArchive(target="T")


def test_the_thresholds_are_where_the_filter_puts_them():
    """Pinned against NOVA-nanobody-filter's search engine, rule by rule."""
    assert is_duplicate(0.95, None)[:1] == (True,)
    assert is_duplicate(0.95, None)[2] is True, "95% whole sequence is a rejection"
    # just under 95 %, with CDR data that triggers nothing: no hit. Without CDR
    # data the looser 90 % fallback below would catch it instead.
    quiet = {"CDR3": 0.50, "CDR1": 0.50, "CDR2": 0.50}
    assert is_duplicate(0.9499, quiet) == (False, "", False)

    hit, _, certain = is_duplicate(0.80, {"CDR3": 0.85})
    assert hit and certain is False, "same clonotype, but far enough to be a warning"
    hit, _, certain = is_duplicate(0.85, {"CDR3": 0.85})
    assert hit and certain is True, "close AND same clonotype is a rejection"

    assert is_duplicate(0.80, {"CDR3": 0.84})[0] is False
    assert is_duplicate(0.80, {"CDR3": 0.80, "CDR1": 0.90, "CDR2": 0.90})[0] is True
    assert is_duplicate(0.80, {"CDR3": 0.80, "CDR1": 0.89, "CDR2": 0.90})[0] is False

    # with no CDR data the rule falls back to a looser whole-sequence bar
    assert is_duplicate(0.90, None)[0] is True
    assert is_duplicate(0.89, None)[0] is False


def test_an_exact_resubmission_is_caught(tmp_path, monkeypatch):
    archive = archive_of([PROBE], tmp_path, monkeypatch)
    hit, reason, certain = archive.duplicate_of(PROBE)
    assert hit and certain and "exact duplicate" in reason


def test_a_point_mutant_is_caught_as_near_identical(tmp_path, monkeypatch):
    archive = archive_of([PROBE], tmp_path, monkeypatch)
    mutant = PROBE[:40] + ("A" if PROBE[40] != "A" else "G") + PROBE[41:]
    hit, reason, certain = archive.duplicate_of(mutant)
    assert hit and certain, f"one substitution should still be >= 95%: {reason}"


def test_an_unrelated_sequence_passes(tmp_path, monkeypatch):
    """A different nanobody must not be called a duplicate of this one."""
    other = ("QVQLQESGGGLVQAGGSLRLSCAASGFTFDDYAIGWFRQAPGKEREGVSCISSSDGSTYYADSVKGRFT"
             "ISSDNAKNTVYLQMNSLKPEDTAVYYCAADPSPYYRGSAYDYWGQGTQVTVSS")
    archive = archive_of([PROBE], tmp_path, monkeypatch)
    assert archive.duplicate_of(other)[0] is False


def test_exclude_replays_an_epoch_whose_submission_is_already_archived(tmp_path, monkeypatch):
    """The validator's exclude_ids: a submission is not a duplicate of itself.

    This is what makes an honest replay possible, and getting it wrong is what
    made the checker refuse sequences the validators had accepted.
    """
    archive = archive_of([PROBE], tmp_path, monkeypatch)
    assert archive.duplicate_of(PROBE)[0] is True
    assert archive.duplicate_of(PROBE, exclude={PROBE})[0] is False


def test_exclusion_is_by_sequence_not_by_neighbourhood(tmp_path, monkeypatch):
    """Excluding one sequence must not excuse a near-identical neighbour.

    Otherwise a replay would quietly pass candidates the validator refuses.
    """
    mutant = PROBE[:40] + ("A" if PROBE[40] != "A" else "G") + PROBE[41:]
    archive = archive_of([PROBE, mutant], tmp_path, monkeypatch)
    assert archive.duplicate_of(PROBE, exclude={mutant})[0] is True
