"""The regions cache, and the failure it used to engrave for good.

A sequence whose IMGT numbering failed was written to the cache with an empty
CDR3, read back as None, and never looked at again -- indistinguishable from a
sequence that genuinely has no CDR3. One transient ANARCI outage was therefore
permanent. It was: on 15 November 6740 of the archive's 17244 rows were marked
unnumberable in one contiguous block, and renumbering them succeeded 6740 times
out of 6740. The generator lost 446 compliant donors binding at 0.72 or better.

    python -m pytest tests/test_regions_cache.py -q
"""
import pandas as pd
import pytest

import nova_lab.nanobodies as mod

PROBE = ("EVQLVESGGGLVQPGGSLRLSCAASGRTFSYNPMGWFRQAPGKGRELVAAISRTGGSTYYPDSVEGRFT"
         "ISRDNAKRMVYLQMNSLRAEDTAVYYCAAAGVRAEDGRVRTLPSEYTFWGQGTQVTVSS")


@pytest.fixture
def archive_dir(tmp_path, monkeypatch):
    monkeypatch.setattr(mod, "DATA_DIR", tmp_path)
    directory = tmp_path / "archive"
    directory.mkdir(parents=True)
    pd.DataFrame({"sequence": [PROBE], "sequence_hash": [mod.seq_hash(PROBE)]}
                 ).to_csv(directory / "T_nanobodies.csv", index=False)
    return directory


def poison(directory):
    """The cache as a failed run leaves it: every region column empty."""
    pd.DataFrame([{"sequence": PROBE, **{k: "" for k in mod.NanobodyArchive.REGION_KEYS}}]
                 ).to_csv(directory / "T_regions.csv", index=False)


def test_a_cached_failure_stays_a_failure_by_default(archive_dir, monkeypatch):
    """The miner must not spend minutes renumbering inside an epoch."""
    poison(archive_dir)
    monkeypatch.setattr(mod, "regions_batch",
                        lambda seqs: pytest.fail("renumbered without being asked"))
    archive = mod.NanobodyArchive(target="T")
    assert archive.regions == [None]


def test_retry_failed_renumbers_what_was_engraved(archive_dir):
    poison(archive_dir)
    archive = mod.NanobodyArchive(target="T", retry_failed=True)
    assert archive.regions[0] is not None, "the probe is numberable"
    assert archive.regions[0]["CDR3"]


def test_the_repair_is_written_back(archive_dir):
    poison(archive_dir)
    mod.NanobodyArchive(target="T", retry_failed=True)
    cached = pd.read_csv(archive_dir / "T_regions.csv", keep_default_na=False)
    assert cached.loc[0, "CDR3"], "the next run must not have to renumber again"
    # and the repaired cache is enough on its own, without retry_failed
    assert mod.NanobodyArchive(target="T").regions[0] is not None


def test_a_successful_cache_is_not_renumbered(archive_dir, monkeypatch):
    mod.NanobodyArchive(target="T", retry_failed=True)      # populate
    monkeypatch.setattr(mod, "regions_batch",
                        lambda seqs: pytest.fail("renumbered an entry it had"))
    assert mod.NanobodyArchive(target="T", retry_failed=True).regions[0] is not None


def test_a_sequence_that_really_cannot_be_numbered_is_not_an_error(archive_dir):
    """It must come back as None, not raise, or one bad row kills the pool."""
    pd.DataFrame({"sequence": ["NOTAPROTEIN"], "sequence_hash": ["x"]}
                 ).to_csv(archive_dir / "T_nanobodies.csv", index=False)
    assert mod.NanobodyArchive(target="T", retry_failed=True).regions == [None]
