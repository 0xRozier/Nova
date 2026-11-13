"""The one stage of the miner that has never run, tested without a GPU.

`boltzgen_rank` is the selector that takes a single entry from 14 % of epochs
won to 48 %, and it is the only stage the production campaign never exercised:
there was no card. What it does itself is small -- write a batch, shell out to
scripts.boltzgen_score, read scores.csv back, rank -- and all of that is
testable with the scorer stubbed out. BoltzGen itself is exercised on Kaggle,
where every scores_*.csv in results/ came from.

The cases below are the ways it can be wrong while looking right, which is the
failure mode this project keeps finding:

  - a leftover scores.csv from an earlier attempt satisfies exists(), and
    BoltzGen reuses an output directory rather than recomputing (the README
    insists on a fresh RUN for exactly this reason), so the miner would rank
    another batch's numbers and submit a sequence chosen from them;
  - a run that scores only some of the batch leaves the rest with no metrics,
    and a sequence with no score is not a sequence that scored badly.

    python -m pytest tests/test_boltzgen_rank.py -q
"""
import subprocess
import time
from pathlib import Path

import pandas as pd
import pytest

import scripts.mine as mine

SEQUENCES = ["EVQLVESGGGAAA", "EVQLVESGGGBBB", "EVQLVESGGGCCC"]
METRICS = {
    "design_iiptm": [0.70, 0.60, 0.50],
    "design_ptm": [0.80, 0.78, 0.76],
    "design_to_target_iptm": [0.65, 0.60, 0.55],
    "min_design_to_target_pae": [2.0, 2.5, 3.0],
    "interaction_pae": [9.0, 10.0, 11.0],
    "plip_hbonds_refolded": [10, 8, 6],
    "plip_saltbridge_refolded": [5, 4, 3],
    "delta_sasa_refolded": [900.0, 850.0, 800.0],
    "liability_score": [200, 210, 220],
    "liability_num_violations": [20, 22, 24],
}


@pytest.fixture
def scores_path(tmp_path, monkeypatch):
    """Redirect DATA_DIR so the stub writes where boltzgen_rank reads."""
    monkeypatch.setattr(mine, "DATA_DIR", tmp_path)
    return tmp_path / "boltzgen" / "run1" / "scores.csv"


def write_scores(path: Path, sequences: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    n = len(sequences)
    pd.DataFrame({"sequence": sequences,
                  **{k: v[:n] for k, v in METRICS.items()}}).to_csv(path, index=False)


def stub(path: Path, sequences: list[str], returncode: int = 0):
    """Stand in for the scorer: write its output, report its exit code."""
    def run(cmd, capture_output=False, text=False, **kwargs):
        if returncode == 0:
            write_scores(path, sequences)
        return subprocess.CompletedProcess(cmd, returncode, stdout="", stderr="boom")
    return run


def test_it_ranks_what_the_scorer_returned(tmp_path, scores_path, monkeypatch):
    monkeypatch.setattr(subprocess, "run", stub(scores_path, SEQUENCES))
    out = mine.boltzgen_rank(SEQUENCES, "run1", tmp_path, tmp_path / "work")
    assert out is not None
    assert list(out["sequence"])[0] == SEQUENCES[0]  # best on every metric
    assert out["rank_sum"].is_monotonic_increasing


def test_it_writes_the_batch_the_scorer_reads(tmp_path, scores_path, monkeypatch):
    seen = {}

    def run(cmd, **kwargs):
        seen["cmd"] = cmd
        write_scores(scores_path, SEQUENCES)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", run)
    mine.boltzgen_rank(SEQUENCES, "run1", tmp_path, tmp_path / "work")
    batch = Path(seen["cmd"][seen["cmd"].index("--file") + 1])
    assert list(pd.read_csv(batch)["sequence"]) == SEQUENCES
    assert seen["cmd"][seen["cmd"].index("--run") + 1] == "run1"


def test_a_failed_run_ranks_nothing(tmp_path, scores_path, monkeypatch):
    monkeypatch.setattr(subprocess, "run", stub(scores_path, SEQUENCES, returncode=1))
    assert mine.boltzgen_rank(SEQUENCES, "run1", tmp_path, tmp_path / "work") is None


def test_it_refuses_a_scores_file_older_than_the_batch(tmp_path, scores_path, monkeypatch):
    """The stale-output trap: BoltzGen reuses a directory instead of recomputing.

    Without this the miner ranks an earlier batch's numbers and submits a
    sequence chosen from them -- the kind of wrong answer that leaves no trace.
    """
    write_scores(scores_path, SEQUENCES)
    old = time.time() - 3600
    import os
    os.utime(scores_path, (old, old))

    def run(cmd, **kwargs):  # succeeds without touching the file
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", run)
    assert mine.boltzgen_rank(SEQUENCES, "run1", tmp_path, tmp_path / "work") is None


def test_it_keeps_only_this_batch_and_says_what_is_missing(tmp_path, scores_path,
                                                           monkeypatch, capsys):
    """A partial run ranks the rest, and an unrelated sequence is not ranked."""
    returned = [SEQUENCES[0], "EVQLVESGGGZZZ"]  # one ours, one from elsewhere
    monkeypatch.setattr(subprocess, "run", stub(scores_path, returned))
    out = mine.boltzgen_rank(SEQUENCES, "run1", tmp_path, tmp_path / "work")
    assert out is not None
    assert list(out["sequence"]) == [SEQUENCES[0]]
    said = capsys.readouterr().out
    assert "no metrics" in said


def test_a_run_that_scored_none_of_ours_ranks_nothing(tmp_path, scores_path, monkeypatch):
    monkeypatch.setattr(subprocess, "run", stub(scores_path, ["EVQLVESGGGZZZ"]))
    assert mine.boltzgen_rank(SEQUENCES, "run1", tmp_path, tmp_path / "work") is None


def test_scores_without_a_sequence_column_rank_nothing(tmp_path, scores_path, monkeypatch):
    def run(cmd, **kwargs):
        scores_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"design_iiptm": [0.7]}).to_csv(scores_path, index=False)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", run)
    assert mine.boltzgen_rank(SEQUENCES, "run1", tmp_path, tmp_path / "work") is None
