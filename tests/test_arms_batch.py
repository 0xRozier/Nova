"""Assembling one batch out of several arms, which is the only comparable way.

BoltzGen is stochastic enough that two draws of the same sequences differ by
about 14.5 places on the rank sum, so an arm scored in its own run cannot be
compared with another. Everything that is to be compared goes in one run, and
the ordering has to survive an interrupted one -- a Kaggle session that dies at
60 % must leave every arm of every epoch measured to some depth, not the first
arm measured completely.

    python -m pytest tests/test_arms_batch.py -q
"""
import pandas as pd
import pytest

import scripts.build_arms_batch as mod


def epoch_batch(directory, target, epoch, candidates, reals):
    rows = []
    for i, s in enumerate(candidates):
        rows.append({"sequence": s, "group": "candidate", "predicted_rank_sum": 100 + i,
                     "true_rank_sum": None, "repo": None, "cdr3": "AAA", "v_identity": 85.0})
    for i, s in enumerate(reals):
        rows.append({"sequence": s, "group": "real", "predicted_rank_sum": None,
                     "true_rank_sum": 40 + i, "repo": "them/repo", "cdr3": None,
                     "v_identity": None})
    directory.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(directory / f"{target}_epoch{epoch}.csv", index=False)


@pytest.fixture
def batches(tmp_path, monkeypatch):
    monkeypatch.setattr(mod, "BATCH_DIR", tmp_path / "batches")
    epoch_batch(tmp_path / "batches", "Q9NZQ7", 26300, ["C1", "C2", "C3"], ["R1", "R2"])
    epoch_batch(tmp_path / "batches", "Q9NZQ7", 26301, ["C4", "C5"], ["R3", "R4"])
    return tmp_path


def humanised_dir(tmp_path, per_epoch):
    directory = tmp_path / "hum"
    directory.mkdir(parents=True, exist_ok=True)
    for epoch, sequences in per_epoch.items():
        pd.DataFrame({"sequence": sequences,
                      "donor_design_iiptm": [0.8] * len(sequences)}
                     ).to_csv(directory / f"e{epoch}.csv", index=False)
    return directory


def test_it_takes_every_arm_of_every_epoch(batches, tmp_path):
    out = tmp_path / "batch.csv"
    import sys
    # --target explicite : la cible configuree a change le 20 novembre (Q9NZQ7
    # vers P05231) et ce test porte sur l'assemblage, pas sur la valeur par
    # defaut. Sans lui il cherchait des lots P05231_epoch*.csv que la fixture
    # n'ecrit pas.
    argv = ["x", "--depth", "2", "--target", "Q9NZQ7", "--out", str(out),
            "--humanised-dir", str(humanised_dir(tmp_path, {26300: ["H1", "H2"],
                                                            26301: ["H3"]}))]
    monkey = pytest.MonkeyPatch()
    monkey.setattr(sys, "argv", argv)
    try:
        assert mod.main() == 0
    finally:
        monkey.undo()
    frame = pd.read_csv(out)
    assert set(frame["arm"]) == {"grafted", "real", "humanised"}
    # depth 2 caps each arm of each epoch
    counts = frame.groupby(["epoch", "arm"]).size().to_dict()
    assert counts[(26300, "grafted")] == 2 and counts[(26300, "real")] == 2
    assert counts[(26301, "humanised")] == 1


def test_an_interrupted_run_still_covers_everything(batches, tmp_path):
    """The first rows must spread across epochs and arms, not finish one arm."""
    frame = mod.interleave(pd.DataFrame({
        "sequence": list("abcdefgh"),
        "epoch": [1, 1, 1, 1, 2, 2, 2, 2],
        "arm": ["grafted", "grafted", "real", "real"] * 2,
    }))
    first_four = frame.head(4)
    assert first_four["epoch"].nunique() == 2
    assert set(first_four["arm"]) == {"grafted", "real"}


def test_duplicate_sequences_are_scored_once(batches, tmp_path):
    """BoltzGen deduplicates anyway, and a duplicate wastes 1.4 min."""
    out = tmp_path / "batch.csv"
    import sys
    argv = ["x", "--depth", "4", "--target", "Q9NZQ7", "--out", str(out),
            "--humanised-dir", str(humanised_dir(tmp_path, {26300: ["C1", "H1"]}))]
    monkey = pytest.MonkeyPatch()
    monkey.setattr(sys, "argv", argv)
    try:
        mod.main()
    finally:
        monkey.undo()
    frame = pd.read_csv(out)
    assert frame["sequence"].is_unique


def test_nothing_to_assemble_is_an_error_not_an_empty_batch(tmp_path, monkeypatch):
    monkeypatch.setattr(mod, "BATCH_DIR", tmp_path / "empty")
    import sys
    monkeypatch.setattr(sys, "argv", ["x", "--out", str(tmp_path / "b.csv")])
    with pytest.raises(SystemExit):
        mod.main()
