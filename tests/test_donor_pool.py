"""The donor pool's freshness, which turns out to decide what the generator can reach.

The epoch fixes the CDR3 length, a donor of the wrong length is no use, and the
archive's good binders are overwhelmingly at length 13. So the pool's ceiling at
the epoch's own length is the ceiling on everything downstream -- and it moves
fast, because it is the leaders' winning sequences that raise it as they get
archived.

Measured across one day (FINDINGS 47):

    lengths with a donor at the field's level (design_iiptm >= 0.72)
      2026-11-13:  2 of 19   (11 % of epochs)
      2026-11-14:  8 of 19   (42 % of epochs)

The 0-for-8 replay of FINDINGS 46 ran on the first pool. At three of its eight
epochs the ceiling has since risen by 0.09 to 0.19.

`scored_but_unarchived` exists because the HuggingFace archive lags the
score-sharing API by about a day: 196 of 230 scored submissions were missing
from it when this was found. A day is ~20 epochs of the freshest and best
donors, so the merge is worth having even though the archive catches up.

    python -m pytest tests/test_donor_pool.py -q
"""
import pandas as pd
import pytest

from scripts.build_nanobody_dataset import METRICS, scored_but_unarchived

SEQ_A = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDRGLGDYWGQGTLVTVSS"
SEQ_B = "QVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDRGLGDYWGQGTLVTVSS"


class FakeArchive:
    """Only `regions_for` is used, and only to number the merged rows."""

    def regions_for(self, sequence):
        return {"CDR3": "AKDRGLGDY"}


def results_csv(tmp_path, target, sequences, iiptm=None):
    rows = []
    for i, s in enumerate(sequences):
        row = {"sequence": s, "hash": f"h{i}", "epoch": 26300 + i, "valid": True}
        row.update({m: 1.0 for m in METRICS})
        row["design_iiptm"] = (iiptm or [0.75] * len(sequences))[i]
        rows.append(row)
    directory = tmp_path / "submissions"
    directory.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(directory / f"{target}_nanobody_results.csv", index=False)


@pytest.fixture
def data_dir(tmp_path, monkeypatch):
    import scripts.build_nanobody_dataset as mod
    monkeypatch.setattr(mod, "DATA_DIR", tmp_path)
    return tmp_path


def test_it_merges_what_the_archive_has_not_got(data_dir):
    results_csv(data_dir, "Q9NZQ7", [SEQ_A, SEQ_B], iiptm=[0.75, 0.79])
    archived = pd.DataFrame({"sequence": [SEQ_A]})
    out = scored_but_unarchived("Q9NZQ7", archived, FakeArchive())
    assert list(out["sequence"]) == [SEQ_B]
    assert out["design_iiptm"].iloc[0] == 0.79
    assert out["epoch"].iloc[0] == 26301
    # The columns the dataset builder goes on to use must all be there.
    for column in ("sequence", "sequence_hash", "epoch", "regions", *METRICS):
        assert column in out.columns


def test_nothing_to_merge_when_the_archive_is_current(data_dir, capsys):
    results_csv(data_dir, "Q9NZQ7", [SEQ_A])
    out = scored_but_unarchived("Q9NZQ7", pd.DataFrame({"sequence": [SEQ_A]}), FakeArchive())
    assert out.empty
    assert "already" in capsys.readouterr().out


def test_a_missing_results_file_says_what_to_run(data_dir, capsys):
    """Silence here would leave the pool on the lagging archive without a word."""
    out = scored_but_unarchived("Q9NZQ7", pd.DataFrame({"sequence": []}), FakeArchive())
    assert out.empty
    said = capsys.readouterr().out
    assert "fetch_submissions" in said and "analyze_nanobodies" in said


def test_rows_without_metrics_are_dropped(data_dir):
    """A submission the API has not scored is not a donor, it is a sequence."""
    results_csv(data_dir, "Q9NZQ7", [SEQ_A, SEQ_B])
    path = data_dir / "submissions" / "Q9NZQ7_nanobody_results.csv"
    frame = pd.read_csv(path)
    frame.loc[frame["sequence"] == SEQ_B, "design_iiptm"] = None
    frame.to_csv(path, index=False)
    out = scored_but_unarchived("Q9NZQ7", pd.DataFrame({"sequence": []}), FakeArchive())
    assert list(out["sequence"]) == [SEQ_A]


def test_a_results_file_missing_a_metric_is_refused(data_dir, capsys):
    results_csv(data_dir, "Q9NZQ7", [SEQ_A])
    path = data_dir / "submissions" / "Q9NZQ7_nanobody_results.csv"
    frame = pd.read_csv(path).drop(columns=["interaction_pae"])
    frame.to_csv(path, index=False)
    out = scored_but_unarchived("Q9NZQ7", pd.DataFrame({"sequence": []}), FakeArchive())
    assert out.empty
    assert "lacks" in capsys.readouterr().out
