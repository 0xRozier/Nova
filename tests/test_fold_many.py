"""Folding many frameworks: batched and resumable, which is what makes §59 askable.

`fold` starts a fresh interpreter per sequence, so NanoBodyBuilder2's weights load
once per framework. For the six scaffolds of FINDINGS 56 that was tolerable; the
question FINDINGS 59 leaves open -- a few designs on many scaffolds rather than
many designs on six -- needs forty, where loading dominates. And a run longer
than the window it runs in has to resume rather than start over.

    python -m pytest tests/test_fold_many.py -q
"""
import subprocess

import pytest

import scripts.design_nanobody as mod


@pytest.fixture
def folder(tmp_path, monkeypatch):
    """A NanoBodyBuilder2 stand-in that records its calls and writes the files."""
    calls = []

    class Done:
        returncode = 0
        stdout = stderr = ""

    def fake(cmd, input=None, **kwargs):
        calls.append(input)
        for line in (input or "").splitlines():
            seq, out = line.split("\t")
            from pathlib import Path
            Path(out).write_text("ATOM      1  CA  ALA H   1       0.0   0.0   0.0\n")
        return Done()

    monkeypatch.setattr(subprocess, "run", fake)
    monkeypatch.setattr(mod, "FOLDER_PYTHON", tmp_path / "python")
    (tmp_path / "python").write_text("")
    monkeypatch.setattr(mod, "renumber_sequentially", lambda p, chain="B": 0)
    return calls


def test_everything_folds_in_one_process(folder, tmp_path):
    jobs = [("SEQ1", tmp_path / "s0.pdb"), ("SEQ2", tmp_path / "s1.pdb"),
            ("SEQ3", tmp_path / "s2.pdb")]
    mod.fold_many(jobs, log=lambda *a: None)
    assert len(folder) == 1, "the weights must load once, not once per framework"
    assert folder[0].count("\n") == 2 and "SEQ3" in folder[0]
    assert all(out.exists() for _, out in jobs)


def test_what_is_already_on_disk_is_not_refolded(folder, tmp_path):
    done = tmp_path / "s0.pdb"
    done.write_text("already")
    jobs = [("SEQ1", done), ("SEQ2", tmp_path / "s1.pdb")]
    mod.fold_many(jobs, log=lambda *a: None)
    assert "SEQ1" not in folder[0] and "SEQ2" in folder[0]
    assert done.read_text() == "already", "an existing structure is left alone"


def test_nothing_to_do_starts_no_process(folder, tmp_path):
    jobs = [("SEQ1", tmp_path / "s0.pdb")]
    jobs[0][1].write_text("already")
    mod.fold_many(jobs, log=lambda *a: None)
    assert folder == [], "a resumed run with nothing left must not load the weights"


def test_a_structure_that_was_not_written_is_an_error(tmp_path, monkeypatch):
    """Silence here would hand BoltzGen a spec pointing at a missing file."""
    class Done:
        returncode = 0
        stdout = stderr = ""
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: Done())
    monkeypatch.setattr(mod, "FOLDER_PYTHON", tmp_path / "python")
    (tmp_path / "python").write_text("")
    with pytest.raises(RuntimeError, match="wrote nothing"):
        mod.fold_many([("SEQ1", tmp_path / "s0.pdb")], log=lambda *a: None)


def test_a_failing_folder_is_an_error(tmp_path, monkeypatch):
    class Failed:
        returncode = 1
        stdout = ""
        stderr = "out of memory"
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: Failed())
    monkeypatch.setattr(mod, "FOLDER_PYTHON", tmp_path / "python")
    (tmp_path / "python").write_text("")
    with pytest.raises(RuntimeError, match="out of memory"):
        mod.fold_many([("SEQ1", tmp_path / "s0.pdb")], log=lambda *a: None)


def test_the_payload_drops_hydrogens_and_keeps_every_residue(tmp_path):
    """Kaggle refuses a notebook over 1 MB, and forty all-atom scaffolds made
    1284 kB. Hydrogens are 48 % of the atoms and are not what a structure model
    reads -- BoltzGen's own reference spec points at a PDB crystal structure,
    which has none. What must NOT change is the residues the spec names.
    """
    from scripts.build_design_notebook import is_hydrogen, without_hydrogens
    pdb = (
        "REMARK  MODELLED\n"
        "ATOM      1  N   GLN B   1     -10.5  -5.6 -18.1  1.00  1.83           N  \n"
        "ATOM      2  CA  GLN B   1     -10.0  -5.0 -17.0  1.00  1.83           C  \n"
        "ATOM      3  H   GLN B   1      -9.0  -4.0 -16.0  1.00  1.83           H  \n"
        "ATOM      4  HA2 GLN B   1      -8.0  -3.0 -15.0  1.00  1.83           H  \n"
        "END\n"
    )
    out = without_hydrogens(pdb)
    assert "  N   GLN" in out and "  CA  GLN" in out, "heavy atoms stay"
    assert " H   GLN" not in out and "HA2" not in out
    assert out.startswith("REMARK") and out.endswith("END\n"), "non-atom lines stay"
    # An atom whose NAME starts with H but whose element column says otherwise:
    # HG is mercury here, and the element column decides. PDB columns are fixed
    # width, so the line has to be built at the right offsets -- element is 77-78.
    line = ("ATOM      5  HG  CYS B   2".ljust(30)
            + f"{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{1.0:6.2f}".ljust(46)
            + "HG")
    assert line[76:78] == "HG"
    assert is_hydrogen(line) is False, "the element column decides, not the atom name"
