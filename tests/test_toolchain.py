"""The install is part of the parity surface, so it gets tests of its own.

Every one of these guards a failure that actually happened and was silent:

  - DSSP absent or 4.x: TNP prints one line and reports a PSH taken on another
    code path (89.50 against 91.21 on the same structure), so eleven
    submissions were refused on a gate we believed we were screening for;
  - the filter's process_pdb patch missing: DSSP then fails on every
    NanoBodyBuilder2 model, which is how the above came about;
  - ImmuneBuilder's set/dict typo: 3 sequences in 37 come back with no model
    and therefore no PSH, which reads like "nothing to refuse".

These run anywhere -- they build fake install trees rather than needing the
real one -- except the two marked `installed`, which check the machine they run
on and skip when the tools are absent.

    python -m pytest tests/test_toolchain.py -q
"""
import pytest

from nova_lab import toolchain
from scripts.tnp_developability import audit, parse

DSSP_FAILURE = "ERROR: Failed to process output PDB. DSSP failed to produce an output"

TNP_LOG = """
TNP is running on sequence c0

Summary Statistics for c0:

    Total IMGT CDR Length: 34 (GREEN flag)
    IMGT CDR3 Length: 18 (GREEN flag)
    CDR Compactness: 1.395982 (GREEN flag)
    Patch CDR Surface Hydrophobicity score is: 91.206400 (GREEN flag)
    Patch CDR Positive Charge score is: 0.251000 (GREEN flag)
    Patch CDR Negative Charge score is: 0.370100 (GREEN flag)
"""


def fake_install(root, *, patched=True, immunebuilder_fixed=True, dssp=True):
    """A tree shaped like the venv TNP lives in, with the defects selectable."""
    bin_dir = root / "bin"
    site = root / "lib" / "python3.12" / "site-packages"
    (site / "scripts").mkdir(parents=True)
    (site / "ImmuneBuilder").mkdir(parents=True)
    bin_dir.mkdir(parents=True)
    for name in ("TNP", "NanoBodyBuilder2"):
        (bin_dir / name).write_text("#!/bin/sh\n")
    (site / "scripts" / "process_pdb.py").write_text(
        "def fix_pdb_for_dssp(f):\n    return f\n" if patched
        else "from Bio import PDB\n")
    threads = ("{'Threads': str(n_threads)}" if immunebuilder_fixed
               else "{'Threads', str(n_threads)}")
    (site / "ImmuneBuilder" / "refine.py").write_text(
        f"simulation = app.Simulation(t, s, i, p, {threads})\n")
    if dssp:
        (bin_dir / "mkdssp").write_text("#!/bin/sh\n")
    return bin_dir / "TNP"


def test_annotated_follows_the_file_dssp_writes(tmp_path):
    """`<name>_annotations.csv` is the per-sequence record that DSSP ran.

    The stdout error line cannot be attributed to a sequence once TNP runs
    several in parallel, so this file is what the guard reads.
    """
    raw = tmp_path / "Raw_Model_Outputs" / "c0"
    raw.mkdir(parents=True)
    assert not toolchain.annotated(tmp_path, "c0")
    (raw / "c0_annotations.csv").write_text("chain,num,ins\n")
    assert toolchain.annotated(tmp_path, "c0")


def test_check_catches_an_unpatched_tnp(tmp_path, monkeypatch):
    monkeypatch.setattr(toolchain, "TNP", fake_install(tmp_path, patched=False))
    rows = dict((what, (ok, detail)) for what, ok, detail in toolchain.check())
    ok, detail = rows["TNP process_pdb patch"]
    assert not ok and "UNPATCHED" in detail


def test_check_catches_the_immunebuilder_typo(tmp_path, monkeypatch):
    monkeypatch.setattr(toolchain, "TNP",
                        fake_install(tmp_path, immunebuilder_fixed=False))
    rows = dict((what, (ok, detail)) for what, ok, detail in toolchain.check())
    ok, detail = rows["ImmuneBuilder refine.py"]
    assert not ok and "SET" in detail


def test_check_passes_a_faithful_tree(tmp_path, monkeypatch):
    """Everything but mkdssp's version, which needs a real binary to report one."""
    monkeypatch.setattr(toolchain, "TNP", fake_install(tmp_path))
    rows = dict((what, ok) for what, ok, _ in toolchain.check())
    assert rows["TNP"] and rows["NanoBodyBuilder2"]
    assert rows["TNP process_pdb patch"] and rows["ImmuneBuilder refine.py"]


def test_tnp_env_puts_both_binaries_on_path(tmp_path, monkeypatch):
    """TNP shells out to a bare NanoBodyBuilder2, Biopython to a bare mkdssp."""
    tnp = fake_install(tmp_path)
    monkeypatch.setattr(toolchain, "TNP", tnp)
    monkeypatch.delenv("NOVA_DSSP", raising=False)
    path = toolchain.tnp_env()["PATH"].split(":")
    assert str(tnp.parent) in path


def test_parse_reads_the_metrics_and_the_flags():
    df = parse(TNP_LOG)
    assert list(df["name"]) == ["c0"]
    row = df.iloc[0]
    assert row["PSH"] == pytest.approx(91.2064)
    assert row["L3"] == 18 and row["reds"] == 0 and not row["rejected"]


def test_parse_marks_a_sequence_dssp_never_annotated(tmp_path):
    (tmp_path / "Raw_Model_Outputs" / "c0").mkdir(parents=True)
    assert not parse(TNP_LOG, tmp_path).iloc[0]["dssp"]


def test_audit_warns_when_dssp_failed(tmp_path, capsys):
    (tmp_path / "Raw_Model_Outputs" / "c0").mkdir(parents=True)
    audit(TNP_LOG + DSSP_FAILURE, tmp_path, wanted=1)
    said = capsys.readouterr().out
    assert "DSSP failed" in said and "NOT the validator's measurement" in said


def test_audit_warns_when_a_sequence_came_back_with_no_model(tmp_path, capsys):
    """No row is not a pass: that is how the gate looked clean while refusing us."""
    audit(TNP_LOG, tmp_path, wanted=3)
    said = capsys.readouterr().out
    assert "2 of 3" in said and "no PSH at all" in said


def test_audit_is_quiet_on_a_clean_run(tmp_path, capsys):
    raw = tmp_path / "Raw_Model_Outputs" / "c0"
    raw.mkdir(parents=True)
    (raw / "c0_annotations.csv").write_text("x\n")
    audit(TNP_LOG, tmp_path, wanted=1)
    assert capsys.readouterr().out == ""


@pytest.mark.installed
def test_this_machine_runs_the_validators_dssp():
    version = toolchain.dssp_version()
    if version is None:
        pytest.skip("no mkdssp on this machine; python -m scripts.doctor says so")
    assert version.startswith("3."), (
        f"mkdssp {version} refuses NanoBodyBuilder2 models (no PDB HEADER), so "
        f"TNP falls back to a PSH that is not the validator's. The filter builds "
        f"{toolchain.DSSP_VERSION} from cmbi/dssp.")


@pytest.mark.installed
def test_this_machine_has_the_whole_gate():
    rows = toolchain.check()
    if not any(ok for what, ok, _ in rows if what == "TNP"):
        pytest.skip("TNP is not installed here; the gate cannot be checked")
    broken = [f"{what}: {detail}" for what, ok, detail in rows if not ok]
    assert not broken, "python -m scripts.doctor — " + "; ".join(broken)
