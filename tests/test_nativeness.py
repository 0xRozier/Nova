"""Calling the validator's own nativeness filter, rather than reconstructing it.

Three validator rules have been reconstructed in this repo and two were wrong:
the novelty gate copied an identity rule from the filter's search engine that
the diversity validator never calls (FINDINGS 54), and MIN_HUMAN_V_IDENTITY /
MIN_HUMAN_J_IDENTITY are ours, not the filter's (FINDINGS 67). So this module
calls the filter's code and only reports.

The trap it exists to close: `NativenessFilter.analyze` catches every exception
and returns None, so a missing IgBLAST database reads exactly like a refused
sequence. Four known-accepted submissions were first scored as refused for that
reason alone.

    python -m pytest tests/test_nativeness.py -q
"""
import pandas as pd
import pytest

from scripts import nativeness


def test_the_nativeness_threshold_is_reported_and_not_enforced():
    """No accepted submission we have reaches 0.80 (median 0.507, max 0.596 over
    60). Enforcing it would refuse everything the validators took."""
    frame = nativeness.verdict(pd.DataFrame([
        {"sequence": "A", "nativeness": 0.507, "humanness": 0.84,
         "hard_reject": False, "reject_reason": None}]))
    assert frame.loc[0, "passes_hard"]
    assert frame.loc[0, "passes_humanness"]
    assert "passes_nativeness" not in frame.columns, "ce seuil ne doit pas etre une porte"
    assert frame.loc[0, "at_least_accepted_nativeness"]


def test_the_humanness_gate_is_the_one_that_discriminates():
    """40 of 40 accepted submissions pass 0.75; only 21 % of the P05231
    archive's binders do."""
    frame = nativeness.verdict(pd.DataFrame([
        {"sequence": "accepted", "nativeness": 0.51, "humanness": 0.84,
         "hard_reject": False, "reject_reason": None},
        {"sequence": "archive", "nativeness": 0.65, "humanness": 0.717,
         "hard_reject": False, "reject_reason": None}]))
    assert list(frame.passes_humanness) == [True, False]
    # and the archive binder scores HIGHER on nativeness while failing humanness
    assert frame.loc[1, "nativeness"] > frame.loc[0, "nativeness"]


def test_a_hard_reject_is_kept_separate_from_the_thresholds():
    frame = nativeness.verdict(pd.DataFrame([
        {"sequence": "junk", "nativeness": None, "humanness": None,
         "hard_reject": True, "reject_reason": "low_v_identity"}]))
    assert not frame.loc[0, "passes_hard"]
    assert not frame.loc[0, "passes_humanness"], "un None ne doit pas passer"


def test_a_missing_piece_is_an_error_not_a_refusal(monkeypatch, tmp_path):
    """The whole point: analyze() returning None must never be read as a verdict."""
    monkeypatch.setattr(nativeness, "IGBLAST_DIR", tmp_path / "nowhere")
    usable, why = nativeness.available()
    assert usable is False
    assert "igblastp" in why
    with pytest.raises(SystemExit, match="n'est pas executable"):
        nativeness.score(["EVQLVESGGG"])


def test_the_human_internal_data_is_checked_too(monkeypatch, tmp_path):
    """igblastp dies on a missing internal_data/human, and the filter swallows it."""
    root = tmp_path / "igblast"
    (root / "bin").mkdir(parents=True)
    (root / "bin" / "igblastp").write_text("")
    (root / "database").mkdir()
    for db in ("camelid_V.phr", "human_V.phr"):
        (root / "database" / db).write_text("")
    (root / "internal_data" / "camelid").mkdir(parents=True)
    monkeypatch.setattr(nativeness, "IGBLAST_DIR", root)
    usable, why = nativeness.available()
    assert usable is False and "internal_data/human" in why
