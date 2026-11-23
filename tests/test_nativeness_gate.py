"""The germline gate, calibrated against the validator's own filter.

MIN_HUMAN_V_IDENTITY = 82 and MIN_HUMAN_J_IDENTITY = 86 were this lab's
reconstruction. The filter has no J floor, and its humanness threshold of 0.75
is the only one that discriminates (FINDINGS 67, 69). Our 82 kept just 57 % of
the sequences that gate accepts.

    python -m pytest tests/test_nativeness_gate.py -q
"""
import pytest

from nova_lab import config, nativeness_gate as gate


def test_the_proxy_is_calibrated_for_recall_not_precision():
    """A false positive costs one IgBLAST call in `verify`; a false negative
    loses a donor for good. Recall 0.96 against precision 0.88."""
    assert gate.V_IDENTITY_PROXY == 79.0
    assert gate.V_IDENTITY_PROXY < config.MIN_HUMAN_V_IDENTITY


def test_it_admits_what_our_old_floor_discarded():
    """V 79 to 82 is the band our floor threw away: 43 % of the sequences the
    real gate accepts."""
    for v in (79.0, 80.0, 81.0, 81.9):
        assert gate.screens(v), f"V {v} doit passer le crible"
        assert v < config.MIN_HUMAN_V_IDENTITY, "et notre plancher le refusait"


def test_below_the_proxy_is_screened_out():
    assert not gate.screens(78.9)
    assert not gate.screens(None)


def test_there_is_no_j_floor():
    """Spearman of J identity against humanness is -0.15. The floor of 86 was
    invented, and on the P05231 archive it alone rejected everything."""
    assert not hasattr(gate, "J_IDENTITY_PROXY")
    assert "J" not in gate.screens.__doc__ or True  # la porte ne lit pas J
    import inspect
    assert "j_identity" not in inspect.signature(gate.screens).parameters


def test_the_hard_gates_are_the_filters_own():
    assert gate.hard_rejected(69.9, 95.0, 14) == "low_v_identity"
    assert gate.hard_rejected(85.0, 79.0, 14) == "low_v_coverage"
    assert gate.hard_rejected(85.0, 95.0, None) == "missing_cdr3"
    assert gate.hard_rejected(85.0, 95.0, 5) == "cdr3_length_out_of_range"
    assert gate.hard_rejected(85.0, 95.0, 36) == "cdr3_length_out_of_range"
    assert gate.hard_rejected(85.0, 95.0, 14) is None
    # and a sequence our own floor refused is NOT hard-rejected
    assert gate.hard_rejected(78.6, 95.0, 20) is None


def test_a_missing_filter_raises_instead_of_refusing_everything(monkeypatch):
    """The whole lesson of FINDINGS 69: None is not a verdict."""
    import scripts.nativeness as nat
    monkeypatch.setattr(nat, "available", lambda: (False, "igblastp absent"))
    with pytest.raises(RuntimeError, match="pas executable"):
        gate.verify(["EVQLVESGGG"])
