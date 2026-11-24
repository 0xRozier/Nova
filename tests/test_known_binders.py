"""The two crystallographic IL-6 binders, as fixtures the filter must reproduce.

FINDINGS 73 measured both against the validator's own nativeness filter. They
are kept here because they anchor three things no synthetic sequence can:

  - that the filter is actually running (a missing pip dependency makes it
    answer hard_reject for everything, which reads as a refusal -- the shape
    of the FINDINGS 69 bug, twice now);
  - the numbers themselves, so a filter or database update that moves them
    shows up as a failed test rather than as a quietly different strategy;
  - that the gate the validators enforce excludes nature's camelid solution
    (VHH6) and admits the humanised one (61H7), which is the whole reason the
    strategy points at humanised frameworks.

The structures: 5FUC (VHH6, a true single-domain VHH, junctional epitope on
IL-6/IL-6R-alpha) and 4O9H (the heavy chain of camelid Fab 61H7, the structure
this epoch's config names).
"""
import pytest

from scripts import nativeness

VHH6 = ("DVQFVESGGGSVHAGGSLRLNCATSGYIYSTYCMGWFRQAPGKEREGVAHIYTNSGRTYYADSVKGRFTISQDNAKNTVY"
        "LQMNSLKPEDTAIYYCAARPSIRCASFSATEYKDWGQGTQVTVSS")
H61H7 = ("EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYRMYWVRQPPGKGLEWVSAISAGGGSTYYGDSVKGRFTISRDNAKNTVY"
         "LQMNSLKPEDTAVYYCANRAGWGMGDYWGQGTQVTVSS")

needs_filter = pytest.mark.skipif(not nativeness.available()[0],
                                  reason=f"filtre indisponible : {nativeness.available()[1]}")


@pytest.fixture(scope="module")
def scored():
    return nativeness.verdict(nativeness.score([VHH6, H61H7])).set_index("sequence")


def test_both_clear_layer_1_length_and_cysteines():
    for seq in (VHH6, H61H7):
        assert 90 <= len(seq) <= 150
        assert seq.count("C") >= 1


@needs_filter
def test_neither_known_binder_is_hard_rejected(scored):
    """A hard reject here would mean the filter is broken, not that they fail."""
    assert not scored.loc[VHH6, "hard_reject"]
    assert not scored.loc[H61H7, "hard_reject"]


@needs_filter
def test_the_true_vhh_is_camelid_but_not_human_enough(scored):
    assert scored.loc[VHH6, "nativeness"] == pytest.approx(0.852, abs=0.02)
    assert scored.loc[VHH6, "humanness"] == pytest.approx(0.414, abs=0.02)
    assert not scored.loc[VHH6, "passes_humanness"]


@needs_filter
def test_the_humanised_camelid_vh_clears_the_enforced_gate(scored):
    assert scored.loc[H61H7, "humanness"] == pytest.approx(0.992, abs=0.02)
    assert scored.loc[H61H7, "passes_humanness"]
    lo, hi = nativeness.ACCEPTED_NATIVENESS_RANGE
    assert lo <= scored.loc[H61H7, "nativeness"] <= hi


@needs_filter
def test_the_gate_separates_them_the_way_the_strategy_assumes(scored):
    """Opposite sides of the only threshold FINDINGS 69 found enforced."""
    assert scored.loc[H61H7, "passes_humanness"]
    assert not scored.loc[VHH6, "passes_humanness"]


def test_available_names_the_pip_dependencies_it_needs():
    """The 69-shaped bug: a missing piece must be reported, never scored as a refusal."""
    assert {"abnumber", "anarci"} <= set(nativeness.PYTHON_DEPS)
