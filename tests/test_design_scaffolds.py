"""Which frameworks a design is allowed to stand on.

Two numbers were conflated for weeks. Grafting mutates CDR1, CDR2 and the CDR3
body, each mutation costing about a point of V identity, so its donors need
2.0 points of headroom above the floor. A design holds the framework and both
other CDRs fixed: measured on the 120 designs of FINDINGS 60, it costs 0.00
points of V identity at the median and 1.02 at worst.

Importing the graft's 2.0 therefore barred exactly the band the field occupies
-- accepted submissions sit at V 83.7 median, 56 % within two points of the
floor, and intra-epoch rho(V identity, design_iiptm) = -0.25, so less germline
binds better. Our 40 scaffolds sat at 84.7-88.8 while 314 compliant frameworks
waited at V 82-84.

    python -m pytest tests/test_design_scaffolds.py -q
"""
import pandas as pd
import pytest

import scripts.design_nanobody as mod


@pytest.fixture
def pool(tmp_path, monkeypatch):
    monkeypatch.setattr(mod, "DATA_DIR", tmp_path)
    # ANARCI cannot number a placeholder, and what is under test is the filter
    # and the ordering, not the numbering. Each sequence gets its own framework
    # signature so none is discarded as a near-duplicate.
    monkeypatch.setattr(mod, "is_vhh", lambda seq: (True, []))
    monkeypatch.setattr(mod, "regions", lambda seq: {
        "FR1": seq[:3], "FR2": "F2", "FR3": "F3", "FR4": "F4",
        "CDR1": "C1", "CDR2": "C2", "CDR3": "ARDPNPLYQFRFDV"})
    (tmp_path / "dataset").mkdir(parents=True)
    frame = pd.DataFrame({
        "sequence":   [f"FW{i}" + "A" * 100 for i in range(6)],
        "v_identity": [81.5, 82.5, 83.2, 84.9, 86.0, 88.5],
        "j_identity": [92.9] * 6,
        "v_genes":    ["IGHV3-74"] * 6,
        "cdr3":       ["ARDPNPLYQFRFDV"] * 6,
        "design_iiptm": [0.5] * 6,
    })
    for col in ("design_ptm", "design_to_target_iptm", "min_design_to_target_pae",
                "interaction_pae", "plip_hbonds_refolded", "plip_saltbridge_refolded",
                "delta_sasa_refolded", "liability_score", "liability_num_violations"):
        frame[col] = 1.0
    frame.to_parquet(tmp_path / "dataset" / "T_nanobodies.parquet")
    return frame


def v_values(frame):
    return sorted(frame["v_identity"])


def test_the_design_margin_is_not_the_graft_margin(pool):
    """1.1 covers the worst measured cost of a graft-free design; 2.0 does not
    describe anything a design spends."""
    assert mod.DESIGN_V_MARGIN == 1.1
    from scripts.propose_nanobody import GERMLINE_MARGIN
    assert GERMLINE_MARGIN == 2.0, "grafting keeps its own, larger headroom"


def test_it_admits_the_band_the_field_occupies(pool):
    got = mod.framework_donors("T", "IGHV3-74", 6)
    assert 83.2 in v_values(got), "V 83.2 is inside the field's own range"
    assert 81.5 not in v_values(got), "below the floor plus margin stays out"


def test_the_old_margin_barred_that_band(pool):
    got = mod.framework_donors("T", "IGHV3-74", 6, v_margin=2.0)
    assert min(v_values(got)) >= 84.0
    assert 83.2 not in v_values(got), "this is what we were doing"


def test_scaffolds_are_drawn_near_the_floor_not_by_archive_rank(pool):
    got = mod.framework_donors("T", "IGHV3-74", 2)
    assert v_values(got)[:2] == [83.2, 84.9], "lowest V first"


def test_ranking_by_archive_score_can_still_be_asked_for(pool):
    got = mod.framework_donors("T", "IGHV3-74", 6, prefer_low_v=False)
    assert len(got) >= 1, "the old behaviour remains reachable"


def test_the_j_floor_is_never_traded_away(pool):
    """A design can cost 7.14 points of J identity at worst, so the floor
    itself is not negotiable even when V headroom is relaxed."""
    frame = pool.copy()
    frame.loc[2, "j_identity"] = 80.0
    frame.to_parquet(mod.DATA_DIR / "dataset" / "T_nanobodies.parquet")
    got = mod.framework_donors("T", "IGHV3-74", 6)
    assert 83.2 not in v_values(got), "J below the floor disqualifies it"
