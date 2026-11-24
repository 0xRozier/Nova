"""The target analysis, locked on 4O9H -- the case FINDINGS 73 did by hand.

Every number asserted here was first obtained by hand and then reproduced by
the script. The point of the file is that the next target gets the same
analysis without a human reading mmCIF columns, and that a silent change in
numbering conventions becomes a red test rather than a batch aimed at the
wrong patch of protein.
"""
from pathlib import Path

import pytest

from nova_lab.nanobodies import regions
from nova_lab.structures import (Structure, best_offset, parse_res_index,
                                 polymers, vh_domain)
from scripts.analyze_target import (classify_chains, epitope_check,
                                    make_submittable, numbering_audit,
                                    paratope_split)

CIF = Path("/home/user/nova-upstream/data/structures/4O9H.cif")
H61H7 = ("EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYRMYWVRQPPGKGLEWVSAISAGGGSTYYGDSVKGRFTIS"
         "RDNAKNTVYLQMNSLKPEDTAVYYCANRAGWGMGDYWGQGTQVTVSS")
LIGHT = ("QTVVTQEPSLSVSPGGTVTLTCGLSSGSVTASNYPGWFQQTPGQAPRALIYSTNDRHSGVPSRFSGSISGN"
         "KAALTITGAQPEDEADYYCALDIGDITEFGGGTHLTVLG")

needs_cif = pytest.mark.skipif(not CIF.exists(), reason=f"{CIF} absent")


def test_parse_res_index_handles_both_forms():
    assert parse_res_index("21..186")[:3] == [21, 22, 23]
    assert len(parse_res_index("21..186")) == 166
    assert parse_res_index("24,77,80,82,131,184..186") == [24, 77, 80, 82, 131, 184, 185, 186]


def test_best_offset_recovers_a_shift():
    reference = "MAPVPPGEDSKDVAAPHRQPLTSSERIDKQIRYILDGISALRKETCNKSNMCESSKEALAENNLN"
    k, frac = best_offset(reference[10:50], reference)
    assert (k, frac) == (10, 1.0)


def test_vh_domain_takes_a_heavy_chain_and_refuses_a_light_one():
    assert vh_domain(H61H7) == H61H7
    assert vh_domain(LIGHT) is None


@needs_cif
def test_the_entity_to_chain_map_is_read_not_guessed():
    """A Fab heavy chain's first residues are often disordered."""
    by_chain = {c: e["sequence"] for e in polymers(CIF) for c in e["chains"]}
    assert set(by_chain) == {"A", "H", "L"}
    assert by_chain["A"].startswith("MAPVPPGEDSKDVAAPHRQ")
    assert vh_domain(by_chain["H"]) == H61H7
    assert vh_domain(by_chain["L"]) is None


@needs_cif
def test_numbering_is_canonical_and_the_window_fits():
    """label - auth is 2 on 4O9H: reading the author numbering misses by two."""
    st = Structure.read(CIF)
    audit = numbering_audit(st, "A", "21..186")
    assert audit["ok"], audit["why"]
    assert audit["label_minus_auth"] == [2]
    assert audit["single_offset"]
    assert audit["observed"] == 152
    assert audit["span"] == (21, 186)
    assert audit["gaps"] == [(50, 63), (106, 109)]


@needs_cif
def test_the_config_binding_site_is_the_real_epitope():
    st = Structure.read(CIF)
    e = epitope_check(st, "A", "24,77,80,82,131,184..186")
    assert len(e["inside"]) == len(e["site"]) == 8
    assert e["missing"] == [] and e["absent"] == []
    assert max(e["epitope"][n] for n in e["inside"]) <= 5.0


@needs_cif
def test_the_site_is_reachable_by_one_domain():
    """6 of the 8 are heavy-chain-only and the light chain reaches 2: the
    light chain's own 7 contacts are excluded from the config site."""
    st = Structure.read(CIF)
    e = epitope_check(st, "A", "24,77,80,82,131,184..186")
    assert e["single_chain_covers"]["H"] == 8
    assert e["single_chain_covers"]["L"] == 2


@needs_cif
def test_61h7_paratope_is_mostly_not_cdr3():
    """25% CDR3 is what makes "keep the framework" coherent for this one."""
    st = Structure.read(CIF)
    kinds = classify_chains(CIF, polymers(CIF)[0]["sequence"])
    assert kinds["target"] == ["A"] and kinds["binder"] == ["H"]
    p = paratope_split(st, "H", ["A"], H61H7)
    assert p["ok"], p["why"]
    assert p["total"] == 16
    assert p["share"]["CDR3"] == pytest.approx(0.25, abs=0.01)
    assert p["share"]["CDR1"] + p["share"]["CDR2"] == pytest.approx(0.69, abs=0.02)


def test_making_61h7_submittable_costs_three_mutations_and_no_cdr():
    fix = make_submittable(H61H7)
    assert fix is not None
    assert fix["n_mutations"] == 3
    assert fix["submittable"] and fix["layer1"] == "ok"
    before, after = regions(H61H7), regions(fix["sequence"])
    for region in ("FR1", "CDR1", "CDR2", "FR3", "CDR3"):
        assert before[region] == after[region], region
    assert after["FR4"] == "WGQGTLVTVSS"


def test_a_sequence_that_already_passes_needs_no_fix():
    """None means "nothing to do", not "could not be fixed"."""
    fix = make_submittable(H61H7)
    assert make_submittable(fix["sequence"]) is None
