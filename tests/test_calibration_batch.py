"""Invariants of the calibration batch, including the one that nearly broke it.

scripts.boltzgen_score names each input by the sequence's SHA-256, so a batch
that repeats a sequence buys one measurement and pays for two. The first draft
of this batch used x3 replicates to beat the sigma of FINDINGS 49 and would
have wasted a third of the GPU hour; verify_boltzgen_notebook caught it. The
duplicate test below is there so it cannot come back.
"""
import pytest

from nova_lab.config import MIN_HUMAN_J_IDENTITY, MIN_HUMAN_V_IDENTITY
from nova_lab.nanobodies import check_sequence, germline_identity, is_vhh, regions
from scripts.build_calibration_batch import BINDERS, CAMELISED, OVERLAP, build


@pytest.fixture(scope="module")
def batch():
    return build()


def test_no_duplicate_sequences(batch):
    """BoltzGen keys by sequence hash: a repeat is paid for and not measured."""
    assert batch.sequence.duplicated().sum() == 0


def test_the_ladder_spans_full_overlap_to_none(batch):
    labels = set(batch[batch.arm == "ref_intact"].label)
    assert labels == set(OVERLAP)
    assert max(OVERLAP.values()) == 8 and min(OVERLAP.values()) == 0


def test_camelising_61h7_touches_no_cdr():
    """The three mutations are at IMGT 49, 50 and in FR4 -- none is a contact."""
    raw, cam = regions(BINDERS["61H7"]), regions(CAMELISED)
    for region in ("FR1", "CDR1", "CDR2", "FR3", "CDR3"):
        assert raw[region] == cam[region], region
    assert len(CAMELISED) == len(BINDERS["61H7"])
    assert sum(1 for a, b in zip(BINDERS["61H7"], CAMELISED) if a != b) == 3


def test_the_camelised_scaffold_is_submittable():
    """Raw 61H7 is refused twice over; this is the whole point of the variant."""
    assert is_vhh(CAMELISED)[0]
    assert not is_vhh(BINDERS["61H7"])[0]
    g = germline_identity(CAMELISED)
    assert g["V"]["identity_pct"] >= MIN_HUMAN_V_IDENTITY
    assert g["J"]["identity_pct"] >= MIN_HUMAN_J_IDENTITY
    assert g["V"]["genes"][0] == "IGHV3-23"
    assert check_sequence(CAMELISED) == []


def test_keep_fw_rows_would_all_pass_the_validator(batch):
    """An arm of sequences we could not submit prices the wrong thing."""
    arm = batch[batch.arm == "keep_fw"]
    assert len(arm) > 0
    assert (arm.layer1 == "ok").all()
    for seq in arm.sequence:
        assert germline_identity(seq)["J"]["identity_pct"] >= MIN_HUMAN_J_IDENTITY


def test_the_scramble_control_keeps_the_germline_encoded_ends(batch):
    """Shuffling the J-encoded tail would fail the control on a rule, not on binding."""
    ref = regions(BINDERS["61H7"])["CDR3"]
    arm = batch[batch.arm == "scramble"]
    assert len(arm) > 0
    for seq in arm.sequence:
        cdr3 = regions(seq)["CDR3"]
        assert cdr3[:3] == ref[:3] and cdr3[-3:] == ref[-3:]
        assert cdr3 != ref
        assert sorted(cdr3) == sorted(ref)


def test_keep_fw_and_scramble_share_the_same_framework(batch):
    """The cross only reads if everything but the loop is held constant."""
    for arm in ("keep_fw", "scramble"):
        for seq in batch[batch.arm == arm].sequence:
            r = regions(seq)
            assert r["FR1"] == regions(CAMELISED)["FR1"]
            assert r["FR2"] == regions(CAMELISED)["FR2"]
            assert r["CDR1"] == regions(CAMELISED)["CDR1"]
            assert r["CDR2"] == regions(CAMELISED)["CDR2"]


def test_the_batch_is_deterministic():
    """A batch nobody can rebuild is a measurement nobody can replay."""
    assert list(build().sequence) == list(build().sequence)
