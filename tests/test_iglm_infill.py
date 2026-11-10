"""IgLM's infill is variable-length; check we read it as such, without a GPU.

Issue #1 asks for these two properties to be verified on the pod before a live
run. They do not need a pod: iglm's own source says what infill returns, and a
fake model standing in for it exercises our side of the contract. From
iglm/model/utils.py:

    def iglm_to_infilled(token_seq, tokenizer):
        infilled_seq = concatenate([
            token_seq[:mask_token_idx],              # prefix
            token_seq[sep_token_idx + 1:cls_token_idx],   # GENERATED infill
            token_seq[mask_token_idx + 1:sep_token_idx],  # suffix
        ])

    def mask_span(seq, start, end, append_span=False):
        masked_seq = seq[:start] + [MASK_TOKEN] + seq[end:] + [SEP_TOKEN]

So the answer is prefix + infill + suffix, the infill runs until the CLS token
(bounded only by max_length=150) and its length is NOT end - start. And
_generate returns `list(decoded_seqs)` built from a `set`, so the order is not
stable between processes.
"""
import random

import pytest

from nova_lab.cdr3_generators import (IgLMGenerator, J_TAIL, _extract_infill,
                                      _in_stable_order)

FRAMEWORK = {
    "FR1": "QVQLVESGGGLVQPGGSLRLSCAAS", "CDR1": "GFTFSSYA", "FR2": "MHWVRQAPGKGLEWVS",
    "CDR2": "AISGSGGST", "FR3": "YYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC",
    "CDR3": "ARDRGLGDY", "FR4": "WGQGTLVTVSS",
}


def parent_of(framework):
    return "".join(framework[k] for k in
                   ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"))


def test_extract_infill_handles_a_longer_answer():
    """A longer infill shifts the suffix right; the slice must follow it."""
    parent = parent_of(FRAMEWORK)
    start = parent.index(FRAMEWORK["CDR3"])
    end = start + 5
    returned = parent[:start] + "WWWWWWWW" + parent[end:]  # 8 residues for 5
    assert _extract_infill(returned, parent, start, end) == "WWWWWWWW"
    # The old fixed-width slice silently mixed infill with what follows it.
    assert returned[start:end] == "WWWWW"


def test_extract_infill_handles_a_shorter_answer():
    parent = parent_of(FRAMEWORK)
    start = parent.index(FRAMEWORK["CDR3"])
    end = start + 5
    returned = parent[:start] + "WW" + parent[end:]
    assert _extract_infill(returned, parent, start, end) == "WW"


def test_extract_infill_rejects_a_mangled_suffix():
    """If FR4 came back changed, the answer is not an infill of this parent."""
    parent = parent_of(FRAMEWORK)
    start = parent.index(FRAMEWORK["CDR3"])
    end = start + 5
    returned = parent[:start] + "WWWWW" + parent[end:-1] + "K"
    assert _extract_infill(returned, parent, start, end) is None


def test_order_is_stable_for_a_given_seed():
    """infill returns list(set(...)): seeded runs must still agree."""
    pool = {"AAA", "BBB", "CCC", "DDD", "EEE"}
    first = _in_stable_order(list(pool), random.Random(7))
    second = _in_stable_order(list(reversed(sorted(pool))), random.Random(7))
    assert first == second


# No D in FR4 ("WGQGTLVTVSS"), so a D-only infill makes it visible when an
# extraction has run past the infill and into the framework.
MARKER = "D"


class FakeIgLM:
    """Stands in for iglm.IgLM: returns infills of the lengths asked for.

    Order is preserved, unlike the real one, so a test can put a decoy first
    and know it is seen first.
    """

    def __init__(self, lengths):
        self.lengths = lengths
        self.calls = []

    def infill(self, sequence, chain_token, species_token, infill_range,
               num_to_generate=1000, **kwargs):
        self.calls.append((infill_range, num_to_generate))
        start, end = infill_range
        return [sequence[:start] + MARKER * length + sequence[end:]
                for length in self.lengths]


def iglm_with(lengths):
    generator = IgLMGenerator.__new__(IgLMGenerator)
    generator._device = "cpu"
    generator._model = FakeIgLM(lengths)
    generator._load = lambda: None
    return generator


@pytest.mark.parametrize("cdr3_length", [9, 12])
def test_generate_cdr3_takes_the_infill_and_not_the_framework(cdr3_length):
    """The body must be the residues IgLM generated, nothing after them.

    Asserting the length alone proves nothing: the fixed-width slice was
    `mutable_end - start` wide by construction, so it always had the right
    length -- while quietly carrying FR4 residues when the infill came back
    short. Hence the check on content.
    """
    donor = "ARDRGLGDYKLM"[:cdr3_length]
    wanted = cdr3_length - J_TAIL
    # A short decoy first, then the answer of the right length. Reading the
    # first answer regardless -- which the old slice forced -- picks up
    # "DD" + "WGQGT..." and fails on the marker check below.
    generator = iglm_with([max(wanted - 2, 1), wanted])
    cdr3 = generator.generate_cdr3(donor, cdr3_length, J_TAIL,
                                   random.Random(0), framework=dict(FRAMEWORK))
    assert len(cdr3) == cdr3_length
    body, tail = cdr3[:wanted], cdr3[wanted:]
    assert set(body) == {MARKER}, f"body {body!r} is not all infill"
    assert tail == donor[cdr3_length - J_TAIL:cdr3_length]


def test_generate_cdr3_falls_back_when_no_answer_fits():
    """No usable infill: return the donor's own CDR3, not a mangled splice."""
    donor = "ARDRGLGDY"
    generator = iglm_with([1, 2, 3])  # none matches 9 - J_TAIL
    cdr3 = generator.generate_cdr3(donor, 9, J_TAIL, random.Random(0),
                                   framework=dict(FRAMEWORK))
    assert len(cdr3) == 9
    assert MARKER * 3 not in cdr3, "a rejected infill leaked into the result"


def test_generate_batch_keeps_only_well_formed_candidates():
    donor_framework = dict(FRAMEWORK)
    cdr3_length = len(donor_framework["CDR3"])
    wanted = cdr3_length - J_TAIL
    generator = iglm_with([wanted - 1, wanted + 4, wanted, wanted])
    batch = generator.generate_batch(donor_framework, cdr3_length, J_TAIL, 2,
                                     random.Random(0))
    assert len(batch) == 2
    for candidate in batch:
        assert len(candidate) == cdr3_length
        assert set(candidate[:wanted]) == {MARKER}
        assert candidate.endswith(donor_framework["CDR3"][-J_TAIL:])
