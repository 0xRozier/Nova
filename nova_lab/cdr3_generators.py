"""CDR3 generation strategies for nanobody design.

The original propose_nanobody.py uses uniform random mutations over 20 amino
acids, but 40% of AAs are hydrophobic while natural nanobody CDR3s are enriched
in polar residues (S=10.5%, Y=10.5%, G=5.8%).  This pushes candidates to a
median PSH of ~98, well above the validator's PSH <= 89 threshold.

Three strategies behind a common interface:

  ChimericGenerator  -- uniform random (baseline, reproduces old behaviour)
  PSSMGenerator      -- position-specific frequencies from archive CDR3s
  IgLMGenerator      -- language model infilling (requires torch + iglm)

Every generator exposes ``generate_cdr3(donor_cdr3, cdr3_length, j_tail, rng)``
and returns a string of exactly *cdr3_length* residues, with the last *j_tail*
residues preserved from the donor.
"""
from __future__ import annotations

import random
from collections import Counter

AAS = "ACDEFGHIKLMNPQRSTVWY"
J_TAIL = 4


# ---------------------------------------------------------------------------
# Chimeric (baseline): uniform random mutations
# ---------------------------------------------------------------------------

class ChimericGenerator:
    """Current behaviour: uniform random point mutations. Baseline."""

    def generate_cdr3(self, donor_cdr3: str, cdr3_length: int,
                      j_tail: int, rng: random.Random) -> str:
        out = list(donor_cdr3[:cdr3_length].ljust(cdr3_length, "G"))
        mutable = range(max(len(out) - j_tail, 0))
        if not mutable:
            return "".join(out)
        n_mut = rng.randint(2, 4)
        for pos in rng.sample(list(mutable), min(n_mut, len(mutable))):
            out[pos] = rng.choice(AAS)
        return "".join(out)


# ---------------------------------------------------------------------------
# PSSM: position-specific weighted mutations from archive CDR3 frequencies
# ---------------------------------------------------------------------------

class PSSMGenerator:
    """Position-specific weighted mutations from archive CDR3 frequencies.

    Builds one PSSM per CDR3 length from the archive sequences passed at
    construction.  Mutations sample from the observed column frequencies
    instead of a uniform distribution over 20 amino acids.
    """

    def __init__(self, cdr3_sequences: list[str], pseudocount: float = 0.5):
        self._pssms: dict[int, list[dict[str, float]]] = {}
        # Group CDR3s by length and build per-position frequency dicts.
        by_length: dict[int, list[str]] = {}
        for seq in cdr3_sequences:
            by_length.setdefault(len(seq), []).append(seq)
        for length, seqs in by_length.items():
            pssm: list[dict[str, float]] = []
            for pos in range(length):
                counts = Counter(s[pos] for s in seqs if pos < len(s))
                total = sum(counts.values()) + pseudocount * len(AAS)
                freq = {a: (counts.get(a, 0) + pseudocount) / total for a in AAS}
                pssm.append(freq)
            self._pssms[length] = pssm
        # Fallback: pool all CDR3s into a global per-position frequency
        # (used when the requested length has no archive examples).
        global_counts: list[Counter] = []
        max_len = max(by_length) if by_length else 0
        for pos in range(max_len):
            c = Counter()
            for seq in cdr3_sequences:
                if pos < len(seq):
                    c[seq[pos]] += 1
            global_counts.append(c)
        self._global: list[dict[str, float]] = []
        for c in global_counts:
            total = sum(c.values()) + pseudocount * len(AAS)
            self._global.append({a: (c.get(a, 0) + pseudocount) / total for a in AAS})

    def _get_pssm(self, length: int) -> list[dict[str, float]]:
        if length in self._pssms:
            return self._pssms[length]
        # No exact-length match: use the global profile, padded or truncated.
        if not self._global:
            return [{a: 1.0 / len(AAS) for a in AAS}] * length
        pssm = self._global[:length]
        while len(pssm) < length:
            pssm.append(pssm[-1])
        return pssm

    def _sample(self, freq: dict[str, float], rng: random.Random) -> str:
        aas = list(freq.keys())
        weights = [freq[a] for a in aas]
        return rng.choices(aas, weights=weights, k=1)[0]

    def generate_cdr3(self, donor_cdr3: str, cdr3_length: int,
                      j_tail: int, rng: random.Random) -> str:
        """Weighted point mutations: same number of mutations as chimeric, but
        drawn from archive column frequencies instead of uniform."""
        out = list(donor_cdr3[:cdr3_length].ljust(cdr3_length, "G"))
        pssm = self._get_pssm(cdr3_length)
        mutable = range(max(len(out) - j_tail, 0))
        if not mutable:
            return "".join(out)
        n_mut = rng.randint(2, 4)
        for pos in rng.sample(list(mutable), min(n_mut, len(mutable))):
            out[pos] = self._sample(pssm[pos], rng)
        return "".join(out)

    def generate_de_novo(self, donor_cdr3: str, cdr3_length: int,
                         j_tail: int, rng: random.Random) -> str:
        """Full PSSM sampling: every mutable position is drawn from archive
        frequencies. The J-tail is preserved from the donor."""
        pssm = self._get_pssm(cdr3_length)
        out = []
        tail_start = max(cdr3_length - j_tail, 0)
        for pos in range(cdr3_length):
            if pos >= tail_start:
                # Preserve J-tail from donor.
                if pos < len(donor_cdr3):
                    out.append(donor_cdr3[pos])
                else:
                    out.append(self._sample(pssm[pos], rng))
            else:
                out.append(self._sample(pssm[pos], rng))
        return "".join(out)


# ---------------------------------------------------------------------------
# IgLM: language model infilling (requires torch + iglm)
# ---------------------------------------------------------------------------

class IgLMGenerator:
    """CDR3 infilling via IgLM camel heavy-chain model.

    Lazy-loads torch and iglm on first use. Raises ImportError with a clear
    message if the dependencies are not installed.
    """

    def __init__(self, device: str = "cpu"):
        self._device = device
        self._model = None
        self._tokenizer = None

    def _load(self):
        if self._model is not None:
            return
        try:
            import torch  # noqa: F401
            from iglm import IgLM
        except ImportError:
            raise ImportError(
                "IgLMGenerator requires torch and iglm. "
                "Install them with: pip install torch iglm"
            )
        model = IgLM()
        if self._device != "cpu":
            model.model.to(self._device)
        self._model = model

    def generate_cdr3(self, donor_cdr3: str, cdr3_length: int,
                      j_tail: int, rng: random.Random) -> str:
        """Infill the mutable part of CDR3 using IgLM, preserving J-tail."""
        self._load()
        tail_start = max(cdr3_length - j_tail, 0)
        tail = donor_cdr3[tail_start:cdr3_length] if tail_start < len(donor_cdr3) else ""
        # IgLM infills a masked span. We mask the mutable region and let the
        # model generate it conditioned on species=camel, chain=heavy.
        parent_seq = donor_cdr3[:cdr3_length].ljust(cdr3_length, "G")
        prefix = parent_seq[:0]  # empty — infill the whole mutable part
        suffix = tail
        num_to_generate = tail_start
        if num_to_generate <= 0:
            return donor_cdr3[:cdr3_length]
        generated, _ = self._model.infill(
            prefix, suffix,
            species="camel",
            chain_type="heavy",
            num_to_generate=1,
            infill_len=num_to_generate,
        )
        if generated:
            infilled = generated[0][:num_to_generate]
            # Pad if the model returned fewer residues than requested.
            while len(infilled) < num_to_generate:
                infilled += rng.choice(AAS)
            return infilled + tail
        # Fallback: return the donor with minimal mutation.
        out = list(parent_seq)
        for pos in rng.sample(list(range(tail_start)), min(2, tail_start)):
            out[pos] = rng.choice(AAS)
        return "".join(out)


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------

GENERATORS = {"chimeric", "pssm", "iglm"}


def make_generator(name: str, cdr3_sequences: list[str] | None = None,
                   device: str = "cpu"):
    """Build a generator by name.  *cdr3_sequences* is required for ``pssm``."""
    if name == "chimeric":
        return ChimericGenerator()
    if name == "pssm":
        if not cdr3_sequences:
            raise ValueError("PSSMGenerator requires cdr3_sequences from the archive")
        return PSSMGenerator(cdr3_sequences)
    if name == "iglm":
        return IgLMGenerator(device=device)
    raise ValueError(f"unknown generator {name!r}; choose from {sorted(GENERATORS)}")
