"""CDR3 generation strategies for nanobody design.

The original propose_nanobody.py uses uniform random mutations over 20 amino
acids, but 40% of AAs are hydrophobic while natural nanobody CDR3s are enriched
in polar residues (S=10.5%, Y=10.5%, G=5.8%).  This pushes candidates to a
median PSH of ~98, well above the validator's PSH <= 89 threshold.

Three strategies behind a common interface:

  ChimericGenerator     -- uniform random (baseline, reproduces old behaviour)
  HydrophilicGenerator  -- the same, biased against the residues PSH counts
  PSSMGenerator      -- position-specific frequencies from archive CDR3s
  IgLMGenerator      -- language model infilling (requires torch + iglm)

Every generator exposes
``generate_cdr3(donor_cdr3, cdr3_length, j_tail, rng, framework=None)``
and returns a string of exactly *cdr3_length* residues, with the last *j_tail*
residues preserved from the donor.  The optional *framework* dict
(``{"FR1":…, "CDR1":…, …}``) provides the full VH context that IgLM needs;
Chimeric and PSSM generators ignore it.
"""
from __future__ import annotations

import random
from collections import Counter

AAS = "ACDEFGHIKLMNPQRSTVWY"
J_TAIL = 4

# Kyte-Doolittle hydropathy: positive is hydrophobic. This is the scale TNP's
# PSH is built on, which is why it is the one to bias against.
HYDROPATHY = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
              "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
              "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
              "Y": -1.3, "V": 4.2}


# ---------------------------------------------------------------------------
# Chimeric (baseline): uniform random mutations
# ---------------------------------------------------------------------------

class ChimericGenerator:
    """Current behaviour: uniform random point mutations. Baseline."""

    def generate_cdr3(self, donor_cdr3: str, cdr3_length: int,
                      j_tail: int, rng: random.Random,
                      framework: dict | None = None) -> str:
        out = list(donor_cdr3[:cdr3_length].ljust(cdr3_length, "G"))
        mutable = range(max(len(out) - j_tail, 0))
        if not mutable:
            return "".join(out)
        n_mut = rng.randint(2, 4)
        for pos in rng.sample(list(mutable), min(n_mut, len(mutable))):
            out[pos] = rng.choice(AAS)
        return "".join(out)


# ---------------------------------------------------------------------------
# Hydrophilic: the chimeric mutation, aimed at the gate that refuses us
# ---------------------------------------------------------------------------

class HydrophilicGenerator:
    """Chimeric mutations, but spent where PSH is made and away from it.

    MEASURED, AND IT DOES NOT WORK. Kept so the idea is not had twice.

    The number it exists to move: at CDR3 18 a fresh chimeric candidate reads a
    median PSH of 100.6 against the accepted field's 83.5, and only 15 % come
    in under the validator's 89
    (results/psh_fresh_candidates_by_length.csv). The README has blamed the
    uniform mutation alphabet for that since this file was written -- "40 % of
    AAs are hydrophobic while natural nanobody CDR3s are enriched in polar
    residues" -- so this biases both halves of the mutation towards polar.

    On the sequence it does exactly that: median CDR3 hydropathy drops by more
    than 2 units (tests/test_generators.py). On the quantity that decides the
    epoch it does almost nothing. 100 candidates per arm at CDR3 18, same
    donors, same seed, profiled on the validator's own install
    (results/psh_generator_ab_cdr3_18.csv):

        chimeric      median PSH 100.2    14/100 at or under 89
        hydrophilic   median PSH  97.4    13/100 at or under 89

    2.8 points of median (Mann-Whitney p = 0.064) and no change at all in the
    rate that matters (Fisher p = 1.0). PSH is a *patch* score on a predicted
    structure: what counts is which residues end up exposed and next to each
    other, not the CDR3's summed hydropathy. The same conclusion comes out of
    trying to predict PSH from the sequence at all -- within a fixed CDR3
    length, a LightGBM on the 147 region features reaches Spearman 0.25-0.31,
    and plain summed CDR hydropathy 0.10-0.22.

    So the lever for a long CDR3 is not a cleverer alphabet, it is measuring
    more candidates and taking the tail: the best of these 100 reads 74.7.
    scripts.mine sizes its pool for that. If this is revisited, the thing to
    bias is exposure, which needs the structure, which is the 3.5 s we are
    trying to avoid.

    Two changes, both on the same axis:

      - the positions to mutate are drawn towards the donor's most hydrophobic
        residues, since those are the ones contributing to the patch;
      - the replacement is drawn away from hydrophobic residues.

    Both are soft, exp(+-beta * hydropathy), because a CDR3 of nothing but
    charged residues is not a binder: `beta` is the one knob and 0 reproduces
    ChimericGenerator exactly. The J tail stays untouched, as everywhere else
    here -- mutating it cost 285 candidates in 800 (README, "the end of the
    CDR3 is coded by the J gene").
    """

    def __init__(self, beta: float = 0.4):
        self.beta = beta
        self._weights = [2.718281828459045 ** (-beta * HYDROPATHY[a]) for a in AAS]

    def generate_cdr3(self, donor_cdr3: str, cdr3_length: int,
                      j_tail: int, rng: random.Random,
                      framework: dict | None = None) -> str:
        out = list(donor_cdr3[:cdr3_length].ljust(cdr3_length, "G"))
        mutable = list(range(max(len(out) - j_tail, 0)))
        if not mutable:
            return "".join(out)
        n_mut = rng.randint(2, 4)
        positions = self._pick_positions(out, mutable, min(n_mut, len(mutable)), rng)
        for pos in positions:
            out[pos] = rng.choices(AAS, weights=self._weights, k=1)[0]
        return "".join(out)

    def _pick_positions(self, residues, mutable, n, rng):
        """Without replacement, favouring the hydrophobic positions."""
        pool = list(mutable)
        chosen = []
        for _ in range(n):
            weights = [2.718281828459045 ** (self.beta * HYDROPATHY.get(residues[p], 0.0))
                       for p in pool]
            pick = rng.choices(range(len(pool)), weights=weights, k=1)[0]
            chosen.append(pool.pop(pick))
        return chosen


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
                      j_tail: int, rng: random.Random,
                      framework: dict | None = None) -> str:
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

def _assemble_regions(r: dict) -> str:
    """Concatenate framework + CDR regions into a full VH sequence."""
    return r["FR1"] + r["CDR1"] + r["FR2"] + r["CDR2"] + r["FR3"] + r["CDR3"] + r["FR4"]


def _cdr3_start(r: dict) -> int:
    """Return the 0-based index where CDR3 begins in the assembled VH."""
    return len(r["FR1"]) + len(r["CDR1"]) + len(r["FR2"]) + len(r["CDR2"]) + len(r["FR3"])


def _extract_infill(returned: str, parent: str, start: int, end: int) -> str | None:
    """The residues IgLM generated in place of ``parent[start:end]``.

    iglm.model.utils.iglm_to_infilled rebuilds the answer as
    ``prefix + generated_infill + suffix``, and the infill's length is whatever
    the model produced -- nothing constrains it to ``end - start`` (generation
    runs to the CLS token, bounded only by max_length=150). So the infill ends
    where the unchanged suffix begins, counted from the RIGHT.

    Slicing ``returned[start:end]`` instead, as this module used to, is wrong
    the moment the infill is a different length: it returns part of the infill
    plus part of FR4. Worse, that slice is ``end - start`` wide by construction,
    so the ``len(cdr3) == cdr3_length`` guard meant to catch variable-length
    infills could never fire, and the first sequence was always accepted.

    Returns None when the suffix does not match, i.e. when the answer is not a
    usable infill of this parent.
    """
    suffix = parent[end:]
    if suffix and not returned.endswith(suffix):
        return None
    infill_end = len(returned) - len(suffix)
    if infill_end < start:
        return None
    return returned[start:infill_end]


def _in_stable_order(sequences, rng: random.Random) -> list[str]:
    """IgLM returns ``list(set(...))``, whose order varies between processes.

    Taking ``seqs[0]`` therefore makes the generator unreproducible from a
    seed, which is how three epochs once submitted the same sequence. Sort,
    then let the seeded rng choose.
    """
    ordered = sorted(sequences)
    return rng.sample(ordered, len(ordered))


class IgLMGenerator:
    """CDR3 infilling via IgLM camel heavy-chain model.

    Lazy-loads torch and iglm on first use.  Raises ImportError with a clear
    message if the dependencies are not installed.

    IgLM's real API::

        model.infill(
            sequence,            # full VH string
            chain_token,         # "[HEAVY]"
            species_token,       # "[CAMEL]"
            infill_range=(start, end),
            num_to_generate=N,
        )

    It returns ``list[str]`` — full VH sequences with the infill region
    replaced by model-generated residues.  The infill can be variable-length,
    so we filter for exact CDR3 length afterwards.
    """

    def __init__(self, device: str = "cpu"):
        self._device = device
        self._model = None

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
                      j_tail: int, rng: random.Random,
                      framework: dict | None = None) -> str:
        """Infill the mutable part of CDR3 using IgLM, preserving J-tail.

        *framework* is a regions dict (FR1, CDR1, …, FR4).  When provided the
        full VH context is given to IgLM so the model can condition on the
        surrounding sequence.  Without it we fall back to chimeric-style
        random mutation (IgLM needs a full VH to work).
        """
        if framework is None:
            # No framework context — fall back to simple random mutation.
            out = list(donor_cdr3[:cdr3_length].ljust(cdr3_length, "G"))
            tail_start = max(cdr3_length - j_tail, 0)
            for pos in rng.sample(list(range(tail_start)), min(2, tail_start)):
                out[pos] = rng.choice(AAS)
            return "".join(out)

        self._load()
        # Build the parent VH with the donor CDR3 in place.
        parent_regions = dict(framework)
        parent_regions["CDR3"] = donor_cdr3[:cdr3_length].ljust(cdr3_length, "G")
        parent_seq = _assemble_regions(parent_regions)
        start = _cdr3_start(parent_regions)
        # Infill just the mutable part (not the J-tail).
        mutable_end = start + max(cdr3_length - j_tail, 0)

        if mutable_end <= start:
            return donor_cdr3[:cdr3_length]

        seqs = self._model.infill(
            parent_seq, "[HEAVY]", "[CAMEL]",
            infill_range=(start, mutable_end),
            num_to_generate=100,
        )

        # Extract CDR3 from returned full VH sequences and filter for length.
        j_tail_seq = donor_cdr3[cdr3_length - j_tail:cdr3_length] if j_tail else ""
        for seq in _in_stable_order(seqs, rng):
            infilled = _extract_infill(seq, parent_seq, start, mutable_end)
            if infilled is None:
                continue
            cdr3 = infilled + j_tail_seq
            if len(cdr3) == cdr3_length:
                return cdr3

        # Fallback: no IgLM result matched the length — random mutation.
        out = list(donor_cdr3[:cdr3_length].ljust(cdr3_length, "G"))
        tail_start = max(cdr3_length - j_tail, 0)
        for pos in rng.sample(list(range(tail_start)), min(2, tail_start)):
            out[pos] = rng.choice(AAS)
        return "".join(out)

    def generate_batch(self, framework: dict, cdr3_length: int,
                       j_tail: int, n: int, rng: random.Random) -> list[str]:
        """Generate up to *n* CDR3 candidates in one IgLM call.

        Much faster than calling generate_cdr3() *n* times because IgLM's
        ``infill()`` supports ``num_to_generate`` to produce many sequences in
        a single forward pass.
        """
        self._load()
        parent_seq = _assemble_regions(framework)
        start = _cdr3_start(framework)
        mutable_end = start + max(cdr3_length - j_tail, 0)

        if mutable_end <= start:
            return [framework["CDR3"][:cdr3_length]] * n

        seqs = self._model.infill(
            parent_seq, "[HEAVY]", "[CAMEL]",
            infill_range=(start, mutable_end),
            num_to_generate=n * 3,
        )

        j_tail_seq = framework["CDR3"][-j_tail:] if j_tail else ""
        results: list[str] = []
        for seq in _in_stable_order(seqs, rng):
            infilled = _extract_infill(seq, parent_seq, start, mutable_end)
            if infilled is None:
                continue
            cdr3 = infilled + j_tail_seq
            if len(cdr3) == cdr3_length:
                results.append(cdr3)
            if len(results) >= n:
                break
        return results


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------

class PSSMDeNovoGenerator(PSSMGenerator):
    """PSSM sampling at every mutable position, J-tail kept from the donor.

    scripts.compare_generators reached this through its own if/elif, so the
    strategy existed but was not usable from propose_nanobody or mine.py, which
    both build generators through make_generator.
    """

    def generate_cdr3(self, donor_cdr3: str, cdr3_length: int,
                      j_tail: int, rng: random.Random,
                      framework: dict | None = None) -> str:
        return self.generate_de_novo(donor_cdr3, cdr3_length, j_tail, rng)


GENERATORS = {"chimeric", "hydrophilic", "pssm", "pssm_denovo", "iglm"}


def make_generator(name: str, cdr3_sequences: list[str] | None = None,
                   device: str = "cpu"):
    """Build a generator by name.  *cdr3_sequences* is required for ``pssm``."""
    if name == "chimeric":
        return ChimericGenerator()
    if name == "hydrophilic":
        return HydrophilicGenerator()
    if name == "pssm":
        if not cdr3_sequences:
            raise ValueError("PSSMGenerator requires cdr3_sequences from the archive")
        return PSSMGenerator(cdr3_sequences)
    if name == "pssm_denovo":
        if not cdr3_sequences:
            raise ValueError("PSSMDeNovoGenerator requires cdr3_sequences from the archive")
        return PSSMDeNovoGenerator(cdr3_sequences)
    if name == "iglm":
        return IgLMGenerator(device=device)
    raise ValueError(f"unknown generator {name!r}; choose from {sorted(GENERATORS)}")
