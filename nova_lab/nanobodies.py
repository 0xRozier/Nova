"""Offline nanobody validity checks, mirroring the NOVA validator.

Reimplements neurons/validator/nanobody_validity.py (and the helpers in
utils/nanobodies.py, utils/germline_identity.py) at UPSTREAM_COMMIT, plus the
sequence-only red flags of NOVA-nanobody-filter's developability filter.

What is covered (CPU, no GPU): the null token, length, allowed amino acids,
homopolymer runs, di-repeats, the signal-peptide heuristic, exact-duplicate and
region-aware similarity against the archive, human V/J germline identity, the
epoch's allowed V gene, VHH classification with IMGT hallmarks, canonical
cysteines, the epoch's allowed CDR3 length, and the CDR-length red flags.

What is NOT covered: the structure-based developability checks (TNP surface
hydrophobic patches, charge patches, CDR3 compactness, Aggrescan3D), which need
structure prediction. `check_sequence` reports them as unknown, so a sequence it
accepts is "valid unless developability rejects it".

Requires `abnumber`, `anarci` and the `hmmscan` binary (HMMER).
"""
import hashlib
from functools import lru_cache

from . import DATA_DIR
from .config import (
    ALLOWED_AAS,
    MAX_DI_REPEAT_PAIRS,
    MAX_HOMOPOLYMER_RUN,
    MAX_SEQUENCE_LENGTH,
    MIN_HUMAN_J_IDENTITY,
    MIN_HUMAN_V_IDENTITY,
    MIN_SEQUENCE_LENGTH,
    SP_HYDRO_MIN_IN_WINDOW,
    SP_SCAN_PREFIX,
    SP_WINDOW,
    TOTAL_CDR_LENGTH_RANGE,
)

HYDROPHOBIC = set("AILMFWV")
# NOVA-nanobody-filter metanano/filters/aggrescan3d.py + config.py: sequence-only
# red flags. The others (compactness, PSH, charge patches) need a structure.
CDR3_RED_FLAG_RANGE = (5, 23)


def normalize(seq: str) -> str:
    return seq.strip().upper()


def seq_hash(seq: str) -> str:
    return hashlib.sha256(seq.encode("ascii")).hexdigest()


def max_run_length(s: str) -> int:
    best = cur = 1 if s else 0
    for i in range(1, len(s)):
        cur = cur + 1 if s[i] == s[i - 1] else 1
        best = max(best, cur)
    return best


def max_di_repeat_pairs(s: str) -> int:
    """Longest ABABAB... run, counted in AB pairs ('GSGSGS' -> 3)."""
    best = 0
    for i in range(len(s) - 3):
        a, b = s[i], s[i + 1]
        if a == b:
            continue  # AA runs are the homopolymer check's business
        pairs, j = 1, i + 2
        while j + 1 < len(s) and s[j] == a and s[j + 1] == b:
            pairs += 1
            j += 2
        best = max(best, pairs)
    return best


def looks_like_signal_peptide(seq: str) -> bool:
    """True if any SP_WINDOW-long window of the N-terminus is strongly hydrophobic."""
    prefix = seq[:SP_SCAN_PREFIX]
    if len(prefix) < SP_WINDOW:
        return False
    count = sum(aa in HYDROPHOBIC for aa in prefix[:SP_WINDOW])
    if count >= SP_HYDRO_MIN_IN_WINDOW:
        return True
    for i in range(SP_WINDOW, len(prefix)):
        count -= prefix[i - SP_WINDOW] in HYDROPHOBIC
        count += prefix[i] in HYDROPHOBIC
        if count >= SP_HYDRO_MIN_IN_WINDOW:
            return True
    return False


@lru_cache(maxsize=4096)
def chain(seq: str):
    """AbNumber IMGT chain, or None when ANARCI does not recognize one domain."""
    from abnumber import Chain
    try:
        return Chain(seq, scheme="imgt", cdr_definition="imgt")
    except Exception:
        return None


def _regions_of(c) -> dict:
    return {"CDR1": c.cdr1_seq, "CDR2": c.cdr2_seq, "CDR3": c.cdr3_seq,
            "FR1": c.fr1_seq, "FR2": c.fr2_seq, "FR3": c.fr3_seq, "FR4": c.fr4_seq}


def regions(seq: str) -> dict | None:
    c = chain(seq)
    return _regions_of(c) if c is not None else None


def regions_batch(sequences: list[str]) -> list[dict | None]:
    """Regions for many sequences at once: AbNumber batches ANARCI, ~100x faster."""
    from abnumber import Chain
    try:
        chains, _ = Chain.batch({str(i): s for i, s in enumerate(sequences)},
                                scheme="imgt", cdr_definition="imgt", use_anarcii=False)
    except Exception:
        return [regions(s) for s in sequences]
    out = []
    for i in range(len(sequences)):
        c = chains.get(str(i))
        out.append(_regions_of(c) if c is not None else None)
    return out


def is_vhh(seq: str) -> tuple[bool, str]:
    """VHH on the IMGT scheme with the FR2 hallmarks (H49 in E/Q, H50 = R)."""
    c = chain(seq)
    if c is None:
        return False, "ANARCI did not return exactly one heavy-chain domain"
    if c.chain_type != "H":
        return False, "not a heavy chain"
    try:
        if c["H49"] not in {"E", "Q"} or c["H50"] != "R":
            return False, f"VHH hallmarks absent (H49={c['H49']}, H50={c['H50']})"
        if c["H23"] != "C" or c["H104"] != "C":
            return False, f"canonical cysteines absent (H23={c['H23']}, H104={c['H104']})"
    except KeyError as e:
        return False, f"missing IMGT position {e}"
    return True, ""


@lru_cache(maxsize=4096)
def germline_identity(seq: str) -> dict | None:
    """Tie-aware V/J identity on ANARCI's 128-state reference alignment.

    Mirrors utils/germline_identity.py: identity is matches over non-gap
    reference positions, and every allele reaching the maximum is reported.
    """
    from fractions import Fraction

    from anarci.anarci import check_for_j, run_hmmer
    from anarci.germlines import all_germlines

    records = [("query", seq)]
    try:
        alignments = run_hmmer(records, ncpu=1, hmmer_species=None)
        check_for_j(records, alignments, "imgt")
    except Exception:
        return None
    _, state_vectors, details = alignments[0]
    if len(state_vectors) != 1 or len(details) != 1 or details[0]["chain_type"] != "H":
        return None

    states = dict(state_vectors[0])
    query = {p: (seq[i] if (i := states.get((p, "m"))) is not None else "-") for p in range(1, 129)}

    out = {}
    for region in ("V", "J"):
        scored = []
        for allele, reference in all_germlines[region]["H"]["human"].items():
            ref = dict(enumerate(reference, start=1))
            positions = [p for p, aa in ref.items() if aa != "-"]
            matches = sum(query.get(p, "-") == ref[p] for p in positions)
            scored.append((Fraction(matches, len(positions)), allele))
        best = max(score for score, _ in scored)
        out[region] = {
            "identity_pct": float(100 * best),
            "genes": sorted({allele.split("*")[0] for score, allele in scored if score == best}),
        }
    return out


def cdr_similarity(a: dict, b: dict) -> dict:
    """Per-CDR similarity, as NOVA-nanobody-filter's search engine computes it:
    position-wise mismatches over the shared length, plus the length difference,
    divided by the longer length."""
    def similarity(x, y):
        if not x and not y:
            return 1.0
        n = min(len(x), len(y))
        distance = sum(p != q for p, q in zip(x[:n], y[:n])) + abs(len(x) - len(y))
        return max(0.0, 1.0 - distance / max(len(x), len(y), 1))
    return {k: similarity(a.get(k, ""), b.get(k, "")) for k in ("CDR1", "CDR2", "CDR3")}


def whole_identity(a: str, b: str) -> float:
    from difflib import SequenceMatcher
    return SequenceMatcher(None, a, b).ratio()


def is_duplicate(identity: float, cdr_sim: dict | None) -> tuple[bool, str, bool]:
    """The archive-duplicate rule of NOVA-nanobody-filter's search engine.

    Returns (hit, reason, certain). `certain` marks the cases a near-identical
    sequence cannot escape: the validator's k-mer prefilter is bound to retrieve
    such a close neighbour, so these are rejections rather than warnings.
    """
    if identity >= 0.95:
        return True, f"whole-sequence identity {identity:.0%}", True
    if cdr_sim is not None:
        if cdr_sim.get("CDR3", 0.0) >= 0.85:
            return True, f"CDR3 identity {cdr_sim['CDR3']:.0%} (same clonotype)", identity >= 0.85
        if (cdr_sim.get("CDR3", 0.0) >= 0.80
                and cdr_sim.get("CDR1", 0.0) >= 0.90 and cdr_sim.get("CDR2", 0.0) >= 0.90):
            return True, "CDR3 >= 80% with CDR1 and CDR2 >= 90%", False
    elif identity >= 0.90:
        return True, f"identity {identity:.0%} without CDR data", False
    return False, "", False


def sequence_red_flags(seq: str) -> list[str]:
    """Developability red flags computable from the sequence alone."""
    r = regions(seq)
    if r is None:
        return ["IMGT regions unavailable"]
    flags = []
    total = len(r["CDR1"]) + len(r["CDR2"]) + len(r["CDR3"])
    lo, hi = TOTAL_CDR_LENGTH_RANGE
    if not lo <= total <= hi:
        flags.append(f"total_cdr_length {total} outside [{lo}, {hi}]")
    lo, hi = CDR3_RED_FLAG_RANGE
    if not lo <= len(r["CDR3"]) <= hi:
        flags.append(f"cdr3_length {len(r['CDR3'])} outside [{lo}, {hi}]")
    return flags


def check_sequence(seq: str, allowed_v: str | None = None, allowed_cdr3_length: int | None = None,
                   archive: "NanobodyArchive | None" = None,
                   archive_exclude: set[str] | None = None,
                   warnings: list | None = None) -> list[str]:
    """Every problem found, in the validator's order. Empty list = valid so far.

    A structure-based developability rejection cannot be ruled out here.
    Near-identical archive matches (>= 95% whole sequence, or an exact
    duplicate) are rejections. The subtler CDR-based rules go to `warnings`: we
    compare against the whole archive while the validator only compares against
    the candidates its k-mer filter retrieves, so we are stricter than it is.
    """
    problems = []
    if "~" in seq:
        return ["null token"]
    seq = normalize(seq)
    if not MIN_SEQUENCE_LENGTH <= len(seq) <= MAX_SEQUENCE_LENGTH:
        problems.append(f"length {len(seq)} outside [{MIN_SEQUENCE_LENGTH}, {MAX_SEQUENCE_LENGTH}]")
    if invalid := set(seq) - ALLOWED_AAS:
        problems.append(f"invalid amino acids {''.join(sorted(invalid))}")
        return problems  # everything below assumes a numberable sequence
    if max_run_length(seq) > MAX_HOMOPOLYMER_RUN:
        problems.append(f"homopolymer run {max_run_length(seq)} > {MAX_HOMOPOLYMER_RUN}")
    if max_di_repeat_pairs(seq) > MAX_DI_REPEAT_PAIRS:
        problems.append(f"di-repeat pairs {max_di_repeat_pairs(seq)} > {MAX_DI_REPEAT_PAIRS}")
    if looks_like_signal_peptide(seq):
        problems.append("looks like a signal peptide")

    germline = germline_identity(seq)
    if germline is None:
        problems.append("ANARCI did not return exactly one heavy-chain domain")
    else:
        if germline["V"]["identity_pct"] < MIN_HUMAN_V_IDENTITY:
            problems.append(f"human V identity {germline['V']['identity_pct']:.1f} < {MIN_HUMAN_V_IDENTITY}")
        if germline["J"]["identity_pct"] < MIN_HUMAN_J_IDENTITY:
            problems.append(f"human J identity {germline['J']['identity_pct']:.1f} < {MIN_HUMAN_J_IDENTITY}")
        if allowed_v and allowed_v not in germline["V"]["genes"]:
            problems.append(f"V gene {germline['V']['genes']} is not the epoch's {allowed_v}")

    ok, reason = is_vhh(seq)
    if not ok:
        problems.append(reason)
    r = regions(seq)
    if r and allowed_cdr3_length and len(r["CDR3"]) != allowed_cdr3_length:
        problems.append(f"CDR3 length {len(r['CDR3'])} is not the epoch's {allowed_cdr3_length}")
    problems += sequence_red_flags(seq)

    if archive is not None:
        hit, reason, certain = archive.duplicate_of(seq, archive_exclude)
        if hit and (certain or warnings is None):
            problems.append(f"too similar to an archived sequence: {reason}")
        elif hit:
            warnings.append(f"too similar to an archived sequence: {reason}")
    return problems


class NanobodyArchive:
    """Archived sequences for a target, for the duplicate and similarity checks."""

    URL = "https://huggingface.co/datasets/Metanova/Submission-Archive/resolve/main/{target}_nanobodies.csv"

    def __init__(self, target: str = "Q9NZQ7", refresh: bool = False, with_regions: bool = True):
        import pandas as pd
        path = DATA_DIR / "archive" / f"{target}_nanobodies.csv"
        if refresh or not path.exists():
            path.parent.mkdir(parents=True, exist_ok=True)
            pd.read_csv(self.URL.format(target=target)).to_csv(path, index=False)
        df = pd.read_csv(path)
        self.sequences = list(df["sequence"])
        self.hashes = set(df["sequence_hash"])
        self.regions = self._cached_regions(target, path) if with_regions else [None] * len(self.sequences)

    def _cached_regions(self, target: str, source) -> list:
        """Numbering the whole archive takes minutes, so keep it on disk."""
        import pandas as pd
        cache = DATA_DIR / "archive" / f"{target}_regions.csv"
        keys = ("CDR1", "CDR2", "CDR3", "FR1", "FR2", "FR3", "FR4")
        if cache.exists() and cache.stat().st_mtime >= source.stat().st_mtime:
            df = pd.read_csv(cache, keep_default_na=False)
            if len(df) == len(self.sequences):
                return [None if not r["CDR3"] else {k: r[k] for k in keys} for _, r in df.iterrows()]
        out = regions_batch(self.sequences)
        pd.DataFrame([r or {k: "" for k in keys} for r in out]).to_csv(cache, index=False)
        return out

    def duplicate_of(self, seq: str, exclude: set[str] | None = None) -> tuple[bool, str, bool]:
        """(hit, reason, certain). `exclude` holds sequences to ignore, to replay
        an epoch whose own submissions are already archived (the validator's
        exclude_ids)."""
        exclude = exclude or set()
        if seq not in exclude and seq_hash(seq) in self.hashes:
            return True, "already submitted (exact duplicate)", True
        mine = regions(seq)
        best = (False, "", False)
        for other, other_regions in zip(self.sequences, self.regions):
            if other in exclude:
                continue
            identity = whole_identity(seq, other)
            if identity < 0.70:  # cheap gate: no rule can trigger this far away
                continue
            sim = cdr_similarity(mine, other_regions) if mine and other_regions else None
            hit, reason, certain = is_duplicate(identity, sim)
            if certain:
                return True, reason, True
            if hit and not best[0]:
                best = (True, reason, False)
        return best
