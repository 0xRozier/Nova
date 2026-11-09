"""Compare CDR3 generation strategies by hydropathy and AA distribution.

Loads the archive, picks CDR3 donors of a given length, generates N candidates
per strategy (chimeric, pssm, pssm_denovo), and reports:

  - median Kyte-Doolittle mean hydropathy per CDR3
  - % of CDR3s with mean hydropathy under a threshold
  - AA frequency table (mutable positions only)
  - Jensen-Shannon divergence vs archive AA frequencies

Runs on CPU, no TNP or GPU needed.

Usage:
    python -m scripts.compare_generators --n 500
    python -m scripts.compare_generators --n 500 --cdr3-length 17
    python -m scripts.compare_generators --generators chimeric pssm pssm_denovo
"""
import argparse
from collections import Counter

import numpy as np
import pandas as pd

from nova_lab import DATA_DIR
from nova_lab.config import NANOBODY_TARGET
from nova_lab.cdr3_generators import (AAS, J_TAIL, ChimericGenerator,
                                      PSSMGenerator)
from nova_lab.nanobody_features import HYDROPATHY

import random


def mean_hydropathy(cdr3: str) -> float:
    """Mean Kyte-Doolittle hydropathy over residues."""
    if not cdr3:
        return 0.0
    return sum(HYDROPATHY.get(a, 0.0) for a in cdr3) / len(cdr3)


def aa_frequencies(cdr3s: list[str], j_tail: int) -> dict[str, float]:
    """AA frequencies over the mutable region (excluding J-tail)."""
    counts: Counter = Counter()
    total = 0
    for seq in cdr3s:
        mutable = seq[:max(len(seq) - j_tail, 0)]
        counts.update(mutable)
        total += len(mutable)
    if total == 0:
        return {a: 0.0 for a in AAS}
    return {a: counts.get(a, 0) / total for a in AAS}


def js_divergence(p: dict[str, float], q: dict[str, float]) -> float:
    """Jensen-Shannon divergence between two AA distributions."""
    keys = sorted(set(p) | set(q))
    pv = np.array([p.get(k, 0.0) for k in keys], dtype=np.float64)
    qv = np.array([q.get(k, 0.0) for k in keys], dtype=np.float64)
    # Normalise to proper distributions.
    pv = pv / pv.sum() if pv.sum() > 0 else pv
    qv = qv / qv.sum() if qv.sum() > 0 else qv
    m = 0.5 * (pv + qv)
    # KL(p||m) and KL(q||m), with 0*log(0) = 0.
    def kl(a, b):
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(a > 0, np.log2(a / np.where(b > 0, b, 1.0)), 0.0)
        return np.sum(a * ratio)
    return float(0.5 * kl(pv, m) + 0.5 * kl(qv, m))


def main():
    parser = argparse.ArgumentParser(description="Compare CDR3 generation strategies")
    parser.add_argument("--n", type=int, default=500, help="candidates per strategy")
    parser.add_argument("--cdr3-length", type=int, default=None,
                        help="CDR3 length (default: most common in archive)")
    parser.add_argument("--target", default=NANOBODY_TARGET[0])
    parser.add_argument("--generators", nargs="+",
                        default=["chimeric", "pssm", "pssm_denovo"],
                        choices=["chimeric", "pssm", "pssm_denovo"])
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--hydro-threshold", type=float, default=-0.5,
                        help="hydropathy threshold for reporting %% under")
    args = parser.parse_args()

    ds = DATA_DIR / "dataset"
    path = ds / f"{args.target}_nanobodies.parquet"
    if not path.exists():
        print(f"archive not found at {path}; run scripts.build_nanobody_dataset first")
        return
    df = pd.read_parquet(path)
    all_cdr3s = df["cdr3"].dropna().tolist()

    # Pick the target CDR3 length.
    lengths = pd.Series([len(c) for c in all_cdr3s])
    cdr3_len = args.cdr3_length or int(lengths.mode().iloc[0])
    donors = [c for c in all_cdr3s if len(c) == cdr3_len]
    print(f"target: {args.target}, CDR3 length: {cdr3_len}, "
          f"archive donors of this length: {len(donors)}")
    if not donors:
        print("no donors of this length in the archive")
        return

    archive_freq = aa_frequencies(donors, J_TAIL)

    # Build generators.
    chimeric = ChimericGenerator()
    pssm = PSSMGenerator(all_cdr3s)
    rng = random.Random(args.seed)

    strategies: dict[str, list[str]] = {}
    for name in args.generators:
        gen_rng = random.Random(rng.randint(0, 2**32))
        cdr3s = []
        for _ in range(args.n):
            donor = gen_rng.choice(donors)
            if name == "chimeric":
                cdr3s.append(chimeric.generate_cdr3(donor, cdr3_len, J_TAIL, gen_rng))
            elif name == "pssm":
                cdr3s.append(pssm.generate_cdr3(donor, cdr3_len, J_TAIL, gen_rng))
            elif name == "pssm_denovo":
                cdr3s.append(pssm.generate_de_novo(donor, cdr3_len, J_TAIL, gen_rng))
        strategies[name] = cdr3s

    # Also add the raw archive donors for comparison.
    strategies["archive"] = donors[:args.n] if len(donors) >= args.n else donors

    # Report.
    print(f"\n{'Strategy':<16} {'Median hydro':>13} {'Mean hydro':>11} "
          f"{'% < {:.1f}'.format(args.hydro_threshold):>10} {'JS vs archive':>14}")
    print("-" * 68)
    for name, cdr3s in strategies.items():
        hydros = [mean_hydropathy(c) for c in cdr3s]
        med = np.median(hydros)
        mn = np.mean(hydros)
        under = 100.0 * sum(1 for h in hydros if h < args.hydro_threshold) / len(hydros)
        freq = aa_frequencies(cdr3s, J_TAIL)
        js = js_divergence(freq, archive_freq)
        print(f"{name:<16} {med:>13.3f} {mn:>11.3f} {under:>9.1f}% {js:>14.4f}")

    # AA frequency table.
    print(f"\nAA frequencies in mutable region (CDR3 positions 1..{cdr3_len - J_TAIL}):")
    header = f"{'AA':<4}" + "".join(f"{name:<14}" for name in strategies)
    print(header)
    print("-" * len(header))
    for aa in AAS:
        row = f"{aa:<4}"
        for name, cdr3s in strategies.items():
            freq = aa_frequencies(cdr3s, J_TAIL)
            row += f"{freq[aa]:>12.3f}  "
        print(row)


if __name__ == "__main__":
    main()
