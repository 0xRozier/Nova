"""What kinds of epoch does the subnet draw, and are we equipped for each?

The epoch's (V gene, CDR3 length) pair comes from the block hash through
`random.Random`, so its distribution follows from the algorithm:
a V neighbourhood category uniformly among four, then a gene uniformly inside
it, then a CDR3 length uniformly in 5..23. IGHV3-23, IGHV3-66 and IGHV3-74 each
carry a whole category, so each is drawn a quarter of the time, while the six
genes sharing category V4 get a twenty-fourth each.

Crossed with the archive, that says how often the generator has material:
a framework donor calling the epoch's gene, and a CDR3 of the required length,
both from sequences that clear today's germline floors.

Usage:
    python -m scripts.epoch_tasks                 # theory + archive coverage
    python -m scripts.epoch_tasks --simulate 200000   # check the draw empirically
    python -m scripts.epoch_tasks --worst 10      # the least covered combinations
"""
import argparse
import random
from collections import Counter

import pandas as pd

from nova_lab import DATA_DIR
from nova_lab.config import (CDR3_LENGTH_RANGE, MIN_HUMAN_J_IDENTITY, MIN_HUMAN_V_IDENTITY,
                             V_NEIGHBORHOOD_CATEGORIES)

GERMLINE_MARGIN = 2.0  # same headroom scripts.propose_nanobody requires of donors
MIN_DONORS = 10        # below this the generator has little to graft onto
MIN_CDR3 = 5


def theoretical_gene_distribution() -> dict:
    """P(gene): 1/4 per category, then uniform inside it."""
    p = {}
    for genes in V_NEIGHBORHOOD_CATEGORIES.values():
        for gene in genes:
            p[gene] = 0.25 / len(genes)
    return p


def simulate(n: int) -> pd.DataFrame:
    """Draw n epochs the way the validator does, from random block hashes."""
    rows = []
    rng_outer = random.Random(0)
    for _ in range(n):
        seed = rng_outer.getrandbits(256)
        rng = random.Random(seed)
        genes = V_NEIGHBORHOOD_CATEGORIES[rng.choice(list(V_NEIGHBORHOOD_CATEGORIES))]
        gene = rng.choice(genes) if len(genes) > 1 else genes[0]
        rows.append((gene, rng.choice(list(CDR3_LENGTH_RANGE))))
    return pd.DataFrame(rows, columns=["gene", "cdr3_length"])


def coverage(target: str = "Q9NZQ7") -> pd.DataFrame:
    """Compliant archived sequences per gene, and CDR3s per length."""
    df = pd.read_parquet(DATA_DIR / "dataset" / f"{target}_nanobodies.parquet")
    compliant = df[(df["v_identity"] >= MIN_HUMAN_V_IDENTITY + GERMLINE_MARGIN)
                   & (df["j_identity"] >= MIN_HUMAN_J_IDENTITY)].copy()
    compliant["len3"] = compliant["cdr3"].str.len()
    rows = []
    for gene, p_gene in theoretical_gene_distribution().items():
        donors = compliant[compliant["v_genes"].str.split(",").apply(lambda g: gene in g)]
        for length in CDR3_LENGTH_RANGE:
            rows.append({
                "gene": gene, "cdr3_length": length,
                "p_epoch": p_gene / len(CDR3_LENGTH_RANGE),
                "framework_donors": len(donors),
                "cdr3_donors": int((compliant["len3"] == length).sum()),
                "exact_match": int(((compliant["len3"] == length)
                                    & compliant["v_genes"].str.split(",").apply(lambda g: gene in g)).sum()),
            })
    out = pd.DataFrame(rows)
    out["equipped"] = (out["framework_donors"] >= MIN_DONORS) & (out["cdr3_donors"] >= MIN_CDR3)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--simulate", type=int, default=0)
    parser.add_argument("--worst", type=int, default=8)
    parser.add_argument("--target", default="Q9NZQ7")
    args = parser.parse_args()

    theory = theoretical_gene_distribution()
    print("How often each V gene is imposed (from the draw itself):")
    for gene, p in sorted(theory.items(), key=lambda kv: -kv[1]):
        print(f"  {gene:<10} {p:6.2%}   1 epoch in {1/p:5.1f}")
    print(f"\nCDR3 length: uniform over {min(CDR3_LENGTH_RANGE)}..{max(CDR3_LENGTH_RANGE)} "
          f"({1/len(CDR3_LENGTH_RANGE):.2%} each)")

    if args.simulate:
        sim = simulate(args.simulate)
        print(f"\nEmpirical check over {args.simulate:,} simulated draws:")
        for gene, share in sim["gene"].value_counts(normalize=True).items():
            print(f"  {gene:<10} {share:6.2%}  (theory {theory[gene]:.2%})")

    cov = coverage(args.target)
    equipped = cov.loc[cov["equipped"], "p_epoch"].sum()
    print(f"\n=== Archive coverage for {args.target} ===")
    print(f"epochs we are equipped for (>= {MIN_DONORS} framework donors "
          f"and >= {MIN_CDR3} CDR3 donors): {equipped:.1%}")

    per_gene = cov.groupby("gene").agg(
        framework_donors=("framework_donors", "first"),
        p_epoch=("p_epoch", "sum"),
        equipped_share=("equipped", "mean"),
    ).sort_values("p_epoch", ascending=False)
    print("\nper gene:")
    print(per_gene.to_string(float_format=lambda v: f"{v:.3f}"))

    print(f"\nleast covered combinations we could still be asked for "
          f"(by framework donors, then exact matches):")
    worst = cov[cov["p_epoch"] > 0].nsmallest(args.worst, ["framework_donors", "exact_match"])
    print(worst[["gene", "cdr3_length", "p_epoch", "framework_donors", "cdr3_donors",
                 "exact_match", "equipped"]].to_string(index=False,
                                                       float_format=lambda v: f"{v:.4f}"))

    out = DATA_DIR / "submissions" / f"{args.target}_epoch_coverage.csv"
    cov.to_csv(out, index=False)
    print(f"\nfull table saved to {out}")


if __name__ == "__main__":
    main()
