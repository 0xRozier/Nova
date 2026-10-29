"""What does it take to win the molecule competition?

Reads data/scores/molecules.parquet (scripts/fetch_scores.py) and prints, per
epoch: the allowed reaction, how many molecules were scored (~20 per valid
submission), and the best achievable 20-molecule total, which is an upper
bound on the winner's final_molecule_score.

Usage:
    python -m scripts.analyze_scores [--target O15427] [--last 30]
"""
import argparse

import pandas as pd

from nova_lab import DATA_DIR
from nova_lab.config import NUM_MOLECULES


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default="O15427")
    parser.add_argument("--last", type=int, default=30)
    args = parser.parse_args()

    df = pd.read_parquet(DATA_DIR / "scores" / "molecules.parquet")
    df = df[df["target"] == args.target].copy()
    df["rxn"] = df["name"].str.split(":").str[1].astype(int)

    print(f"{len(df)} scored molecules on {args.target}, epochs {df.epoch.min()} -> {df.epoch.max()}")
    print(f"score: mean {df.score.mean():.4f}  p90 {df.score.quantile(.9):.4f}  "
          f"p99 {df.score.quantile(.99):.4f}  max {df.score.max():.4f}\n")

    rows = []
    for epoch, g in df.groupby("epoch"):
        top = g.nlargest(NUM_MOLECULES, "score")["score"]
        rows.append({
            "epoch": epoch,
            "rxn": g["rxn"].mode().iat[0],
            "n_mols": len(g),
            "~subs": round(len(g) / NUM_MOLECULES, 1),
            "best_mol": g["score"].max(),
            "top20_sum": top.sum() if len(top) == NUM_MOLECULES else float("nan"),
            "median_mol": g["score"].median(),
        })
    per_epoch = pd.DataFrame(rows).set_index("epoch")
    with pd.option_context("display.float_format", "{:.4f}".format, "display.width", 120):
        print(per_epoch.tail(args.last))

        print("\nBy reaction (all epochs):")
        by_rxn = df.groupby("rxn")["score"].agg(["count", "mean", "median", "max"])
        by_rxn["p95"] = df.groupby("rxn")["score"].quantile(0.95)
        print(by_rxn.sort_values("p95", ascending=False))


if __name__ == "__main__":
    main()
