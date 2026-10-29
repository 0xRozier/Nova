"""Turn the scraped scores into a training set: resolve SMILES, featurize.

Output: data/dataset/<target>.parquet (name, epoch, rxn, smiles, score, ...)
and data/dataset/<target>_X.npy (features, same row order).

Usage:
    python -m scripts.build_dataset [--target O15427] [--workers 4]
"""
import argparse
from multiprocessing import Pool

import numpy as np
import pandas as pd

from nova_lab import DATA_DIR
from nova_lab.chem import boltz_score
from nova_lab.features import FEATURE_NAMES, featurize
from nova_lab.reactions import name_to_smiles


def _resolve(name):
    smi = name_to_smiles(name)
    x = featurize(smi) if smi else None
    return smi, x


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default="O15427")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    df = pd.read_parquet(DATA_DIR / "scores" / "molecules.parquet")
    df = df[df["target"] == args.target]
    # A molecule scored in several epochs: keep the mean (Boltz is noisy).
    cols = ["score", "affinity_probability_binary", "affinity_pred_value", "confidence_score", "iptm"]
    df = df.groupby("name", as_index=False).agg({**{c: "mean" for c in cols}, "epoch": "min"})
    df["rxn"] = df["name"].str.split(":").str[1].astype(int)

    with Pool(args.workers) as pool:
        out = pool.map(_resolve, df["name"], chunksize=500)
    df["smiles"] = [s for s, _ in out]
    ok = np.array([x is not None for _, x in out])
    print(f"{ok.sum()}/{len(df)} molecules resolved")

    df = df[ok].reset_index(drop=True)
    X = np.vstack([x for _, x in out if x is not None])
    # Older targets were scored with other formulas: recompute with today's.
    heavy = X[:, FEATURE_NAMES.index("heavy_atoms")]
    df["score"] = [boltz_score(p, v, n) for p, v, n in
                   zip(df["affinity_probability_binary"], df["affinity_pred_value"], heavy)]

    out_dir = DATA_DIR / "dataset"
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_dir / f"{args.target}.parquet", index=False)
    np.save(out_dir / f"{args.target}_X.npy", X)
    print(f"saved {X.shape} to {out_dir}")


if __name__ == "__main__":
    main()
