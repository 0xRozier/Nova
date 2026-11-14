"""Build the nanobody training set: archived sequences, metrics and features.

Sources:
  - the public archive (every scored sequence for the target, with its
    10 BoltzGen metrics) -> the bulk of the training data, no epoch attached;
  - data/scores/nanobodies.jsonl (scripts.analyze_nanobodies) -> the epoch each
    submission was scored in, needed for a temporal backtest.

Output: data/dataset/<target>_nanobodies.parquet and _X.npy (same row order).

Usage:
    python -m scripts.build_nanobody_dataset [--target Q9NZQ7]
"""
import argparse
import json

import numpy as np
import pandas as pd

from nova_lab import DATA_DIR
from nova_lab.nanobodies import NanobodyArchive, germline_identity, regions, seq_hash
from nova_lab.nanobody_features import featurize

METRICS = ["design_iiptm", "design_ptm", "design_to_target_iptm", "min_design_to_target_pae",
           "interaction_pae", "plip_hbonds_refolded", "plip_saltbridge_refolded",
           "delta_sasa_refolded", "liability_score", "liability_num_violations"]


def scored_but_unarchived(target: str, archived: pd.DataFrame, archive) -> pd.DataFrame:
    """Submissions the validators scored that the HuggingFace archive has not got yet.

    The archive lags, and the generator's donor pool was built from it alone --
    so the freshest and best donors were invisible. Measured on 2026-11-14:
    196 of the 230 scored submissions of epochs 26299-26306 were missing, at a
    median design_iiptm of 0.694 and a maximum of 0.779, and they sat at
    exactly the CDR3 lengths where the pool had nothing. Merging them takes the
    lengths with a donor at the field's level (0.72) from **2 of 19 to 7 of
    19** -- from 11 % of epochs to 37 %, on eight epochs of catch-up alone.

    This is why it matters: the epoch fixes the CDR3 length, a donor of the
    wrong length is no use, and the archive's good binders are almost all at
    length 13. The leaders' own winning sequences, at every length, go through
    this path and nowhere else.

    Needs scripts.fetch_submissions then scripts.analyze_nanobodies to have run.
    """
    path = DATA_DIR / "submissions" / f"{target}_nanobody_results.csv"
    if not path.exists():
        print(f"no {path.name}: the donor pool is the HuggingFace archive alone, "
              f"which lags. Run scripts.fetch_submissions then "
              f"scripts.analyze_nanobodies to add what the validators scored since.")
        return pd.DataFrame()
    scored = pd.read_csv(path)
    fresh = scored[~scored["sequence"].isin(set(archived["sequence"]))].copy()
    missing = [m for m in METRICS if m not in fresh.columns]
    if missing:
        print(f"{path.name} lacks {missing}; not merging")
        return pd.DataFrame()
    fresh = fresh.dropna(subset=METRICS)
    if fresh.empty:
        print("the scored submissions are all in the archive already")
        return pd.DataFrame()
    out = pd.DataFrame({
        "sequence": fresh["sequence"].to_numpy(),
        "sequence_hash": fresh["hash"].to_numpy(),
        "epoch": fresh["epoch"].to_numpy(),
        **{m: fresh[m].to_numpy() for m in METRICS},
    })
    out["regions"] = [archive.regions_for(s) if hasattr(archive, "regions_for")
                      else regions(s) for s in out["sequence"]]
    numberable = sum(r is not None for r in out["regions"])
    print(f"+{len(out)} scored submissions the archive has not got yet "
          f"({numberable} numberable), design_iiptm median "
          f"{fresh['design_iiptm'].median():.3f}, max {fresh['design_iiptm'].max():.3f}")
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default="Q9NZQ7")
    parser.add_argument("--no-refresh", action="store_true",
                        help="rebuild from the archive already on disk instead of "
                             "re-downloading it (for a replay at a fixed snapshot)")
    args = parser.parse_args()

    # Re-download by default. NanobodyArchive only fetches when the file is
    # absent, so rebuilding from the cached copy reproduces the stale pool it
    # was run to fix -- and reports the same lag, which is how the staleness
    # survived a "rebuild" on 2026-11-09: the archive on disk was three days
    # old and the epoch labels stopped 63 epochs back.
    archive = NanobodyArchive(args.target, refresh=not args.no_refresh)
    df = pd.read_csv(DATA_DIR / "archive" / f"{args.target}_nanobodies.csv")
    df["regions"] = archive.regions
    print(f"{len(df)} archived sequences, {sum(r is not None for r in archive.regions)} numberable")

    # epoch per sequence hash, where we know it
    epochs = {}
    cache = DATA_DIR / "scores" / "nanobodies.jsonl"
    if cache.exists():
        for line in cache.open():
            row = json.loads(line)
            if row["metrics"]:
                epochs.setdefault(row["hash"], row["epoch"])
    df["epoch"] = df["sequence_hash"].map(epochs)
    print(f"{df['epoch'].notna().sum()} sequences with a known epoch "
          f"({df['epoch'].min():.0f}-{df['epoch'].max():.0f})")
    df = pd.concat([df, scored_but_unarchived(args.target, df, archive)],
                   ignore_index=True)

    rows, features = [], []
    for row in df.itertuples():
        if row.regions is None:
            continue
        germline = germline_identity(row.sequence)
        x = featurize(row.sequence, row.regions, germline)
        if x is None:
            continue
        rows.append({"sequence": row.sequence, "hash": row.sequence_hash, "epoch": row.epoch,
                     "v_genes": ",".join(germline["V"]["genes"]) if germline else "",
                     "v_identity": germline["V"]["identity_pct"] if germline else np.nan,
                     "j_identity": germline["J"]["identity_pct"] if germline else np.nan,
                     **{m: getattr(row, m) for m in METRICS},
                     **{f"cdr{i}": row.regions[f"CDR{i}"] for i in (1, 2, 3)}})
        features.append(x)
        if len(rows) % 1000 == 0:
            print(f"  featurized {len(rows)}")

    out = pd.DataFrame(rows)
    X = np.vstack(features)
    ok = np.isfinite(X).all(axis=1)
    out, X = out[ok].reset_index(drop=True), X[ok]

    out_dir = DATA_DIR / "dataset"
    out_dir.mkdir(parents=True, exist_ok=True)
    out.to_parquet(out_dir / f"{args.target}_nanobodies.parquet", index=False)
    np.save(out_dir / f"{args.target}_nanobodies_X.npy", X)
    print(f"saved {X.shape} to {out_dir}")


if __name__ == "__main__":
    main()
