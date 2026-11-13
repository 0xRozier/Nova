"""Can a CPU model predict a nanobody's BoltzGen metrics from its sequence?

One LightGBM model per metric, then two honest evaluations:

  1. **Random split** (optimistic): the test sequences have close relatives in
     training, as they would when we mutate around known-good sequences.
  2. **Temporal split** (strict): train on the epochs before a cutoff, test on
     the ones after, then reproduce the validator's ranking per epoch and ask
     where our best-predicted sequence would actually have ranked.

The second is the number that matters: it says whether picking by prediction
would have won epochs. The molecule surrogate looked fine on a backtest and
collapsed out of distribution, so believe the strict split only.

Usage:
    python -m scripts.train_nanobody_surrogate [--target Q9NZQ7]
"""
import argparse

import lightgbm as lgb
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from nova_lab import DATA_DIR
from nova_lab.nanobody_features import feature_names
from nova_lab.surrogate import Surrogate

# config/boltzgen_config.yaml: metric -> rank direction ("min" = lower is better)
METRIC_MODES = {
    "design_iiptm": "max", "design_ptm": "max", "design_to_target_iptm": "max",
    "min_design_to_target_pae": "min", "interaction_pae": "min",
    "plip_hbonds_refolded": "max", "plip_saltbridge_refolded": "max",
    "delta_sasa_refolded": "max", "liability_score": "min", "liability_num_violations": "min",
}
PARAMS = dict(objective="regression", learning_rate=0.05, num_leaves=63,
              min_child_samples=20, feature_fraction=0.5, bagging_fraction=0.8,
              bagging_freq=1, lambda_l2=1.0, verbose=-1)


def train_all(X, df, train_idx, names, rounds=400) -> dict:
    models = {}
    for metric in METRIC_MODES:
        y = df[metric].to_numpy()
        models[metric] = lgb.train(PARAMS, lgb.Dataset(X[train_idx], y[train_idx], feature_name=names), rounds)
    return models


def predict_all(models, X) -> pd.DataFrame:
    return pd.DataFrame({m: model.predict(X) for m, model in models.items()})


def rank_sum(frame: pd.DataFrame) -> pd.Series:
    """The validator's score: dense rank per metric, summed. Lower wins."""
    ranks = pd.DataFrame({m: frame[m].rank(method="dense", ascending=(mode == "min"))
                          for m, mode in METRIC_MODES.items()})
    return ranks.sum(axis=1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default="Q9NZQ7")
    parser.add_argument("--test-epochs", type=int, default=15)
    args = parser.parse_args()

    ds = DATA_DIR / "dataset"
    df = pd.read_parquet(ds / f"{args.target}_nanobodies.parquet")
    X = np.load(ds / f"{args.target}_nanobodies_X.npy")
    names = feature_names()
    print(f"{len(df)} sequences, {X.shape[1]} features, "
          f"{df['epoch'].notna().sum()} with an epoch\n")

    # ---- 1. random split: how well do metrics follow from the sequence at all?
    rng = np.random.default_rng(0)
    perm = rng.permutation(len(df))
    cut = int(0.8 * len(df))
    models = train_all(X, df, perm[:cut], names)
    pred = predict_all(models, X[perm[cut:]])
    print("random split (optimistic) — Spearman per metric:")
    for metric in METRIC_MODES:
        rho = spearmanr(pred[metric], df[metric].to_numpy()[perm[cut:]]).statistic
        print(f"  {metric:28s} {rho:5.2f}")
    truth = df.iloc[perm[cut:]].reset_index(drop=True)
    print(f"  {'rank_sum':28s} {spearmanr(rank_sum(pred), rank_sum(truth)).statistic:5.2f}")

    # ---- the model the miner loads: trained on everything there is.
    # This used to be the backtest's model, trained only on what preceded the
    # held-out epochs, and it was saved after the backtest -- so a machine with
    # fewer than --test-epochs epochs took the early return and ended up with
    # NO model at all, which is what a fresh install looks like. The deployed
    # model also has no reason to discard the fifteen most recent epochs, which
    # are the donors the competition is actually at.
    Surrogate(main=train_all(X, df, np.arange(len(df)), names)).save(
        f"{args.target}_nanobody")
    print(f"\nsaved per-metric models trained on all {len(df)} sequences to {ds}")

    # ---- 2. temporal split: would picking by prediction have won epochs?
    known = df[df["epoch"].notna()].copy()
    epochs = np.sort(known["epoch"].unique())
    if len(epochs) <= args.test_epochs:
        print("\nnot enough epochs for a temporal split (needs more than "
              f"{args.test_epochs}; {len(epochs)} known). Run "
              "scripts.fetch_submissions --last 60 then scripts.analyze_nanobodies "
              "to label more, then rebuild. The model above is already saved.")
        return
    split = epochs[-args.test_epochs]
    train_mask = (df["epoch"] < split) | df["epoch"].isna()
    # Sequences without an epoch come from earlier rounds of this same target,
    # so they are safe to train on; anything from a test epoch is held out.
    train_idx = np.flatnonzero(train_mask & ~df["hash"].isin(known[known["epoch"] >= split]["hash"]))
    print(f"\ntemporal split (strict): train on {len(train_idx)} sequences from before epoch {split:.0f}")

    models = train_all(X, df, train_idx, names)
    rows = []
    for epoch, g in known[known["epoch"] >= split].groupby("epoch"):
        idx = g.index.to_numpy()
        if len(idx) < 5:
            continue
        true_ranks = rank_sum(g.reset_index(drop=True))
        pred_ranks = rank_sum(predict_all(models, X[idx]))
        our_pick = int(np.argmin(pred_ranks.to_numpy()))
        order = np.argsort(true_ranks.to_numpy())  # best first
        rows.append({
            "epoch": int(epoch), "n": len(idx),
            "our_pick_true_place": int(np.where(order == our_pick)[0][0]) + 1,
            "best_place": 1, "spearman": spearmanr(pred_ranks, true_ranks).statistic,
        })
    bt = pd.DataFrame(rows).set_index("epoch")
    with pd.option_context("display.width", 120, "display.max_rows", 100):
        print(bt.round(3))
    print(f"\nmean Spearman within epoch: {bt['spearman'].mean():.2f}")
    print(f"our pick's true place: median {bt['our_pick_true_place'].median():.0f} "
          f"of {bt['n'].median():.0f} sequences; "
          f"won {(bt['our_pick_true_place'] == 1).sum()}/{len(bt)} epochs")
    print(f"a random pick would average place {(bt['n'] + 1).mean() / 2:.1f}")
    print("  (the backtest's model is not the one on disk: that one is trained "
          "on every sequence, including these test epochs)")


if __name__ == "__main__":
    main()
