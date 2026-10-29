"""Train and backtest a CPU surrogate of the validator's Boltz-2 molecule score.

Temporal split: the model only sees epochs before the test window, like it
would live. Reported:
  - Spearman correlation on the test molecules (overall and within epoch);
  - a selection backtest: in each test epoch, pick 20 of the scored molecules
    by predicted score and compare their true total with the best possible
    20 (oracle) and with a random 20.

Usage:
    python -m scripts.train_surrogate [--target O15427] [--test-epochs 30]
    python -m scripts.train_surrogate --aux P40261   # add a past target's model as a feature
"""
import argparse

import lightgbm as lgb
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from nova_lab import DATA_DIR
from nova_lab.config import NUM_MOLECULES
from nova_lab.features import FEATURE_NAMES
from nova_lab.surrogate import Surrogate

PARAMS = dict(
    objective="regression",
    learning_rate=0.05,
    num_leaves=127,
    min_child_samples=20,
    feature_fraction=0.3,
    bagging_fraction=0.8,
    bagging_freq=1,
    lambda_l2=1.0,
    verbose=-1,
)


def backtest(test: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    for epoch, g in test.groupby("epoch"):
        if len(g) < 2 * NUM_MOLECULES:
            continue
        rows.append({
            "epoch": epoch,
            "pool": len(g),
            "model": g.nlargest(NUM_MOLECULES, "pred")["score"].sum(),
            "oracle": g.nlargest(NUM_MOLECULES, "score")["score"].sum(),
            "random": g.sample(NUM_MOLECULES, random_state=int(rng.integers(1 << 31)))["score"].sum(),
            "spearman": spearmanr(g["pred"], g["score"]).statistic,
        })
    return pd.DataFrame(rows).set_index("epoch")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default="O15427")
    parser.add_argument("--test-epochs", type=int, default=30)
    parser.add_argument("--rounds", type=int, default=1500)
    parser.add_argument("--aux", nargs="*", default=[], help="past targets built with scripts.build_dataset")
    args = parser.parse_args()

    ds_dir = DATA_DIR / "dataset"
    aux_models = {}
    for target in args.aux:
        aux_df = pd.read_parquet(ds_dir / f"{target}.parquet")
        aux_X = np.load(ds_dir / f"{target}_X.npy", mmap_mode="r")
        print(f"training auxiliary model on {target} ({len(aux_df)} mols)")
        aux_models[target] = lgb.train(PARAMS, lgb.Dataset(aux_X[:], aux_df["score"].to_numpy()), 600)
    bundle = Surrogate(main=None, aux=aux_models)
    feature_names = FEATURE_NAMES + [f"aux_{t}" for t in aux_models]

    df = pd.read_parquet(ds_dir / f"{args.target}.parquet")
    X = bundle.augment(np.load(ds_dir / f"{args.target}_X.npy"))

    epochs = np.sort(df["epoch"].unique())
    split = epochs[-args.test_epochs]
    train_idx = np.flatnonzero(df["epoch"] < split)
    test_idx = np.flatnonzero(df["epoch"] >= split)
    print(f"train {len(train_idx)} mols (epochs < {split}), test {len(test_idx)} mols")

    # Last 10% of the training epochs serve as early-stopping validation.
    val_split = np.quantile(df["epoch"].iloc[train_idx], 0.9)
    fit_idx = train_idx[df["epoch"].iloc[train_idx] < val_split]
    val_idx = train_idx[df["epoch"].iloc[train_idx] >= val_split]
    y = df["score"].to_numpy()

    train_set = lgb.Dataset(X[fit_idx], y[fit_idx], feature_name=feature_names, free_raw_data=True)
    val_set = lgb.Dataset(X[val_idx], y[val_idx], reference=train_set)
    model = lgb.train(PARAMS, train_set, args.rounds, valid_sets=[val_set],
                      callbacks=[lgb.early_stopping(100), lgb.log_evaluation(200)])

    test = df.iloc[test_idx].copy()
    test["pred"] = model.predict(X[test_idx], num_iteration=model.best_iteration)
    print(f"\nSpearman (all test molecules): {spearmanr(test['pred'], test['score']).statistic:.3f}")

    bt = backtest(test, np.random.default_rng(0))
    with pd.option_context("display.float_format", "{:.3f}".format, "display.width", 120):
        print(bt)
        print("\nmean:", bt.mean().round(3).to_dict())
    gap = (bt["model"] - bt["random"]) / (bt["oracle"] - bt["random"])
    print(f"model closes {gap.mean():.0%} of the random->oracle gap on average")

    importance = pd.Series(model.feature_importance("gain"), index=feature_names).nlargest(15)
    print("\ntop features:\n", (importance / importance.sum()).round(3).to_string())

    # Final model on everything, for use by the candidate generator.
    bundle.main = lgb.train(PARAMS, lgb.Dataset(X, y, feature_name=feature_names), model.best_iteration)
    bundle.save(args.target)
    print(f"\nsaved final model ({model.best_iteration} rounds) to {ds_dir}")


if __name__ == "__main__":
    main()
