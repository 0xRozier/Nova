"""Predict design_iiptm from the sequence, to select candidates without a GPU.

Why this exists, and what it replaces.

Our generator already produces competitive candidates: over 84 grafted
candidates scored by BoltzGen, the best read design_iiptm 0.692, above the
field's median of 0.664. But the eight we actually send to BoltzGen are chosen
on PSH -- and FINDINGS 3ter measured that PSH does not predict the validator's
verdict below 96. So we were discarding the good candidates before looking at
them: the eight of the epoch-26238 dry run topped out at 0.375 while the pool
they came from reached 0.692.

Selecting on the real metric is out of reach: BoltzGen costs 150 s per sequence
on 2x T4, so 84 candidates is 3.5 hours against a 67-minute epoch. Predicting
it costs milliseconds, and there are now 5956 labelled sequences -- the whole
scored field plus our own scored candidates.

What it is worth, measured three ways, because the first number was too good to
take at face value:

  random 5-fold over everything       rho = 0.913   <- inflated
  5-fold grouped BY EPOCH             rho = 0.882   near-duplicates across
                                                    epochs explain only 0.03
  trained on the FIELD, ranking OUR   rho = 0.569   the task that matters
  candidates only                     (p = 6e-12, n = 124)

The third is the one to quote. A model that merely told November field
sequences (0.66) from grafts (0.34) would score well on the first two and be
useless for ranking our own pool, which is the mistake this project has made
three times in other guises -- a real signal that explains none of the part of
the result you care about.

The gain on selection, on our 124 scored candidates:

  our 8 best by the surrogate   real design_iiptm median 0.572, max 0.649
  8 at random                                       median 0.332
  the pool                                          median 0.335, max 0.667

So the surrogate finds nearly the best candidate the pool contains. What it
does NOT do is close the gap: 0.649 still sits under the field's per-epoch best
of 0.757. This takes us from the bottom 12% of the field to roughly its median.
Necessary, not sufficient.

Usage:
    python -m scripts.train_iiptm_surrogate              # train, report, save
    python -m scripts.train_iiptm_surrogate --no-save    # report only
"""
import argparse
import glob
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from nova_lab import DATA_DIR, REPO_ROOT
from nova_lab.config import NANOBODY_TARGET
from nova_lab.nanobodies import regions_batch
from nova_lab.nanobody_features import featurize

TARGET_METRIC = "design_iiptm"


def labelled(target: str) -> pd.DataFrame:
    """Every sequence we know a design_iiptm for, field and ours."""
    field = pd.read_csv(DATA_DIR / "submissions" / f"{target}_nanobody_results.csv")
    field = field[["sequence", TARGET_METRIC, "epoch"]].dropna()
    ours = [pd.read_csv(f)[["sequence", TARGET_METRIC]]
            for f in sorted(glob.glob(str(REPO_ROOT / "results" /
                                          "scores_nanobody_observation_*.csv")))]
    ours = (pd.concat(ours).dropna().assign(epoch=-1) if ours
            else pd.DataFrame(columns=["sequence", TARGET_METRIC, "epoch"]))
    both = pd.concat([field.assign(source="field"), ours.assign(source="ours")])
    return both.drop_duplicates("sequence").reset_index(drop=True)


def features(sequences: list[str]) -> tuple[np.ndarray, list[int]]:
    """Feature matrix, and which input rows made it in."""
    regions = regions_batch(sequences)
    rows, kept = [], []
    for i, (sequence, region) in enumerate(zip(sequences, regions)):
        if not region:
            continue
        vector = featurize(sequence, region)
        if vector is None:
            continue
        rows.append(vector)
        kept.append(i)
    return np.vstack(rows), kept


def model():
    from sklearn.ensemble import HistGradientBoostingRegressor
    return HistGradientBoostingRegressor(max_iter=400, learning_rate=0.06,
                                         random_state=0)


def report(X: np.ndarray, y: np.ndarray, source: np.ndarray,
           epoch: np.ndarray) -> None:
    """The three numbers, in order of how much they should be trusted."""
    from sklearn.model_selection import GroupKFold, KFold
    from scipy import stats

    predicted = np.zeros_like(y)
    for train, test in KFold(5, shuffle=True, random_state=0).split(X):
        predicted[test] = model().fit(X[train], y[train]).predict(X[test])
    print(f"  random 5-fold over everything      rho = "
          f"{stats.spearmanr(predicted, y).statistic:.3f}  (inflated)")

    is_field = source == "field"
    Xf, yf, ef = X[is_field], y[is_field], epoch[is_field]
    predicted = np.zeros_like(yf)
    for train, test in GroupKFold(5).split(Xf, yf, groups=ef):
        predicted[test] = model().fit(Xf[train], yf[train]).predict(Xf[test])
    print(f"  5-fold grouped by epoch            rho = "
          f"{stats.spearmanr(predicted, yf).statistic:.3f}")

    is_ours = source == "ours"
    if is_ours.sum() < 20:
        print("  (too few scored candidates of our own to run the third test)")
        return
    fitted = model().fit(Xf, yf)
    p, yt = fitted.predict(X[is_ours]), y[is_ours]
    rho = stats.spearmanr(p, yt)
    print(f"  field -> ranking OUR candidates    rho = {rho.statistic:.3f}  "
          f"(p = {rho.pvalue:.1g}, n = {is_ours.sum()})  <- the one to quote")
    best8 = np.argsort(-p)[:8]
    rng = np.random.default_rng(0)
    print(f"\n  selecting 8 of our own candidates:")
    print(f"    by the surrogate   median {np.median(yt[best8]):.3f}, "
          f"max {yt[best8].max():.3f}")
    print(f"    at random          median {np.median(rng.choice(yt, 8)):.3f}")
    print(f"    the whole pool     median {np.median(yt):.3f}, max {yt.max():.3f}")
    print(f"\n  field per-epoch best is 0.757: the surrogate's pick closes much of")
    print(f"  the gap to the pool's own ceiling, not the gap to the field.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default=NANOBODY_TARGET[0])
    parser.add_argument("--no-save", action="store_true")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    data = labelled(args.target)
    X, kept = features(data["sequence"].tolist())
    data = data.iloc[kept].reset_index(drop=True)
    y = data[TARGET_METRIC].to_numpy()
    print(f"{len(y)} labelled sequences, {X.shape[1]} features "
          f"({(data['source'] == 'ours').sum()} of our own)\n")
    report(X, y, data["source"].to_numpy(), data["epoch"].to_numpy())

    if args.no_save:
        return
    # Fit on everything for use: the held-out numbers above are what the model
    # is worth, and there is no reason to withhold data from the shipped one.
    fitted = model().fit(X, y)
    out = args.out or DATA_DIR / "dataset" / f"{args.target}_iiptm_surrogate.joblib"
    out.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump({"model": fitted, "metric": TARGET_METRIC, "n_train": len(y)}, out)
    print(f"\nsaved to {out}")


if __name__ == "__main__":
    main()
