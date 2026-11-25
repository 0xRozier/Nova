"""A surrogate for the validators' score, trained on what the score actually is.

FINDINGS 77 showed the rank-sum is one measurement in four collinear copies and
that it does not track real binding. So the useful surrogate does not predict
biology -- it predicts `design_iiptm`, from which the rank-sum follows at
rho -0.92. Every earlier surrogate in this repo was aimed at a quantity we
believed meant something; this one is aimed at the oracle itself.

Three things make this honest rather than flattering:

  1. ONE TARGET AT A TIME. A score is a function of the (sequence, target) pair.
     2 132 of our 2 217 scored rows are against Q9NZQ7/8AOK and 85 against
     P05231/4O9H, and exactly ONE sequence was scored against both -- so
     cross-target transfer is not measurable from what we have, and pooling the
     two would be inventing it. Trained on the old target, where the data is.

  2. REPEATS AVERAGED, AND THE NOISE CEILING STATED. 248 sequences were scored
     in more than one run; their within-sequence spread is 0.0121 (median).
     That is the between-run measurement noise, and it caps how well ANY model
     can correlate with a single-run label. Averaging repeats raises the
     ceiling, and the ceiling is printed next to the score so a result cannot
     be read as better than the data allows.

  3. GROUPED SPLITS. Our batches contain arms of 32 sequences sharing a
     framework and differing by one loop. A random split puts near-twins on
     both sides and the model recites the framework instead of predicting the
     loop. Sequences are clustered by 5-mer Jaccard and whole clusters are held
     out, which is the only split that answers the question we care about:
     given a sequence we have never queried, what will the oracle say?

Usage:
    python -m scripts.train_score_surrogate
    python -m scripts.train_score_surrogate --target 4O9H   # the 85, for contrast
                                                            # (the live structure)
"""
from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from nova_lab import REPO_ROOT
from nova_lab.config import NANOBODY_STRUCTURE, PREVIOUS_NANOBODY
from nova_lab.nanobodies import regions

# Which structure each scored file was produced against. Derived from the
# config rather than written out, so a target switch moves this with it instead
# of silently mislabelling two incomparable datasets as one (FINDINGS 76).
CURRENT = NANOBODY_STRUCTURE
PREVIOUS = PREVIOUS_NANOBODY["structure"]

# Which result files were produced against the LIVE target. A registry, not a
# guess: FINDINGS 80 measured the transfer between targets at rho +0.003, so a
# file filed under the wrong target does not degrade the model, it trains it on
# an unrelated function. Every new live-target run must be added here.
LIVE_TARGET_FILES = {
    "scores_nanobody_calibration.csv",   # FINDINGS 77, 85 sequences
    "scores_nanobody_transfer.csv",      # FINDINGS 80, 56 sequences
    "scores_nanobody_learning.csv",      # FINDINGS 83, 70 sequences
}
# Two defensible severities, and they answer different questions.
#
#   0.80  groups sequences that share a framework and differ by one loop
#         (measured: ~0.81 for our keep_fw arm), so a held-out cluster is a
#         framework the model has NEVER queried. This is the conservative read
#         and the default.
#   0.90  is the validator's own diversity rule, under which those same
#         sequences count as distinct (two scrambles of one CDR3 sit at 0.854).
#         It answers the easier question -- choosing a loop inside a framework
#         we have already queried -- and is worth about +0.05 of rho.
#
# Reporting only the second would be flattering: rho runs 0.63 at 0.90 and
# 0.58 at 0.80, so part of the higher number is near-twins on both sides of a
# split. Both are printed.
JACCARD_LINK = 0.80
VALIDATOR_LINK = 0.90
KMER_K = 5
AAS = "ACDEFGHIKLMNPQRSTVWY"
HYDRO = set("AVILMFWCY")
AROM = set("FWY")


def load(target: str) -> pd.DataFrame:
    """One row per distinct sequence, label = mean iiptm over its runs."""
    frames = []
    for path in sorted(glob.glob("results/*.csv")):
        is_live = Path(path).name in LIVE_TARGET_FILES
        if (target == CURRENT) != is_live:
            continue
        try:
            f = pd.read_csv(path)
        except Exception:
            continue
        if {"sequence", "design_iiptm"} <= set(f.columns):
            frames.append(f[["sequence", "design_iiptm"]])
    if not frames:
        raise SystemExit(f"aucun score trouve pour {target}")
    raw = pd.concat(frames, ignore_index=True).dropna()
    raw = raw[raw.sequence.str.fullmatch(f"[{AAS}]+")]
    out = raw.groupby("sequence").design_iiptm.agg(["mean", "std", "count"])
    return out.rename(columns={"mean": "y", "std": "spread", "count": "runs"}).reset_index()


def noise_ceiling(frame: pd.DataFrame) -> tuple[float, float]:
    """(sigma between runs, the best Spearman any model could reach).

    With y = signal + noise, a perfect predictor of the signal correlates with
    the observed label at most sd(signal) / sd(observed). Averaging r repeats
    divides the noise variance by r, which is why `runs` matters.
    """
    repeats = frame[frame.runs > 1]
    sigma = float(repeats.spread.median()) if len(repeats) else float("nan")
    if not np.isfinite(sigma):
        return sigma, float("nan")
    mean_runs = float(frame.runs.mean())
    var_total = float(frame.y.var(ddof=1))
    var_noise = sigma ** 2 / mean_runs
    var_signal = max(var_total - var_noise, 1e-12)
    return sigma, float(np.sqrt(var_signal / var_total))


def descriptors(sequence: str) -> dict:
    r = regions(sequence) or {}
    cdr3 = r.get("CDR3", "")
    d = {"length": len(sequence), "cdr3_length": len(cdr3)}
    for name, s in (("all", sequence), ("cdr3", cdr3)):
        n = max(len(s), 1)
        d[f"{name}_hydro"] = sum(c in HYDRO for c in s) / n
        d[f"{name}_arom"] = sum(c in AROM for c in s) / n
        d[f"{name}_charge"] = (s.count("K") + s.count("R") - s.count("D") - s.count("E")) / n
        d[f"{name}_basic"] = (s.count("K") + s.count("R") + s.count("H")) / n
        d[f"{name}_gly"] = s.count("G") / n
        d[f"{name}_pro"] = s.count("P") / n
        for aa in AAS:
            d[f"{name}_{aa}"] = s.count(aa) / n
    return d


def kmers(sequence: str, k: int = 3) -> set[str]:
    return {sequence[i:i + k] for i in range(len(sequence) - k + 1)}


def cluster(sequences: list[str], link: float = JACCARD_LINK) -> np.ndarray:
    """Single-linkage clusters on 5-mer Jaccard, as union-find.

    The same rule the validator's diversity filter uses (k=5 Jaccard), so a
    held-out cluster is also what the subnet would call a distinct sequence.
    """
    sets = [kmers(s, KMER_K) for s in sequences]
    parent = list(range(len(sequences)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            a, b = sets[i], sets[j]
            if not a or not b:
                continue
            inter = len(a & b)
            if inter / (len(a) + len(b) - inter) >= link:
                ri, rj = find(i), find(j)
                if ri != rj:
                    parent[ri] = rj
    return np.array([find(i) for i in range(len(sequences))])


def grouped_spearman(X: np.ndarray, y: np.ndarray, groups: np.ndarray,
                     model, folds: int = 5, seed: int = 77) -> tuple[float, float]:
    """Out-of-fold Spearman, whole clusters held out. (rho, sd over folds)."""
    from sklearn.base import clone
    rng = np.random.default_rng(seed)
    unique = np.unique(groups)
    rng.shuffle(unique)
    assignment = {g: i % folds for i, g in enumerate(unique)}
    fold_of = np.array([assignment[g] for g in groups])
    predicted = np.full(len(y), np.nan)
    per_fold = []
    for f in range(folds):
        train, test = fold_of != f, fold_of == f
        if test.sum() < 5 or train.sum() < 20:
            continue
        m = clone(model).fit(X[train], y[train])
        predicted[test] = m.predict(X[test])
        if len(np.unique(y[test])) > 2:
            per_fold.append(spearmanr(y[test], predicted[test]).statistic)
    ok = ~np.isnan(predicted)
    return float(spearmanr(y[ok], predicted[ok]).statistic), float(np.std(per_fold))


def build_features(frame: pd.DataFrame) -> tuple[np.ndarray, list[str], np.ndarray, list[str]]:
    desc = pd.DataFrame([descriptors(s) for s in frame.sequence]).fillna(0.0)
    from sklearn.feature_extraction.text import CountVectorizer
    vec = CountVectorizer(analyzer="char", ngram_range=(3, 3), min_df=5)
    counts = vec.fit_transform(frame.sequence).toarray().astype(float)
    counts /= np.maximum(counts.sum(axis=1, keepdims=True), 1)
    return desc.to_numpy(), list(desc.columns), counts, list(vec.get_feature_names_out())


def decision_value(y: np.ndarray, predicted: np.ndarray, bands: dict,
                   pool: int = 30, draws: int = 4000, seed: int = 77) -> list[dict]:
    """What the surrogate is worth where it will be used, not where it looks best.

    A Spearman over the whole set is flattered by dynamic range: our scores run
    0.11 to 0.79 because they pool de novo designs with archive donors, while
    the pool an epoch actually chooses from spans a tenth of that. So the
    surrogate is also scored inside narrow bands, on the only question that
    matters operationally: drawing `pool` candidates at random, how much of the
    gap between a random pick and the oracle's best does picking by surrogate
    recover?
    """
    rng = np.random.default_rng(seed)
    rows = []
    for label, mask in bands.items():
        sel = mask & ~np.isnan(predicted)
        if sel.sum() < 2 * pool:
            continue
        yy, pp = y[sel], predicted[sel]
        got = np.empty((draws, 3))
        for i in range(draws):
            idx = rng.choice(len(yy), size=pool, replace=False)
            got[i] = (yy[idx][np.argmax(pp[idx])], yy[idx].mean(), yy[idx].max())
        chosen, chance, best = got.mean(axis=0)
        rows.append({"bande": label, "n": int(sel.sum()), "sd": float(yy.std()),
                     "rho": float(spearmanr(yy, pp).statistic),
                     "substitut": chosen, "hasard": chance, "oracle": best,
                     "part_captee": (chosen - chance) / (best - chance),
                     "gain_en_sd": (chosen - chance) / float(yy.std())})
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default=PREVIOUS, choices=[PREVIOUS, CURRENT],
                        help=f"{PREVIOUS} carries the data; {CURRENT} is the live "
                             f"target and has only the calibration run")
    parser.add_argument("--folds", type=int, default=5)
    args = parser.parse_args()

    frame = load(args.target)
    sigma, ceiling = noise_ceiling(frame)
    print(f"cible {args.target} : {len(frame)} sequences distinctes, "
          f"{int(frame.runs.sum())} requetes GPU payees")
    print(f"  label = iiptm moyen ; {int((frame.runs > 1).sum())} sequences notees plusieurs fois")
    print(f"  bruit entre runs (median intra-sequence) : {sigma:.4f}")
    print(f"  PLAFOND de Spearman atteignable avec ce bruit : {ceiling:.3f}")
    print(f"  dispersion du label : sd {frame.y.std():.4f}, "
          f"etendue [{frame.y.min():.4f}, {frame.y.max():.4f}]")

    groups = cluster(list(frame.sequence))
    sizes = pd.Series(groups).value_counts()
    print(f"\n{len(sizes)} grappes a Jaccard 5-mer >= {JACCARD_LINK} "
          f"(assez severe pour regrouper une charpente partagee)")
    print(f"  plus grosse grappe : {sizes.iloc[0]} sequences ; "
          f"grappes de taille 1 : {(sizes == 1).sum()}")

    Xd, names_d, Xk, names_k = build_features(frame)
    y = frame.y.to_numpy()
    from sklearn.dummy import DummyRegressor
    from sklearn.ensemble import HistGradientBoostingRegressor
    from sklearn.linear_model import RidgeCV
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    trials = [
        ("temoin (predit la moyenne)", DummyRegressor(), Xd),
        ("ridge sur descripteurs", make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-2, 4, 25))), Xd),
        ("gradient boosting sur descripteurs", HistGradientBoostingRegressor(max_iter=400, learning_rate=0.05, random_state=77), Xd),
        ("ridge sur 3-mers", make_pipeline(StandardScaler(with_mean=False), RidgeCV(alphas=np.logspace(-1, 5, 25))), Xk),
        ("gradient boosting sur 3-mers", HistGradientBoostingRegressor(max_iter=400, learning_rate=0.05, random_state=77), Xk),
        ("ridge sur descripteurs + 3-mers", make_pipeline(StandardScaler(with_mean=False), RidgeCV(alphas=np.logspace(-1, 5, 25))), np.hstack([Xd, Xk])),
        ("gradient boosting sur tout", HistGradientBoostingRegressor(max_iter=500, learning_rate=0.05, random_state=77), np.hstack([Xd, Xk])),
    ]
    print(f"\nvalidation croisee par grappes, {args.folds} plis "
          f"(une grappe entiere est toujours d'un seul cote)")
    print(f"{'modele':<36} {'rho':>7} {'sd plis':>9} {'% du plafond':>13}")
    results = []
    for label, model, X in trials:
        rho, sd = grouped_spearman(X, y, groups, model, args.folds)
        results.append((label, rho, sd))
        share = f"{100 * rho / ceiling:.0f}%" if np.isfinite(ceiling) and rho > 0 else "-"
        print(f"{label:<36} {rho:+7.3f} {sd:9.3f} {share:>13}")

    # The decision-relevant read, with the best model's out-of-fold predictions.
    from sklearn.base import clone
    rng = np.random.default_rng(77)
    uniq = np.unique(groups); rng.shuffle(uniq)
    fold_of = np.array([{g: i % args.folds for i, g in enumerate(uniq)}[g] for g in groups])
    Xall = np.hstack([Xd, Xk])
    predicted = np.full(len(y), np.nan)
    for f in range(args.folds):
        tr, te = fold_of != f, fold_of == f
        if te.sum() < 5 or tr.sum() < 20:
            continue
        m = clone(HistGradientBoostingRegressor(max_iter=500, learning_rate=0.05,
                                                random_state=77)).fit(Xall[tr], y[tr])
        predicted[te] = m.predict(Xall[te])
    bands = {"tout le jeu": np.ones(len(y), bool),
             "bande realiste 0,15-0,45": (y >= 0.15) & (y <= 0.45),
             "haut du vivier 0,55-0,80": (y >= 0.55) & (y <= 0.80)}
    rows = decision_value(y, predicted, bands)
    if rows:
        print("\nce que le substitut vaut la ou il servira : choisir 1 parmi 30")
        d = pd.DataFrame(rows)
        d["part_captee"] = (100 * d.part_captee).round(0).astype(int).astype(str) + "%"
        d["gain_en_sd"] = d.gain_en_sd.round(2)
        print(d.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
        print("  part_captee : fraction de l'ecart hasard -> oracle que le substitut recupere")
        print("  gain_en_sd  : le meme gain en ecarts-types du vivier, comparable au GPU")

    best = max(results, key=lambda r: r[1])
    print(f"\nmeilleur : {best[0]} a rho {best[1]:+.3f} "
          f"(plafond {ceiling:.3f}, temoin {results[0][1]:+.3f})")
    loose = cluster(list(frame.sequence), link=VALIDATOR_LINK)
    rho_loose, _ = grouped_spearman(
        Xall, y, loose,
        HistGradientBoostingRegressor(max_iter=500, learning_rate=0.05, random_state=77),
        args.folds)
    print(f"  charpente jamais vue (grappes a {JACCARD_LINK}) : rho {best[1]:+.3f}")
    print(f"  boucle choisie dans une charpente deja interrogee "
          f"(grappes a {VALIDATOR_LINK}) : rho {rho_loose:+.3f}")
    out = REPO_ROOT / "results" / f"surrogate_{args.target}.csv"
    pd.DataFrame(results, columns=["modele", "rho", "sd_plis"]).to_csv(out, index=False)
    print(f"-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
