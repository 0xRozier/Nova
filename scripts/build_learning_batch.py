"""The 70 queries that take the live target's surrogate to its plateau.

FINDINGS 80 bis measured the learning curve on the previous target: 200 queries
reach rho 0.42, three quarters of what 1 157 reach, and the return collapses
after. The live target has 141 queries and sits at rho 0.235 +/- 0.14. About 70
more put it at the plateau, for 2.5 h of GPU.

Which 70 is not a detail. The live set is lopsided: 39 of its 132 sequences --
30% -- share the single camelised 61H7 framework, because the calibration batch
was built to answer a question about that framework (FINDINGS 74). A surrogate
trained on it extrapolates everywhere else. So:

  cover   45 sequences chosen by farthest-point sampling in feature space,
          seeded with the 132 already queried. Greedy max-min distance: each
          pick is the pool sequence furthest from everything chosen so far.
          This buys coverage, which is what reduces extrapolation; it is not a
          random sample, and a random sample would re-draw the archive's own
          lopsidedness instead of correcting ours.
  exploit 25 sequences the current surrogate predicts highest. These serve
          twice: they populate the top of the range, where every operational
          decision is taken, and they are a PROSPECTIVE test of the surrogate
          -- the only kind cross-validation cannot fake.

Drawn from the 17 634 archive sequences never queried against the live target.
The archive is the previous target's submissions, which is the right pool
anyway: it is what our generator selects from, so the surrogate should be
trained on the distribution it will be asked to rank.

READINGS FIXED BEFORE THE RUN:

  1. Retrained on 202 queries, rho should approach the plateau the previous
     target showed at 200 (+0.42 +/- 0.10). If it stays near 0.235, the live
     target is harder than the previous one and the learning curve does not
     transfer either -- which would be worth knowing before buying GPU.
  2. The 25 `exploit` picks should score above the median of the 45 `cover`
     picks. At rho 0.235 the expected lift is small but measurable, and this is
     a real held-out test of whether the surrogate is usable at all. If the
     exploit arm does NOT beat the cover arm, the surrogate's cross-validated
     rho is not translating into picks, and the whole plan in FINDINGS 80 bis
     needs rethinking.
  3. Watch for range: if the 70 come back inside 0.12-0.39 like the existing
     132, the live target's score range really is this narrow, and the 0.0133
     noise is then a fifth of the usable spread rather than a thirtieth.

Usage:
    python -m scripts.build_learning_batch --out inputs/kaggle_learning.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from nova_lab import DATA_DIR, REPO_ROOT
from nova_lab.config import NANOBODY_STRUCTURE, PREVIOUS_NANOBODY
from nova_lab.nanobodies import check_sequence
from scripts.train_score_surrogate import build_features, load

ARCHIVE = DATA_DIR / "archive" / f"{PREVIOUS_NANOBODY['target']}_regions.csv"
N_COVER = 45
N_EXPLOIT = 25
POOL = 4000          # subsample of the archive, kept for a tractable distance pass
COMPONENTS = 40      # PCA dimensions the distances are taken in
SEED = 81
MINUTES_PER_SEQUENCE = 2.2


def candidate_pool(exclude: set[str], seed: int = SEED,
                   pool: int = POOL) -> pd.DataFrame:
    """Archive sequences never queried against the live target.

    `pool` subsamples because the farthest-point pass is O(pool x known) in a
    40-dimensional space; 4 000 is plenty to find the far corners of 17 634 and
    keeps the tests from costing minutes.
    """
    if not ARCHIVE.exists():
        raise SystemExit(f"{ARCHIVE} absent")
    arch = pd.read_csv(ARCHIVE).dropna(subset=["CDR3"]).drop_duplicates("sequence")
    fresh = arch[~arch.sequence.isin(exclude)][["sequence"]].reset_index(drop=True)
    if len(fresh) > pool:
        fresh = fresh.sample(pool, random_state=seed).reset_index(drop=True)
    return fresh


def project(sequences: list[str], components: int = COMPONENTS) -> np.ndarray:
    """Features reduced to a space where euclidean distance is meaningful.

    The raw features mix counts, fractions and lengths on wildly different
    scales, so an unstandardised distance would be dominated by `length` alone
    and 'coverage' would mean 'a spread of lengths'.
    """
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    desc, _, kmers, _ = build_features(pd.DataFrame({"sequence": sequences}))
    X = StandardScaler().fit_transform(np.hstack([desc, kmers]))
    n = min(components, X.shape[0], X.shape[1])
    return PCA(n_components=n, random_state=SEED).fit_transform(X)


def farthest_points(space: np.ndarray, n_seeds: int, take: int) -> list[int]:
    """Greedy max-min: each pick is furthest from everything already held.

    Indices below `n_seeds` are the already-queried sequences and are never
    returned; they only push the picks away from what we have.
    """
    distance = np.linalg.norm(space[n_seeds:, None, :] - space[None, :n_seeds, :],
                              axis=2).min(axis=1)
    chosen: list[int] = []
    for _ in range(min(take, len(distance))):
        pick = int(np.argmax(distance))
        chosen.append(pick)
        extra = np.linalg.norm(space[n_seeds:] - space[n_seeds + pick], axis=1)
        distance = np.minimum(distance, extra)
    return chosen


def surrogate_predictions(known: pd.DataFrame, pool: list[str]) -> np.ndarray:
    """The current live-target surrogate, applied to the pool."""
    from sklearn.ensemble import HistGradientBoostingRegressor
    desc_k, _, kmer_k, _ = build_features(known[["sequence"]])
    both = pd.DataFrame({"sequence": list(known.sequence) + list(pool)})
    desc_a, _, kmer_a, _ = build_features(both)
    X = np.hstack([desc_a, kmer_a])
    split = len(known)
    model = HistGradientBoostingRegressor(max_iter=500, learning_rate=0.05,
                                          random_state=SEED)
    model.fit(X[:split], known.y.to_numpy())
    return model.predict(X[split:])


def build(seed: int = SEED, pool_size: int = POOL,
          n_cover: int = N_COVER, n_exploit: int = N_EXPLOIT) -> pd.DataFrame:
    known = load(NANOBODY_STRUCTURE)
    pool = candidate_pool(set(known.sequence), seed, pool_size)
    sequences = list(known.sequence) + list(pool.sequence)
    space = project(sequences)

    cover_idx = farthest_points(space, len(known), n_cover)
    predicted = surrogate_predictions(known, list(pool.sequence))
    order = np.argsort(-predicted)
    exploit_idx: list[int] = []
    for i in order:
        if len(exploit_idx) == n_exploit:
            break
        if i not in cover_idx:
            exploit_idx.append(int(i))

    rows = []
    for rank, i in enumerate(cover_idx, 1):
        rows.append({"sequence": pool.sequence[i], "arm": "cover",
                     "label": f"cover{rank:02d}",
                     "note": f"couverture, prediction {predicted[i]:.4f}",
                     "predicted": float(predicted[i])})
    for rank, i in enumerate(exploit_idx, 1):
        rows.append({"sequence": pool.sequence[i], "arm": "exploit",
                     "label": f"exploit{rank:02d}",
                     "note": f"predit haut, prediction {predicted[i]:.4f}",
                     "predicted": float(predicted[i])})
    frame = pd.DataFrame(rows).drop_duplicates("sequence").reset_index(drop=True)
    frame["group"] = "candidate"
    frame["epoch"] = pd.NA
    frame["true_rank_sum"] = pd.NA
    # Recorded, never used to drop a row: this batch trains a surrogate on the
    # oracle, so a sequence the validator would refuse is a valid data point.
    frame["layer1"] = [", ".join(check_sequence(s)) or "ok" for s in frame.sequence]
    frame["depth"] = frame.groupby("arm").cumcount()
    return frame.sort_values(["depth", "arm"]).drop(columns="depth").reset_index(drop=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path,
                        default=REPO_ROOT / "inputs" / "kaggle_learning.csv")
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()

    frame = build(args.seed)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.out, index=False)

    known = load(NANOBODY_STRUCTURE)
    print(f"{len(frame)} sequences -> {args.out}")
    for arm, group in frame.groupby("arm"):
        p = group.predicted
        print(f"  {arm:<8} {len(group):3d}  prediction {p.min():.4f} a {p.max():.4f} "
              f"(mediane {p.median():.4f})")
    print(f"  ok a la couche 1 : {(frame.layer1 == 'ok').sum()}/{len(frame)}")
    print(f"  deja interrogees sur {NANOBODY_STRUCTURE} : {len(known)} "
          f"-> {len(known) + len(frame)} apres ce lot")
    print(f"  {len(frame) * MINUTES_PER_SEQUENCE / 60:.1f} h sur 2x T4 "
          f"a {MINUTES_PER_SEQUENCE} min/sequence")
    print("  les trois lectures sont fixees dans le docstring, avant le run")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
