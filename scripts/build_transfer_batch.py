"""The batch that says whether 1 490 paid queries are reusable.

FINDINGS 78 trained a surrogate on the previous target (1 490 queries against
8AOK) that reaches rho 0.580 on unseen frameworks. It has to serve on the live
target, 4O9H -- and exactly ONE sequence in the repo was ever scored against
both, so the transfer is unmeasured. Pooling the two datasets without measuring
it would be inventing the thing that makes the surrogate worth anything.

This batch measures that, and nothing else:

  transfer  48 sequences already scored against the previous target, scored now
            against the live one. Stratified over the full previous-target
            range in eight bands, six per band, and drawn only from the 248
            sequences with more than one run, so the x-axis of the comparison
            carries less noise than a single query would.
  repeat     8 sequences from the calibration run, scored a second time. We
            have never measured the between-run noise ON THE LIVE TARGET;
            without it a weak correlation is ambiguous between "the transfer
            fails" and "the labels are noisy". These eight separate the two.

Stratification is the point, not a nicety: the surrogate's value was measured
inside a band of sd 0.055, and a sample drawn from one band would attenuate the
correlation by range restriction and understate a transfer that is real.

What the result decides, fixed here before the run:

  rho >= 0.6   the queries carry over. The surrogate is usable on the live
               target now, and the GPU spent on the previous one was not wasted.
  0.3 to 0.6   partial transfer: usable as a prior, retrained as live queries
               accumulate.
  rho < 0.3    the score is target-specific. The surrogate must be rebuilt for
               each target, and the 1 490 queries are history, not an asset.

And one confound to settle afterwards from the same data, with no extra
sequences: `design_ptm` is the binder's own fold confidence and does not
involve the target at all. If the transfer correlation survives partialling
`design_ptm` out, it is genuinely about the interface; if it does not, we are
only rediscovering that some sequences fold confidently anywhere, which would
not help us pick a winner for a given epoch.

Usage:
    python -m scripts.build_transfer_batch --out inputs/kaggle_transfer.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from nova_lab import REPO_ROOT
from nova_lab.config import NANOBODY_STRUCTURE, PREVIOUS_NANOBODY
from nova_lab.nanobodies import check_sequence
from scripts.train_score_surrogate import load

CALIBRATION = REPO_ROOT / "results" / "scores_nanobody_calibration.csv"
BANDS = 8
PER_BAND = 6
REPEATS = 8
SEED = 78
MINUTES_PER_SEQUENCE = 2.0


def transfer_rows(seed: int = SEED) -> pd.DataFrame:
    """Sequences scored against the previous target, spread over its full range.

    Only sequences with more than one run are eligible: their label is a mean
    over runs, so the comparison's x-axis is measured better than its y-axis
    will be, which is the right way round for a correlation we intend to act on.
    """
    previous = load(PREVIOUS_NANOBODY["structure"])
    pool = previous[previous.runs > 1].copy()
    if pool.empty:
        raise SystemExit("aucune sequence notee plusieurs fois : rien a comparer")
    edges = np.linspace(pool.y.min(), pool.y.max(), BANDS + 1)
    pool["band"] = np.clip(np.digitize(pool.y, edges[1:-1]), 0, BANDS - 1)
    picked = []
    for band, group in pool.groupby("band"):
        take = min(PER_BAND, len(group))
        picked.append(group.sample(take, random_state=seed + int(band)))
        if take < PER_BAND:
            print(f"  note : bande {band} n'offre que {take} sequences sur {PER_BAND}")
    frame = pd.concat(picked, ignore_index=True)
    frame["arm"] = "transfer"
    frame["label"] = [f"band{int(b)}" for b in frame.band]
    frame["note"] = [f"iiptm {y:.4f} contre {PREVIOUS_NANOBODY['structure']} "
                     f"sur {int(n)} run(s)" for y, n in zip(frame.y, frame.runs)]
    return frame.rename(columns={"y": "previous_iiptm", "runs": "previous_runs"})


def repeat_rows(seed: int = SEED) -> pd.DataFrame:
    """Calibration sequences to score again: the live target's own noise."""
    if not CALIBRATION.exists():
        raise SystemExit(f"{CALIBRATION} absent : pas de run de calibration a repeter")
    cal = pd.read_csv(CALIBRATION).dropna(subset=["sequence", "design_iiptm"])
    cal = cal.drop_duplicates("sequence").sort_values("design_iiptm")
    # Spread over the calibration range, so the noise estimate is not taken
    # from one corner of it.
    index = np.linspace(0, len(cal) - 1, REPEATS).round().astype(int)
    frame = cal.iloc[index][["sequence", "design_iiptm", "arm", "label"]].copy()
    frame["note"] = [f"deja note {v:.4f} dans le run calibration, bras {a}"
                     for v, a in zip(frame.design_iiptm, frame.arm)]
    frame = frame.rename(columns={"design_iiptm": "calibration_iiptm",
                                  "arm": "calibration_arm",
                                  "label": "calibration_label"})
    frame["arm"] = "repeat"
    frame["label"] = frame.calibration_label
    frame["previous_iiptm"] = pd.NA
    frame["previous_runs"] = pd.NA
    return frame.reset_index(drop=True)


def build(seed: int = SEED) -> pd.DataFrame:
    transfer = transfer_rows(seed)
    repeats = repeat_rows(seed)
    overlap = set(transfer.sequence) & set(repeats.sequence)
    if overlap:
        # A sequence in both arms would be scored once and serve neither.
        repeats = repeats[~repeats.sequence.isin(overlap)]
    keep = ["sequence", "arm", "label", "note", "previous_iiptm", "previous_runs"]
    frame = pd.concat([transfer[keep], repeats[keep]], ignore_index=True)
    frame = frame.drop_duplicates("sequence").reset_index(drop=True)
    frame["group"] = "candidate"
    frame["epoch"] = pd.NA
    frame["true_rank_sum"] = pd.NA
    # Recorded, never used to drop a row: this batch measures the oracle, so a
    # sequence the validator would refuse is still a valid measurement.
    frame["layer1"] = [", ".join(check_sequence(s)) or "ok" for s in frame.sequence]
    frame["depth"] = frame.groupby("arm").cumcount()
    return frame.sort_values(["depth", "arm"]).drop(columns="depth").reset_index(drop=True)


def analyse(scores: Path, batch: Path) -> dict:
    """Read the run against the thresholds in this module's docstring.

    Written before the data exists, so the verdict cannot be chosen after
    seeing it. Reports the transfer correlation, the live target's own noise
    from the repeat arm, and the correlation with the binder's own fold
    confidence partialled out -- the confound named in the docstring.
    """
    from scipy.stats import spearmanr
    got = pd.read_csv(scores).dropna(subset=["sequence", "design_iiptm"])
    want = pd.read_csv(batch)
    joined = want.merge(got[["sequence", "design_iiptm", "design_ptm"]],
                        on="sequence", how="left")
    moved = joined[joined.arm == "transfer"].dropna(subset=["design_iiptm",
                                                            "previous_iiptm"])
    out: dict = {"n_transfer": len(moved)}
    if len(moved) >= 8:
        x = moved.previous_iiptm.astype(float).to_numpy()
        y = moved.design_iiptm.astype(float).to_numpy()
        rho = spearmanr(x, y)
        out |= {"rho": float(rho.statistic), "p": float(rho.pvalue),
                "se": float(1 / np.sqrt(len(moved) - 1))}
        # The confound: design_ptm is the binder's own fold confidence and
        # never involves the target. Partial it out of both sides by ranks.
        ptm = moved.design_ptm.astype(float).to_numpy()
        ranks = lambda v: pd.Series(v).rank().to_numpy()
        rx = ranks(x) - np.poly1d(np.polyfit(ranks(ptm), ranks(x), 1))(ranks(ptm))
        ry = ranks(y) - np.poly1d(np.polyfit(ranks(ptm), ranks(y), 1))(ranks(ptm))
        out["rho_sans_ptm"] = float(spearmanr(rx, ry).statistic)
    again = joined[joined.arm == "repeat"].dropna(subset=["design_iiptm"])
    if len(again) >= 3:
        delta = (again.design_iiptm.astype(float)
                 - again.note.str.extract(r"deja note ([0-9.]+)")[0].astype(float))
        out |= {"n_repeat": len(again), "biais_run": float(delta.mean()),
                "bruit_cible_vivante": float(delta.std(ddof=1) / np.sqrt(2))}
    return out


def verdict(result: dict) -> str:
    rho = result.get("rho")
    if rho is None:
        return "pas assez de points pour conclure"
    if rho >= 0.6:
        return ("les requetes payees se transferent : le substitut est utilisable "
                "sur la cible vivante des maintenant")
    if rho >= 0.3:
        return ("transfert partiel : utilisable comme a priori, a reentrainer "
                "au fur et a mesure des requetes sur la cible vivante")
    return ("le score est propre a la cible : le substitut doit etre reconstruit "
            "par cible, et les 1 490 requetes sont de l'histoire, pas un actif")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path,
                        default=REPO_ROOT / "inputs" / "kaggle_transfer.csv")
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--analyse", type=Path,
                        help="scores CSV returned by the run; reads it against "
                             "the thresholds instead of building a batch")
    args = parser.parse_args()

    if args.analyse:
        result = analyse(args.analyse, args.out)
        print(f"{result['n_transfer']} sequences comparables")
        if "rho" in result:
            print(f"  Spearman(ancienne cible, cible vivante) = "
                  f"{result['rho']:+.3f} +/- {result['se']:.2f}  p={result['p']:.2g}")
            print(f"  avec la confiance de repliement retiree : "
                  f"{result['rho_sans_ptm']:+.3f}")
        if "n_repeat" in result:
            print(f"  bruit de la cible vivante ({result['n_repeat']} repetitions) : "
                  f"{result['bruit_cible_vivante']:.4f}, biais de run "
                  f"{result['biais_run']:+.4f}")
        print(f"\n-> {verdict(result)}")
        return 0

    frame = build(args.seed)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.out, index=False)

    n_transfer = int((frame.arm == "transfer").sum())
    print(f"{len(frame)} sequences -> {args.out}")
    for arm, group in frame.groupby("arm"):
        print(f"  {arm:<9} {len(group):3d}")
    previous = frame[frame.arm == "transfer"].previous_iiptm.astype(float)
    print(f"  etendue couverte contre {PREVIOUS_NANOBODY['structure']} : "
          f"{previous.min():.3f} a {previous.max():.3f}")
    print(f"  notation contre {NANOBODY_STRUCTURE} (la cible vivante)")
    print(f"  {len(frame) * MINUTES_PER_SEQUENCE / 60:.1f} h sur 2x T4 "
          f"a {MINUTES_PER_SEQUENCE:.0f} min/sequence")
    print(f"  erreur-type attendue du Spearman : "
          f"~{1 / np.sqrt(max(n_transfer - 1, 1)):.2f} sur {n_transfer} points")
    print("  seuils de decision fixes dans le docstring du module, avant le run")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
