"""Does selecting on predicted `design_ptm` raise the score, or only the ptm?

FINDINGS 87 found the one metric that behaves differently from the other nine.
`design_ptm` is the nanobody's own fold confidence, so it never involves the
target, and it has three properties nothing else has: it TRANSFERS between
targets (rho +0.487 against +0.043 for the interface), it PREDICTS from
sequence (rho +0.661 on 1 701 pooled queries against +0.260 for iiptm), and it
carries a THRESHOLD -- below ptm 0.749 the iiptm sits flat at 0.205 whatever
else happens, above it the iiptm climbs, and our median is parked on the
threshold at 0.748 while the field is at 0.830.

That makes it the only lever that would act on the DISTRIBUTION rather than on
the search inside it (FINDINGS 85 measured our distribution topping out at
0.726 of iiptm against the field's 0.75-0.79), and it would cost no GPU:
generate many, keep only what is predicted above 0.78.

None of which is proof. The relation is correlational, and FINDINGS 78 is the
standing lesson on what those are worth when pushed without a pre-registered
reading. So this batch is a dose-response, not a high-versus-control:

  four bands of PREDICTED ptm, 15 sequences each, spanning 0.70 to 0.82 --
  across the 0.749 threshold and into the field's band

Dose-response rather than two arms because two arms cannot distinguish "the
selection works" from "the top band happens to be special". A monotone ladder
can, and it also prices the threshold rather than assuming it.

CONFOUND, and how it is handled. A high-ptm candidate may differ from a low-ptm
one in ways that drive the iiptm by themselves -- CDR3 length above all. The
bands are therefore matched on CDR3 length as closely as the pool allows, and
the match quality is printed so a result can be read against it rather than
taken on trust. Length, V identity and composition are recorded per band for
the same reason.

READINGS FIXED BEFORE THE RUN:

  1. THE PREDICTOR GENERALISES. rho(predicted ptm, measured ptm) over the 60
     must reach 0.4. The 0.661 of FINDINGS 87 was cross-validated; this is
     prospective, and prospective is what failed in FINDINGS 83 (rho +0.145
     where cross-validation said 0.235). Under 0.3 the lever closes here, and
     nothing below matters.
  2. THE THRESHOLD HOLDS ON THE LIVE TARGET. Measured iiptm in the top band
     must beat the bottom band by more than twice the measurement noise
     (0.0133), so by more than 0.027. This is the operational claim: that
     pushing ptm pushes the score.
  3. DOSE-RESPONSE. Spearman(measured ptm, measured iiptm) over the 60 must
     reach 0.3. Our existing live-target data gave -0.066, but over a ptm range
     of sd 0.0148 -- too narrow to see anything. Spanning the range resolves
     whether that null was the relation's absence or the range's.
  4. THE CLEAN NEGATIVE. If 1 holds and 2 and 3 fail, `design_ptm` is
     predictable and useless as a proxy: the FINDINGS 87 lever is closed, and
     that is worth knowing in 2.2 h rather than after a season of generation
     work built on it.

Usage:
    python -m scripts.build_ptm_batch --out inputs/kaggle_ptm.csv
"""
from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd

from nova_lab import DATA_DIR, REPO_ROOT
from nova_lab.config import NANOBODY_STRUCTURE, PREVIOUS_NANOBODY
from nova_lab.nanobodies import check_sequence, germline_identity, regions
from scripts.train_score_surrogate import build_features, load

ARCHIVE = DATA_DIR / "archive" / f"{PREVIOUS_NANOBODY['target']}_regions.csv"
# Bands of predicted ptm. The first sits under the 0.749 threshold, the last in
# the field's own band, and the two middles price the climb between them.
BANDS = [(0.000, 0.750, "sous le seuil"),
         (0.750, 0.772, "juste au-dessus"),
         (0.772, 0.792, "milieu"),
         (0.792, 1.000, "bande du champ")]
PER_BAND = 15
POOL = 6000
SEED = 88
MINUTES_PER_SEQUENCE = 2.2
NOISE = 0.0133          # between-run noise on the live target, FINDINGS 80


def training_set() -> pd.DataFrame:
    """Every sequence we have a measured `design_ptm` for, both targets pooled.

    Pooled deliberately: FINDINGS 87 measured that this metric transfers
    (rho +0.487), which is the whole reason it can be learned once on the
    history instead of being re-bought on each week's target.
    """
    frames = []
    for path in sorted(glob.glob("results/*.csv")):
        try:
            f = pd.read_csv(path)
        except Exception:
            continue
        if {"sequence", "design_ptm"} <= set(f.columns):
            frames.append(f[["sequence", "design_ptm"]].dropna())
    if not frames:
        raise SystemExit("aucun design_ptm mesure dans results/")
    d = pd.concat(frames).groupby("sequence").design_ptm.mean().reset_index()
    d = d[d.sequence.str.fullmatch("[ACDEFGHIKLMNPQRSTVWY]+")]
    return d.rename(columns={"design_ptm": "y"}).reset_index(drop=True)


def predict_pool(seed: int = SEED, pool: int = POOL) -> pd.DataFrame:
    """Archive sequences never queried, with a predicted ptm each."""
    from sklearn.ensemble import HistGradientBoostingRegressor
    known = training_set()
    live = set(load(NANOBODY_STRUCTURE).sequence)
    arch = pd.read_csv(ARCHIVE).dropna(subset=["CDR3"]).drop_duplicates("sequence")
    fresh = arch[~arch.sequence.isin(set(known.sequence) | live)][["sequence", "CDR3"]]
    fresh = fresh.sample(min(pool, len(fresh)), random_state=seed).reset_index(drop=True)
    both = pd.DataFrame({"sequence": list(known.sequence) + list(fresh.sequence)})
    desc, _, kmer, _ = build_features(both)
    X = np.hstack([desc, kmer])
    model = HistGradientBoostingRegressor(max_iter=500, learning_rate=0.05,
                                          random_state=seed)
    model.fit(X[:len(known)], known.y.to_numpy())
    fresh["predicted_ptm"] = model.predict(X[len(known):])
    fresh["cdr3_length"] = fresh.CDR3.str.len()
    return fresh


def matched_bands(pool: pd.DataFrame, per_band: int = PER_BAND,
                  seed: int = SEED) -> pd.DataFrame:
    """One slice per band, CDR3 lengths matched to the pool's own distribution.

    Unmatched bands would let CDR3 length ride along with predicted ptm, and a
    rise in iiptm could then be the loop's length rather than the fold
    confidence. The target distribution is the whole pool's, so no band is
    favoured.
    """
    rng = np.random.default_rng(seed)
    wanted = (pool.cdr3_length.value_counts(normalize=True) * per_band).round().astype(int)
    picked = []
    for lo, hi, label in BANDS:
        band = pool[(pool.predicted_ptm >= lo) & (pool.predicted_ptm < hi)]
        if band.empty:
            print(f"  note : la bande {label} est vide")
            continue
        taken = []
        for length, count in wanted.items():
            if count <= 0:
                continue
            same = band[band.cdr3_length == length]
            if same.empty:
                continue
            taken.append(same.sample(min(count, len(same)),
                                     random_state=int(rng.integers(1 << 30))))
        chosen = pd.concat(taken) if taken else band.head(0)
        if len(chosen) < per_band:          # top up from the band at large
            rest = band[~band.index.isin(chosen.index)]
            if not rest.empty:
                chosen = pd.concat([chosen, rest.sample(
                    min(per_band - len(chosen), len(rest)),
                    random_state=int(rng.integers(1 << 30)))])
        chosen = chosen.head(per_band).copy()
        chosen["band"] = label
        chosen["band_low"], chosen["band_high"] = lo, hi
        picked.append(chosen)
    return pd.concat(picked, ignore_index=True)


def build(seed: int = SEED, pool_size: int = POOL,
          per_band: int = PER_BAND) -> pd.DataFrame:
    frame = matched_bands(predict_pool(seed, pool_size), per_band, seed)
    frame = frame.drop_duplicates("sequence").reset_index(drop=True)
    frame["arm"] = "ptm_" + frame.band.str.replace(" ", "_")
    frame["label"] = [f"{b}#{i:02d}" for b, i in
                      zip(frame.band.str.replace(" ", "_"),
                          frame.groupby("band").cumcount() + 1)]
    frame["note"] = [f"ptm predit {p:.4f}, CDR3 L={l}"
                     for p, l in zip(frame.predicted_ptm, frame.cdr3_length)]
    frame["group"] = "candidate"
    frame["epoch"] = pd.NA
    frame["true_rank_sum"] = pd.NA
    germline = [germline_identity(s) for s in frame.sequence]
    frame["v_identity"] = [round(g["V"]["identity_pct"], 1) if g else None for g in germline]
    # Recorded, never used to drop a row: this batch measures the oracle, so a
    # sequence the validator would refuse is still a valid data point.
    frame["layer1"] = [", ".join(check_sequence(s)) or "ok" for s in frame.sequence]
    frame["depth"] = frame.groupby("arm").cumcount()
    return frame.sort_values(["depth", "arm"]).drop(columns=["depth", "CDR3"]).reset_index(drop=True)


def analyse(scores: Path, batch: Path) -> dict:
    """The three readings of this module's docstring, and nothing else."""
    from scipy.stats import spearmanr
    got = pd.read_csv(scores).dropna(subset=["sequence", "design_iiptm"])
    want = pd.read_csv(batch)
    j = want.merge(got[["sequence", "design_iiptm", "design_ptm"]], on="sequence")
    out = {"n": len(j)}
    if len(j) < 12:
        return out
    out["rho_predictor"] = float(spearmanr(j.predicted_ptm, j.design_ptm).statistic)
    out["rho_dose"] = float(spearmanr(j.design_ptm, j.design_iiptm).statistic)
    order = [label for _, _, label in BANDS]
    by = {b: g for b, g in j.groupby("band")}
    present = [b for b in order if b in by]
    if len(present) >= 2:
        low, high = by[present[0]], by[present[-1]]
        out |= {"iiptm_bas": float(low.design_iiptm.median()),
                "iiptm_haut": float(high.design_iiptm.median()),
                "ecart": float(high.design_iiptm.median() - low.design_iiptm.median()),
                "seuil_ecart": 2 * NOISE}
    out["bandes"] = {b: {"n": len(g), "ptm_predit": round(g.predicted_ptm.median(), 4),
                         "ptm_mesure": round(g.design_ptm.median(), 4),
                         "iiptm": round(g.design_iiptm.median(), 4),
                         "cdr3_median": float(g.cdr3_length.median())}
                     for b, g in by.items()}
    return out


def verdict(result: dict) -> str:
    if "rho_predictor" not in result:
        return "pas assez de points pour conclure"
    if result["rho_predictor"] < 0.3:
        return ("le predicteur ne generalise pas : le rho 0,661 du §87 etait "
                "de la validation croisee, et le levier se ferme ici")
    climbs = result.get("ecart", 0) > result.get("seuil_ecart", 1)
    dose = result.get("rho_dose", 0) >= 0.3
    if climbs and dose:
        return ("le seuil tient sur la cible vivante : pousser le ptm pousse le "
                "score, et c'est une contrainte de generation utilisable")
    if climbs or dose:
        return ("un seul des deux criteres passe : signal reel mais insuffisant "
                "pour en faire une contrainte, a reprendre avec plus de points")
    return ("design_ptm est predictible et inutile comme procuration : le levier "
            "du §87 est ferme, pour 2,2 h de GPU au lieu d'une saison")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path,
                        default=REPO_ROOT / "inputs" / "kaggle_ptm.csv")
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--analyse", type=Path, help="scores CSV from the run")
    args = parser.parse_args()

    if args.analyse:
        r = analyse(args.analyse, args.out)
        print(f"{r['n']} sequences comparables")
        if "rho_predictor" in r:
            print(f"  1. le predicteur : rho(predit, mesure) = {r['rho_predictor']:+.3f} "
                  f"(seuil 0,40 ; sous 0,30 le levier se ferme)")
            print(f"  2. le seuil : iiptm bas {r['iiptm_bas']:.4f} -> haut "
                  f"{r['iiptm_haut']:.4f}, ecart {r['ecart']:+.4f} "
                  f"(il faut > {r['seuil_ecart']:.4f})")
            print(f"  3. dose-reponse : rho(ptm mesure, iiptm) = {r['rho_dose']:+.3f} "
                  f"(seuil 0,30)")
            print("\n  par bande :")
            for b, v in r["bandes"].items():
                print(f"    {b:<18} n={v['n']:>2}  ptm predit {v['ptm_predit']:.4f} "
                      f"-> mesure {v['ptm_mesure']:.4f}  iiptm {v['iiptm']:.4f}  "
                      f"CDR3 median {v['cdr3_median']:.0f}")
        print(f"\n-> {verdict(r)}")
        return 0

    frame = build(args.seed)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.out, index=False)
    print(f"{len(frame)} sequences -> {args.out}")
    print(f"{'bande':<18} {'n':>3} {'ptm predit':>22} {'CDR3 median':>12} "
          f"{'V id median':>12} {'couche 1 ok':>12}")
    for (lo, hi, label) in BANDS:
        g = frame[frame.band == label]
        if g.empty:
            continue
        print(f"{label:<18} {len(g):>3} {g.predicted_ptm.min():>10.4f}-"
              f"{g.predicted_ptm.max():<11.4f} {g.cdr3_length.median():>12.0f} "
              f"{g.v_identity.median():>12.1f} {(g.layer1 == 'ok').sum():>12}")
    lengths = frame.groupby("band").cdr3_length.median()
    print(f"\n  appariement des longueurs de CDR3 entre bandes : "
          f"median {lengths.min():.0f} a {lengths.max():.0f} "
          f"({'bon' if lengths.max() - lengths.min() <= 2 else 'A LIRE AVEC PRUDENCE'})")
    print(f"  {len(frame) * MINUTES_PER_SEQUENCE / 60:.1f} h sur 2x T4")
    print("  les quatre lectures sont fixees dans le docstring, avant le run")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
