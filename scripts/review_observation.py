"""What the observed candidates were actually worth, once BoltzGen has spoken.

scripts.observe says what it would have submitted; it cannot say whether that
would have won, because the surrogate only orders candidates. This script takes
the BoltzGen scores of an observation batch — our shortlist and the epoch's real
submissions, scored in the same run — ranks them the way a validator does and
reports the place our candidate would have taken.

Scoring both sides in one run matters: BoltzGen is stochastic, and the measured
spread between two runs of the same sequences was wide enough (around 14 places
on 46 sequences) that metrics from different runs cannot be compared.

Usage:
    python -m scripts.review_observation results/boltzgen_batch4_epoch25989.csv
    python -m scripts.review_observation --summary       # every reviewed epoch
"""
import argparse
import json
from pathlib import Path

import pandas as pd

from nova_lab import DATA_DIR
from scripts.analyze_nanobodies import CATEGORIES, METRICS
from scripts.observe import JOURNAL, OBS_DIR

REVIEWS = OBS_DIR / "reviews.jsonl"


def rank_pool(df: pd.DataFrame) -> pd.DataFrame:
    """Dense rank per metric, summed — the validator's ranking, applied to the
    pool actually scored (our candidates plus the real submissions)."""
    out = df.drop_duplicates("sequence").copy()
    for metric, (_, mode) in METRICS.items():
        out[f"{metric}_rank"] = out[metric].rank(method="dense", ascending=(mode == "min"))
    out["pool_rank_sum"] = out[[f"{m}_rank" for m in METRICS]].sum(axis=1)
    for cat in CATEGORIES:
        out[f"{cat}_rank_sum"] = out[[f"{m}_rank" for m, (c, _) in METRICS.items() if c == cat]].sum(axis=1)
    order = ["pool_rank_sum", *[f"{c}_rank_sum" for c in CATEGORIES]]
    return out.sort_values(order).reset_index(drop=True)


def epoch_of(path: Path, journal: list[dict]) -> int | None:
    for entry in journal:
        if entry.get("sequence") and entry["sequence"] in set(pd.read_csv(path)["sequence"]):
            return entry["epoch"]
    return None


def review_all(path: Path, journal: list[dict], epoch: int | None = None) -> list[dict]:
    """A combined batch carries its own epoch column: rank each epoch on its own,
    the way a validator does, instead of pooling epochs together."""
    scored = pd.read_csv(path)
    if "epoch" not in scored.columns or scored["epoch"].nunique() <= 1:
        return [review(path, journal, epoch, scored)]
    out = []
    for ep, part in scored.groupby("epoch"):
        try:
            out.append(review(path, journal, int(ep), part))
        except SystemExit as err:
            print(f"epoch {int(ep)}: skipped ({err})\n")
        print()
    return out


def review(path: Path, journal: list[dict], epoch: int | None = None,
           scored: pd.DataFrame | None = None) -> dict:
    scored = pd.read_csv(path) if scored is None else scored
    missing = [m for m in METRICS if m not in scored.columns]
    if missing:
        raise SystemExit(f"{path} lacks the metrics {missing}: score it with the BoltzGen notebook")
    scored = scored.dropna(subset=list(METRICS))
    if "group" not in scored.columns:
        raise SystemExit(f"{path} has no 'group' column; it is not an observation batch")

    epoch = epoch if epoch is not None else epoch_of(path, journal)
    entry = next((e for e in journal if e["epoch"] == epoch), {})
    ranked = rank_pool(scored)
    ours = ranked[ranked["group"] == "candidate"]
    theirs = ranked[ranked["group"] != "candidate"]
    if ours.empty or theirs.empty:
        raise SystemExit("needs both our candidates and real submissions scored together")

    best_ours = ours.iloc[0]
    beaten = int((theirs["pool_rank_sum"] > best_ours["pool_rank_sum"]).sum())
    place = int((theirs["pool_rank_sum"] < best_ours["pool_rank_sum"]).sum()) + 1
    result = {
        "epoch": epoch, "scores": str(path),
        "allowed_v": entry.get("allowed_v"), "allowed_cdr3_length": entry.get("allowed_cdr3_length"),
        "candidates_scored": int(len(ours)), "rivals_scored": int(len(theirs)),
        "our_best_rank_sum": float(best_ours["pool_rank_sum"]),
        "best_rival_rank_sum": float(theirs["pool_rank_sum"].min()),
        "place_among_rivals": place, "rivals_beaten": beaten,
        "would_have_won": bool(place == 1),
        "sequence": best_ours["sequence"],
    }

    print(f"epoch {epoch}: V {result['allowed_v']}, CDR3 {result['allowed_cdr3_length']}")
    print(f"pool: {len(ours)} of ours + {len(theirs)} real submissions, scored in one BoltzGen run")
    cols = ["group", "repo", "pool_rank_sum", "design_iiptm", "interaction_pae", "liability_score"]
    cols = [c for c in cols if c in ranked.columns]
    with pd.option_context("display.width", 180):
        print(ranked[cols].head(12).to_string())
    print(f"\nour best scores {result['our_best_rank_sum']:.0f}, best rival "
          f"{result['best_rival_rank_sum']:.0f}")
    print(f"  -> place {place} of {len(theirs) + 1}; beats {beaten} of {len(theirs)} real submissions")
    print("  caution: only the real submissions carried into this batch are in the pool, and only "
          "those the validator itself accepted; developability rejections stay invisible offline.")
    return result


def summary() -> None:
    if not REVIEWS.exists():
        raise SystemExit("no review yet: run this script on a scored batch first")
    rows = [json.loads(l) for l in REVIEWS.open() if l.strip()]
    df = pd.DataFrame(rows).drop_duplicates("epoch", keep="last").sort_values("epoch")
    with pd.option_context("display.width", 180):
        print(df[["epoch", "allowed_v", "allowed_cdr3_length", "rivals_scored",
                  "our_best_rank_sum", "best_rival_rank_sum", "place_among_rivals",
                  "would_have_won"]].to_string(index=False))
    n = len(df)
    print(f"\n{int(df['would_have_won'].sum())} of {n} reviewed epochs would have been won "
          f"({df['would_have_won'].mean():.0%})")
    print(f"median place among the scored rivals: {df['place_among_rivals'].median():.0f}")
    if n < 20:
        print(f"{n} epochs is too few to call this a win rate: the spread between two BoltzGen "
              f"runs of the same sequences is itself several places.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("scores", nargs="*", type=Path, help="BoltzGen score CSVs for observation batches")
    parser.add_argument("--epoch", type=int, help="label the batch with this epoch "
                        "when the journal cannot identify it (older batches)")
    parser.add_argument("--summary", action="store_true")
    args = parser.parse_args()

    if args.scores:
        journal = [json.loads(l) for l in JOURNAL.open() if l.strip()] if JOURNAL.exists() else []
        OBS_DIR.mkdir(parents=True, exist_ok=True)
        for path in args.scores:
            for result in review_all(path, journal, args.epoch):
                with REVIEWS.open("a") as f:
                    f.write(json.dumps(result) + "\n")
    if args.summary or not args.scores:
        summary()


if __name__ == "__main__":
    main()
