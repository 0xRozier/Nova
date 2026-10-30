"""Exact nanobody-competition results per epoch.

For every decrypted submission (scripts.fetch_submissions) fetch the BoltzGen
metrics validators published for (sequence hash, epoch); a sequence has
metrics only if it was valid. Then rank like upstream boltzgen_wrapper +
ranking.py: dense rank per metric across the epoch's unique sequences, sum of
ranks (lower wins), ties broken by the confidence / physical interaction /
developability rank sums, then commitment block.

Usage:
    python -m scripts.analyze_nanobodies [--target Q9NZQ7]
"""
import argparse
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
import requests

from nova_lab import DATA_DIR

API = "https://vali-score-share-api.metanova-labs.ai/api/v1/nanobodies/{hash}/{epoch}"
CACHE = DATA_DIR / "scores" / "nanobodies.jsonl"
# config/boltzgen_config.yaml: metric -> (category, mode)
METRICS = {
    "design_iiptm": ("confidence", "max"),
    "design_ptm": ("confidence", "max"),
    "design_to_target_iptm": ("confidence", "max"),
    "min_design_to_target_pae": ("confidence", "min"),
    "interaction_pae": ("confidence", "min"),
    "plip_hbonds_refolded": ("physical_interaction", "max"),
    "plip_saltbridge_refolded": ("physical_interaction", "max"),
    "delta_sasa_refolded": ("physical_interaction", "max"),
    "liability_score": ("developability", "min"),
    "liability_num_violations": ("developability", "min"),
}
CATEGORIES = ("confidence", "physical_interaction", "developability")


def seq_hash(seq: str) -> str:
    return hashlib.sha256(seq.encode("ascii")).hexdigest()


def load_cache() -> dict:
    cache = {}
    if CACHE.exists():
        for line in CACHE.open():
            row = json.loads(line)
            cache[(row["hash"], row["epoch"])] = row["metrics"]
    return cache


def fetch(key):
    h, epoch = key
    try:
        resp = requests.get(API.format(hash=h, epoch=epoch), timeout=20)
    except requests.RequestException:
        return key, None, False
    if resp.status_code == 404:
        return key, None, True
    if resp.status_code != 200:
        return key, None, False
    averages = resp.json().get("target_averages") or {}
    return key, averages, True


def rank_epoch(g: pd.DataFrame) -> pd.DataFrame:
    """Upstream ranking over the epoch's unique valid sequences."""
    uniq = g.drop_duplicates("hash").copy()
    for metric, (_, mode) in METRICS.items():
        uniq[f"{metric}_rank"] = uniq[metric].rank(method="dense", ascending=(mode == "min"))
    uniq["rank_sum"] = uniq[[f"{m}_rank" for m in METRICS]].sum(axis=1)
    for cat in CATEGORIES:
        uniq[f"{cat}_rank_sum"] = uniq[[f"{m}_rank" for m, (c, _) in METRICS.items() if c == cat]].sum(axis=1)
    cols = ["hash", "rank_sum", *[f"{c}_rank_sum" for c in CATEGORIES]]
    return g.merge(uniq[cols], on="hash").sort_values(
        ["rank_sum", *[f"{c}_rank_sum" for c in CATEGORIES], "block"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default="Q9NZQ7")
    args = parser.parse_args()

    subs = [json.loads(line) for line in (DATA_DIR / "submissions" / "submissions.jsonl").open()]
    rows = []
    for s in subs:
        if not s["sequences"]:
            continue
        seq = s["sequences"][0]  # num_sequences = 1
        rows.append({"epoch": s["epoch"], "uid": s["uid"], "hotkey": s["hotkey"], "block": s["block"],
                     "repo": "/".join(s["path"].split("/")[:2]), "sequence": seq, "hash": seq_hash(seq)})
    df = pd.DataFrame(rows)

    cache = load_cache()
    todo = sorted({(h, e) for h, e in zip(df["hash"], df["epoch"])} - set(cache))
    with ThreadPoolExecutor(8) as pool, CACHE.open("a") as f:
        for key, averages, final in pool.map(fetch, todo):
            if final:  # do not cache transient errors
                cache[key] = averages
                f.write(json.dumps({"hash": key[0], "epoch": key[1], "metrics": averages}) + "\n")
    print(f"metrics fetched for {len(todo)} new (hash, epoch) pairs")

    def metrics_for(h, e):
        m = (cache.get((h, e)) or {}).get(args.target)
        return m if m and all(m.get(k) is not None for k in METRICS) else None

    df["metrics"] = [metrics_for(h, e) for h, e in zip(df["hash"], df["epoch"])]
    df["valid"] = df["metrics"].notna()
    valid = df[df["valid"]].copy()
    for k in METRICS:
        valid[k] = valid["metrics"].map(lambda m, k=k: m[k])
    valid = valid.drop(columns="metrics")

    ranked = pd.concat([rank_epoch(g) for _, g in valid.groupby("epoch")])
    table = []
    for epoch, g in ranked.groupby("epoch"):
        w = g.iloc[0]
        table.append({"epoch": epoch, "submitted": int((df["epoch"] == epoch).sum()), "valid": len(g),
                      "unique": g["hash"].nunique(), "winner_uid": w["uid"], "winner_repo": w["repo"],
                      "rank_sum": w["rank_sum"], "second": g["rank_sum"].iloc[1] if len(g) > 1 else None})
    table = pd.DataFrame(table).set_index("epoch")
    with pd.option_context("display.width", 160, "display.max_rows", 200):
        print(table)
        print(f"\ninvalid submissions: {(~df['valid']).sum()} / {len(df)}")
        print("\nwins by repo:")
        print(table["winner_repo"].value_counts().head(10).to_string())

    out = DATA_DIR / "submissions" / f"{args.target}_nanobody_results.csv"
    ranked.to_csv(out, index=False)
    print(f"\nranked submissions saved to {out}")


if __name__ == "__main__":
    main()
