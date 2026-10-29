"""Exact molecule-competition results per epoch.

Joins the decrypted submissions (scripts.fetch_submissions) with the published
per-molecule Boltz-2 scores (scripts.fetch_scores). Validators only score valid
submissions, so a submission whose 20 molecules all have a score for that
epoch was valid; its final_molecule_score is their sum. The winner is the best
total, ties going to the earliest commitment block (then push time, which we
do not have).

Usage:
    python -m scripts.analyze_winners [--target O15427]
"""
import argparse
import json

import pandas as pd

from nova_lab import DATA_DIR
from nova_lab.config import NUM_MOLECULES


def load_results(target: str) -> pd.DataFrame:
    subs = [json.loads(line) for line in (DATA_DIR / "submissions" / "submissions.jsonl").open()]
    scores = pd.read_parquet(DATA_DIR / "scores" / "molecules.parquet")
    scores = scores[scores["target"] == target]
    lookup = dict(zip(zip(scores["name"], scores["epoch"]), scores["score"]))

    rows = []
    for s in subs:
        if not s["molecules"]:
            continue
        found = [lookup.get((m, s["epoch"])) for m in s["molecules"]]
        found = [x for x in found if x is not None]
        valid = len(s["molecules"]) == NUM_MOLECULES and len(found) == NUM_MOLECULES
        rows.append({
            "epoch": s["epoch"], "uid": s["uid"], "hotkey": s["hotkey"], "block": s["block"],
            "repo": "/".join(s["path"].split("/")[:2]),
            "rxn": s["molecules"][0].split(":")[1] if s["molecules"][0].startswith("rxn:") else "?",
            "valid": valid,
            "total": sum(found) if valid else None,
            "uses_k": any(m.count(":") >= 4 and m.split(":")[-1].startswith("k") for m in s["molecules"]),
        })
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default="O15427")
    args = parser.parse_args()

    res = load_results(args.target)
    scored_epochs = set(pd.read_parquet(DATA_DIR / "scores" / "molecules.parquet", columns=["epoch"])["epoch"])
    res = res[res["epoch"].isin(scored_epochs)]

    per_epoch = []
    for epoch, g in res.groupby("epoch"):
        v = g[g["valid"]].sort_values(["total", "block"], ascending=[False, True])
        per_epoch.append({
            "epoch": epoch,
            "rxn": g["rxn"].mode().iat[0],
            "submitted": len(g),
            "valid": len(v),
            "winner_uid": v["uid"].iat[0] if len(v) else None,
            "winner": v["total"].iat[0] if len(v) else None,
            "second": v["total"].iat[1] if len(v) > 1 else None,
            "median": v["total"].median() if len(v) else None,
        })
    table = pd.DataFrame(per_epoch).set_index("epoch")
    with pd.option_context("display.float_format", "{:.3f}".format, "display.width", 140, "display.max_rows", 200):
        print(table)
        print(f"\nwinning total: mean {table['winner'].mean():.3f}, median {table['winner'].median():.3f}, "
              f"min {table['winner'].min():.3f}, max {table['winner'].max():.3f}")
        print(f"winner - second: median gap {(table['winner'] - table['second']).median():.3f}")
        print(f"invalid submissions: {(~res['valid']).sum()} / {len(res)}")

        wins = res.merge(table[["winner_uid"]], left_on="epoch", right_index=True)
        wins = wins[wins["uid"] == wins["winner_uid"]]
        print("\nwins by hotkey (miner identity):")
        board = wins.groupby(["hotkey", "repo"]).size().sort_values(ascending=False).rename("wins")
        print(board.head(15).to_string())
        print(f"\nwinners using :k product indices: {wins['uses_k'].mean():.0%}")

    out = DATA_DIR / "submissions" / f"{args.target}_results.csv"
    res.to_csv(out, index=False)
    print(f"\nper-submission results saved to {out}")


if __name__ == "__main__":
    main()
