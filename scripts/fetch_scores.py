"""Download the Boltz-2 molecule scores published by NOVA validators.

Source: the public score-sharing API used by validators to average their
(non-deterministic) Boltz runs. Every scored molecule appears there with its
per-target metrics and the epoch it was scored in.

Output: data/scores/molecules.jsonl (raw, one item per line, resumable) and
data/scores/molecules.parquet (flattened, deduplicated on name+epoch+target).

Usage:
    python scripts/fetch_scores.py                 # full download
    python scripts/fetch_scores.py --max-pages 50  # only the newest pages
"""
import argparse
import json
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
import requests

API = "https://vali-score-share-api.metanova-labs.ai/api/v1/molecules"
PAGE_SIZE = 100
OUT_DIR = Path(__file__).resolve().parent.parent / "data" / "scores"


def fetch_page(page: int, retries: int = 5) -> list[dict]:
    for attempt in range(retries):
        try:
            resp = requests.get(API, params={"page": page, "page_size": PAGE_SIZE}, timeout=30)
            if resp.status_code == 200:
                return resp.json()["items"]
            if resp.status_code == 429:
                time.sleep(5 * (attempt + 1))
                continue
        except requests.RequestException:
            pass
        time.sleep(2 ** attempt)
    print(f"page {page}: giving up")
    return []


def flatten(items: list[dict]) -> pd.DataFrame:
    rows = []
    for item in items:
        for target, metrics in (item.get("target_averages") or {}).items():
            row = {"name": item["name"], "epoch": item["epoch"], "target": target,
                   "created_at": item.get("created_at")}
            row.update({k: v for k, v in metrics.items() if v is not None})
            rows.append(row)
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.drop_duplicates(subset=["name", "epoch", "target"], keep="last")
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-pages", type=int, default=None)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    first = requests.get(API, params={"page": 1, "page_size": PAGE_SIZE}, timeout=30).json()
    total_pages = first["total_pages"]
    if args.max_pages:
        total_pages = min(total_pages, args.max_pages)
    print(f"{first['total']} molecules, fetching {total_pages} pages")

    # New results are prepended while we download, so pages drift slightly;
    # dedup on (name, epoch, target) absorbs the overlap.
    raw_path = OUT_DIR / "molecules.jsonl"
    items: list[dict] = []
    with ThreadPoolExecutor(args.workers) as pool, raw_path.open("w") as raw:
        for i, page_items in enumerate(pool.map(fetch_page, range(1, total_pages + 1)), 1):
            for item in page_items:
                raw.write(json.dumps(item) + "\n")
            items.extend(page_items)
            if i % 100 == 0:
                print(f"  {i}/{total_pages} pages")

    df = flatten(items)
    df.to_parquet(OUT_DIR / "molecules.parquet", index=False)
    print(f"saved {len(df)} rows ({df['name'].nunique()} unique molecules) to {OUT_DIR}")


if __name__ == "__main__":
    main()
