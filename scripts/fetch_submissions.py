"""Rebuild what every miner submitted in past epochs, like a validator would.

For each epoch: read the miners' on-chain commitments at the epoch's last block
(archive node), fetch the GitHub file each one points to, and decrypt it now
that its drand round has passed.

Output: data/submissions/submissions.jsonl, one line per (epoch, uid). Epochs
already in the file are skipped, so re-running only adds new ones.

Usage:
    python -m scripts.fetch_submissions --last 50
    python -m scripts.fetch_submissions --from-epoch 25956 --to-epoch 26005
"""
import argparse
import json
from concurrent.futures import ThreadPoolExecutor

import bittensor as bt
import requests
from bittensor.core.chain_data.utils import decode_metadata

from nova_lab import DATA_DIR
from nova_lab import timelock
from nova_lab.config import EPOCH_LENGTH, NO_SUBMISSION_BLOCKS

NETUID = 68
ARCHIVE_NODE = "wss://archive.chain.opentensor.ai:443"
OUT = DATA_DIR / "submissions" / "submissions.jsonl"
FILE_CACHE = DATA_DIR / "submissions" / "files"


def fetch_file(path: str) -> bytes | None:
    cached = FILE_CACHE / path
    if cached.exists():
        return cached.read_bytes()
    try:
        resp = requests.get(f"https://raw.githubusercontent.com/{path}", timeout=15,
                            headers={"Range": f"bytes=0-{timelock.MAX_RESPONSE_SIZE}"})
    except requests.RequestException:
        return None
    if resp.status_code not in (200, 206):
        return None
    cached.parent.mkdir(parents=True, exist_ok=True)
    cached.write_bytes(resp.content)
    return resp.content


def epoch_submissions(st, epoch: int) -> list[dict]:
    start, end = epoch * EPOCH_LENGTH, (epoch + 1) * EPOCH_LENGTH
    block_hash = st.determine_block_hash(end)
    hotkeys = st.metagraph(NETUID, block=end).hotkeys
    rows = []
    for uid, hotkey in enumerate(hotkeys):
        commit = st.substrate.query("Commitments", "CommitmentOf", [NETUID, hotkey], block_hash=block_hash).value
        if not commit or not (start < commit["block"] < end - NO_SUBMISSION_BLOCKS):
            continue
        path = decode_metadata(commit)
        rows.append({"epoch": epoch, "uid": uid, "hotkey": hotkey, "block": commit["block"], "path": path})

    def resolve(row):
        content = fetch_file(row["path"]) if "/" in row["path"] else None
        parsed = timelock.parse_file(row["path"], content) if content else None
        payload = timelock.decrypt(row["uid"], *parsed) if parsed else None
        row["status"] = "ok" if payload else ("undecryptable" if parsed else "no_file")
        mols, seqs = [], []
        if payload and payload.count("|") == 1:
            mol_part, seq_part = payload.strip().split("|")
            mols = [] if mol_part == "~" else [m for m in mol_part.split(",") if m]
            seqs = [] if seq_part == "~" else [s for s in seq_part.split(",") if s]
        row.update(molecules=mols, sequences=seqs)
        return row

    with ThreadPoolExecutor(8) as pool:
        return list(pool.map(resolve, rows))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--last", type=int, help="number of most recent finished epochs")
    parser.add_argument("--from-epoch", type=int)
    parser.add_argument("--to-epoch", type=int)
    args = parser.parse_args()

    st = bt.Subtensor(ARCHIVE_NODE)
    last_done = st.get_current_block() // EPOCH_LENGTH - 1
    if args.last:
        epochs = range(last_done - args.last + 1, last_done + 1)
    else:
        epochs = range(args.from_epoch, min(args.to_epoch, last_done) + 1)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if OUT.exists():
        done = {json.loads(line)["epoch"] for line in OUT.open()}

    with OUT.open("a") as f:
        for epoch in epochs:
            if epoch in done:
                continue
            rows = epoch_submissions(st, epoch)
            for row in rows:
                f.write(json.dumps(row) + "\n")
            f.flush()
            n_mol = sum(1 for r in rows if r["molecules"])
            n_seq = sum(1 for r in rows if r["sequences"])
            n_bad = sum(1 for r in rows if r["status"] != "ok")
            print(f"epoch {epoch}: {len(rows)} commitments, {n_mol} with molecules, "
                  f"{n_seq} with nanobodies, {n_bad} unreadable")


if __name__ == "__main__":
    main()
