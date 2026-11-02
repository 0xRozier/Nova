"""Observation-mode miner: generate, check and journal what we *would* submit.

This script never touches a wallet, never writes to the chain and never pushes
to GitHub. It only reads: the epoch parameters from the chain, the archive from
disk. For every epoch it

  1. resolves the imposed V gene and CDR3 length,
  2. runs the generator (scripts.propose_nanobody.propose) on CPU,
  3. journals what it would have committed, in data/observation/journal.jsonl,
  4. writes a GPU batch so BoltzGen can later say what the candidate was worth.

For a replayed epoch the batch also carries that epoch's real submissions, once
they are decrypted (scripts.fetch_submissions): scoring them in the same run is
the only honest comparison, since BoltzGen metrics are not comparable between
runs. For a live epoch the competitors are still encrypted, so the batch holds
only our shortlist and the review has to wait for the decryption.

Usage:
    python -m scripts.observe --once                  # the current epoch
    python -m scripts.observe --epochs 25980,25989    # replay named epochs
    python -m scripts.observe --epochs 25980-26000    # replay a range
    python -m scripts.observe --watch                 # one proposal per epoch, forever
"""
import argparse
import json
import time
from datetime import datetime, timezone

import pandas as pd

from nova_lab import DATA_DIR
from nova_lab.config import EPOCH_LENGTH, NANOBODY_TARGET
from scripts.analyze_nanobodies import seq_hash
from scripts.propose_nanobody import propose

# This miner observes. Submitting means encrypting, pushing to GitHub and
# committing the pointer on chain; none of that lives in this file on purpose.
SUBMIT = False

OBS_DIR = DATA_DIR / "observation"
JOURNAL = OBS_DIR / "journal.jsonl"
BATCH_DIR = OBS_DIR / "batches"
CALIBRATION = 4  # real submissions carried into the batch, best-ranked first


def parse_epochs(spec: str) -> list[int]:
    out = []
    for part in spec.split(","):
        part = part.strip()
        if "-" in part:
            lo, hi = part.split("-")
            out += list(range(int(lo), int(hi) + 1))
        elif part:
            out.append(int(part))
    return out


def journalled() -> set[int]:
    if not JOURNAL.exists():
        return set()
    return {json.loads(line)["epoch"] for line in JOURNAL.open() if line.strip()}


def real_submissions(epoch: int, target: str) -> pd.DataFrame:
    """That epoch's decrypted nanobody submissions, with their real rank sums
    when scripts.analyze_nanobodies has already computed them."""
    path = DATA_DIR / "submissions" / "submissions.jsonl"
    if not path.exists():
        return pd.DataFrame()
    rows = []
    for line in path.open():
        s = json.loads(line)
        if s["epoch"] != epoch or not s.get("sequences"):
            continue
        rows.append({"sequence": s["sequences"][0], "uid": s["uid"],
                     "repo": "/".join(s["path"].split("/")[:2])})
    df = pd.DataFrame(rows).drop_duplicates("sequence")
    if df.empty:
        return df
    results = DATA_DIR / "submissions" / f"{target}_nanobody_results.csv"
    if results.exists():
        ranked = pd.read_csv(results)
        ranked = ranked[ranked["epoch"] == epoch][["sequence", "rank_sum"]].drop_duplicates("sequence")
        df = df.merge(ranked, on="sequence", how="left")
    else:
        df["rank_sum"] = pd.NA
    return df.sort_values("rank_sum", na_position="last")


def write_batch(epoch: int, target: str, shortlist: pd.DataFrame, rivals: pd.DataFrame) -> str:
    """Our shortlist plus the epoch's real submissions, in the format the Kaggle
    BoltzGen notebook reads (scripts.build_boltzgen_notebook)."""
    BATCH_DIR.mkdir(parents=True, exist_ok=True)
    ours = shortlist.assign(group="candidate", true_rank_sum=pd.NA, repo="", stratum=0.0)
    cols = ["sequence", "predicted_rank_sum", "stratum", "cdr3", "v_identity", "group",
            "true_rank_sum", "repo"]
    frames = [ours.reindex(columns=cols)]
    if not rivals.empty:
        theirs = rivals.head(CALIBRATION).assign(
            group="real", predicted_rank_sum=pd.NA, stratum=1.0, cdr3="", v_identity=pd.NA)
        theirs = theirs.rename(columns={"rank_sum": "true_rank_sum"})
        frames.append(theirs.reindex(columns=cols))
    # Interleave ours and theirs: a run interrupted halfway still compares the
    # two sides instead of scoring only our candidates.
    frames = [f.reset_index(drop=True).rename_axis("slot").reset_index() for f in frames]
    batch = pd.concat(frames, ignore_index=True).sort_values(["slot", "stratum"])
    batch = batch.drop(columns="slot").reset_index(drop=True)
    path = BATCH_DIR / f"{target}_epoch{epoch}.csv"
    batch.to_csv(path, index=False)
    return str(path.relative_to(DATA_DIR.parent))


def report(target: str) -> None:
    """What the journal holds so far. Offline only: the place shown is the
    surrogate's guess against a reference pool, not a result."""
    if not JOURNAL.exists():
        raise SystemExit("nothing journalled yet")
    rows = [json.loads(l) for l in JOURNAL.open() if l.strip()]
    df = pd.DataFrame([r for r in rows if r.get("target") == target])
    df = df.drop_duplicates("epoch", keep="last").sort_values("epoch")
    cols = ["epoch", "allowed_v", "allowed_cdr3_length", "status", "generated", "valid",
            "framework_donors", "cdr3_donors", "surrogate_place_vs_reference",
            "real_submissions_known"]
    with pd.option_context("display.width", 200, "display.max_rows", 300):
        print(df[cols].to_string(index=False))
    ok = df[df["status"] == "ok"]
    print(f"\n{len(ok)} of {len(df)} observed epochs produced a submittable candidate")
    if len(ok):
        print(f"median valid candidates per epoch: {ok['valid'].median():.0f} "
              f"of {ok['generated'].median():.0f} generated")
        print("surrogate place against the 200-sequence reference pool: "
              f"median {ok['surrogate_place_vs_reference'].median():.0f}")
    print("none of this was submitted; the surrogate cannot rank candidates "
          "(Spearman ~0.2), so only BoltzGen can say what these were worth: "
          "score the batches, then scripts.review_observation.")


def combined_batch(target: str, depth: int = 2) -> None:
    """One batch for every journalled epoch, since the target never changes.

    BoltzGen scores against the same structure at every epoch, so 21 separate
    Kaggle runs are a waste: one batch covers them all. Rows are ordered by
    depth, not by epoch — round 1 holds each epoch's best candidate and best
    rival, round 2 the second of each, and so on. An interrupted run therefore
    gives a shallow verdict on every epoch instead of a deep one on the first.
    """
    if not JOURNAL.exists():
        raise SystemExit("nothing journalled yet")
    rows = [json.loads(l) for l in JOURNAL.open() if l.strip()]
    seen, entries = set(), []
    for entry in reversed(rows):  # last observation of each epoch wins
        if entry.get("target") == target and "batch" in entry and entry["epoch"] not in seen:
            seen.add(entry["epoch"])
            entries.append(entry)
    entries.sort(key=lambda e: e["epoch"])

    frames = []
    for entry in entries:
        path = BATCH_DIR / f"{target}_epoch{entry['epoch']}.csv"
        if not path.exists():
            continue
        part = pd.read_csv(path)
        part["epoch"] = entry["epoch"]
        part["allowed_v"] = entry["allowed_v"]
        # slot = depth within the epoch, counted separately for both sides
        part["slot"] = part.groupby("group").cumcount()
        frames.append(part[part["slot"] < depth])
    if not frames:
        raise SystemExit("no per-epoch batch found; run scripts.observe first")
    batch = pd.concat(frames, ignore_index=True).sort_values(["slot", "epoch", "stratum"])
    batch = batch.drop(columns="slot").reset_index(drop=True)
    out = BATCH_DIR / f"{target}_combined.csv"
    batch.to_csv(out, index=False)
    per_epoch = batch.groupby("group").size().to_dict()
    print(f"{len(batch)} sequences over {batch['epoch'].nunique()} epochs -> {out}")
    print(f"  {per_epoch} ; about {len(batch) * 3.7 / 60:.1f} h on 2x T4 at 3.7 min/sequence")
    print("  rows are ordered by depth: an interrupted run still covers every epoch")


def rebuild_batches(target: str) -> None:
    """Re-write journalled batches now that more rivals are decrypted.

    A live epoch's competitors are still timelocked when we observe it, so its
    batch leaves with our shortlist alone. Once scripts.fetch_submissions and
    scripts.analyze_nanobodies have run, this puts the epoch's real submissions
    back into the batch, which is what makes the review comparable.
    """
    if not JOURNAL.exists():
        raise SystemExit("nothing journalled yet")
    for line in JOURNAL.open():
        entry = json.loads(line)
        if entry.get("target") != target or "batch" not in entry:
            continue
        path = BATCH_DIR / f"{target}_epoch{entry['epoch']}.csv"
        if not path.exists():
            continue
        old = pd.read_csv(path)
        ours = old[old["group"] == "candidate"].drop(columns=["group", "true_rank_sum", "repo"])
        rivals = real_submissions(entry["epoch"], target)
        n_before = int((old["group"] != "candidate").sum())
        write_batch(entry["epoch"], target, ours, rivals)
        print(f"epoch {entry['epoch']}: {n_before} -> {min(len(rivals), CALIBRATION)} "
              f"real submissions in the batch ({len(rivals)} decrypted)")


def observe_epoch(params: dict, target: str, args) -> dict:
    epoch = params["epoch"]
    # The seed must carry the epoch. Two epochs drawing the same (V gene, CDR3
    # length) are otherwise generated identically: the 21-epoch sweep proposed
    # one sequence three times for (IGHV3-23, 19), and in production the second
    # submission would be an exact archive duplicate, rejected outright.
    result = propose(params, target, candidates=args.candidates, n_donors=args.donors,
                     n_cdr3=args.cdr3_pool, shortlist=args.shortlist,
                     seed=args.seed + epoch)
    shortlist = result["shortlist"]
    rivals = real_submissions(epoch, target)

    entry = {
        "epoch": epoch,
        "observed_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "target": target,
        "allowed_v": result["allowed_v"],
        "allowed_cdr3_length": result["allowed_cdr3_length"],
        "status": result["status"],
        "generated": result["generated"],
        "valid": result["valid"],
        "framework_donors": result.get("framework_donors"),
        "cdr3_donors": result.get("cdr3_donors"),
        "rejections": dict(sorted(result["reasons"].items(), key=lambda kv: -kv[1])[:6]),
        "surrogate_place_vs_reference": result["place"],
        "reference_size": result["reference_size"],
        "predicted_rank_sum": result["best_rank_sum"],
        "real_submissions_known": int(len(rivals)),
        "submitted": SUBMIT,  # always false: this miner only observes
    }
    if not shortlist.empty:
        top = shortlist.iloc[0]
        entry.update({
            "sequence": top["sequence"],
            "sequence_hash": seq_hash(top["sequence"]),
            "cdr3": top["cdr3"],
            "v_identity": float(top["v_identity"]),
            "shortlist_size": int(len(shortlist)),
            "batch": write_batch(epoch, target, shortlist, rivals),
        })
    OBS_DIR.mkdir(parents=True, exist_ok=True)
    with JOURNAL.open("a") as f:
        f.write(json.dumps(entry) + "\n")
    print(f"\njournalled epoch {epoch} ({entry['status']})"
          + (f"; batch -> {entry['batch']}" if "batch" in entry else ""))
    return entry


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--once", action="store_true", help="observe the current epoch and stop")
    parser.add_argument("--epochs", help="replay these epochs, e.g. 25980,25989 or 25980-26000")
    parser.add_argument("--watch", action="store_true", help="keep observing each new epoch")
    parser.add_argument("--target", default=NANOBODY_TARGET[0])
    parser.add_argument("--candidates", type=int, default=3000)
    parser.add_argument("--donors", type=int, default=40)
    parser.add_argument("--cdr3-pool", type=int, default=60)
    parser.add_argument("--shortlist", type=int, default=10)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--redo", action="store_true", help="re-observe epochs already journalled")
    parser.add_argument("--report", action="store_true", help="print the journal and stop")
    parser.add_argument("--combined-batch", action="store_true",
                        help="one GPU batch covering every journalled epoch, and stop")
    parser.add_argument("--depth", type=int, default=2,
                        help="candidates and rivals per epoch in the combined batch")
    parser.add_argument("--rebuild-batches", action="store_true",
                        help="put newly decrypted rivals into the journalled batches and stop")
    args = parser.parse_args()
    if args.report:
        report(args.target)
        return
    if args.combined_batch:
        combined_batch(args.target, args.depth)
        return
    if args.rebuild_batches:
        rebuild_batches(args.target)
        return
    if not (args.once or args.epochs or args.watch):
        parser.error("choose --once, --epochs or --watch")

    from scripts.check_nanobody import resolve_epoch  # imports bittensor lazily

    done = set() if args.redo else journalled()

    if args.epochs:
        for epoch in parse_epochs(args.epochs):
            if epoch in done:
                print(f"epoch {epoch} already journalled; --redo to repeat")
                continue
            observe_epoch(resolve_epoch(epoch), args.target, args)
        if not args.watch:
            return

    if args.once or args.watch:
        while True:
            params = resolve_epoch(None)
            if params["epoch"] in done:
                wait = 60
                print(f"epoch {params['epoch']} already journalled; sleeping {wait}s")
            else:
                observe_epoch(params, args.target, args)
                done.add(params["epoch"])
                wait = 60
            if args.once:
                return
            time.sleep(wait)


if __name__ == "__main__":
    main()
