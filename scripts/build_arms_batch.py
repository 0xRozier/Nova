"""One GPU batch, several arms, the same rivals: the only way to compare them.

BoltzGen is stochastic -- two draws of the same 46 sequences differ by about
14.5 places on the rank sum -- so metrics from different runs are not
comparable. Every arm we want to compare, and the real submissions we want to
beat, therefore go into the **same** run.

The arms this assembles:

  grafted     the generator as the miner runs it, from the fresh donor pool
              (scripts.observe, whose batches already interleave the rivals)
  humanised   the archive's good binders moved onto a compliant framework
              (scripts.humanise) -- the piste of FINDINGS 47
  real        the best submissions of each epoch, with their official rank

Rows come out ordered by arm within each epoch's depth, so an interrupted run
still covers every epoch and every arm rather than finishing one of them.

Usage:
    python -m scripts.build_arms_batch --humanised-dir <dir> --depth 4 --out batch.csv
"""
import argparse
from pathlib import Path

import pandas as pd

from nova_lab import DATA_DIR
from nova_lab.config import NANOBODY_TARGET

BATCH_DIR = DATA_DIR / "observation" / "batches"


def grafted_and_real(target: str, depth: int,
                     epochs: set[int] | None = None) -> pd.DataFrame:
    """Each epoch's journalled batch: our grafted candidates and its rivals.

    `epochs` restricts which journalled batches are taken. Without it this
    picks up every file in the directory, including ones left by earlier runs
    -- and those were produced by a different generator, a different selector,
    or, before 2026-11-16, the archive cutoff that leaked the future
    (docs/REPLAY_AUDIT.md). Mixing them into one batch spends GPU hours
    measuring a mixture nobody can interpret.
    """
    frames = []
    for path in sorted(BATCH_DIR.glob(f"{target}_epoch*.csv")):
        epoch = int(path.stem.split("epoch")[1])
        if epochs is not None and epoch not in epochs:
            continue
        frame = pd.read_csv(path)
        frame["epoch"] = frame.get("epoch", epoch)
        frame["arm"] = frame["group"].map({"candidate": "grafted", "real": "real"})
        for arm in ("grafted", "real"):
            part = frame[frame["arm"] == arm].head(depth)
            frames.append(part)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def humanised(directory: Path, depth: int) -> pd.DataFrame:
    """scripts.humanise output, one CSV per epoch named e<epoch>.csv."""
    frames = []
    for path in sorted(directory.glob("e*.csv")):
        epoch = int(path.stem.lstrip("e"))
        frame = pd.read_csv(path).head(depth)
        frame["epoch"] = epoch
        frame["arm"] = "humanised"
        frame["group"] = "candidate"
        frames.append(frame)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def parse_epochs(spec: str | None) -> set[int] | None:
    """"26304,26345-26347" -> {26304, 26345, 26346, 26347}; None stays None."""
    if not spec:
        return None
    out: set[int] = set()
    for piece in spec.split(","):
        piece = piece.strip()
        if "-" in piece:
            lo, hi = (int(x) for x in piece.split("-", 1))
            out.update(range(lo, hi + 1))
        elif piece:
            out.add(int(piece))
    return out or None


def interleave(frame: pd.DataFrame) -> pd.DataFrame:
    """Depth-major, then epoch, then arm: a short run still covers everything."""
    frame = frame.copy()
    frame["depth"] = frame.groupby(["epoch", "arm"]).cumcount()
    return frame.sort_values(["depth", "epoch", "arm"]).drop(columns="depth")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default=NANOBODY_TARGET[0])
    parser.add_argument("--humanised-dir", type=Path,
                        help="directory of scripts.humanise CSVs (e<epoch>.csv)")
    parser.add_argument("--depth", type=int, default=4,
                        help="sequences per arm per epoch")
    parser.add_argument("--epochs", help="only these epochs, e.g. 26304,26345-26347. "
                                        "Without it every journalled batch on disk is "
                                        "taken, including stale ones from earlier runs.")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    epochs = parse_epochs(args.epochs)
    parts = [grafted_and_real(args.target, args.depth, epochs)]
    if args.humanised_dir:
        parts.append(humanised(args.humanised_dir, args.depth))
    parts = [p for p in parts if not p.empty]
    if not parts:
        raise SystemExit("nothing to assemble: run scripts.observe first")

    keep = ["sequence", "epoch", "arm", "group", "true_rank_sum", "repo",
            "predicted_rank_sum", "cdr3", "v_identity", "donor_design_iiptm"]
    frame = pd.concat(parts, ignore_index=True)
    for column in keep:
        if column not in frame.columns:
            frame[column] = pd.NA
    frame = frame[keep].dropna(subset=["sequence"]).drop_duplicates("sequence")
    frame = interleave(frame)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.out, index=False)
    counts = frame["arm"].value_counts().to_dict()
    minutes = len(frame) * 1.4  # measured: 92.2 min for 64 on 2x T4
    print(f"{len(frame)} sequences over {frame['epoch'].nunique()} epochs -> {args.out}")
    print(f"  {counts}")
    print(f"  about {minutes / 60:.1f} h on 2x T4 at 1.4 min/sequence")
    print("  rows are ordered by depth then epoch then arm, so an interrupted "
          "run still covers every arm of every epoch")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
