"""Bring a fresh machine to the point where the miner can run.

What each piece costs, and whether it can be cached:

| piece                    | source                  | cost        | cacheable |
|--------------------------|-------------------------|-------------|-----------|
| submission archive       | HuggingFace             | seconds     | no, it grows every epoch |
| IMGT regions of the archive | computed             | ~10 min cold, 0.6 s with the cache | yes |
| training dataset         | computed from the archive | ~10 min   | yes, rebuild daily |
| surrogate model          | trained                 | minutes     | yes, 23 MB |
| MSA and target structure | the upstream clone      | with the clone | yes |
| BoltzGen weights         | HuggingFace             | ~6 GB       | yes |

Only the archive must be fresh at every epoch, and with the region cache in
place that refresh costs well under a second. Everything else is derived and
can be copied from a machine that already has it, or rebuilt off the critical
path.

Usage:
    python -m scripts.bootstrap            # report what is present and what it costs
    python -m scripts.bootstrap --build    # build what is missing
"""
import argparse
import time
from pathlib import Path

from nova_lab import DATA_DIR
from nova_lab.config import NANOBODY_TARGET


def state(target: str) -> list[tuple[str, Path, str]]:
    ds = DATA_DIR / "dataset"
    return [
        ("submission archive", DATA_DIR / "archive" / f"{target}_nanobodies.csv",
         "refreshed from HuggingFace every epoch"),
        ("IMGT region cache", DATA_DIR / "archive" / f"{target}_regions.csv",
         "~10 min to build cold, then 0.6 s per refresh — copy it rather than rebuild"),
        ("training dataset", ds / f"{target}_nanobodies.parquet",
         "~10 min; feeds the donor pool, rebuild daily not hourly"),
        ("surrogate model", ds / f"{target}_nanobody_surrogate.joblib",
         "23 MB; only orders the shortlist, and is worse than chance at it"),
    ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default=NANOBODY_TARGET[0])
    parser.add_argument("--build", action="store_true", help="build what is missing")
    parser.add_argument("--install-prebuilt", action="store_true",
                        help="copy data/prebuilt into place first (the 20 minutes that travel)")
    args = parser.parse_args()

    if args.install_prebuilt:
        import shutil
        src = DATA_DIR / "prebuilt"
        for name, dest in ((f"{args.target}_regions.csv", DATA_DIR / "archive"),
                           (f"{args.target}_nanobodies.parquet", DATA_DIR / "dataset")):
            if (src / name).exists():
                dest.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(src / name, dest / name)
                print(f"installed {name} -> {dest / name}")

    rows = state(args.target)
    # The dataset is two files written together, and --install-prebuilt only
    # ships the parquet: without this, the parquet alone looked complete, the
    # rebuild was skipped and train_nanobody_surrogate died on the missing
    # feature matrix.
    def present(name: str, path: Path) -> bool:
        if not path.exists():
            return False
        if "dataset" in name:
            return (DATA_DIR / "dataset" / f"{args.target}_nanobodies_X.npy").exists()
        return True

    missing = [(name, path, note) for name, path, note in rows
               if not present(name, path)]
    for name, path, note in rows:
        here = present(name, path)
        size = f"{path.stat().st_size / 1e6:.1f} MB" if here else "MISSING"
        mark = "ok " if here else "-- "
        if here and "dataset" in name:
            import pandas as pd
            dated = pd.read_parquet(path, columns=["epoch"])["epoch"].dropna()
            if not dated.empty:
                note = f"covers up to epoch {dated.max():.0f} — REBUILD if that is days old"
                mark = "ok?"
        print(f"  {mark}{name:<22}{size:>10}   {note}")
    if not missing:
        print("\neverything the miner needs is present")
    elif not args.build:
        print(f"\n{len(missing)} missing; pass --build, or copy them from a machine that has them")
        return

    if not args.build:
        return
    import subprocess
    import sys
    for name, path, _ in missing:
        if "archive" in name or "region" in name:
            print(f"\nbuilding {name}...")
            t = time.time()
            from nova_lab.nanobodies import NanobodyArchive
            NanobodyArchive(args.target, refresh=True)
            print(f"  {time.time() - t:.1f}s")
        elif "dataset" in name:
            print(f"\nbuilding {name}...")
            subprocess.run([sys.executable, "-m", "scripts.build_nanobody_dataset",
                            "--target", args.target], check=True)
        elif "surrogate" in name:
            print(f"\ntraining {name}...")
            subprocess.run([sys.executable, "-m", "scripts.train_nanobody_surrogate",
                            "--target", args.target], check=True)


if __name__ == "__main__":
    main()
