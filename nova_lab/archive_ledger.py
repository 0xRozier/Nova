"""When each archived sequence was first seen, because the archive has no dates.

The public submission archive is a flat list of sequences with no submission
date. That is harmless for live mining -- everything in it was submitted before
today -- but it makes an honest replay impossible: re-downloading the archive
injects today's sequences, and a replay cannot tell them from six-week-old
ones. On 15 November this was measured to matter a great deal: of the 197
undated donors that bind at 0.72 or better and clear the germline floors, none
existed at the epochs FINDINGS 48, 49 and 51 replayed. All 197 arrived after
the only dated view we had.

Until now that dated view was `data/prebuilt/<target>_regions.csv`, frozen when
the archive held 9466 entries, so every replay after epoch 26025 saw a
six-week-old archive and nothing newer. This ledger fixes that going forward:
each time the archive is downloaded, the sequences that are new to us are
stamped with the epoch we saw them at. A sequence's stamp is an upper bound on
its real submission epoch -- it was submitted at or before the moment we first
saw it -- which is the safe direction: a replay may miss a donor that existed,
never invent one that did not.

It is append-only. A sequence already in the ledger keeps its first stamp,
because seeing it again says nothing new.
"""
import csv
from pathlib import Path

from nova_lab import DATA_DIR

# The archive held 9466 entries when data/prebuilt/<target>_regions.csv was
# built, which covers everything submitted up to this epoch.
SNAPSHOT_EPOCH = 26025


def path(target: str) -> Path:
    return DATA_DIR / "archive" / f"{target}_first_seen.csv"


def load(target: str) -> dict[str, int]:
    """sequence -> the epoch we first saw it in the archive."""
    file = path(target)
    if not file.exists():
        return {}
    with file.open() as handle:
        return {row["sequence"]: int(row["first_seen_epoch"])
                for row in csv.DictReader(handle)}


def bootstrap(target: str) -> dict[str, int]:
    """The frozen regions file, stamped at the epoch it covers.

    Used when the ledger does not exist yet, so the first run is no worse than
    what the repo did before.
    """
    source = DATA_DIR / "prebuilt" / f"{target}_regions.csv"
    if not source.exists():
        return {}
    with source.open() as handle:
        return {row["sequence"]: SNAPSHOT_EPOCH for row in csv.DictReader(handle)
                if row.get("sequence")}


def record(target: str, sequences, epoch: int) -> int:
    """Stamp sequences we have not seen before. Returns how many were new."""
    known = load(target) or bootstrap(target)
    new = {s for s in sequences if s and s not in known}
    if not new:
        return 0
    known.update({s: int(epoch) for s in new})
    file = path(target)
    file.parent.mkdir(parents=True, exist_ok=True)
    tmp = file.with_suffix(".csv.tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sequence", "first_seen_epoch"])
        writer.writerows(sorted(known.items()))
    tmp.replace(file)          # atomic: a killed run leaves the old ledger
    return len(new)


def seen_before(target: str, epoch: int) -> set[str]:
    """Archived sequences we had already seen before this epoch.

    This is what a replay of `epoch` is allowed to graft from, on the archive
    side. Dated submissions are filtered separately, on their own epoch.
    """
    known = load(target) or bootstrap(target)
    return {s for s, e in known.items() if e < epoch}
