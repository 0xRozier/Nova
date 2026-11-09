"""The chain calls the miner depends on, in one place, with a version guard.

bittensor's Subtensor is not a stable surface. requirements.txt pins 10.4.1,
but this container was found running 11.3.0, where the three calls the miner
makes are all gone: get_current_block, get_uid_for_hotkey_on_subnet and
set_commitment. The miner did not notice until it read the chain, which is
the first thing it does inside an epoch -- so the whole pipeline died with a
bare AttributeError, and the epoch was lost.

So: read the block head through current_block(), which works on both APIs,
and assert the write path exists before an epoch starts rather than after
forty minutes of CPU.
"""
from __future__ import annotations

# What the miner calls. current_block is not in the list: it is shimmed below.
SUBMISSION_API = ("get_uid_for_hotkey_on_subnet", "set_commitment")
PINNED = "10.4.1"


def current_block(st) -> int:
    """The chain head. 10.x exposes a method, 11.x a property."""
    if hasattr(st, "get_current_block"):
        return int(st.get_current_block())
    return int(st.block)


def missing_api(st, needed: tuple[str, ...] = SUBMISSION_API) -> list[str]:
    """Which of the calls we need this Subtensor does not have."""
    return [name for name in needed if not hasattr(st, name)]


def require_api(st, needed: tuple[str, ...] = SUBMISSION_API) -> None:
    """Refuse to start when the installed bittensor cannot submit."""
    missing = missing_api(st, needed)
    if not missing:
        return
    import bittensor as bt
    raise RuntimeError(
        f"bittensor {getattr(bt, '__version__', '?')} has no "
        f"{', '.join(missing)}; the miner cannot submit with it. "
        f"requirements.txt pins {PINNED}: pip install bittensor=={PINNED}"
    )
