"""Check nanobody sequences against the validator's rules, offline, before submitting.

Resolves the current epoch's allowed V gene and CDR3 length from the chain
(or takes them as flags), then reports every problem per sequence. A sequence
reported OK can still be refused by the structure-based developability checks,
which need a GPU pipeline (see nova_lab.nanobodies).

Usage:
    python -m scripts.check_nanobody EVQLVESGGG...                  # current epoch
    python -m scripts.check_nanobody --file candidates.txt          # one sequence per line
    python -m scripts.check_nanobody --epoch 26005 --file seqs.txt  # replay a past epoch
    python -m scripts.check_nanobody --allowed-v IGHV3-23 --cdr3 17 SEQ   # no chain access
"""
import argparse
import sys
from pathlib import Path

from nova_lab.chain import current_block
from nova_lab.config import EPOCH_LENGTH, NANOBODY_TARGET
from nova_lab.nanobodies import NanobodyArchive, check_sequence, germline_identity, regions


def resolve_epoch(epoch: int | None, max_retries: int = 5, base_delay: float = 30) -> dict:
    import time
    import bittensor as bt
    from nova_lab.reactions import reaction_ids
    from scripts.epoch_info import challenge_params
    url = "wss://archive.chain.opentensor.ai:443" if epoch else "finney"
    for attempt in range(max_retries + 1):
        try:
            st = bt.Subtensor(url)
            epoch_val = epoch if epoch is not None else current_block(st) // EPOCH_LENGTH
            params = challenge_params(st.determine_block_hash(epoch_val * EPOCH_LENGTH), len(reaction_ids()))
            return {"epoch": epoch_val, **params}
        except Exception as exc:
            msg = f"{type(exc).__name__}: {exc}"
            transient = any(tok in msg for tok in ("429", "ConnectionClosed",
                                                   "TimeoutError", "WebSocket",
                                                   "Timeout", "timed out"))
            if not transient or attempt == max_retries:
                raise
            delay = base_delay * (2 ** attempt)
            print(f"[resolve_epoch] transient error (attempt {attempt + 1}/"
                  f"{max_retries + 1}): {msg}; retrying in {delay:.0f}s")
            time.sleep(delay)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("sequences", nargs="*")
    parser.add_argument("--file", type=Path, help="file with one sequence per line")
    parser.add_argument("--epoch", type=int, help="replay this epoch instead of the current one")
    parser.add_argument("--allowed-v", help="skip the chain lookup and use this V gene")
    parser.add_argument("--cdr3", type=int, help="skip the chain lookup and use this CDR3 length")
    parser.add_argument("--target", default=NANOBODY_TARGET[0])
    parser.add_argument("--no-archive", action="store_true", help="skip the similarity checks (faster)")
    args = parser.parse_args()

    sequences = list(args.sequences)
    if args.file:
        sequences += [l.strip() for l in args.file.read_text().splitlines() if l.strip()]
    if not sequences:
        parser.error("give sequences as arguments or with --file")

    if args.allowed_v and args.cdr3:
        allowed_v, cdr3, label = args.allowed_v, args.cdr3, "given"
    else:
        p = resolve_epoch(args.epoch)
        allowed_v, cdr3, label = p["allowed_v"], p["allowed_cdr3_length"], f"epoch {p['epoch']}"
    print(f"{label}: V gene {allowed_v}, CDR3 length {cdr3}\n")

    archive = None if args.no_archive else NanobodyArchive(args.target)
    if archive:
        print(f"archive: {len(archive.sequences)} sequences already submitted for {args.target}\n")

    n_ok = 0
    for seq in sequences:
        warnings = []
        problems = check_sequence(seq, allowed_v, cdr3, archive, warnings=warnings)
        r = regions(seq)
        g = germline_identity(seq)
        head = f"{seq[:40]}... ({len(seq)} aa)"
        if r and g:
            head += (f"  V {g['V']['identity_pct']:.1f}% {'/'.join(g['V']['genes'])}"
                     f"  J {g['J']['identity_pct']:.1f}%  CDR3 {r['CDR3']} ({len(r['CDR3'])})")
        print(head)
        for p in problems:
            print(f"  REJECTED: {p}")
        for w in warnings:
            print(f"  WARNING: {w} (we check the whole archive, the validator only "
                  f"its k-mer candidates: it may still accept this)")
        if not problems:
            n_ok += 1
            print("  OK (unless structure-based developability rejects it)")
        print()

    print(f"{n_ok}/{len(sequences)} sequences pass the offline checks")
    return 0 if n_ok else 1


if __name__ == "__main__":
    sys.exit(main())
