"""Timelock encryption and decryption for NOVA submissions (drand quicknet).

Miners encrypt with ideal-lab's `timelock` (drand quicknet). bittensor_drand
wraps the same ciphertext format in a SCALE-encoded `UserData {encrypted_data:
Vec<u8>, reveal_round: u64}`, so we add that envelope and let it decrypt.
"""
import ast
import hashlib
import struct
from functools import lru_cache

import bittensor_drand
import requests

QUICKNET_CHAIN = "52db9ba70e0cc0f6eaf7803dd07447a1f5477735fd3f661792ba94600c84e971"
DRAND_URLS = ["https://api.drand.sh", "https://drand.cloudflare.com", "https://api2.drand.sh"]
MAX_RESPONSE_SIZE = 20 * 1024  # same cap as the validator


def _get(path: str) -> dict:
    last_error = None
    for base in DRAND_URLS:
        try:
            resp = requests.get(f"{base}/{QUICKNET_CHAIN}/{path}", timeout=10)
            if resp.status_code == 200:
                return resp.json()
            last_error = f"HTTP {resp.status_code}"
        except requests.RequestException as e:
            last_error = e
    raise RuntimeError(f"drand {path}: {last_error}")


def latest_round() -> int:
    return _get("public/latest")["round"]


@lru_cache(maxsize=4096)
def signature(round_number: int) -> str:
    return _get(f"public/{round_number}")["signature"]


def _compact(n: int) -> bytes:
    if n < 1 << 6:
        return bytes([n << 2])
    if n < 1 << 14:
        return struct.pack("<H", (n << 2) | 1)
    if n < 1 << 30:
        return struct.pack("<I", (n << 2) | 2)
    raise ValueError("ciphertext too large")


def parse_file(path: str, content: bytes) -> tuple[int, bytes] | None:
    """Validate a submission file like the validator: hash-named, literal (round, bytes)."""
    text = content.decode("utf-8", errors="replace")
    if not path.endswith(f"/{hashlib.sha256(text.encode('utf-8')).hexdigest()[:20]}.txt"):
        return None
    if len(text) > MAX_RESPONSE_SIZE:
        return None
    try:
        value = ast.literal_eval(text)
    except (SyntaxError, ValueError, MemoryError, RecursionError, TypeError):
        return None
    if isinstance(value, tuple) and len(value) == 2 and isinstance(value[0], int) and isinstance(value[1], bytes):
        return value
    return None


# utils/btdr.py: quicknet genesis and period, and the block/round conversion the
# validator's own encrypt uses (one 12 s block is four 3 s drand rounds).
NET_START = 1692803367
ROUND_DURATION = 3
ROUNDS_PER_BLOCK = 4


def current_round() -> int:
    """The round drand is on now, from the clock rather than the network."""
    import time
    return int(time.time() - NET_START) // ROUND_DURATION


def reveal_round(current_block: int, epoch_length: int, no_submission_blocks: int) -> int:
    """The round that unlocks where the epoch stops accepting submissions.

    Mirrors utils/btdr.py: the target block is the next epoch boundary minus the
    no-submission window, and the round is scaled from the current one.
    """
    next_boundary = ((current_block // epoch_length) + 1) * epoch_length
    target_block = next_boundary - no_submission_blocks
    return current_round() + (target_block - current_block) * ROUNDS_PER_BLOCK


def encrypt(uid: int, round_number: int, message: str) -> bytes:
    """The raw ciphertext a submission file carries, for `{uid}:{message}`.

    bittensor_drand returns the SCALE `UserData` envelope, while the file holds
    the bare ciphertext — the same thing `decrypt` below takes and wraps again.
    Stripping the envelope here is what keeps the two sides symmetric, and what
    makes our files parse like the real submissions we decrypted.
    """
    import bittensor_drand

    blob, got = bittensor_drand.encrypt_at_round(f"{uid}:{message}".encode("ascii"), round_number)
    if got != round_number:
        raise RuntimeError(f"asked for round {round_number}, got {got}")
    blob = bytes(blob)
    body = blob[:-8]  # trailing u64 reveal_round
    first = body[0] & 0b11
    header = {0: 1, 1: 2, 2: 4}.get(first)
    if header is None:
        raise RuntimeError("unexpected SCALE length prefix in the ciphertext envelope")
    length = int.from_bytes(body[:header], "little") >> 2
    ciphertext = body[header:]
    if len(ciphertext) != length:
        raise RuntimeError(f"envelope says {length} bytes, found {len(ciphertext)}")
    return ciphertext


def submission_file(uid: int, round_number: int, message: str) -> tuple[str, str]:
    """`(filename, content)` for the GitHub file the commitment points at.

    The validator reads the file as a literal `(round, bytes)` tuple and demands
    the name be the first 20 hex characters of the content's SHA-256.
    """
    content = repr((round_number, encrypt(uid, round_number, message)))
    return f"{hashlib.sha256(content.encode('utf-8')).hexdigest()[:20]}.txt", content


def decrypt(uid: int, round_number: int, ciphertext: bytes) -> str | None:
    """Plaintext payload ('mols|seqs') if it decrypts and carries this UID's prefix."""
    envelope = _compact(len(ciphertext)) + ciphertext + struct.pack("<Q", round_number)
    try:
        plaintext = bittensor_drand.decrypt_with_signature(envelope, signature(round_number))
    except ValueError:
        return None
    text = plaintext.decode("utf-8", errors="replace") if isinstance(plaintext, bytes) else str(plaintext)
    prefix = f"{uid}:"
    return text[len(prefix):] if text.startswith(prefix) else None
