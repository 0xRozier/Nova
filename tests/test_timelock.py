"""The submission file we build must be the one a validator can read.

`nova_lab.timelock.decrypt` is the function that decrypted real miners'
submissions, and `parse_file` reimplements the validator's own validation, so
round-tripping through both is the check that our encryption produces the real
format rather than merely a self-consistent one.
"""
import ast

from nova_lab import timelock


def test_submission_file_round_trips():
    # A round already in the past, so its drand signature exists and the
    # ciphertext can be opened straight away.
    reveal = timelock.latest_round() - 20
    uid, payload = 42, "~|EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQ"

    filename, content = timelock.submission_file(uid, reveal, payload)

    # The validator reads the file as a literal (round, bytes) tuple.
    round_number, ciphertext = ast.literal_eval(content)
    assert round_number == reveal
    assert isinstance(ciphertext, bytes)

    # And demands the name be the content hash, which parse_file enforces.
    assert timelock.parse_file(f"owner/repo/main/{filename}", content.encode()) is not None

    # The payload comes back, UID prefix and all.
    assert timelock.decrypt(uid, round_number, ciphertext) == payload
    # A different UID must not claim it.
    assert timelock.decrypt(uid + 1, round_number, ciphertext) is None


def test_reveal_round_lands_on_the_submission_deadline():
    epoch_length, no_submission = 361, 10
    block = 100 * epoch_length + 5          # early in an epoch
    deadline = 101 * epoch_length - no_submission
    expected = timelock.current_round() + (deadline - block) * timelock.ROUNDS_PER_BLOCK
    assert timelock.reveal_round(block, epoch_length, no_submission) == expected
