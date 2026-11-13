"""Which gate refuses the field's submissions, and which would refuse ours?

PSH > 89 explains 36 % of the rejections in epochs 26299-26306 and nothing
about the other 64 % (FINDINGS 43). That matters directly: a gate we do not
reproduce is a gate we cannot screen for, and the production run was refused by
exactly such a gate for eleven epochs running.

So this takes the submissions whose verdict we know -- the validator either
published BoltzGen metrics for them or it did not -- and runs every gate we can
reproduce, to see which one fires on the refused and stays quiet on the
accepted. Three outcomes are possible for each gate, and they mean different
things:

  fires on refused, quiet on accepted   -> we reproduce a real rejection reason
  fires on both                         -> our reproduction is wrong, or stricter
  quiet on both                         -> the reason is somewhere we cannot see

The structures are reused rather than refolded: TNP already built them when
PSH was measured, and Aggrescan3D scores the same models. Pass --tnp-dir to
point at that run.

Usage:
    python -m scripts.why_refused --tnp-dir <dir from the PSH run>
    python -m scripts.why_refused --tnp-dir <dir> --out results/why_refused.csv
"""
import argparse
import json
import time
from collections import Counter
from pathlib import Path

import pandas as pd

from nova_lab import DATA_DIR
from nova_lab.config import EPOCH_LENGTH, NANOBODY_TARGET
from nova_lab.nanobodies import NanobodyArchive, check_sequence, seq_hash
from nova_lab.psh_calibration import THRESHOLD, refusal_probability


def epoch_params(epochs: list[int]) -> dict[int, dict]:
    """The allowed V gene and CDR3 length per epoch, from the chain."""
    import bittensor as bt

    from nova_lab.reactions import reaction_ids
    from scripts.epoch_info import challenge_params
    n_reactions = len(reaction_ids())
    out = {}
    for epoch in epochs:
        for attempt in range(5):
            try:
                st = bt.Subtensor("wss://archive.chain.opentensor.ai:443")
                block_hash = st.determine_block_hash(int(epoch) * EPOCH_LENGTH)
                out[epoch] = challenge_params(block_hash, n_reactions)
                break
            except Exception as error:  # noqa: BLE001
                if attempt == 4:
                    raise
                print(f"  epoch {epoch}: {type(error).__name__}, retrying", flush=True)
                time.sleep(5 * 2 ** attempt)
    return out


def labelled(target: str) -> pd.DataFrame:
    """One row per decrypted submission, with the validator's verdict."""
    scored = set(pd.read_csv(
        DATA_DIR / "submissions" / f"{target}_nanobody_results.csv")["hash"])
    rows = []
    for line in (DATA_DIR / "submissions" / "submissions.jsonl").open():
        s = json.loads(line)
        if s["sequences"]:
            sequence = s["sequences"][0]
            rows.append({"sequence": sequence, "epoch": s["epoch"],
                         "accepted": seq_hash(sequence) in scored})
    return pd.DataFrame(rows).drop_duplicates("sequence")


def structural(frame: pd.DataFrame, tnp_dir: Path) -> pd.DataFrame:
    """TNP's six flags and the Aggrescan3D caps, from models already built."""
    from scripts.aggrescan_gate import MAX_AVERAGE, MAX_PEAK
    from scripts.aggrescan_gate import screen as a3d_screen
    from scripts.tnp_developability import parse

    log = tnp_dir / "tnp_stdout.log"
    if not log.exists():
        raise SystemExit(f"{log} is missing; point --tnp-dir at a finished TNP run")
    fasta = (tnp_dir / "s.fasta").read_text().split()
    names = dict(zip([n.lstrip(">") for n in fasta[0::2]], fasta[1::2]))
    output = tnp_dir / "TNP_output"
    tnp = parse(log.read_text(), output)
    tnp = tnp[tnp["dssp"]] if "dssp" in tnp else tnp
    by_sequence = {names[n]: row for n, row in zip(tnp["name"], tnp.to_dict("records"))
                   if n in names}

    wanted = {n: s for n, s in names.items() if s in set(frame["sequence"])}
    _, values = a3d_screen(wanted, output / "Final_Models",
                           tnp_dir / "a3d_why", with_values=True)
    out = []
    for sequence in frame["sequence"]:
        row = by_sequence.get(sequence, {})
        average, peak = values.get(sequence, (float("nan"), float("nan")))
        out.append({
            "psh": row.get("PSH", float("nan")),
            "tnp_reds": row.get("reds", float("nan")),
            "tnp_red_on": row.get("red_on", ""),
            "a3d_average": average, "a3d_peak": peak,
            "a3d_over": (average > MAX_AVERAGE) or (peak > MAX_PEAK)
            if average == average else False,
        })
    return pd.DataFrame(out, index=frame.index)


def report(frame: pd.DataFrame) -> None:
    accepted = frame[frame["accepted"]]
    refused = frame[~frame["accepted"]]
    print(f"\n{len(frame)} submissions whose verdict we know: "
          f"{len(accepted)} accepted, {len(refused)} refused\n")

    gates = {
        "deterministic (offline checker)": frame["deterministic"].ne(""),
        f"PSH > {THRESHOLD:.0f}": frame["psh"] > THRESHOLD,
        "TNP red flag": frame["tnp_reds"].fillna(0) > 0,
        "Aggrescan3D past a cap": frame["a3d_over"].fillna(False),
    }
    print(f"{'gate':34} {'on refused':>12} {'on accepted':>12}   {'reading':<30}")
    for name, fires in gates.items():
        on_refused = fires[~frame["accepted"]]
        on_accepted = fires[frame["accepted"]]
        if on_accepted.mean() > 0.15:
            verdict = "FIRES ON BOTH: ours is stricter"
        elif on_refused.sum() == 0:
            verdict = "silent: not where they refuse"
        elif on_refused.mean() > 2 * max(on_accepted.mean(), 0.01):
            verdict = "a real rejection reason"
        else:
            verdict = "no separation"
        print(f"{name:34} {on_refused.sum():>5}/{len(on_refused):<6} "
              f"{on_accepted.sum():>5}/{len(on_accepted):<6}   {verdict:<30}")

    any_gate = pd.concat(gates.values(), axis=1).any(axis=1)
    explained = any_gate[~frame["accepted"]]
    print(f"\nrefusals ANY of our gates explains: {explained.sum()}/{len(explained)} "
          f"({explained.mean():.0%})")
    print(f"accepted submissions we would wrongly refuse: "
          f"{any_gate[frame['accepted']].sum()}/{len(accepted)} "
          f"({any_gate[frame['accepted']].mean():.0%})")

    unexplained = refused[~any_gate[~frame["accepted"]]]
    if len(unexplained):
        print(f"\nthe {len(unexplained)} refusals nothing we have explains:")
        print(f"  PSH median {unexplained['psh'].median():.1f} against "
              f"{refused['psh'].median():.1f} for all refused and "
              f"{accepted['psh'].median():.1f} for accepted")
        print(f"  per epoch: {dict(Counter(unexplained['epoch']).most_common())}")
        print("  Those are gates we do not reproduce at all. The validator also "
              "drops a UID when ANY of its sequences fails, and submits are "
              "per-hotkey, so some of these may be collateral.")

    print("\nwhat our own candidate would face, by PSH reading:")
    for value in (84.0, 87.0, 88.5, 91.0):
        print(f"  PSH {value:.1f} over 3 folds -> "
              f"{refusal_probability(value, 3):.0%} chance of refusal on PSH alone")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default=NANOBODY_TARGET[0])
    parser.add_argument("--tnp-dir", type=Path, required=True,
                        help="a finished TNP run (tnp_stdout.log + s.fasta + TNP_output)")
    parser.add_argument("--out", type=Path, help="write the full table here")
    args = parser.parse_args()

    frame = labelled(args.target)
    params = epoch_params([int(e) for e in sorted(frame["epoch"].unique())])
    archive = NanobodyArchive(args.target, refresh=False)

    # Sequences submitted at or after the epoch under test are in today's
    # archive but were not in it then, so they cannot count as duplicates.
    by_epoch: dict[int, set] = {}
    for row in frame.itertuples():
        by_epoch.setdefault(row.epoch, set()).add(row.sequence)

    problems = []
    for row in frame.itertuples():
        later = {s for e, seqs in by_epoch.items() if e >= row.epoch for s in seqs}
        found = check_sequence(row.sequence, params[row.epoch]["allowed_v"],
                               params[row.epoch]["allowed_cdr3_length"],
                               archive, later, warnings=[])
        problems.append("; ".join(found))
    frame["deterministic"] = problems
    frame = pd.concat([frame, structural(frame, args.tnp_dir)], axis=1)

    report(frame)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(args.out, index=False)
        print(f"\nfull table saved to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
