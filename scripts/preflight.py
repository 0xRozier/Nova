"""Would the validator accept this sequence? Every gate, in its order.

Each gate below is reproduced from the validator's own code, and the point of
this script is that the reproduction is *checked against ground truth* rather
than against intuition: `--calibrate` runs it on the 36 submissions the
validators accepted in epoch 26174 plus the one of ours they refused, and
reports whether the verdicts come out right. Every wrong conclusion in this
project came from a number nobody had checked that way.

The gates split in two:

  deterministic — length, amino acids, homopolymers, di-repeats, signal
    peptide, germline V/J identity, the epoch's V gene, VHH hallmarks,
    canonical cysteines, CDR3 length, archive uniqueness and similarity, TNP's
    six flags and the Aggrescan3D caps. These we reproduce exactly, and a
    failure here is a refusal.

  PSH — the validator refuses PSH > 89, and that is what rejected all 11 of our
    production submissions. It cannot be reproduced exactly: PSH is computed on
    a predicted structure, NanoBodyBuilder2 is not deterministic, and one
    measurement deviates by 4.01 (ten sequences measured twice). On our scale
    16 of those 36 accepted submissions read above 89, so a hard threshold here
    would refuse sequences the validator takes. This script therefore reports a
    REJECTION PROBABILITY rather than a verdict.

Usage:
    python -m scripts.preflight --calibrate
    python -m scripts.preflight --file candidates.txt --epoch 26174
    python -m scripts.preflight --file c.txt --allowed-v IGHV3-74 --cdr3 8 --folds 3
"""
import argparse
import math
import sys
from pathlib import Path

import pandas as pd

from nova_lab import DATA_DIR
from nova_lab.config import NANOBODY_TARGET
from nova_lab.config import PSH_STRICT_REJECTION_THRESHOLD as PSH_STRICT
from nova_lab.nanobodies import NanobodyArchive, check_sequence

PSH_SIGMA = 4.01  # one measurement, from 10 sequences measured twice
A3D_MAX_AVERAGE = -0.68002724866200
A3D_MAX_PEAK = 1.58023066762925


def rejection_probability(mean: float, folds: int) -> float:
    """P(the validator's single PSH draw exceeds 89), given our k-fold mean.

    The validator takes one draw with deviation PSH_SIGMA; our estimate of the
    true value carries PSH_SIGMA/sqrt(k). The two combine in quadrature.
    """
    spread = PSH_SIGMA * math.sqrt(1 + 1 / max(folds, 1))
    return 0.5 * (1 - math.erf((PSH_STRICT - mean) / (spread * math.sqrt(2))))


def structural(sequences: list[str], work: Path, folds: int) -> pd.DataFrame:
    """TNP flags, PSH over k folds, and the Aggrescan3D caps."""
    from scripts.aggrescan_gate import screen as a3d_screen
    from scripts.tnp_developability import TNP, parse, run

    if not TNP.exists():
        print(f"TNP missing at {TNP}: structural gates not checked", file=sys.stderr)
        return pd.DataFrame({"sequence": sequences})

    names = {f"s{i}": s for i, s in enumerate(sequences)}
    rows = {s: {"sequence": s, "psh_folds": []} for s in sequences}
    first = None
    for fold in range(1, folds + 1):
        sub = work / f"fold{fold}"
        sub.mkdir(parents=True, exist_ok=True)
        (sub / "seqs.fasta").write_text("".join(f">{n}\n{s}\n" for n, s in names.items()))
        df = parse(run(sub / "seqs.fasta", 4, sub))
        if first is None:
            first = df
        for n, psh in zip(df["name"], df["PSH"]):
            if n in names:
                rows[names[n]]["psh_folds"].append(float(psh))

    for n, red, greys in zip(first["name"], first["reds"], first["greys"]):
        if n in names:
            rows[names[n]]["tnp_reds"] = int(red)
            rows[names[n]]["tnp_greys"] = int(greys)
            rows[names[n]]["tnp_red_on"] = first.loc[first["name"] == n, "red_on"].iat[0]

    a3d = a3d_screen(names, work / "fold1" / "TNP_output" / "Final_Models", work / "a3d")
    for s in sequences:
        rows[s]["a3d_pass"] = a3d.get(s)
        folds_seen = rows[s]["psh_folds"]
        rows[s]["psh_mean"] = sum(folds_seen) / len(folds_seen) if folds_seen else float("nan")
        rows[s]["psh_n"] = len(folds_seen)
    return pd.DataFrame(list(rows.values()))


def run_preflight(sequences: list[str], allowed_v: str | None, cdr3_len: int | None,
                  folds: int, target: str, work: Path,
                  archive: NanobodyArchive | None) -> pd.DataFrame:
    rows = []
    for seq in sequences:
        warnings: list[str] = []
        problems = check_sequence(seq, allowed_v, cdr3_len, archive, warnings=warnings)
        rows.append({"sequence": seq, "deterministic": "; ".join(problems) or "",
                     "archive_warning": "; ".join(warnings)})
    out = pd.DataFrame(rows)
    struct = structural(sequences, work, folds)
    out = out.merge(struct, on="sequence", how="left")
    if "psh_mean" in out:
        out["psh_reject_p"] = [rejection_probability(m, int(n or 1)) if m == m else float("nan")
                               for m, n in zip(out["psh_mean"], out.get("psh_n", 1))]
    out["refused"] = out["deterministic"].ne("") | out.get("tnp_reds", 0).fillna(0).gt(0) \
        | out.get("a3d_pass", True).fillna(True).eq(False)
    return out


def report(out: pd.DataFrame) -> None:
    for r in out.itertuples():
        head = f"{r.sequence[:34]}... "
        if r.deterministic:
            print(f"{head}REFUSED — {r.deterministic}")
            continue
        bits = []
        if getattr(r, "tnp_reds", 0):
            bits.append(f"TNP red on {r.tnp_red_on}")
        if getattr(r, "a3d_pass", None) is False:
            bits.append("Aggrescan3D over cap")
        if bits:
            print(f"{head}REFUSED — {', '.join(bits)}")
            continue
        p = getattr(r, "psh_reject_p", float("nan"))
        m = getattr(r, "psh_mean", float("nan"))
        if m != m:
            print(f"{head}deterministic gates pass; PSH not measured")
        else:
            risk = "low" if p < 0.1 else ("even" if p < 0.5 else "high")
            print(f"{head}deterministic gates pass; PSH {m:.1f} over {int(r.psh_n)} folds "
                  f"-> {p:.0%} chance of refusal ({risk})")
        if r.archive_warning:
            print(f"    note: {r.archive_warning}")


def calibrate(folds: int, target: str) -> int:
    """Check the reproduction against known verdicts. Returns an exit code."""
    truth = DATA_DIR.parent / "results" / "tnp_epoch26174.csv"
    if not truth.exists():
        sys.exit(f"{truth} is missing; it carries the ground truth this calibrates on")
    t = pd.read_csv(truth)
    accepted = list(t[t["group"] == "accepte"]["sequence"])
    refused = list(t[t["group"] != "accepte"]["sequence"])
    print(f"calibrating on epoch 26174: {len(accepted)} submissions the validators "
          f"accepted, {len(refused)} they refused\n")

    work = DATA_DIR / "preflight" / "calibration"
    out = run_preflight(accepted + refused, "IGHV3-74", 8, folds, target, work, None)
    out["truth"] = ["accepted"] * len(accepted) + ["refused"] * len(refused)

    ok = out[out["truth"] == "accepted"]
    bad = int(ok["refused"].sum())
    print(f"deterministic gates on the accepted set: {len(ok) - bad}/{len(ok)} pass")
    if bad:
        print("  WRONGLY REFUSED — the reproduction is broken, not the sequences:")
        for r in ok[ok["refused"]].itertuples():
            why = r.deterministic or (f"TNP red on {r.tnp_red_on}" if getattr(r, "tnp_reds", 0)
                                      else "Aggrescan3D over cap")
            print(f"    {r.sequence[:34]}... {why}")

    if "psh_mean" in out:
        a = out[out["truth"] == "accepted"]["psh_mean"].dropna()
        print(f"\nPSH on the accepted set (all of them cleared the validator's 89): "
              f"median {a.median():.1f}, {(a > PSH_STRICT).mean():.0%} read above 89 here")
        print("  that share is why PSH is reported as a probability, not a verdict")
        for r in out[out["truth"] == "refused"].itertuples():
            if getattr(r, "psh_mean", float("nan")) == getattr(r, "psh_mean", float("nan")):
                print(f"\nthe sequence they refused: PSH {r.psh_mean:.1f} "
                      f"-> {r.psh_reject_p:.0%} chance of refusal")

    out.to_csv(DATA_DIR / "preflight" / "calibration.csv", index=False)
    print(f"\nsaved to {DATA_DIR / 'preflight' / 'calibration.csv'}")
    return 1 if bad else 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("sequences", nargs="*")
    parser.add_argument("--file", type=Path, help="one sequence per line")
    parser.add_argument("--epoch", type=int, help="resolve the gene and CDR3 length from the chain")
    parser.add_argument("--allowed-v")
    parser.add_argument("--cdr3", type=int)
    parser.add_argument("--target", default=NANOBODY_TARGET[0])
    parser.add_argument("--folds", type=int, default=3, help="structure folds averaged for PSH")
    parser.add_argument("--no-archive", action="store_true")
    parser.add_argument("--calibrate", action="store_true",
                        help="check the reproduction against known verdicts and stop")
    args = parser.parse_args()

    if args.calibrate:
        sys.exit(calibrate(args.folds, args.target))

    seqs = list(args.sequences)
    if args.file:
        seqs += [l.strip() for l in args.file.read_text().splitlines() if l.strip()]
    if not seqs:
        parser.error("give sequences as arguments, with --file, or use --calibrate")

    if args.allowed_v and args.cdr3:
        allowed_v, cdr3_len = args.allowed_v, args.cdr3
    elif args.epoch is not None:
        from scripts.check_nanobody import resolve_epoch
        p = resolve_epoch(args.epoch)
        allowed_v, cdr3_len = p["allowed_v"], p["allowed_cdr3_length"]
    else:
        parser.error("give --epoch, or --allowed-v with --cdr3")
    print(f"epoch gates: V {allowed_v}, CDR3 {cdr3_len}\n")

    archive = None if args.no_archive else NanobodyArchive(args.target)
    out = run_preflight(seqs, allowed_v, cdr3_len, args.folds, args.target,
                        DATA_DIR / "preflight" / "run", archive)
    report(out)
    path = DATA_DIR / "preflight" / "preflight.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)
    print(f"\nsaved to {path}")


if __name__ == "__main__":
    main()
