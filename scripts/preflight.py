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
    production submissions. On the install it runs (nova_lab.toolchain checks
    it: TNP at 29dcac7, the filter's process_pdb patch, mkdssp 3.1.4) our
    reading IS theirs: of those 36 accepted submissions, folded three times
    each, none reads above 89 and the median is 83.5. What remains is the
    stochastic refinement, so one measurement deviates by 1.74 and the
    validator takes a single draw of its own; this script therefore reports a
    REJECTION PROBABILITY rather than a verdict. That probability is a cost to
    weigh against a candidate's quality, not a line to refuse at: see
    nova_lab.psh_calibration.

    Earlier versions of this file said the reproduction was impossible and put
    the deviation at 4.01. Both came from an install whose DSSP step was
    failing silently, where 44-62 % of accepted submissions read above 89.

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
from nova_lab import developability, psh_calibration
from nova_lab.nanobodies import NanobodyArchive, check_sequence

PSH_SIGMA = psh_calibration.SIGMA  # one measurement, pooled over 37 x 3 folds
A3D_MAX_AVERAGE = -0.68002724866200
A3D_MAX_PEAK = 1.58023066762925
# Aggrescan3D rests on the same predicted structure as PSH, so it carries the
# same kind of noise. Measured on 37 sequences folded three times each:
#
#   average  deviation 0.0087
#   peak     deviation 0.0211 median, 0.1248 at the 90th percentile, 0.2959 max
#
# The peak is a maximum over residues, so it is far noisier and heavy-tailed; a
# first estimate from three folds of ONE sequence gave 0.0019, eleven times too
# small, and the calibration duly refused an accepted submission whose peak
# breached by 0.084 — well inside the real noise. Both statistics are therefore
# averaged over the folds, which divides the deviation by root k, and a breach
# counts as a refusal only beyond two deviations of that mean.
A3D_SIGMA_AVERAGE = 0.0087
A3D_SIGMA_PEAK = 0.0211


def rejection_probability(mean: float, folds: int) -> float:
    """P(the validator's single PSH draw exceeds 89), given our k-fold mean.

    One line, kept for the callers already written against it; the model and
    the deviation it uses live in nova_lab.psh_calibration, which is where the
    calibration evidence is documented.
    """
    return psh_calibration.refusal_probability(mean, folds)


def structural(sequences: list[str], work: Path, folds: int) -> pd.DataFrame:
    """TNP flags, PSH over k folds, and the Aggrescan3D caps."""
    from scripts.aggrescan_gate import screen as a3d_screen
    from scripts.tnp_developability import METRICS as TNP_METRICS
    from scripts.tnp_developability import TNP, parse, run

    if not TNP.exists():
        print(f"TNP missing at {TNP}: structural gates not checked", file=sys.stderr)
        return pd.DataFrame({"sequence": sequences})

    names = {f"s{i}": s for i, s in enumerate(sequences)}
    rows = {s: {"sequence": s, "psh_folds": []} for s in sequences}
    # Every TNP metric is averaged over the folds, not just PSH. Taking the
    # others from fold 1 refused an accepted submission on a PPC of 1.2858
    # whose three folds mean 1.1687 and pass -- PPC deviates by 0.037 and up to
    # 0.23, so one reading next to a bound says nothing (FINDINGS 45).
    per_metric: dict[str, dict[str, list[float]]] = {s: {} for s in sequences}
    for fold in range(1, folds + 1):
        sub = work / f"fold{fold}"
        sub.mkdir(parents=True, exist_ok=True)
        (sub / "seqs.fasta").write_text("".join(f">{n}\n{s}\n" for n, s in names.items()))
        df = parse(run(sub / "seqs.fasta", 4, sub), sub)
        if "dssp" in df.columns:
            df = df[df["dssp"]]
        for row in df.to_dict("records"):
            sequence = names.get(row["name"])
            if sequence is None:
                continue
            rows[sequence]["psh_folds"].append(float(row["PSH"]))
            for metric in TNP_METRICS:
                if metric in row and row[metric] == row[metric]:
                    per_metric[sequence].setdefault(metric, []).append(float(row[metric]))

    for sequence, samples in per_metric.items():
        means = {m: sum(v) / len(v) for m, v in samples.items() if v}
        n_folds = min((len(v) for v in samples.values()), default=0) or 1
        flags = developability.breaches(means, n_folds)
        rows[sequence]["tnp_reds"] = len(flags)
        rows[sequence]["tnp_red_on"] = "; ".join(flags)
        rows[sequence]["tnp_borderline"] = "; ".join(
            developability.borderline(means, n_folds))
        for metric, value in means.items():
            rows[sequence][f"tnp_{metric}"] = value

    # Score every fold's structures and average, so the peak's noise shrinks.
    per_fold: dict[str, list[tuple[float, float]]] = {s: [] for s in sequences}
    for fold in range(1, folds + 1):
        models = work / f"fold{fold}" / "TNP_output" / "Final_Models"
        if not models.exists():
            continue
        _, values = a3d_screen(names, models, work / f"a3d{fold}", with_values=True)
        for seq, pair in values.items():
            per_fold[seq].append(pair)

    band_avg = 2 * A3D_SIGMA_AVERAGE / math.sqrt(max(folds, 1))
    band_peak = 2 * A3D_SIGMA_PEAK / math.sqrt(max(folds, 1))
    for s in sequences:
        seen = per_fold.get(s) or []
        if not seen:
            rows[s]["a3d_average"] = rows[s]["a3d_peak"] = float("nan")
            rows[s]["a3d_refused"] = rows[s]["a3d_borderline"] = False
            rows[s]["a3d_folds"] = 0
            continue
        avg = sum(a for a, _ in seen) / len(seen)
        peak = sum(p for _, p in seen) / len(seen)
        rows[s]["a3d_average"], rows[s]["a3d_peak"] = avg, peak
        rows[s]["a3d_folds"] = len(seen)
        over_avg = avg - A3D_MAX_AVERAGE
        over_peak = peak - A3D_MAX_PEAK
        rows[s]["a3d_refused"] = (over_avg > band_avg) or (over_peak > band_peak)
        rows[s]["a3d_borderline"] = (not rows[s]["a3d_refused"]) and (
            over_avg > 0 or over_peak > 0)
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
    # Without TNP, structural() returns the sequence column alone and said so,
    # so these columns are absent rather than NaN -- and out.get(name, default)
    # then hands back the scalar default, which has no .fillna. The script used
    # to die here instead of reporting the deterministic gates it can still
    # check on a machine with no TNP installed.
    def flags(name: str, default):
        if name in out:
            return out[name].fillna(default)
        return pd.Series(default, index=out.index)

    out["refused"] = (out["deterministic"].ne("")
                      | flags("tnp_reds", 0).gt(0)
                      | flags("a3d_refused", False).astype(bool))
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
        if getattr(r, "a3d_refused", False):
            bits.append(f"Aggrescan3D past the cap by more than two deviations "
                        f"(avg {r.a3d_average:.4f}, peak {r.a3d_peak:.4f})")
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
        if getattr(r, "a3d_borderline", False):
            print(f"    note: Aggrescan3D sits on its cap (avg {r.a3d_average:.4f}, "
                  f"peak {r.a3d_peak:.4f}) — within the noise of the limit")
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
                                      else "Aggrescan3D past the cap")
            print(f"    {r.sequence[:34]}... {why}")

    if "psh_mean" in out:
        a = out[out["truth"] == "accepted"]["psh_mean"].dropna()
        over = int((a > PSH_STRICT).sum())
        print(f"\nPSH on the accepted set (the validators took all of them, so "
              f"every one of them read at or under {PSH_STRICT} on THEIR "
              f"install): {len(a)} measured, median {a.median():.1f}, "
              f"max {a.max():.1f}, {over} above {PSH_STRICT} here")
        # What a faithful install reads is "a median near 83.5 and a handful at
        # the top", not "none above 89". The deviation of a 3-fold mean is 1.0,
        # so a submission truly sitting at 88.9 -- which is where the top of
        # this set sits -- reads above 89 half the time. An earlier version of
        # this message demanded 0 and called anything else a drift; that was
        # calibrated on one lucky run, and the next run of the same 36
        # sequences duly produced a 1 and an alarm that meant nothing.
        share = over / max(len(a), 1)
        if share > 0.14:
            print(f"  {over} of {len(a)} is too many: on this scale the "
                  f"accepted field reads a median of "
                  f"{psh_calibration.ACCEPTED_MEDIAN}, so our scale has drifted "
                  f"off theirs and the gate is screening something else. Run "
                  f"python -m scripts.doctor.")
        else:
            print(f"  {over} of {len(a)} is what the fold noise predicts "
                  f"(deviation of a 3-fold mean is "
                  f"{psh_calibration.SIGMA / 3 ** 0.5:.1f}, and the top of this "
                  f"set sits at 88.9). The median is the statistic to read: "
                  f"{a.median():.1f} against "
                  f"{psh_calibration.ACCEPTED_MEDIAN} measured.")
        for r in out[out["truth"] == "refused"].itertuples():
            if getattr(r, "psh_mean", float("nan")) == getattr(r, "psh_mean", float("nan")):
                print(f"\nthe sequence they refused: PSH {r.psh_mean:.1f} "
                      f"-> {r.psh_reject_p:.0%} chance of refusal "
                      f"(their log for it says 89.57)")

    path = DATA_DIR / "preflight" / "calibration.csv"
    path.parent.mkdir(parents=True, exist_ok=True)  # absent on a fresh machine
    out.to_csv(path, index=False)
    print(f"\nsaved to {path}")
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
