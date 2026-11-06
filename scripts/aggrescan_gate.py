"""The last validator gate we had not tested: Aggrescan3D.

The validator's developability stage requires BOTH the TNP criteria and the
Aggrescan3D caps (`aggrescan_rejection` in utils/aggrescan3d.py). A3D scores the
very structure TNP already built, so no refolding is needed:

    a3d_average <= -0.68002724866200
    a3d_maximum <=  1.58023066762925

Those cutoffs are NOVA policy, not physics (docs/usage/aggrescan3d-rust.md), and
the scorer is the Rust one vendored in the nanobody-filter repository, built with
its pinned FreeSASA 2.0.1.

Run scripts.tnp_developability first: this reuses its NanoBodyBuilder2 models.
As there, the real submissions are the positive control — the validators
accepted them, so red flags among them mean the setup is wrong.

Usage:
    python -m scripts.aggrescan_gate
"""
import argparse
import csv
import hashlib
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

from nova_lab import DATA_DIR

MAX_AVERAGE = -0.68002724866200
MAX_PEAK = 1.58023066762925
TOOLS = Path("/home/user/nanofilter/.tools/bin")


def screen(sequences: dict[str, str], models: Path, work: Path) -> dict[str, bool]:
    """True per sequence if Aggrescan3D passes both caps. Empty if it cannot run.

    `sequences` maps the TNP model name to the sequence, so the models
    scripts.tnp_developability already built are reused rather than refolded —
    which is also what the validator's adapter does.
    """
    scorer, freesasa = TOOLS / "aggrescan3d-rs", TOOLS / "freesasa"
    if not (scorer.exists() and freesasa.exists()):
        return {}
    structures = work / "structures"
    structures.mkdir(parents=True, exist_ok=True)
    records = {}
    for name, sequence in sequences.items():
        model = models / f"{name}_NanoBodyBuilder2_Model.pdb"
        if not model.exists():
            continue
        digest = hashlib.sha256(sequence.encode()).hexdigest()
        shutil.copyfile(model, structures / f"{digest}.pdb")
        records[digest] = sequence
    if not records:
        return {}
    fasta = work / "sequences.fasta"
    fasta.write_text("".join(f">{d}\n{s}\n" for d, s in records.items()))
    summary = work / "summary"
    proc = subprocess.run([str(scorer), "--fasta", str(fasta), "--structures", str(structures),
                           "--out", str(summary), "--residue-out", str(work / "residues"),
                           "--freesasa", str(freesasa)], capture_output=True, text=True)
    out = summary / "aggrescan3d.csv"
    if proc.returncode != 0 or not out.is_file():
        return {}
    verdict = {}
    for row in csv.DictReader(out.open(newline="")):
        try:
            avg, peak = float(row["a3d_average"]), float(row["a3d_maximum"])
        except (KeyError, ValueError):
            continue
        sequence = records.get(row.get("sequence_hash", ""))
        if sequence is not None:
            verdict[sequence] = avg <= MAX_AVERAGE and peak <= MAX_PEAK
    return verdict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--models", type=Path,
                        default=DATA_DIR / "tnp" / "TNP_output" / "Final_Models")
    parser.add_argument("--meta", type=Path, default=DATA_DIR / "tnp" / "meta.csv")
    parser.add_argument("--work", type=Path, default=DATA_DIR / "aggrescan")
    args = parser.parse_args()

    scorer, freesasa = TOOLS / "aggrescan3d-rs", TOOLS / "freesasa"
    for tool in (scorer, freesasa):
        if not tool.exists():
            sys.exit(f"{tool} is missing; build it with install/install-aggrescan3d.sh")

    meta = pd.read_csv(args.meta)
    structures = args.work / "structures"
    structures.mkdir(parents=True, exist_ok=True)
    records, missing = {}, 0
    for row in meta.itertuples():
        model = args.models / f"{row.name}_NanoBodyBuilder2_Model.pdb"
        if not model.exists():
            missing += 1
            continue
        # The wrapper names each structure after the sequence's SHA-256 and
        # gives the scorer a FASTA keyed the same way; that pairing is what
        # carries provenance, so reproduce it exactly.
        digest = hashlib.sha256(row.sequence.encode()).hexdigest()
        shutil.copyfile(model, structures / f"{digest}.pdb")
        records[digest] = row.sequence
    print(f"{len(records)} structures ready" + (f", {missing} models missing" if missing else ""))

    fasta = args.work / "sequences.fasta"
    fasta.write_text("".join(f">{d}\n{s}\n" for d, s in records.items()))
    summary = args.work / "summary"
    cmd = [str(scorer), "--fasta", str(fasta), "--structures", str(structures),
           "--out", str(summary), "--residue-out", str(args.work / "residues"),
           "--freesasa", str(freesasa)]
    print(" ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        print((proc.stdout + proc.stderr)[-1500:])
        sys.exit(f"aggrescan3d-rs failed with {proc.returncode}")

    out = summary / "aggrescan3d.csv"
    if not out.is_file():
        sys.exit(f"no summary at {out}")
    rows = list(csv.DictReader(out.open(newline="")))
    print(f"scored {len(rows)} structures; columns: {', '.join(rows[0])}")

    def number(row, *names):
        for n in names:
            if n in row and row[n] not in ("", None):
                try:
                    return float(row[n])
                except ValueError:
                    pass
        return float("nan")

    df = pd.DataFrame([{
        "digest": r.get("sequence_hash") or r.get("name") or r.get("id"),
        "a3d_average": number(r, "a3d_average", "average", "avg_score"),
        "a3d_maximum": number(r, "a3d_maximum", "maximum", "max_score"),
    } for r in rows])
    df["sequence"] = df["digest"].map(records)
    df["passed"] = (df["a3d_average"] <= MAX_AVERAGE) & (df["a3d_maximum"] <= MAX_PEAK)
    merged = meta.merge(df, on="sequence", how="inner")

    out_csv = args.work / "aggrescan_results.csv"
    merged.to_csv(out_csv, index=False)
    print(f"\ncaps: average <= {MAX_AVERAGE}, peak <= {MAX_PEAK}")
    for group, sub in merged.groupby("group"):
        rejected = int((~sub["passed"]).sum())
        label = "positive control (accepted by the validators)" if group == "real" else "our candidates"
        print(f"\n{label}: {rejected} of {len(sub)} rejected ({rejected / len(sub):.0%})")
        print(f"  average: median {sub['a3d_average'].median():.3f}, worst {sub['a3d_average'].max():.3f}")
        print(f"  peak:    median {sub['a3d_maximum'].median():.3f}, worst {sub['a3d_maximum'].max():.3f}")
        print(f"  fails on average {int((sub['a3d_average'] > MAX_AVERAGE).sum())}, "
              f"on peak {int((sub['a3d_maximum'] > MAX_PEAK).sum())}")
    ctrl = merged[merged["group"] == "real"]
    if len(ctrl) and (~ctrl["passed"]).mean() > 0.2:
        print("\nWARNING: the control fails, so read nothing into our candidates' numbers.")
    print(f"\nfull table saved to {out_csv}")


if __name__ == "__main__":
    main()
