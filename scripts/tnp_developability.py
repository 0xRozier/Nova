"""Does our output pass the developability gate the validator applies first?

The validator's nanobody filter calls TNP (the Therapeutic Nanobody Profiler,
OPIG Oxford, BSD-3-Clause) and rejects a sequence on any red flag. Its six
metrics are the total CDR length, the CDR3 length, the CDR3 compactness, and
the CDR patches of surface hydrophobicity, positive charge and negative charge.
Our offline checker only covers the first two; the other four need a predicted
structure, which TNP builds on CPU with NanoBodyBuilder2.

Checked against the validator's own config: TNP's red thresholds match it
exactly (C 0.5634/1.6132, PSH 73.4039/155.4699, PPC 1.1853, PNC 1.8838), so
this is the real gate rather than an approximation.

The batch carries the real submissions as a positive control. Those were
accepted by the validators, so if TNP flags them red, this setup is wrong and
the run says nothing about our candidates.

Usage:
    python -m scripts.tnp_developability --fasta data/tnp/sequences.fasta
    python -m scripts.tnp_developability --parse data/tnp/tnp_stdout.log   # re-read a finished run
"""
import argparse
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

from nova_lab import DATA_DIR, toolchain
from nova_lab.config import PSH_STRICT_REJECTION_THRESHOLD as PSH_STRICT

TNP = toolchain.TNP
METRICS = ["L", "L3", "C", "PSH", "PPC", "PNC"]
LABELS = {
    "Total IMGT CDR Length": "L",
    "IMGT CDR3 Length": "L3",
    "CDR Compactness": "C",
    "Patch CDR Surface Hydrophobicity score is": "PSH",
    "Patch CDR Positive Charge score is": "PPC",
    "Patch CDR Negative Charge score is": "PNC",
}
HEAD = re.compile(r"Summary Statistics for (\S+?):")
LINE = re.compile(r"^\s{4}(.+?): ([-\d.eE+]+) \((\w+) flag\)$")


def parse(text: str, out_dir: Path | None = None) -> pd.DataFrame:
    """One row per sequence, with each metric's value and flag.

    With `out_dir` (the directory given to --output), each row also carries
    `dssp`: whether TNP got the DSSP annotation for that sequence. A row with
    dssp False holds a PSH computed on a different code path from the
    validator's -- 89.50 against 91.21 on the one candidate measured both ways
    -- so it is a reading of something else, not a noisy reading of theirs.
    """
    rows, current = [], None
    for line in text.splitlines():
        if m := HEAD.search(line):
            current = {"name": m.group(1)}
            rows.append(current)
            continue
        if current is None:
            continue
        if m := LINE.match(line.rstrip()):
            label, value, flag = m.group(1), m.group(2), m.group(3).lower()
            key = LABELS.get(label.strip())
            if key:
                current[key] = float(value)
                current[f"{key}_flag"] = flag
    if not rows:
        # A run can profile nothing at all -- TNP exiting on an import error
        # profiles 0 of 100 -- and callers then got a KeyError on 'name' three
        # frames away from the cause. An empty run is a legitimate answer
        # ("nothing was screened"), so it comes back shaped like a result.
        empty = pd.DataFrame(columns=["name", *METRICS, *(f"{m}_flag" for m in METRICS),
                                      "reds", "greys", "rejected", "red_on"])
        return empty.assign(dssp=pd.Series(dtype=bool)) if out_dir is not None else empty
    df = pd.DataFrame(rows).drop_duplicates("name", keep="last")
    flags = [f"{m}_flag" for m in METRICS if f"{m}_flag" in df.columns]
    df["reds"] = df[flags].eq("red").sum(axis=1)
    df["greys"] = df[flags].eq("grey").sum(axis=1)
    df["rejected"] = df["reds"] > 0
    df["red_on"] = [", ".join(m for m in METRICS
                              if r.get(f"{m}_flag") == "red") for _, r in df.iterrows()]
    if out_dir is not None:
        root = out_dir / "TNP_output" if (out_dir / "TNP_output").is_dir() else out_dir
        df["dssp"] = [toolchain.annotated(root, n) for n in df["name"]]
    return df


def run(fasta: Path, ncores: int, out_dir: Path) -> str:
    if not TNP.exists():
        sys.exit(f"TNP not installed at {TNP}")
    out_dir.mkdir(parents=True, exist_ok=True)
    # TNP resolves --file at parse time and then chdirs into --output, so a
    # relative path would be looked up from the wrong directory — and its bare
    # `except` reports that as "not FASTA format".
    fasta = fasta.resolve()
    cmd = [str(TNP), "--file", str(fasta), "--ncores", str(ncores),
           "--output", str(out_dir / "TNP_output")]
    # TNP shells out to a bare `NanoBodyBuilder2` and Biopython to a bare
    # `mkdssp`, so both have to be on PATH; toolchain.tnp_env puts them there.
    print(" ".join(cmd), flush=True)
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=out_dir,
                          env=toolchain.tnp_env())
    text = proc.stdout + "\n" + proc.stderr
    (out_dir / "tnp_stdout.log").write_text(text)
    print(f"TNP exit {proc.returncode}; log -> {out_dir / 'tnp_stdout.log'}")
    audit(text, out_dir / "TNP_output", wanted=fasta.read_text().count(">"))
    return text


def audit(text: str, output: Path, wanted: int) -> None:
    """Say out loud what TNP swallowed, because it carries on either way.

    TNP reports a PSH whether or not DSSP annotated the model, and reports
    nothing at all when NanoBodyBuilder2 fails -- a sequence then simply has no
    row, which reads as "not screened" and must never read as "passed". Both
    happened here: DSSP 4.x refuses NanoBodyBuilder2 models outright, and
    ImmuneBuilder 1.2's set/dict typo lost 3 of 37 sequences. Eleven production
    submissions were refused on PSH while the screen looked like it was working.
    """
    profiled = text.count("Summary Statistics for ")
    if profiled < wanted:
        print(f"  WARNING: {wanted - profiled} of {wanted} sequences came back with "
              f"no model, so they have no PSH at all. Run scripts.doctor: "
              f"ImmuneBuilder's refine.py typo loses sequences this way.")
    if "DSSP failed to produce an output" in text:
        names = [p.name for p in (output / "Raw_Model_Outputs").glob("*")
                 if p.is_dir()] if (output / "Raw_Model_Outputs").is_dir() else []
        blind = [n for n in names if not toolchain.annotated(output, n)]
        print(f"  WARNING: DSSP failed on {len(blind) or 'some'} of {wanted} "
              f"sequences. Their PSH is NOT the validator's measurement "
              f"(89.50 against 91.21 on the one candidate measured both ways). "
              f"Run scripts.doctor.")


def report(df: pd.DataFrame, meta: pd.DataFrame) -> None:
    m = meta.merge(df, on="name", how="left")
    scored = m[m["reds"].notna()]
    print(f"\n{len(scored)} of {len(m)} sequences profiled by TNP")
    if scored.empty:
        return

    ctrl = scored[scored["group"] == "real"]
    ours = scored[scored["group"] == "candidate"]
    print("\n=== Positive control: real submissions the validators accepted ===")
    print(f"{len(ctrl)} profiled, {int(ctrl['rejected'].sum())} would be rejected "
          f"({ctrl['rejected'].mean():.0%})")
    if ctrl["rejected"].mean() > 0.2:
        print("  WARNING: the control fails. Either this setup is wrong, or the validator")
        print("  does not gate on TNP red flags the way the config suggests. Read nothing")
        print("  into our candidates' numbers until this is resolved.")

    print("\n=== Our candidates ===")
    print(f"{len(ours)} profiled, {int(ours['rejected'].sum())} would be rejected "
          f"({ours['rejected'].mean():.0%})")

    print(f"\nstrict PSH rejection (validator: PSH > {PSH_STRICT} is refused, "
          f"whatever TNP's own flag says):")
    for group, sub in scored.groupby("group"):
        over = sub["PSH"] > PSH_STRICT
        print(f"  {group:<10} {int(over.sum())}/{len(sub)} over {PSH_STRICT} "
              f"({over.mean():.0%}), median {sub['PSH'].median():.1f}")
    print("  caution: our PSH is not calibrated to the validator's — 62% of the real")
    print("  submissions it accepted are over 89 here, so read the comparison between")
    print("  groups, not the absolute value.")

    print("\nred flags by metric (share of each group):")
    table = pd.DataFrame({
        g: [sub[f"{k}_flag"].eq("red").mean() for k in METRICS]
        for g, sub in scored.groupby("group")}, index=METRICS)
    print(table.to_string(float_format=lambda v: f"{v:.1%}"))

    print("\nmedian value by metric:")
    print(pd.DataFrame({g: [sub[k].median() for k in METRICS]
                        for g, sub in scored.groupby("group")},
                       index=METRICS).to_string(float_format=lambda v: f"{v:.2f}"))

    out = DATA_DIR / "tnp" / "tnp_results.csv"
    m.to_csv(out, index=False)
    print(f"\nfull table saved to {out}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=Path, default=DATA_DIR / "tnp" / "sequences.fasta")
    parser.add_argument("--meta", type=Path, default=DATA_DIR / "tnp" / "meta.csv")
    parser.add_argument("--ncores", type=int, default=4)
    parser.add_argument("--parse", type=Path, help="parse an existing TNP log instead of running")
    args = parser.parse_args()

    text = args.parse.read_text() if args.parse else run(args.fasta, args.ncores,
                                                         DATA_DIR / "tnp")
    report(parse(text), pd.read_csv(args.meta))


if __name__ == "__main__":
    main()
