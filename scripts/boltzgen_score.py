"""Score nanobody sequences with BoltzGen exactly like the NOVA validator.

Reproduces external_tools/boltzgen/src/boltzgen/boltzgen_wrapper.py and
config/boltzgen_config.yaml at UPSTREAM_COMMIT: the same input YAML (target
structure chain, residue range, binding site, MSA; the design as entity B),
the same protocol and steps, the same step/noise scales.

BoltzGen is the nova-vendored build, not the PyPI one. It is not installed here:
run this where it is (the Kaggle notebook from scripts.build_boltzgen_notebook
installs it), and point --boltzgen at that interpreter's entry point.

Usage:
    python -m scripts.boltzgen_score --file inputs/kaggle_nanobody_batch.csv
    python -m scripts.boltzgen_score --file seqs.csv --run batch1 --nova-dir /kaggle/temp/nova
"""
import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

from nova_lab import DATA_DIR

# config/config.yaml protein_selection.nanobody + config/boltzgen_config.yaml
TARGET = "Q9NZQ7"
STRUCTURE = "8AOK"
STRUCTURE_CHAIN = "A"
STRUCTURE_RES_INDEX = "1..115"
BINDING_SITE = "37,39,41..44,46,49,51,56,59,94..96,98,100,104..108"
PROTOCOL = "nanobody-anything"
NUM_DESIGNS = 1
STEP_SCALE = 2.0
NOISE_SCALE = 0.88
EXECUTE_STEPS = ["design", "folding", "design_folding", "analysis"]
METRICS = ["design_iiptm", "design_ptm", "design_to_target_iptm", "min_design_to_target_pae",
           "interaction_pae", "plip_hbonds_refolded", "plip_saltbridge_refolded",
           "delta_sasa_refolded", "liability_score", "liability_num_violations"]
METRIC_MODES = {"design_iiptm": "max", "design_ptm": "max", "design_to_target_iptm": "max",
                "min_design_to_target_pae": "min", "interaction_pae": "min",
                "plip_hbonds_refolded": "max", "plip_saltbridge_refolded": "max",
                "delta_sasa_refolded": "max", "liability_score": "min",
                "liability_num_violations": "min"}


def subprocess_env() -> dict:
    """BoltzGen imports matplotlib through biotite. A notebook host exports
    MPLBACKEND=module://matplotlib_inline.backend_inline, which our separate
    interpreter cannot provide, so matplotlib raises on import: pin a headless
    backend instead."""
    return {**os.environ, "MPLBACKEND": "Agg"}


def seq_hash(seq: str) -> str:
    return hashlib.sha256(seq.encode("ascii")).hexdigest()


def input_yaml(sequence: str, structure: Path, msa: Path) -> str:
    return f"""entities:
- file:
    path: {structure}
    include:
        - chain:
            id: {STRUCTURE_CHAIN}
            res_index: {STRUCTURE_RES_INDEX}
            msa: {msa}
    binding_types:
        - chain:
            id: {STRUCTURE_CHAIN}
            binding: {BINDING_SITE}
- protein:
    id: B
    sequence: "{sequence}"
    msa: empty
"""


def write_inputs(sequences: list[str], in_dir: Path, structure: Path, msa: Path) -> dict:
    in_dir.mkdir(parents=True, exist_ok=True)
    records = {}
    for seq in sequences:
        rid = seq_hash(seq)
        (in_dir / f"{rid}_{TARGET}_input.yaml").write_text(input_yaml(seq, structure, msa))
        records[rid] = seq
    return records


def rank_sum(frame: pd.DataFrame) -> pd.Series:
    ranks = pd.DataFrame({m: frame[m].rank(method="dense", ascending=(mode == "min"))
                          for m, mode in METRIC_MODES.items()})
    return ranks.sum(axis=1)


def read_metrics(out_dir: Path, records: dict) -> pd.DataFrame:
    path = out_dir / "intermediate_designs" / "aggregate_metrics_analyze.csv"
    if not path.exists():
        print(f"no metrics file at {path}")
        return pd.DataFrame()
    df = pd.read_csv(path)
    missing = [m for m in METRICS if m not in df.columns]
    if missing:
        print(f"metrics missing from BoltzGen output: {missing}")
    df["record"] = df["id"].astype(str).str.split("_").str[0]
    df["sequence"] = df["record"].map(records)
    keep = ["sequence", "record"] + [m for m in METRICS if m in df.columns]
    return df[keep].dropna(subset=["sequence"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", required=True, help="CSV with a 'sequence' column")
    parser.add_argument("--run", default=None)
    parser.add_argument("--nova-dir", type=Path, default=Path("/kaggle/temp/nova"),
                        help="where the vendored boltzgen, 8AOK.cif and the MSA live")
    parser.add_argument("--boltzgen", default=shutil.which("boltzgen") or "boltzgen")
    parser.add_argument("--top", type=int)
    args = parser.parse_args()

    meta = pd.read_csv(args.file)
    if args.top:
        meta = meta.head(args.top)
    sequences = list(dict.fromkeys(meta["sequence"]))
    run = args.run or time.strftime("run_%Y%m%d_%H%M%S")
    base = DATA_DIR / "boltzgen" / run

    structure = args.nova_dir / "data" / "structures" / f"{STRUCTURE}.cif"
    msa = args.nova_dir / "data" / "msa_files" / f"{TARGET}.a3m"
    for path in (structure, msa):
        if not path.exists():
            sys.exit(f"missing {path} — see scripts.build_boltzgen_notebook")

    records = write_inputs(sequences, base / "inputs", structure, msa)
    print(f"{len(records)} sequences to score, inputs in {base / 'inputs'}")

    t = time.time()
    configure = [args.boltzgen, "configure", str(base / "inputs"), "--output", str(base / "out"),
                 "--protocol", PROTOCOL, "--num_designs", str(NUM_DESIGNS),
                 "--skip_inverse_folding", "--step_scale", str(STEP_SCALE),
                 "--noise_scale", str(NOISE_SCALE)]
    env = subprocess_env()
    subprocess.run(configure, check=True, env=env)
    subprocess.run([args.boltzgen, "execute", str(base / "out"), "--steps", *EXECUTE_STEPS],
                   check=True, env=env)
    elapsed = time.time() - t

    res = read_metrics(base / "out", records)
    if res.empty:
        sys.exit("BoltzGen produced no metrics")
    res = res.merge(meta.drop_duplicates("sequence"), on="sequence", how="left")
    res["rank_sum"] = rank_sum(res)
    res = res.sort_values("rank_sum")
    res.to_csv(base / "scores.csv", index=False)

    cols = [c for c in ["group", "rank_sum", "true_rank_sum", "predicted_rank_sum",
                        "design_iiptm", "interaction_pae"] if c in res]
    with pd.option_context("display.width", 160, "display.max_rows", 200):
        print(res[cols].to_string(index=False))
    print(f"\n{len(res)} sequences in {elapsed / 60:.1f} min "
          f"({elapsed / max(len(res), 1) / 60:.1f} min/sequence)")

    if "group" in res and {"candidate", "calibration"} <= set(res["group"].dropna()):
        ours = res[res["group"] == "candidate"]["rank_sum"]
        theirs = res[res["group"] == "calibration"]["rank_sum"]
        print(f"our candidates rank {ours.min():.0f}-{ours.max():.0f}; "
              f"the real submissions {theirs.min():.0f}-{theirs.max():.0f} "
              f"(lower is better, ranked in this pool)")
    if "true_rank_sum" in res:
        known = res.dropna(subset=["true_rank_sum"])
        if len(known) >= 3:
            print(f"calibration: Spearman {known['rank_sum'].corr(known['true_rank_sum'], method='spearman'):.2f} "
                  f"against the validators' own ranking of these sequences")
    print(f"saved to {base / 'scores.csv'}")


if __name__ == "__main__":
    main()
