"""Score molecules with Boltz-2 exactly like the NOVA validator (CPU or GPU).

Reproduces upstream external_tools/boltz/boltz_wrapper.py + config/boltz_config.yaml:
same input YAML (target clipped from UniProt, validator's MSA, ligand as
binder), same predict settings and seed. Needs Boltz 2.2.1 (Python <= 3.12):

    python3.12 -m venv .venv-boltz && .venv-boltz/bin/pip install boltz==2.2.1
    (CPU only: install torch from https://download.pytorch.org/whl/cpu first)

Usage:
    python -m scripts.boltz_score rxn:33:2661472:131680 rxn:41:2633028:2588908
    python -m scripts.boltz_score --file data/proposals/O15427_rxn41_shortlist.csv --top 5
    python -m scripts.boltz_score --file inputs/kaggle_batch.csv --accelerator gpu --extra=--no_kernels

Extra columns of --file (group, pred, published_score...) are kept in the output.
Re-running with the same --run resumes: Boltz skips records already predicted.
"""
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
import requests
import yaml

from nova_lab import DATA_DIR
from nova_lab.chem import boltz_score, heavy_atoms
from nova_lab.reactions import name_to_smiles

# protein_selection.small_molecule in upstream config.yaml
TARGETS = {"O15427": (1, 465)}
MSA_URL = "https://raw.githubusercontent.com/metanova-labs/nova/main/data/msa_files/{target}.a3m"
BASE_SEED = 68
PREDICT_ARGS = [  # config/boltz_config.yaml
    "--recycling_steps", "3",
    "--diffusion_samples", "1",
    "--sampling_steps", "100",
    "--sampling_steps_affinity", "100",
    "--diffusion_samples_affinity", "3",
    "--affinity_mw_correction",
    "--output_format", "mmcif",
    "--seed", str(BASE_SEED),
    "--num_workers", "0",
]


def record_id(smiles: str) -> int:
    h = hashlib.sha256(smiles.encode()).digest()
    return (int.from_bytes(h[:8], "little") ^ BASE_SEED) % (2**31 - 1)


def protein_sequence(target: str) -> str:
    """Clipped target sequence, cached in data/msa/<target>.seq (the Kaggle notebook ships it)."""
    cached = DATA_DIR / "msa" / f"{target}.seq"
    if cached.exists():
        return cached.read_text().strip()
    fasta = requests.get(f"https://rest.uniprot.org/uniprotkb/{target}.fasta", timeout=30).text
    seq = "".join(line.strip() for line in fasta.splitlines() if not line.startswith(">"))
    start, end = TARGETS[target]
    cached.parent.mkdir(parents=True, exist_ok=True)
    cached.write_text(seq[start:end])
    return seq[start:end]


def msa_path(target: str) -> Path:
    path = DATA_DIR / "msa" / f"{target}.a3m"
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        resp = requests.get(MSA_URL.format(target=target), timeout=60)
        resp.raise_for_status()
        path.write_bytes(resp.content)
    return path


def write_inputs(names: list[str], target: str, in_dir: Path) -> dict:
    in_dir.mkdir(parents=True, exist_ok=True)
    seq, msa = protein_sequence(target), str(msa_path(target))
    records = {}
    for name in names:
        smiles = name_to_smiles(name)
        if not smiles:
            print(f"skipping {name}: cannot resolve SMILES")
            continue
        rid = f"{record_id(smiles)}_{target}"
        doc = {
            "version": 1,
            "sequences": [
                {"protein": {"id": "A", "sequence": seq, "msa": msa}},
                {"ligand": {"id": "B", "smiles": smiles}},
            ],
            "properties": [{"affinity": {"binder": "B"}}],
        }
        (in_dir / f"{rid}.yaml").write_text(yaml.safe_dump(doc, sort_keys=False, default_flow_style=False))
        records[rid] = {"name": name, "smiles": smiles}
    return records


def read_outputs(results_dir: Path, records: dict) -> pd.DataFrame:
    rows = []
    for rid, info in records.items():
        pred_dir = results_dir / "predictions" / rid
        metrics = {}
        for f in pred_dir.glob("*.json"):
            if f.name.startswith(("affinity", "confidence")):
                metrics.update(json.loads(f.read_text()))
        if "affinity_probability_binary" not in metrics:
            print(f"no affinity output for {info['name']}")
            continue
        p, v = metrics["affinity_probability_binary"], metrics["affinity_pred_value"]
        rows.append({**info, "affinity_probability_binary": p, "affinity_pred_value": v,
                     "score": boltz_score(p, v, heavy_atoms(info["smiles"]))})
    return pd.DataFrame(rows, columns=["name", "smiles", "affinity_probability_binary", "affinity_pred_value", "score"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("names", nargs="*")
    parser.add_argument("--file", help="CSV with a 'name' column (e.g. a proposer shortlist)")
    parser.add_argument("--top", type=int, default=None)
    parser.add_argument("--target", default="O15427")
    parser.add_argument("--boltz", default=shutil.which("boltz") or str(Path(sys.executable).parent / "boltz"),
                        help="boltz CLI")
    parser.add_argument("--accelerator", default="cpu", choices=["cpu", "gpu"])
    parser.add_argument("--run", default=None, help="run name (default: timestamp)")
    parser.add_argument("--extra", default="", help="extra boltz predict flags, e.g. --extra=--no_kernels")
    args = parser.parse_args()

    names = list(args.names)
    meta = pd.DataFrame({"name": names})
    if args.file:
        meta = pd.concat([meta, pd.read_csv(args.file).head(args.top)], ignore_index=True)
        names = list(meta["name"])
    run = args.run or time.strftime("run_%Y%m%d_%H%M%S")
    base = DATA_DIR / "boltz" / run
    records = write_inputs(names, args.target, base / "inputs")

    t = time.time()
    subprocess.run([args.boltz, "predict", str(base / "inputs"), "--out_dir", str(base / "out"),
                    "--accelerator", args.accelerator, *PREDICT_ARGS, *args.extra.split()], check=True)
    elapsed = time.time() - t

    # boltz names its results folder after the input folder
    res = read_outputs(base / "out" / "boltz_results_inputs", records)

    res = res.merge(meta.drop_duplicates("name"), on="name", how="left")
    scores_path = DATA_DIR / "scores" / "molecules.parquet"
    if "published_score" not in res and scores_path.exists():
        published = pd.read_parquet(scores_path)
        published = published[published["target"] == args.target].groupby("name")["score"].mean()
        res["published_score"] = res["name"].map(published)
    res.to_csv(base / "scores.csv", index=False)
    cols = [c for c in ["name", "group", "pred", "score", "published_score"] if c in res]
    with pd.option_context("display.float_format", "{:.4f}".format, "display.width", 140, "display.max_rows", 500):
        print(res[cols])
    if "published_score" in res and res["published_score"].notna().sum() >= 3:
        known = res.dropna(subset=["published_score"])
        print(f"\ncalibration: corr {known['score'].corr(known['published_score']):.3f}, "
              f"mean abs diff {(known['score'] - known['published_score']).abs().mean():.4f}")
    if "group" in res and (res["group"] == "submission").any():
        print(f"proposed submission true total: {res.loc[res['group'] == 'submission', 'score'].sum():.3f} "
              f"({(res['group'] == 'submission').sum()} molecules scored)")
    print(f"\n{len(res)} molecules in {elapsed / 60:.1f} min ({elapsed / max(len(res), 1) / 60:.1f} min/molecule)")
    print(f"saved to {base / 'scores.csv'}")


if __name__ == "__main__":
    main()
