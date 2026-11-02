"""Has the subnet changed the rules under us?

Everything in nova_lab.config was copied from the validators' config.yaml at a
pinned commit. Validators run an auto-updater that pulls main hourly, so a
commit upstream changes the competition within the hour — most importantly the
target protein, which turns every archive, MSA, dataset and GPU batch we hold
into material for the wrong protein.

Measured cadence of nanobody target changes (git history of config/config.yaml,
April to October 2026): 7, 14, 21, 28, 36, 38 and 53 days apart, median 28.
Targets are announced on Discord about 24 h ahead, so this check is the floor,
not a substitute for watching the announcements.

Usage:
    python -m scripts.check_config
"""
import argparse
import sys

import requests
import yaml

from nova_lab import UPSTREAM_COMMIT
from nova_lab.config import (EPOCH_LENGTH, MIN_HUMAN_J_IDENTITY, MIN_HUMAN_V_IDENTITY,
                             NANOBODY_TARGET, NO_SUBMISSION_BLOCKS, NUM_MOLECULES,
                             SMALL_MOLECULE_TARGET, CDR3_LENGTH_RANGE,
                             V_NEIGHBORHOOD_CATEGORIES)

RAW = "https://raw.githubusercontent.com/metanova-labs/nova/{ref}/config/config.yaml"
API = "https://api.github.com/repos/metanova-labs/nova/commits/main"


def fetch(ref: str) -> dict:
    resp = requests.get(RAW.format(ref=ref), timeout=30)
    resp.raise_for_status()
    return yaml.safe_load(resp.text)


def head_commit() -> tuple[str, str]:
    try:
        data = requests.get(API, timeout=30).json()
        return data["sha"], data["commit"]["committer"]["date"]
    except Exception:
        return "", ""


def compare(live: dict) -> list[str]:
    """Differences on the values this lab depends on."""
    sel = live["protein_selection"]
    comp = live["competition"]
    nano = live["nanobody_requirements"]
    germ = live["germline_filtering"]
    checks = [
        ("small molecule target", SMALL_MOLECULE_TARGET, sel["small_molecule"]["target"].split(",")),
        ("nanobody target", NANOBODY_TARGET, sel["nanobody"]["target"].split(",")),
        ("nanobody structure", None, sel["nanobody"].get("structure")),
        ("no_submission_blocks", NO_SUBMISSION_BLOCKS, comp["no_submission_blocks"]),
        ("num_molecules", NUM_MOLECULES, live["molecule_requirements"]["num_molecules"]),
        ("min human V identity", MIN_HUMAN_V_IDENTITY, nano.get("min_human_v_identity")),
        ("min human J identity", MIN_HUMAN_J_IDENTITY, nano.get("min_human_j_identity")),
        ("CDR3 length range", [min(CDR3_LENGTH_RANGE), max(CDR3_LENGTH_RANGE)],
         [germ["cdr3_length_range"]["min"], germ["cdr3_length_range"]["max"]]),
        ("V neighbourhood categories", {k: list(v) for k, v in V_NEIGHBORHOOD_CATEGORIES.items()},
         {k: list(v) for k, v in germ["v_neighborhood_categories"].items()}),
        ("nanobody weight", None, comp["nanobody_weight"]),
    ]
    problems = []
    for name, ours, theirs in checks:
        if ours is None:
            print(f"  {name}: {theirs} (not mirrored locally)")
            continue
        same = ours == theirs
        print(f"  {'ok  ' if same else 'DIFF'} {name}: ours {ours} / upstream {theirs}")
        if not same:
            problems.append(f"{name}: ours {ours}, upstream {theirs}")
    return problems


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", default="main", help="upstream ref to compare against")
    args = parser.parse_args()

    sha, when = head_commit()
    if sha:
        print(f"upstream main: {sha[:7]} ({when})")
        print(f"our pin:       {UPSTREAM_COMMIT[:7]}"
              + ("  — same commit" if sha.startswith(UPSTREAM_COMMIT[:7]) else "  — BEHIND"))
    print(f"\nconfig/config.yaml at {args.ref}:")
    problems = compare(fetch(args.ref))

    print(f"\nepoch length {EPOCH_LENGTH} blocks")
    if problems:
        print("\nThe rules moved. Validators auto-update hourly, so this is already live:")
        for p in problems:
            print(f"  - {p}")
        print("\nIf the nanobody target changed, every archive, MSA, dataset and GPU batch "
              "here is for the wrong protein: rebuild before scoring anything.")
        sys.exit(1)
    print("\nNothing we depend on has moved.")


if __name__ == "__main__":
    main()
