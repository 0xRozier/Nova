"""Check a built BoltzGen notebook before anyone spends a GPU session on it.

Three Kaggle sessions were lost to defects a local check would have caught: a
PAYLOAD cell referenced before it ran, a subprocess that inherited Kaggle's
MPLBACKEND, and a missing accelerator the notebook merely printed. Each of
those is an assertion here. Run it under Kaggle's own MPLBACKEND, which is how
the failures reproduce:

    MPLBACKEND="module://matplotlib_inline.backend_inline" \
        python scripts/verify_boltzgen_notebook.py kaggle_minimal2_gpu.ipynb

It does not launch anything and needs no GPU.
"""

from __future__ import annotations

import argparse
import ast
import base64
import csv
import io
import json
import os
import subprocess
import re
import sys
import tarfile
import tempfile
from pathlib import Path

from nova_lab.config import NANOBODY_STRUCTURE, NANOBODY_TARGET

ROOT = Path(__file__).resolve().parent.parent
# The structure is whichever one nova_lab.config names, because that is the one
# scripts.boltzgen_score writes into its yaml. Pinning a PDB id here would pass
# a notebook built for the previous target.
REQUIRED_IN_PAYLOAD = (
    "scripts/boltzgen_score.py",
    "nova/boltzgen/pyproject.toml",
    f"nova/data/structures/{NANOBODY_STRUCTURE}.cif",
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("notebook")
    parser.add_argument("--batch", help="batch path as embedded, e.g. inputs/minimal2_gpu.csv")
    parser.add_argument("--run", help="expected run name")
    parser.add_argument("--controls", action="store_true",
                        help="the batch pairs variants with their parents: check every parent rides along")
    args = parser.parse_args()

    cells = json.loads(Path(args.notebook).read_text())["cells"]
    code = [(i, "".join(c["source"])) for i, c in enumerate(cells) if c["cell_type"] == "code"]
    for index, source in code:
        try:
            ast.parse(source)
        except SyntaxError as exc:
            raise SystemExit(f"cellule {index} ne s analyse pas : {exc}")
    print(f"{len(code)} cellules de code, toutes analysables")

    defined = next((i for i, s in code if "PAYLOAD = " in s), None)
    if defined is None:
        raise SystemExit("aucune cellule ne definit PAYLOAD")
    used = [i for i, s in code if "PAYLOAD" in s and i != defined]
    if not used or any(u < defined for u in used):
        raise SystemExit(f"PAYLOAD defini en {defined}, utilise en {used} : ordre invalide")
    print(f"PAYLOAD defini en cellule {defined}, utilise en {used}")

    blob = next(s for i, s in code if i == defined).split('PAYLOAD = "')[1].strip().rstrip('"')
    with tempfile.TemporaryDirectory() as tmp:
        unpacked = Path(tmp)
        with tarfile.open(fileobj=io.BytesIO(base64.b64decode(blob)), mode="r:xz") as tar:
            tar.extractall(unpacked, filter="data")
            names = set(tar.getnames())
        for need in REQUIRED_IN_PAYLOAD + ((args.batch,) if args.batch else ()):
            if need not in names:
                raise SystemExit(f"la charge utile n embarque pas {need}")
        print(f"charge utile deballee : {len(names)} entrees")

        # The failure that cost 6 GPU minutes of a 2 h run: the notebook had
        # downloaded P05231.a3m while the embedded scorer still carried a
        # hardcoded TARGET = "Q9NZQ7" and reopened the MSA under that name.
        # Nothing disagreed locally, because the name is only ever used to
        # build a path. Checked here against the embedded copy of the scorer,
        # not against the local one, so a stale payload cannot pass.
        embedded = (unpacked / "scripts" / "boltzgen_score.py").read_text()
        found = re.search(r'^TARGET = (.+)$', embedded, re.M)
        if found is None:
            raise SystemExit("la charge utile ne definit pas TARGET dans boltzgen_score")
        literal = re.fullmatch(r'"([A-Z0-9]+)"', found.group(1).strip())
        scorer_target = literal.group(1) if literal else NANOBODY_TARGET[0]
        msa = re.search(r'^MSA_NAME = "(.+?)"$', "\n".join(
            s2 for _, s2 in code), re.M)
        if msa is None:
            raise SystemExit("le carnet ne definit pas MSA_NAME")
        if msa.group(1) != f"{scorer_target}.a3m":
            raise SystemExit(
                f"le carnet telecharge {msa.group(1)} mais le scoreur embarque "
                f"cherchera {scorer_target}.a3m : regenerer le carnet")
        print(f"MSA et cible du scoreur concordent : {msa.group(1)}")

        if args.batch:
            rows = list(csv.DictReader((unpacked / args.batch).open()))
            sequences = [r["sequence"] for r in rows]
            if len(set(sequences)) != len(sequences):
                raise SystemExit("sequence en double : BoltzGen ne la noterait qu une fois")
            print(f"lot : {len(rows)} sequences, sans doublon")
            if args.controls:
                roles = [r.get("role") for r in rows]
                parents = {r["sequence"] for r in rows if r.get("role") == "parent"}
                orphans = [r for r in rows if r.get("role") == "variant"
                           and r.get("parent_sequence") not in parents]
                if orphans:
                    raise SystemExit(
                        f"{len(orphans)} variantes dont le parent n est pas dans le lot : "
                        "sans controle interne la mesure ne vaut rien"
                    )
                print(f"  {roles.count('variant')} variantes, {roles.count('parent')} parents "
                      "comme controles internes")

    whole = "".join("".join(c["source"]) for c in cells)
    for marker, what in (("existe déjà", "garde anti-reprise de run"),
                         ("Aucun GPU", "garde accelerateur absent")):
        if marker not in whole:
            raise SystemExit(f"{what} absente du carnet")
    print("gardes presentes : nom de run neuf, refus si pas de GPU")
    if args.run and f'RUN = "{args.run}"' not in whole:
        raise SystemExit(f'le carnet ne porte pas RUN = "{args.run}"')

    inherited = os.environ.get("MPLBACKEND", "(absent)")
    probe = subprocess.run(
        [sys.executable, "-c",
         "from scripts.boltzgen_score import subprocess_env; print(subprocess_env()['MPLBACKEND'])"],
        cwd=ROOT, capture_output=True, text=True,
    )
    if probe.returncode != 0 or probe.stdout.strip() != "Agg":
        raise SystemExit(f"MPLBACKEND des enfants non force a Agg : {probe.stdout}{probe.stderr[-500:]}")
    print(f"MPLBACKEND herite {inherited} -> enfants forces a Agg")
    print("\ncarnet verifiable en local : rien d autre ne se teste sans GPU")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
