"""Build a self-contained Kaggle notebook that scores nanobodies with BoltzGen.

BoltzGen is scored with the build vendored in metanova-labs/nova, not the PyPI
one, and Kaggle cannot reliably fetch GitHub, so its source (~300 kB packed),
the target structure and our scoring code are embedded in the notebook. The
MSA comes from jsDelivr checked against its hash, and the model weights (~6 GB)
from Hugging Face, which works there.

Needs the upstream clone: --nova-dir, default ../nova-upstream.

Usage:
    python -m scripts.build_boltzgen_notebook [--batch inputs/kaggle_nanobody_batch.csv]
"""
import argparse
import base64
import hashlib
import io
import json
import tarfile
from pathlib import Path

from nova_lab import REPO_ROOT, UPSTREAM_COMMIT

KAGGLE_IMPORT_LIMIT = 1_000_000
MSA_URLS = [
    "https://cdn.jsdelivr.net/gh/metanova-labs/nova@{commit}/data/msa_files/{target}.a3m",
    "https://raw.githubusercontent.com/metanova-labs/nova/{commit}/data/msa_files/{target}.a3m",
]
OUR_FILES = [
    "nova_lab/__init__.py",
    "nova_lab/config.py",
    "scripts/__init__.py",
    "scripts/boltzgen_score.py",
    "kaggle/bittensor.py",
]


def payload(batch: str, nova_dir: Path, target: str) -> str:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:xz") as tar:
        for rel in OUR_FILES + [batch]:
            tar.add(REPO_ROOT / rel, arcname=rel)
        # The vendored BoltzGen: source and packaging only, no weights or docs.
        tar.add(nova_dir / "external_tools" / "boltzgen" / "src", arcname="nova/boltzgen/src",
                filter=lambda info: None if "__pycache__" in info.name else info)
        tar.add(nova_dir / "external_tools" / "boltzgen" / "pyproject.toml",
                arcname="nova/boltzgen/pyproject.toml")
        tar.add(nova_dir / "data" / "structures" / "8AOK.cif", arcname="nova/data/structures/8AOK.cif")
    return base64.b64encode(buf.getvalue()).decode()


def md(text):
    return {"cell_type": "markdown", "metadata": {}, "source": text.strip("\n").splitlines(True)}


def code(text):
    return {"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [],
            "source": text.strip("\n").splitlines(True)}


def build(batch: str, nova_dir: Path, target: str, run: str | None = None) -> dict:
    # One run name per batch: BoltzGen resumes an existing run directory, so a
    # stale name silently returns the previous results instead of scoring.
    run = run or Path(batch).stem
    n_seq = sum(1 for _ in (REPO_ROOT / batch).open()) - 1
    msa_path = nova_dir / "data" / "msa_files" / f"{target}.a3m"
    msa_sha = hashlib.sha256(msa_path.read_bytes()).hexdigest()
    msa_urls = [u.format(commit=UPSTREAM_COMMIT, target=target) for u in MSA_URLS]

    cells = [
        md(f"""
# NOVA SN68 — scoring BoltzGen de nanobodies sur GPU Kaggle

Note des séquences avec **la version de BoltzGen embarquée par les validateurs NOVA**, avec leur structure cible ({target} / 8AOK), leur site de liaison, leur protocole et leurs échelles.

Notebook **autonome** : la source BoltzGen, la structure, le lot (`{batch}`) et le code de scoring sont intégrés. La MSA est téléchargée via jsDelivr et vérifiée par empreinte SHA-256 ; les poids (~6 Go) viennent de Hugging Face.

**Réglages Kaggle :** *Accelerator* **GPU T4 x2**, *Internet* **On**.
- Lot complet ({n_seq} séquences, ~{n_seq * 3.7 / 60:.1f} h à 3,7 min/séquence sur 2×T4) : **Save Version → Save & Run All (Commit)**.
- Test rapide : `TOP = 2` puis **Run All** dans l'éditeur.

Le lot alterne candidats et séquences de calibration, donc même un run interrompu compare les deux.

Les noyaux cuEquivariance se désactivent d'eux-mêmes sur T4 (capacité 7.5 < 8). Mesuré : **3,7 min/séquence** sur 2×T4 (design 162 s, repliement 213 s, analyse 23 s pour 2 séquences).
"""),
        code("""
# Paramètres
TOP = None          # None = tout le lot ; 2 pour un test rapide
RUN = "__RUN__"
REPLICATE = 1       # incrémenter pour une réplique : BoltzGen est stochastique, mais il
                    # REPREND les résultats déjà présents, donc un nom de run neuf est
                    # indispensable pour obtenir un vrai second tirage.

import base64, io, os, pathlib, subprocess, sys, tarfile, time
T0 = time.time()
TEMP = pathlib.Path("/kaggle/temp")
WORK = TEMP / "Nova"
NOVA = TEMP / "nova"
VENV = TEMP / "bgvenv"
OUT = pathlib.Path("/kaggle/working")

RUN = RUN if REPLICATE == 1 else f"{RUN}_rep{REPLICATE}"

def step(msg):
    print(f"[{time.strftime('%H:%M:%S')}] +{(time.time() - T0) / 60:5.1f} min  {msg}", flush=True)

def sh(cmd, check=True, **kw):
    # a string goes through the shell, a list does not
    res = subprocess.run(cmd, shell=isinstance(cmd, str), text=True, capture_output=True, **kw)
    out = (res.stdout + res.stderr).strip()
    if out:
        print(out[-3000:])
    if check and res.returncode != 0:
        raise RuntimeError(f"échec ({res.returncode}) : {cmd}")
    return res.stdout

step("démarrage")
# Un notebook réimporté repart avec l'accélérateur par défaut : sans GPU, BoltzGen
# serait inutilisable, autant le dire tout de suite.
import shutil
if shutil.which("nvidia-smi") is None:
    raise SystemExit("Aucun GPU : mettre 'Accelerator' sur 'GPU T4 x2' dans les "
                     "paramètres du notebook (panneau de droite), puis relancer.")
sh("nvidia-smi --query-gpu=name,memory.total --format=csv")
TEMP.mkdir(parents=True, exist_ok=True)
sh("df -h /kaggle/temp | tail -1", check=False)
""".replace("__RUN__", run)),
        code("""
step("test de l'accès Internet")
for url in ["https://huggingface.co", "https://pypi.org", "https://cdn.jsdelivr.net"]:
    rc = subprocess.run(["curl", "-s", "-o", "/dev/null", "-m", "15", "-w", "%{http_code}", url],
                        capture_output=True, text=True).stdout or "000"
    print(f"  {url}: {rc}")
    if rc == "000":
        raise SystemExit("Pas d'accès Internet : activer 'Internet' dans les paramètres du notebook.")
"""),
        code("""
step("extraction du code, de la source BoltzGen et de la structure")
for d in (WORK, OUT):
    d.mkdir(parents=True, exist_ok=True)
with tarfile.open(fileobj=io.BytesIO(base64.b64decode(PAYLOAD)), mode="r:xz") as tar:
    tar.extractall(WORK)
(TEMP / "nova").mkdir(parents=True, exist_ok=True)
sh(f"cp -r {WORK}/nova/. {NOVA}/ && rm -rf {WORK}/nova")
step("téléchargement de la MSA des validateurs")
import hashlib
msa = NOVA / "data/msa_files" / MSA_NAME
msa.parent.mkdir(parents=True, exist_ok=True)
for url in MSA_URLS:
    subprocess.run(["curl", "-sSL", "-m", "180", "--retry", "2", "-o", str(msa), url])
    if msa.exists() and hashlib.sha256(msa.read_bytes()).hexdigest() == MSA_SHA256:
        print(f"  MSA OK depuis {url}")
        break
    print(f"  échec depuis {url}")
else:
    raise RuntimeError("MSA introuvable ou corrompue")
sh(f"find {NOVA} -maxdepth 3 -type d | sort && ls {WORK}")
step("étape 1 terminée")
"""),
        md("""
## Environnement BoltzGen
La version des validateurs épingle `numpy==2.0.2`, `torch==2.7.1` et une trentaine d'autres paquets : elle est installée dans un venv séparé (compter ~10 min).
"""),
        code("""
step("installation de uv + Python 3.12")
sh(f"{sys.executable} -m pip install -q uv")
sh(f"uv venv -q --python 3.12 {VENV}")
py = str(VENV / "bin/python")
step("installation de torch 2.7.1 CUDA")
sh(f"uv pip install -q --python {py} torch==2.7.1 --index-url https://download.pytorch.org/whl/cu126")
step("installation de la BoltzGen des validateurs (+ pandas, requests)")
sh(f"uv pip install -q --python {py} {NOVA}/boltzgen pandas requests")
# La BoltzGen de NOVA importe bittensor pour neuf appels de journalisation ; un
# substitut minimal evite d'installer ce paquet lourd dans l'environnement epingle.
site = pathlib.Path(sh([py, "-c", "import sysconfig; print(sysconfig.get_paths()['purelib'])"]).strip())
sh(["cp", str(WORK / "kaggle" / "bittensor.py"), str(site / "bittensor.py")])
sh(f'{py} -c "import torch, boltzgen, bittensor; print(torch.__version__, torch.cuda.is_available(), '
   f'torch.cuda.get_device_capability(0))"')
step("étape 2 terminée")
"""),
        md("""
## Scoring
Le premier lancement télécharge les poids BoltzGen (~6 Go). La progression s'affiche toutes les 5 minutes avec la dernière ligne du journal ; le journal complet est dans `/kaggle/working/boltzgen_<RUN>.log`.
"""),
        code(f"""
import pandas as pd
BATCH = "{batch}"
os.chdir(WORK)
n_total = len(pd.read_csv(BATCH).head(TOP) if TOP else pd.read_csv(BATCH))
log_path = OUT / f"boltzgen_{{RUN}}.log"
cmd = [py, "-m", "scripts.boltzgen_score", "--file", BATCH, "--run", RUN,
       "--nova-dir", str(NOVA), "--boltzgen", str(VENV / "bin/boltzgen")]
if TOP:
    cmd += ["--top", str(TOP)]

def last_log_line():
    lines = [l.strip() for l in log_path.read_text(errors="replace").replace("\\r", "\\n").splitlines() if l.strip()]
    return lines[-1][:200] if lines else "(journal vide)"

out_dir = WORK / f"data/boltzgen/{{RUN}}/out"
if out_dir.exists():
    raise SystemExit(f"{{out_dir}} existe déjà : BoltzGen reprendrait les résultats au lieu de "
                     f"recalculer. Incrémenter REPLICATE pour un nouveau tirage.")
step(f"lancement de BoltzGen sur {{n_total}} séquences (run {{RUN}})")
with open(log_path, "w") as log:
    proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)
    while proc.poll() is None:
        time.sleep(300)
        step(f"en cours — {{last_log_line()}}")

print(log_path.read_text(errors="replace").replace("\\r", "\\n")[-5000:])
step(f"BoltzGen terminé (code {{proc.returncode}})")
assert proc.returncode == 0, "BoltzGen a échoué : voir le journal ci-dessus"
"""),
        code("""
import shutil
src = WORK / f"data/boltzgen/{RUN}/scores.csv"
out = OUT / f"scores_nanobody_{RUN}.csv"
shutil.copy(src, out)
df = pd.read_csv(out)
step(f"résultats : {out} ({len(df)} séquences)")
df
"""),
    ]
    cells.insert(2, code(
        f"# Données intégrées : BoltzGen des validateurs, structure 8AOK, {batch}, code de scoring (tar.xz base64)\n"
        f"MSA_NAME = \"{target}.a3m\"\n"
        f"MSA_SHA256 = \"{msa_sha}\"\n"
        f"MSA_URLS = {json.dumps(msa_urls)}\n"
        f"PAYLOAD = \"{payload(batch, nova_dir, target)}\""))
    return {"cells": cells,
            "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                         "language_info": {"name": "python"}},
            "nbformat": 4, "nbformat_minor": 5}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch", default="inputs/kaggle_nanobody_batch4.csv")
    parser.add_argument("--run", help="run name (default: the batch file's stem)")
    parser.add_argument("--target", default="Q9NZQ7")
    parser.add_argument("--nova-dir", type=Path, default=REPO_ROOT.parent / "nova-upstream")
    parser.add_argument("--out", default="notebooks/kaggle_boltzgen.ipynb")
    args = parser.parse_args()
    if not (args.nova_dir / "external_tools" / "boltzgen" / "src").is_dir():
        raise SystemExit(f"upstream clone not found at {args.nova_dir}; pass --nova-dir")

    nb = build(args.batch, args.nova_dir, args.target, args.run)
    path = REPO_ROOT / args.out
    path.write_text(json.dumps(nb, indent=1, ensure_ascii=False))
    size = path.stat().st_size
    print(f"wrote {path} ({size / 1e3:.0f} kB)")
    if size >= KAGGLE_IMPORT_LIMIT:
        raise SystemExit(f"too big for Kaggle's {KAGGLE_IMPORT_LIMIT / 1e6:.0f} MB import limit")


if __name__ == "__main__":
    main()
