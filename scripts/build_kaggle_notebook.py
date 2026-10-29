"""Build a self-contained Kaggle notebook that scores a batch with Boltz-2.

The notebook carries its own code (nova_lab + scripts/boltz_score.py), the
batch CSV and the target's MSA as an embedded archive: fetching anything from
GitHub proved unreliable on Kaggle (git clone hung, codeload returned non-gzip).
The MSA is downloaded (jsDelivr, then GitHub) and checked against its hash;
otherwise only Hugging Face (molecule DB) and PyPI/PyTorch are used.
Run scripts.boltz_score once locally first so data/msa/<target>.{a3m,seq} exist.

Usage:
    python -m scripts.build_kaggle_notebook [--batch inputs/kaggle_batch.csv] [--target O15427]
"""
import argparse
import base64
import hashlib
import io
import json
import tarfile

from nova_lab import REPO_ROOT, UPSTREAM_COMMIT

KAGGLE_IMPORT_LIMIT = 1_000_000
# The MSA is too big to embed under Kaggle's 1 MB import limit: download it,
# jsDelivr first (GitHub itself is unreliable from Kaggle), checked by hash.
MSA_URLS = [
    "https://cdn.jsdelivr.net/gh/metanova-labs/nova@{commit}/data/msa_files/{target}.a3m",
    "https://raw.githubusercontent.com/metanova-labs/nova/{commit}/data/msa_files/{target}.a3m",
]

FILES = [
    "nova_lab/__init__.py",
    "nova_lab/config.py",
    "nova_lab/reactions.py",
    "nova_lab/filters.py",
    "nova_lab/chem.py",
    "scripts/__init__.py",
    "scripts/boltz_score.py",
]


def payload(batch: str, target: str) -> str:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:xz") as tar:
        for rel in FILES + [batch]:
            tar.add(REPO_ROOT / rel, arcname=rel)
        tar.add(REPO_ROOT / "data" / "msa" / f"{target}.seq", arcname=f"data/msa/{target}.seq")
    return base64.b64encode(buf.getvalue()).decode()


def md(text):
    return {"cell_type": "markdown", "metadata": {}, "source": text.strip("\n").splitlines(True)}


def code(text):
    return {"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [],
            "source": text.strip("\n").splitlines(True)}


def build(batch: str, target: str) -> dict:
    msa_sha = hashlib.sha256((REPO_ROOT / "data" / "msa" / f"{target}.a3m").read_bytes()).hexdigest()
    msa_urls = [u.format(commit=UPSTREAM_COMMIT, target=target) for u in MSA_URLS]
    cells = [
        md(f"""
# NOVA SN68 — scoring Boltz-2 sur GPU Kaggle

Note un lot de molécules avec Boltz-2 **exactement comme les validateurs NOVA** (même entrée, même MSA, mêmes paramètres, même graine).

Notebook **autonome** : le code, le lot (`{batch}`) et la séquence de {target} sont intégrés (cellule de données). La MSA des validateurs est téléchargée via jsDelivr (GitHub en secours) et vérifiée par empreinte SHA-256.

**Réglages Kaggle :** *Accelerator* **GPU T4 x2**, *Internet* **On**.
- Test rapide : mettre `TOP = 2` ci-dessous puis **Run All** dans l'éditeur.
- Lot complet : `TOP = None`, puis **Save Version → Save & Run All (Commit)** ; récupérer `scores_<RUN>.csv` dans l'onglet *Output*.

Seuls `scores_<RUN>.csv` et `boltz_<RUN>.log` sont écrits dans `/kaggle/working`. Usage de recherche uniquement : rien n'est soumis au réseau Bittensor.
"""),
        code("""
# Paramètres
TOP = None          # None = tout le lot ; ex. 2 pour un test de durée rapide
RUN = "kaggle_batch1"
MULTI_GPU = False   # True pour essayer les 2 T4 (plus rapide, moins éprouvé)

import base64, io, os, pathlib, subprocess, sys, tarfile, time
T0 = time.time()
WORK = pathlib.Path("/kaggle/temp/Nova")
VENV = pathlib.Path("/kaggle/temp/bvenv")
OUT = pathlib.Path("/kaggle/working")

def step(msg):
    print(f"[{time.strftime('%H:%M:%S')}] +{(time.time() - T0) / 60:5.1f} min  {msg}", flush=True)

def sh(cmd, **kw):
    \"\"\"Run a shell command, show its output, stop the notebook if it fails.\"\"\"
    res = subprocess.run(cmd, shell=True, text=True, capture_output=True, **kw)
    out = (res.stdout + res.stderr).strip()
    if out:
        print(out[-3000:])
    if res.returncode != 0:
        raise RuntimeError(f"échec ({res.returncode}) : {cmd}")
    return res.stdout

step("démarrage")
sh("nvidia-smi --query-gpu=name,memory.total --format=csv")
"""),
        code("""
step("test de l'accès Internet")
for url in ["https://huggingface.co", "https://pypi.org", "https://download.pytorch.org/whl/cu126/"]:
    code = subprocess.run(["curl", "-s", "-o", "/dev/null", "-m", "15", "-w", "%{http_code}", url],
                          capture_output=True, text=True).stdout or "000"
    print(f"  {url}: {code}")
    if code == "000":
        raise SystemExit("Pas d'accès Internet : activer 'Internet' dans les paramètres du notebook, puis relancer.")
"""),
        code("""
step("extraction du code, du lot et de la MSA intégrés")
WORK.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)
with tarfile.open(fileobj=io.BytesIO(base64.b64decode(PAYLOAD)), mode="r:xz") as tar:
    tar.extractall(WORK)
os.chdir(WORK)
step("téléchargement de la MSA des validateurs")
import hashlib
msa = WORK / "data/msa/" / MSA_NAME
for url in MSA_URLS:
    subprocess.run(["curl", "-sSL", "-m", "120", "--retry", "2", "-o", str(msa), url])
    if msa.exists() and hashlib.sha256(msa.read_bytes()).hexdigest() == MSA_SHA256:
        print(f"  MSA OK depuis {url}")
        break
    print(f"  échec depuis {url}")
else:
    raise RuntimeError("MSA introuvable ou corrompue")
sh("find . -type f | sort")
step("téléchargement de la base combinatoire (~120 Mo)")
sh('wget -nv --timeout=60 --tries=3 -O data/molecules.sqlite '
   '"https://huggingface.co/datasets/Metanova/Mol-Rxn-DB/resolve/main/molecules.sqlite"')
assert (WORK / "data/molecules.sqlite").stat().st_size > 100_000_000, "base incomplète"
step("étape 1 terminée")
"""),
        md("""
## Environnement isolé pour Boltz 2.2.1
Boltz 2.2.1 exige `numpy<2` et Python ≤ 3.12 : il est installé dans un venv séparé (~2 min).
"""),
        code("""
step("installation de uv + Python 3.12")
sh(f"{sys.executable} -m pip install -q uv")
sh(f"uv venv -q --python 3.12 {VENV}")
py = str(VENV / "bin/python")
step("installation de torch 2.7.1 CUDA")
sh(f"uv pip install -q --python {py} torch==2.7.1 --index-url https://download.pytorch.org/whl/cu126")
step("installation de boltz 2.2.1")
sh(f'uv pip install -q --python {py} boltz==2.2.1 "pyarrow<17" pandas requests pyyaml')
sh(f'{py} -c "import torch, boltz; print(torch.__version__, torch.cuda.is_available(), torch.cuda.device_count())"')
step("étape 2 terminée")
"""),
        md("""
## Précision de calcul
Les validateurs tournent en `bf16-mixed` sur des GPU récents. Les T4 (architecture 7.5) ne gèrent pas le bf16 : on passe Boltz en 32 bits. Léger écart, mesuré par les molécules de calibration du lot.
"""),
        code("""
cap = sh(f"{py} -c \\"import torch; print('%d.%d' % torch.cuda.get_device_capability(0))\\"").strip()
n_gpus = int(sh(f'{py} -c "import torch; print(torch.cuda.device_count())"'))
main_py = pathlib.Path(sh(f'{py} -c "import boltz.main as m; print(m.__file__)"').strip())
if float(cap) < 8.0:
    main_py.write_text(main_py.read_text().replace('"bf16-mixed"', '32'))
    print("bf16 non supporté -> précision 32 bits")
EXTRA = f"--no_kernels --devices {n_gpus if MULTI_GPU else 1}"
step(f"réglages Boltz : capacité {cap}, {EXTRA}")
"""),
        md("""
## Scoring
Le premier lancement télécharge les poids Boltz-2 (~6 Go). Une ligne s'affiche à chaque molécule terminée (avec le temps restant estimé) et toutes les 10 minutes sinon. Relancer avec le même `RUN` reprend là où le calcul s'était arrêté.
"""),
        code(f"""
import pandas as pd
BATCH = "{batch}"
n_total = len(pd.read_csv(BATCH).head(TOP) if TOP else pd.read_csv(BATCH))
cmd = [py, "-m", "scripts.boltz_score", "--file", BATCH, "--accelerator", "gpu", "--run", RUN,
       "--boltz", str(VENV / "bin/boltz"), f"--extra={{EXTRA}}"] + (["--top", str(TOP)] if TOP else [])
log_path = OUT / f"boltz_{{RUN}}.log"
pred_dir = WORK / f"data/boltz/{{RUN}}/out/boltz_results_inputs/predictions"

def last_log_line():
    lines = [l.strip() for l in log_path.read_text(errors="replace").replace("\\r", "\\n").splitlines() if l.strip()]
    return lines[-1][:200] if lines else "(journal vide)"

step(f"lancement de Boltz sur {{n_total}} molécules")
with open(log_path, "w") as log:
    proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)
    done, first_done_at, last_report = 0, None, time.time()
    while proc.poll() is None:
        time.sleep(30)
        now_done = len(list(pred_dir.glob("*/affinity_*.json")))
        if now_done > done:
            first_done_at = first_done_at or time.time()
            rate = (time.time() - first_done_at) / 60 / (now_done - 1) if now_done > 1 else None
            eta = f" | {{rate:.1f}} min/molécule, ~{{rate * (n_total - now_done):.0f}} min restantes" if rate else ""
            step(f"{{now_done}}/{{n_total}} molécules notées{{eta}}")
            done, last_report = now_done, time.time()
        elif time.time() - last_report > 600:
            step(f"en cours ({{done}}/{{n_total}}) — {{last_log_line()}}")
            last_report = time.time()

print(log_path.read_text(errors="replace").replace("\\r", "\\n")[-4000:])
step(f"Boltz terminé (code {{proc.returncode}})")
assert proc.returncode == 0, "Boltz a échoué : voir le journal ci-dessus"
"""),
        code("""
import shutil
out = OUT / f"scores_{RUN}.csv"
shutil.copy(WORK / f"data/boltz/{RUN}/scores.csv", out)
df = pd.read_csv(out)
step(f"résultats : {out} ({len(df)} molécules)")
df
"""),
    ]
    # The payload must exist before the extraction cell runs: define it right after the parameters.
    cells.insert(2, code(f"# Données intégrées : code nova_lab, scripts/boltz_score.py, {batch}, séquence {target} (tar.xz en base64)\n"
                         f"MSA_NAME = \"{target}.a3m\"\n"
                         f"MSA_SHA256 = \"{msa_sha}\"\n"
                         f"MSA_URLS = {json.dumps(msa_urls)}\n"
                         f"PAYLOAD = \"{payload(batch, target)}\""))
    return {"cells": cells,
            "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                         "language_info": {"name": "python"}},
            "nbformat": 4, "nbformat_minor": 5}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch", default="inputs/kaggle_batch.csv")
    parser.add_argument("--target", default="O15427")
    parser.add_argument("--out", default="notebooks/kaggle_boltz.ipynb")
    args = parser.parse_args()
    nb = build(args.batch, args.target)
    path = REPO_ROOT / args.out
    path.write_text(json.dumps(nb, indent=1, ensure_ascii=False))
    size = path.stat().st_size
    print(f"wrote {path} ({size / 1e3:.0f} kB)")
    if size >= KAGGLE_IMPORT_LIMIT:
        raise SystemExit(f"too big for Kaggle's {KAGGLE_IMPORT_LIMIT / 1e6:.0f} MB import limit")


if __name__ == "__main__":
    main()
