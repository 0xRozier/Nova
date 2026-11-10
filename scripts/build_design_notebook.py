"""Build a Kaggle notebook that DESIGNS CDR3s with BoltzGen instead of grafting.

scripts/build_boltzgen_notebook.py scores sequences we already have. This one
generates them: same vendored BoltzGen, same target and binding site, but run
as the binder design model it is, with the framework held fixed and only the
CDR3 designed (see scripts/design_nanobody.py for the spec and why each region
is fixed or free).

The notebook is deliberately ordered so the cheap failure comes first:

  1. `boltzgen check` on the spec, which needs no GPU time worth speaking of
     and writes back the structure with the design mask applied. The cell
     asserts on it and prints the designed span, so the one assumption
     design_nanobody.py cannot verify offline -- that res_index is the
     sequential mmcif index and not IMGT -- is settled before anything
     expensive runs.
  2. the design itself, then the sequences and their metrics.

What to look for in the output: design_iiptm. Grafting gives us 0.375; the
field's median is 0.664 and the per-epoch best 0.757 (FINDINGS 10). If the
designed CDR3s land near 0.7, the method transfers and the generator should be
replaced. If they land near 0.4, it does not, and we have spent one Kaggle
session finding that out instead of a month of tuning.

The validator's gates are NOT applied here: ANARCI and hmmer would have to be
installed on Kaggle for that. The designs come back as a CSV and are screened
locally with scripts.preflight, which is where that code already lives.

Usage:
    python -m scripts.design_nanobody --allowed-v IGHV3-66 --cdr3 16 --frameworks 4 --fold
    python -m scripts.build_design_notebook --spec inputs/design
"""
import argparse
import base64
import hashlib
import io
import json
import tarfile
from pathlib import Path

from nova_lab import REPO_ROOT, UPSTREAM_COMMIT

MSA_URLS = [
    "https://cdn.jsdelivr.net/gh/metanova-labs/nova@{commit}/data/msa_files/{target}.a3m",
    "https://raw.githubusercontent.com/metanova-labs/nova@{commit}/data/msa_files/{target}.a3m",
]
OUR_FILES = ["nova_lab/__init__.py", "nova_lab/config.py", "kaggle/bittensor.py"]


def payload(spec_dir: Path, nova_dir: Path, target: str) -> str:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:xz") as tar:
        for rel in OUR_FILES:
            tar.add(REPO_ROOT / rel, arcname=rel)
        # The spec, the scaffold yamls and their structures.
        tar.add(spec_dir, arcname="design",
                filter=lambda info: None if "__pycache__" in info.name else info)
        tar.add(nova_dir / "external_tools" / "boltzgen" / "src",
                arcname="nova/boltzgen/src",
                filter=lambda info: None if "__pycache__" in info.name else info)
        tar.add(nova_dir / "external_tools" / "boltzgen" / "pyproject.toml",
                arcname="nova/boltzgen/pyproject.toml")
        tar.add(nova_dir / "data" / "structures" / "8AOK.cif",
                arcname="design/8AOK.cif")
    return base64.b64encode(buf.getvalue()).decode()


def md(text):
    return {"cell_type": "markdown", "metadata": {},
            "source": text.strip("\n").splitlines(True)}


def code(text):
    return {"cell_type": "code", "metadata": {}, "execution_count": None,
            "outputs": [], "source": text.strip("\n").splitlines(True)}


def build(spec_dir: Path, nova_dir: Path, target: str, run: str,
          num_designs: int) -> dict:
    report = json.loads((spec_dir / "report.json").read_text())
    msa_path = nova_dir / "data" / "msa_files" / f"{target}.a3m"
    msa_sha = hashlib.sha256(msa_path.read_bytes()).hexdigest()
    msa_urls = [u.format(commit=UPSTREAM_COMMIT, target=target) for u in MSA_URLS]
    spans = ", ".join(f"{r['scaffold']} {r['res_index_designed']}" for r in report)

    cells = [
        md(f"""
# NOVA SN68 — **conception** de CDR3 par BoltzGen (GPU Kaggle)

Ne note pas des séquences : elle en **fabrique**. Même BoltzGen embarqué par les validateurs, même cible ({target} / 8AOK) et même site de liaison, mais lancé comme le modèle de conception de binders qu'il est — charpente figée, CDR3 seul conçu.

**Pourquoi.** Le greffage de CDR3 archivés nous donne `design_iiptm` = 0,375. La médiane du champ est à 0,664 et le meilleur par époque à 0,757. Si les CDR3 conçus ici atterrissent vers 0,7, la méthode se transfère et il faut remplacer le générateur. Si c'est 0,4, elle ne se transfère pas — et on l'aura su en une session Kaggle au lieu d'un mois de réglages.

**Ce qui est figé et pourquoi** (mesuré, voir `docs/FINDINGS.md` §7ter) : CDR1 et CDR2 (notés contre la référence V), les 2 premiers résidus du CDR3 (codés par le gène V : `…YYCAR` pour IGHV3-66) et les 3 derniers (segment J). Les laisser libres coûtait 58 % de refus sur l'identité J et 38 % de pertes de l'appel V.

Spans conçus : {spans}

**Réglages Kaggle :** *Accelerator* **GPU T4 x2**, *Internet* **On**.

La cellule `check` s'exécute d'abord et **lève** si la spécification est mauvaise : c'est la seule hypothèse invérifiable hors GPU (`res_index` est l'index mmcif séquentiel, pas la numérotation IMGT qu'écrit NanoBodyBuilder2). Elle coûte quelques secondes ; la conception coûte des dizaines de minutes.
"""),
        code(f"""
# Paramètres
NUM_DESIGNS = {num_designs}   # designs generes au total, puis filtres par le pipeline
RUN = "{run}"
PROTOCOL = "nanobody-anything"

import base64, io, os, pathlib, subprocess, sys, tarfile, time
T0 = time.time()
TEMP = pathlib.Path("/kaggle/temp")
WORK = TEMP / "Nova"
NOVA = TEMP / "nova"
VENV = TEMP / "bgvenv"
OUT = pathlib.Path("/kaggle/working")
DESIGN = WORK / "design"

def step(msg):
    print(f"[{{time.strftime('%H:%M:%S')}}] +{{(time.time() - T0) / 60:5.1f}} min  {{msg}}", flush=True)

def sh(cmd, **kw):
    res = subprocess.run(cmd, shell=isinstance(cmd, str), capture_output=True, text=True, **kw)
    if res.stdout:
        print(res.stdout[-4000:])
    if res.returncode != 0:
        print(res.stderr[-4000:])
        raise RuntimeError(f"echec ({{res.returncode}}) : {{cmd}}")
    return res.stdout

for d in (WORK, NOVA):
    d.mkdir(parents=True, exist_ok=True)
with tarfile.open(fileobj=io.BytesIO(base64.b64decode(PAYLOAD)), mode="r:xz") as tar:
    for m in tar.getmembers():
        dest = NOVA if m.name.startswith("nova/") else WORK
        m.name = m.name[len("nova/"):] if m.name.startswith("nova/") else m.name
        tar.extract(m, dest)
step(f"charge utile depliee : {{len(list(DESIGN.iterdir()))}} fichiers dans {{DESIGN}}")
print(sorted(p.name for p in DESIGN.iterdir()))
"""),
        code(f"""
# La MSA de la cible, verifiee par empreinte
import hashlib, urllib.request
msa = DESIGN / "{target}.a3m"
for url in {msa_urls!r}:
    try:
        data = urllib.request.urlopen(url, timeout=120).read()
    except Exception as exc:
        step(f"{{url}} : {{exc}}")
        continue
    if hashlib.sha256(data).hexdigest() == "{msa_sha}":
        msa.write_bytes(data)
        step(f"MSA verifiee ({{len(data)}} octets)")
        break
else:
    raise RuntimeError("MSA introuvable ou corrompue")
"""),
        md("""
## Environnement BoltzGen
La version des validateurs épingle `numpy==2.0.2`, `torch==2.7.1` et une trentaine d'autres paquets : venv séparé, ~10 min.
"""),
        code("""
step("installation de uv + Python 3.12")
sh(f"{sys.executable} -m pip install -q uv")
sh(f"uv venv -q --python 3.12 {VENV}")
py = str(VENV / "bin/python")
bg = str(VENV / "bin/boltzgen")
step("installation de torch 2.7.1 CUDA")
sh(f"uv pip install -q --python {py} torch==2.7.1 --index-url https://download.pytorch.org/whl/cu126")
step("installation de la BoltzGen des validateurs")
sh(f"uv pip install -q --python {py} {NOVA}/boltzgen pandas requests")
site = pathlib.Path(sh([py, "-c", "import sysconfig; print(sysconfig.get_paths()['purelib'])"]).strip())
sh(["cp", str(WORK / "kaggle" / "bittensor.py"), str(site / "bittensor.py")])
sh(f'{py} -c "import torch, boltzgen, bittensor; print(torch.__version__, torch.cuda.is_available())"')
step("environnement pret")
"""),
        md("""
## 1. `check` — l'indexation est-elle bonne ?

**Cette cellule lève si la spécification est invalide**, avant toute dépense GPU. Elle écrit aussi la structure avec le masque de conception appliqué : `check_out/` est téléchargeable depuis *Output* et s'ouvre dans molstar.org/viewer pour confirmer visuellement que le span conçu est bien le CDR3 et pas une portion de charpente.
"""),
        code("""
os.chdir(DESIGN)
check_out = OUT / "check_out"
step("boltzgen check sur la specification")
try:
    out = sh([bg, "check", "design_spec.yaml", "--output", str(check_out)])
except RuntimeError:
    print()
    print("=== SI CE CHECK ECHOUE, LIRE CECI AVANT DE BRICOLER ===")
    print("La cause la plus probable est le couple exclude / design_insertions.")
    print("Les scaffolds excluent TOUT le milieu concu du CDR3 puis reinserent")
    print("exactement le nombre de residus voulu. L exemple amont")
    print("example/nanobody_scaffolds/7eow.yaml n exclut qu une PARTIE du span")
    print("concu, ce qui est peut-etre la seule forme supportee.")
    print()
    print("Correctif a essayer, dans scripts/design_nanobody.py, scaffold_spec() :")
    print("  garder quelques residus du donneur dans le span concu, c est-a-dire")
    print("  exclure start+V_HEAD+k .. donor_middle_end et inserer")
    print("  designable_count-k residus, avec k petit (1 ou 2).")
    print("L arithmetique de longueur reste : 2 (tete) + designable_count + 3 (queue).")
    print()
    print("L indexation, elle, est deja verifiee hors ligne : design_nanobody.py")
    print("relit les residus que le span nomme dans la structure et les compare")
    print("au CDR3 du donneur prive de sa tete et de sa queue (verify_span).")
    raise
step("specification valide")
print("structures ecrites :", sorted(p.name for p in check_out.iterdir()) if check_out.exists() else "(aucune)")
"""),
        md(f"""
## 2. Conception

`--protocol nanobody-anything` avec `--num_designs NUM_DESIGNS`. Le premier lancement télécharge les poids (~6 Go). Mesuré en mode notation sur 2×T4 : design 49 s/séquence, repliement 68 s, analyse 4,5 s. La conception génère davantage de candidats par séquence retenue, donc compter large.
"""),
        code("""
log_path = OUT / f"design_{RUN}.log"
out_dir = WORK / "designs" / RUN
if out_dir.exists():
    raise SystemExit(f"{out_dir} existe deja : BoltzGen REPREND les resultats au lieu "
                     f"de recalculer. Changer RUN pour un nouveau tirage.")
configure = [bg, "configure", "design_spec.yaml", "--output", str(out_dir),
             "--protocol", PROTOCOL, "--num_designs", str(NUM_DESIGNS)]
step(" ".join(configure))
sh(configure)

t_design = time.time()
step(f"conception de {NUM_DESIGNS} designs (run {RUN})")
with open(log_path, "w") as log:
    proc = subprocess.Popen([bg, "execute", str(out_dir)], stdout=log, stderr=subprocess.STDOUT)
    while proc.poll() is None:
        time.sleep(300)
        lines = [l.strip() for l in log_path.read_text(errors="replace").replace("\\r","\\n").splitlines() if l.strip()]
        step(f"en cours — {lines[-1][:160] if lines else '(journal vide)'}")
elapsed = time.time() - t_design
print(log_path.read_text(errors="replace").replace("\\r", "\\n")[-6000:])
step(f"conception terminee (code {proc.returncode}) en {elapsed/60:.1f} min")
assert proc.returncode == 0, "BoltzGen a echoue : voir le journal ci-dessus"
"""),
        md("""
## 3. Les séquences conçues et leurs métriques

**Le chiffre qui décide : `design_iiptm`.** Greffage = 0,375. Médiane du champ = 0,664. Meilleur par époque = 0,757.

Les verrous du validateur (identité germinale, marqueurs VHH, PSH, Aggrescan3D, unicité d'archive) ne sont **pas** appliqués ici : ANARCI et hmmer ne sont pas installés. Le CSV se criblera localement avec `scripts.preflight`.
"""),
        code("""
import pandas as pd, glob
rows = []
for f in sorted(glob.glob(str(out_dir / "**" / "*.csv"), recursive=True)):
    try:
        d = pd.read_csv(f)
    except Exception:
        continue
    if "design_iiptm" in d.columns or "sequence" in d.columns:
        print(f"{f}  ({len(d)} lignes, colonnes : {list(d.columns)[:10]})")
        rows.append(d.assign(_source=pathlib.Path(f).name))
if not rows:
    raise SystemExit("aucun CSV de resultats : inspecter l arborescence ci-dessous")
res = pd.concat(rows, ignore_index=True)
dest = OUT / f"designs_{RUN}.csv"
res.to_csv(dest, index=False)
print(f"\\n{len(res)} lignes -> {dest}")
if "design_iiptm" in res.columns:
    s = res.design_iiptm.dropna()
    print(f"\\n=== design_iiptm des CDR3 CONCUS ===")
    print(f"  n={len(s)}, mediane {s.median():.3f}, max {s.max():.3f}")
    print(f"  greffage (notre generateur actuel) : 0.375")
    print(f"  mediane du champ reel              : 0.664")
    print(f"  meilleur par epoque du champ       : 0.757")
    verdict = ("la methode se transfere, il faut remplacer le generateur"
               if s.max() >= 0.60 else
               "elle ne se transfere pas telle quelle a ce cadre")
    print(f"  -> {verdict}")
res.head(20)
"""),
        code("""
sh(f"find {out_dir} -maxdepth 3 | head -60")
"""),
    ]
    # Before the cell that unpacks it, not after: cells[1] both sets the
    # parameters and opens the tarball, so PAYLOAD has to exist by then. The
    # scoring notebook inserts at 2 because there those are two separate cells,
    # and copying that index shipped a notebook that died on NameError at the
    # first cell.
    cells.insert(1, code(
        f"# Charge utile : BoltzGen des validateurs, 8AOK, la specification de conception "
        f"et ses charpentes (tar.xz base64)\nPAYLOAD = \"{payload(spec_dir, nova_dir, target)}\""))
    sources = ["".join(c["source"]) for c in cells]
    defines = next(i for i, s in enumerate(sources) if s.startswith("# Charge utile"))
    uses = next(i for i, s in enumerate(sources) if "b64decode(PAYLOAD)" in s)
    if defines >= uses:
        raise AssertionError(
            f"PAYLOAD is defined in cell {defines} and used in cell {uses}: the "
            f"notebook would fail on NameError at the first run.")
    return {"cells": cells,
            "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python",
                                        "name": "python3"},
                         "language_info": {"name": "python"}},
            "nbformat": 4, "nbformat_minor": 5}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--spec", type=Path, default=Path("inputs/design"),
                        help="directory scripts.design_nanobody wrote")
    parser.add_argument("--target", default="Q9NZQ7")
    parser.add_argument("--run", default="design1")
    parser.add_argument("--num-designs", type=int, default=20)
    parser.add_argument("--nova-dir", type=Path,
                        default=REPO_ROOT.parent / "nova-upstream")
    parser.add_argument("--out", default="notebooks/kaggle_design.ipynb")
    args = parser.parse_args()

    spec = args.spec / "design_spec.yaml"
    if not spec.exists():
        raise SystemExit(f"{spec} absent : lancer scripts.design_nanobody --fold d'abord")
    missing = [r["scaffold"].replace(".yaml", ".pdb") for r in
               json.loads((args.spec / "report.json").read_text())
               if not (args.spec / r["scaffold"].replace(".yaml", ".pdb")).exists()]
    if missing:
        raise SystemExit(f"structures manquantes : {missing}. Relancer "
                         f"scripts.design_nanobody avec --fold")

    nb = build(args.spec, args.nova_dir, args.target, args.run, args.num_designs)
    out = REPO_ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(nb))
    print(f"wrote {out} ({out.stat().st_size // 1024} kB)")


if __name__ == "__main__":
    main()
