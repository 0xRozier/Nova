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


def target_structure(target: str) -> str:
    from scripts.design_nanobody import TARGET_STRUCTURE
    return TARGET_STRUCTURE


def is_hydrogen(line: str) -> bool:
    element = line[76:78].strip()
    return element == "H" if element else line[12:16].strip().startswith("H")


def without_hydrogens(text: str) -> str:
    """Drop hydrogens from a PDB, which halves it.

    NanoBodyBuilder2 writes all-atom models refined with OpenMM, so 48 % of the
    atoms are hydrogens -- 150 kB per scaffold, 655 kB compressed for forty.
    Kaggle refuses a notebook over 1 MB on import, and the forty-scaffold batch
    hit 1284 kB.

    Hydrogens are not what a structure model reads: BoltzGen's own reference
    spec, example/nanobody_scaffolds/7eow.yaml, points at a crystal structure
    from the PDB, which has none. Stripping them makes our input the ordinary
    case rather than the exotic one. `verify_span` re-reads every residue the
    spec names out of the stripped file, so a parser that disagreed would be
    caught here and not after an hour of Kaggle GPU.
    """
    return "".join(l for l in text.splitlines(True)
                   if not (l.startswith(("ATOM", "HETATM")) and is_hydrogen(l)))


def payload(spec_dir: Path, nova_dir: Path, target: str) -> str:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:xz") as tar:
        for rel in OUR_FILES:
            tar.add(REPO_ROOT / rel, arcname=rel)
        # The spec, the scaffold yamls and their structures. The PDBs go in
        # without their hydrogens; everything else verbatim.
        for item in sorted(spec_dir.iterdir()):
            if "__pycache__" in item.name or not item.is_file():
                continue
            if item.suffix == ".pdb":
                data = without_hydrogens(item.read_text()).encode()
                info = tarfile.TarInfo(f"design/{item.name}")
                info.size = len(data)
                tar.addfile(info, io.BytesIO(data))
            else:
                tar.add(item, arcname=f"design/{item.name}")
        # The validators' BoltzGen is NOT embedded: this notebook installs the
        # original 0.3.2 from PyPI (FINDINGS 55), so the vendored source was
        # 219 kB of payload that nothing referenced.
        tar.add(nova_dir / "data" / "structures" / f"{target_structure(target)}.cif",
                arcname=f"design/{target_structure(target)}.cif")
    return base64.b64encode(buf.getvalue()).decode()


def md(text):
    return {"cell_type": "markdown", "metadata": {},
            "source": text.strip("\n").splitlines(True)}


def code(text):
    return {"cell_type": "code", "metadata": {}, "execution_count": None,
            "outputs": [], "source": text.strip("\n").splitlines(True)}


def build(spec_dir: Path, nova_dir: Path, target: str, run: str,
          num_designs: int, steps: list[str] | None = None) -> dict:
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

> **À régler AVANT de lancer, dans l'éditeur : *Settings → Accelerator → GPU T4 x2*, et *Internet → On*.**
> Un carnet neuf démarre **sans** GPU, et BoltzGen en réclame un dès `configure`. Le carnet s'arrête maintenant tout de suite avec ce message plutôt que de mourir deux minutes plus tard dans une trace CUDA.

La cellule `check` s'exécute d'abord et **lève** si la spécification est mauvaise : c'est la seule hypothèse invérifiable hors GPU (`res_index` est l'index mmcif séquentiel, pas la numérotation IMGT qu'écrit NanoBodyBuilder2). Elle coûte quelques secondes ; la conception coûte des dizaines de minutes.
"""),
        code(f"""
# Paramètres
NUM_DESIGNS = {num_designs}   # designs generes au total, puis filtres par le pipeline
RUN = "{run}"
PROTOCOL = "nanobody-anything"
STEPS = {steps!r}   # None = les cinq etapes

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

# Kaggle exports MPLBACKEND=module://matplotlib_inline.backend_inline. BoltzGen
# imports biotite, which imports matplotlib, which reads that at import time and
# refuses it because matplotlib_inline lives in Kaggle's python and not in this
# venv. scripts/boltzgen_score.py already carried this fix in subprocess_env();
# writing a bare subprocess.run here lost it and cost a Kaggle run.
BG_ENV = {{**os.environ, "MPLBACKEND": "Agg"}}

def sh(cmd, **kw):
    kw.setdefault("env", BG_ENV)
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
# LA BoltzGen ORIGINALE, pas celle des validateurs. Le guide amont le dit:
# "l implementation incluse dans le depot NOVA est adaptee pour servir de
# moteur de NOTATION et non de pipeline generatif. Vous utiliserez
# l implementation originale pour creer vos designs."
# Mesure a l appui (FINDINGS 55): la 0.2.0 embarquee manque 173 lignes de
# inverse_fold.py, dont build_constraint_logit_mask, donc l evitement de la
# cysteine impose par filter_cysteine=true s applique AUSSI aux positions que
# la specification a gelees. Resultat des deux runs precedents: 210 designs
# sur 210 a zero cysteine, Cys23 du pont disulfure canonique remplacee par un
# glutamate, et des sequences que min_cysteines aurait refusees.
step("installation de la BoltzGen ORIGINALE 0.3.2 (conception)")
sh(f"uv pip install -q --python {py} boltzgen==0.3.2 pandas requests")
site = pathlib.Path(sh([py, "-c", "import sysconfig; print(sysconfig.get_paths()['purelib'])"]).strip())
sh(["cp", str(WORK / "kaggle" / "bittensor.py"), str(site / "bittensor.py")])
# ASSERTION, pas affichage. Le lancement precedent a imprime
# "2.7.1+cu126 False" a 80 s puis a continue, pour mourir 20 s plus tard dans
# configure sur "Found no NVIDIA driver". L accelerateur est un reglage PAR
# CARNET : un carnet neuf demarre sans GPU, et aucun code ne peut y remedier.
probe = sh([py, "-c", "import torch; print(torch.__version__, torch.cuda.is_available(), "
                      "torch.cuda.device_count())"]).split()
if probe[1] != "True":
    print("AUCUN GPU VISIBLE : torch.cuda.is_available() = False")
    print("Dans l editeur Kaggle : Settings -> Accelerator -> GPU T4 x2, puis relancer.")
    print("BoltzGen reclame le GPU des l etape configure, qui lit la capacite de")
    print("calcul du materiel : il n y a pas de repli CPU possible.")
    raise SystemExit("pas de GPU : regler l accelerateur sur GPU T4 x2")
sh(f'{py} -c "import boltzgen, bittensor"')
step(f"environnement pret : torch {probe[0]}, {probe[2]} GPU")
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
    print("=== LIRE CECI AVANT DE BRICOLER LA SPECIFICATION ===")
    print("Cette meme specification a ete validee hors ligne : boltzgen check")
    print("sort en 0 et ecrit design_spec.cif, chaine A 115 residus et chaine B")
    print("la charpente moins le milieu exclu. L indexation est verifiee contre")
    print("la structure elle-meme (verify_span), et le comportement est identique")
    print("a celui de l exemple amont 7eow.yaml.")
    print()
    print("Donc une erreur ici est d abord une erreur d ENVIRONNEMENT, pas de")
    print("specification. Lire la trace ci-dessus : si elle ne parle pas de")
    print("chain id, de res_index ou de design spec, ce n est pas le yaml.")
    print()
    print("Si la trace EST une erreur de specification, la piste restante est le")
    print("couple exclude / design_insertions : garder k residus du donneur dans")
    print("le span concu, exclure le reste, inserer designable_count-k, k petit.")
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
# Les cinq etapes du pipeline, ou un sous-ensemble. Nous notons avec la build
# des validateurs, donc le repliement et l'analyse de BoltzGen sont du cout pur :
# mesure sur les runs a 120 designs, design 2096 s et inverse_folding 67 s contre
# folding 5466 s -- trois quarts du GPU pour des metriques qu'on jette ensuite.
configure = [bg, "configure", "design_spec.yaml", "--output", str(out_dir),
             "--protocol", PROTOCOL, "--num_designs", str(NUM_DESIGNS)]
if STEPS:
    configure += ["--steps", *STEPS]
step(" ".join(configure))
sh(configure)

t_design = time.time()
step(f"conception de {NUM_DESIGNS} designs (run {RUN})")
with open(log_path, "w") as log:
    # env=BG_ENV, like every other call: this one does not go through sh(), and
    # omitting it cost a Kaggle run to the same MPLBACKEND error check had
    # already been fixed for.
    proc = subprocess.Popen([bg, "execute", str(out_dir)], stdout=log,
                            stderr=subprocess.STDOUT, env=BG_ENV)
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
# L arborescence AVANT la lecture : si le parsing echoue, papermill arrete le
# carnet et la cellule d inspection qui suivait ne tournait jamais -- on perdait
# exactement l information necessaire pour corriger la lecture.
print("=== arborescence des resultats ===")
sh(f"find {out_dir} -maxdepth 4 | head -80")

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
    # Avec --steps design inverse_folding il n'y a ni analyse ni filtrage, donc
    # aucun CSV : les sequences ne vivent que dans les structures. On les en sort.
    # Sans ce repli le carnet mourait ici apres des heures de GPU, sur un lot bon.
    print("aucun CSV : extraction des sequences depuis les structures")
    AA3 = dict(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E",
               GLY="G", HIS="H", ILE="I", LEU="L", LYS="K", MET="M", PHE="F",
               PRO="P", SER="S", THR="T", TRP="W", TYR="Y", VAL="V")
    def chains_of(path):
        out, seen = dict(), set()
        text = pathlib.Path(path).read_text(errors="ignore")
        if path.endswith(".cif"):
            cols, atoms, inloop = dict(), [], False
            for line in text.splitlines():
                if line.startswith("_atom_site."):
                    cols[line.strip().split(".")[1]] = len(cols); inloop = True
                elif inloop and line.startswith(("ATOM", "HETATM")):
                    atoms.append(line.split())
                elif inloop and line.startswith("#"):
                    inloop = False
            need = ("label_atom_id", "label_comp_id", "label_asym_id", "label_seq_id")
            if not all(k in cols for k in need):
                return out
            top = max(cols[k] for k in need)
            for f in atoms:
                if len(f) <= top or f[cols["label_atom_id"]].strip(chr(34)) != "CA":
                    continue
                ch, rid = f[cols["label_asym_id"]], f[cols["label_seq_id"]]
                if (ch, rid) in seen:
                    continue
                seen.add((ch, rid))
                out[ch] = out.get(ch, "") + AA3.get(f[cols["label_comp_id"]].upper(), "X")
        else:
            for line in text.splitlines():
                if not line.startswith("ATOM") or line[12:16].strip() != "CA":
                    continue
                ch, rid = line[21], line[22:27]
                if (ch, rid) in seen:
                    continue
                seen.add((ch, rid))
                out[ch] = out.get(ch, "") + AA3.get(line[17:20].strip().upper(), "X")
        return out
    # Choisir la bonne chaine. La cible fait 115 residus et nos nanobodies
    # 123-125, donc "la plus courte" prendrait la CIBLE. Deux discriminants,
    # dans l'ordre : la chaine B, que notre specification nomme pour la
    # charpente ; sinon la chaine qui VARIE d'un fichier a l'autre, la cible
    # etant identique dans tous.
    per_file = []
    for pat in ("**/intermediate_designs_inverse_folded/**/*", "**/*"):
        for f in sorted(glob.glob(str(out_dir / pat), recursive=True)):
            if f.endswith((".cif", ".pdb")):
                ch = chains_of(f)
                if ch:
                    per_file.append((f, ch))
        if per_file:
            break
    counts = dict()
    for _, ch in per_file:
        for seq in ch.values():
            counts[seq] = counts.get(seq, 0) + 1
    constant = set(seq for seq, n in counts.items()
                   if n > 1 and n >= len(per_file) * 0.5)
    found = []
    for f, ch in per_file:
        seq = ch.get("B")
        if seq is None:
            usable = [v for v in ch.values() if v not in constant]
            seq = max(usable, key=len) if usable else None
        if seq and 90 <= len(seq) <= 150 and "X" not in seq:
            found.append(dict(designed_chain_sequence=seq,
                              id=pathlib.Path(f).stem,
                              _source=pathlib.Path(f).name))
    print(f"{len(per_file)} structures lues, {len(constant)} sequence(s) constante(s) "
          f"ecartee(s) comme cible")
    if not found:
        raise SystemExit("aucun CSV et aucune structure exploitable")
    rows = [pd.DataFrame(found).drop_duplicates("designed_chain_sequence")]
    print(f"{len(rows[0])} sequences distinctes extraites des structures")
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

    # Every code cell must parse. An escaping slip in a generated cell -- a \\n
    # that the enclosing literal turned into a real newline -- produced a cell
    # with an unterminated string, which only showed up when the cell ran.
    import ast
    for i, (cell, source) in enumerate(zip(cells, sources)):
        if cell["cell_type"] != "code":
            continue
        try:
            ast.parse(source)
        except SyntaxError as error:
            raise AssertionError(
                f"cell {i} is not valid Python: {error}") from None

    # Every subprocess call must carry BG_ENV. Fixing sh() was not enough: the
    # design step uses Popen directly, so it kept Kaggle's MPLBACKEND and died
    # on the error check had already been fixed for. A claim of coverage is
    # worth nothing unless something checks it.
    import re
    for i, source in enumerate(sources):
        for line in source.splitlines():
            if "subprocess.Popen(" in line or "subprocess.run(" in line:
                statement = source[source.index(line):]
                statement = statement[:statement.index(")\n") + 1
                                      if ")\n" in statement else len(statement)]
                through_sh = "kw.setdefault" in source and "def sh(" in source
                if "env=" not in statement and not through_sh:
                    raise AssertionError(
                        f"cell {i} calls subprocess without env=BG_ENV:\n  "
                        f"{line.strip()}\nKaggle's MPLBACKEND breaks every "
                        f"BoltzGen import that inherits it.")
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
    parser.add_argument("--steps", nargs="+",
                        choices=["design", "inverse_folding", "design_folding",
                                 "folding", "affinity", "analysis", "filtering"],
                        help="etapes BoltzGen ; omettre pour les cinq. "
                             "'design inverse_folding' rend les sequences seules, "
                             "un quart du temps GPU.")
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

    nb = build(args.spec, args.nova_dir, args.target, args.run,
               args.num_designs, args.steps)
    out = REPO_ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(nb))
    print(f"wrote {out} ({out.stat().st_size // 1024} kB)")


if __name__ == "__main__":
    main()
