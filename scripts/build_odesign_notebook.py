"""A Kaggle notebook that designs our CDR3 with ODesign instead of BoltzGen.

Why. FINDINGS 64 measured BoltzGen's design ceiling at design_iiptm 0.698, and
showed it does not move with scale: eight times more designs and four times more
scaffolds left the bound unchanged to a hundredth. That leaves one question the
repo cannot answer with BoltzGen alone, and FINDINGS 65 states it:

  - is 0.698 **BoltzGen's** ceiling, so another model under the same constraint
    would do better;
  - or is it the **constraint's** -- a frozen VHH framework with 9 to 20 free
    residues cannot pass 0.70, whatever designs it?

ODesign (The-Institute-for-AI-Molecular-Design) separates the two. It is a
different model, from a different lab, and its `partial_diff` mode does exactly
what our BoltzGen spec does: keep a structure, redesign named residues.

What this uses that nothing else has. FINDINGS 66 found that `8AOK.cif`, which
this repo has always passed as "the target", is a **complex**: chain A is PD-L1
and chain B is a crystallographic VHH bound to it, CDR3 `AAEGTLYGSGGRTHQSAYDY`
at residues B/98-117. So the reference complex ODesign's partial diffusion needs
has been on disk all along, with the pose measured rather than predicted --
where BoltzGen had to invent it.

What comes back, and what it is not. The designed chain keeps chain B's
framework, whose germline identity is V 78.6 and J 85.7: **below the validator's
floors**, so no output of this notebook is submittable as it stands. What it
produces is CDR3 loops designed in the true pose against the true epitope, to be
grafted onto compliant frameworks afterwards. Whether that graft keeps the
binding is what the `xtal_26301` batch measures on the natural CDR3, and this
notebook is worth running only if that one says yes.

Usage:
    python -m scripts.build_odesign_notebook --n-sample 40 --seeds 42,43
"""
from __future__ import annotations

import argparse
import base64
import io
import json
import tarfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# The complex, and the CDR3 inside its nanobody chain (FINDINGS 66).
TARGET_CHAIN, TARGET_SPAN = "A", "1-115"
DESIGN_CHAIN, DESIGN_SPAN = "B", "2-128"
CDR3_SPAN = "98-117"
# The epitope, verified against the complex: all 21 of these residues are within
# 5 A of chain B, and they cover 91 % of the real interface (FINDINGS 66).
BINDING_SITE = "37,39,41,42,43,44,46,49,51,56,59,94,95,96,98,100,104,105,106,107,108"

ODESIGN_REPO = "https://github.com/The-Institute-for-AI-Molecular-Design/ODesign.git"
# Protenix publishes the CCD files ODesign requires, at a direct URL. The
# project's own README points at a Google Drive FOLDER, which a notebook cannot
# fetch reliably: no per-file id. Both were verified by range request --
# 408 MB and 120 MB, and the CIF starts with `data_000 / _chem_comp.id`.
CCD_BASE = "https://af3-dev.tos-cn-beijing.volces.com/release_data"
CCD_FILES = ["components.v20240608.cif", "components.v20240608.cif.rdkit_mol.pkl"]
HF = "https://huggingface.co"
CKPTS = {
    "odesign_base_prot_rigid.pt":
        f"{HF}/The-Institute-for-AI-Molecular-Design/ODesign/resolve/main/ckpt/odesign_base_prot_rigid.pt",
    "oinvfold_protein.ckpt":
        f"{HF}/The-Institute-for-AI-Molecular-Design/OInvFold/resolve/main/oinvfold_protein.ckpt",
    # src/model/modules/invfold/proteinmpnn/wrapper.py reads this name verbatim.
    "v_48_020.pt":
        "https://raw.githubusercontent.com/dauparas/ProteinMPNN/main/vanilla_proteinmpnn/vanilla_model_weights/v_48_020.pt",
}
# torch 2.3.1 + cu121 is what requirements.txt pins, and torch_scatter,
# torch_cluster and torch_geometric are imported by 12 files of src/. Their
# wheels exist for pt23cu121 on cp310, cp311 and cp312 -- checked against
# data.pyg.org -- so the venv is built on 3.12 and the wheels come from there
# rather than being compiled.
TORCH = "2.3.1"
CUDA = "cu121"
PYG_INDEX = f"https://data.pyg.org/whl/torch-{TORCH}+{CUDA}.html"


def payload(nova_dir: Path) -> str:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:xz") as tar:
        tar.add(nova_dir / "data" / "structures" / "8AOK.cif", arcname="8AOK.cif")
    return base64.b64encode(buf.getvalue()).decode()


def md(text: str) -> dict:
    return {"cell_type": "markdown", "metadata": {},
            "source": text.strip("\n").splitlines(True)}


def code(text: str) -> dict:
    return {"cell_type": "code", "metadata": {}, "execution_count": None,
            "outputs": [], "source": text.strip("\n").splitlines(True)}


def build(nova_dir: Path, run: str, n_sample: int, seeds: list[int],
          snr: float) -> dict:
    cells = [
        md(f"""
# NOVA SN68 — conception du CDR3 par **ODesign** (GPU Kaggle)

Ne note rien : fabrique des boucles CDR3. Et pas avec BoltzGen — avec un autre
modèle, pour répondre à la question que le §64 laisse ouverte : le plafond de
**0,698** est-il celui de BoltzGen, ou celui de la contrainte ?

**Ce que ce carnet utilise et que rien d'autre n'a utilisé.** `8AOK.cif` est un
*complexe* : chaîne A = PD-L1, chaîne B = un VHH cristallographique lié à elle
(§66). La pose est donc **mesurée**, là où BoltzGen devait l'inventer. ODesign
redessine le CDR3 de la chaîne B *in situ*, par `partial_diff` sur
`B/{CDR3_SPAN}`, contre l'épitope vérifié (21 résidus, tous à moins de 5 Å de la
chaîne B, couvrant 91 % de l'interface réelle).

**Ce qui en sort n'est pas soumissible tel quel** : la charpente de la chaîne B
est à V 78,6 et J 85,7, sous les planchers du validateur. Le produit utile est
la **boucle**, à greffer ensuite sur une charpente conforme — ce que le lot
`xtal_26301` mesure sur le CDR3 naturel. Ne lancez ce carnet que si celui-là dit
que la greffe conserve la liaison.

Run : `{run}` · {n_sample} échantillons par graine · graines {seeds} · SNR {snr}
"""),
        code(f'''
RUN = "{run}"
N_SAMPLE = {n_sample}
SEEDS = "{json.dumps(seeds).replace(' ', '')}"
SNR = {snr}
MODEL = "odesign_base_prot_rigid"
PAYLOAD = "{payload(nova_dir)}"

import base64, io, json, os, pathlib, subprocess, sys, tarfile, time
TEMP = pathlib.Path("/kaggle/temp"); OUT = pathlib.Path("/kaggle/working")
WORK = TEMP / "odesign_work"; CKPT = TEMP / "ckpt"; DATA = TEMP / "ccd"
for d in (WORK, CKPT, DATA):
    d.mkdir(parents=True, exist_ok=True)
with tarfile.open(fileobj=io.BytesIO(base64.b64decode(PAYLOAD))) as t:
    t.extractall(WORK)
T0 = time.time()
def step(msg):
    print(f"\\n[{{time.time()-T0:7.0f}}s] === {{msg}}", flush=True)
def sh(cmd, **kw):
    if isinstance(cmd, str): cmd = cmd.split()
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        print(r.stdout[-3000:]); print(r.stderr[-3000:])
        raise SystemExit(f"echec: {{' '.join(map(str, cmd))}}")
    return r.stdout
step(f"charge utile depliee dans {{WORK}}")
print(sh(f"ls -la {{WORK}}"))
'''),
        md("""
## 1. L'environnement

`torch 2.3.1+cu121`, puis les roues PyG **préconstruites** pour `pt23cu121`
(elles existent en cp312, vérifié sur `data.pyg.org` : sans elles
`torch_scatter` se compile pendant une heure), puis les dépendances d'ODesign.
L'accélérateur est une **assertion**, pas un affichage : un carnet neuf démarre
sans GPU et ODesign meurt au chargement du modèle.
"""),
        code(f'''
step("uv + venv python 3.12")
sh(f"{{sys.executable}} -m pip install -q uv")
VENV = TEMP / "odvenv"
sh(f"uv venv --python 3.12 {{VENV}}")
py = str(VENV / "bin/python")

step("torch {TORCH}+{CUDA}")
sh(f"uv pip install -q --python {{py}} torch=={TORCH} torchvision torchaudio "
   f"--index-url https://download.pytorch.org/whl/{CUDA}")

step("roues PyG preconstruites (pas de compilation)")
sh(f"uv pip install -q --python {{py}} torch_scatter torch_cluster -f {PYG_INDEX}")
sh(f"uv pip install -q --python {{py}} torch_geometric")

step("ODesign")
sh(f"git clone --depth 1 {ODESIGN_REPO} {{WORK / 'ODesign'}}")
REQ = WORK / "ODesign" / "requirements.txt"
# torch, les paquets PyG et les roues nvidia sont deja poses, et jupyter n'a
# rien a faire ici : les reinstaller casserait la version CUDA en place.
skip = ("torch", "torchvision", "torchaudio", "triton", "nvidia-", "pyg-lib",
        "torch-scatter", "torch-sparse", "torch-cluster", "torch-spline-conv",
        "torch-geometric", "jupyter", "notebook", "ipykernel", "ipywidgets",
        "jupyterlab")
keep = [l.strip() for l in REQ.read_text().splitlines()
        if l.strip() and not l.startswith("#")
        and not any(l.lower().startswith(s) for s in skip)]
(TEMP / "req.txt").write_text("\\n".join(keep) + "\\n")
sh(f"uv pip install -q --python {{py}} -r {{TEMP / 'req.txt'}}")

probe = sh([py, "-c", "import torch; print(torch.__version__, "
                      "torch.cuda.is_available(), torch.cuda.device_count())"]).split()
print("torch:", probe)
if probe[1] != "True":
    print("AUCUN GPU VISIBLE. Settings -> Accelerator -> GPU, puis relancer.")
    raise SystemExit("pas de GPU")
sh([py, "-c", "import torch_scatter, torch_cluster, torch_geometric, protenix"])
step(f"environnement pret : torch {{probe[0]}}, {{probe[2]}} GPU")
'''),
        md(f"""
## 2. Les poids et le dictionnaire de composants

Trois fichiers de poids, aux noms que le code lit **verbatim**
(`src/utils/inference/infer_runner.py` et le wrapper ProteinMPNN) : le modèle
ODesign ({'4,32 Go'}), `oinvfold_protein.ckpt`, et `v_48_020.pt` de ProteinMPNN.

Et les deux fichiers CCD, que `src/utils/data/ccd.py` charge **sans condition** :
le README d'ODesign les met sur un *dossier* Google Drive, qu'un carnet ne peut
pas récupérer faute d'identifiant par fichier. Ils viennent donc de Protenix en
HTTP direct (408 Mo et 120 Mo, vérifiés par requête de plage).
"""),
        code('''
CKPT_URLS = ''' + json.dumps(CKPTS, indent=4) + '''
def fetch(url, dest):
    if dest.exists() and dest.stat().st_size > 1_000_000:
        print(f"  deja la : {dest.name} ({dest.stat().st_size/1e9:.2f} Go)"); return
    sh(["curl", "-sSL", "--retry", "4", "--retry-delay", "5", "-C", "-",
        "-o", str(dest), url])
    print(f"  {dest.name} : {dest.stat().st_size/1e9:.2f} Go")

step("poids")
for name, url in CKPT_URLS.items():
    fetch(url, CKPT / name)
step("dictionnaire de composants (CCD)")
for name in ''' + json.dumps(CCD_FILES) + ''':
    fetch(f"{''' + f'"{CCD_BASE}"' + '''}/{name}", DATA / name)
head = (DATA / "components.v20240608.cif").read_text(errors="ignore")[:40]
assert head.startswith("data_"), f"CCD inattendu : {head!r}"
print("CCD verifie :", head.splitlines()[0])
'''),
        md(f"""
## 3. La spécification

`partial_diff` sur `B/{CDR3_SPAN}` — les 20 résidus du CDR3 — avec la chaîne A
et le reste de la chaîne B conservés depuis la structure, et le `hotspot` sur
l'épitope vérifié. C'est la traduction exacte de notre spécification BoltzGen,
à une différence près qui est tout l'intérêt : **la pose vient de la structure**.
"""),
        code(f'''
spec = [{{
    "name": "pdl1_nb_cdr3",
    "ref_file": str(WORK / "8AOK.cif"),
    "chains": [
        {{"chain_type": "proteinChain", "sequence": "{TARGET_CHAIN}/{TARGET_SPAN}"}},
        {{"chain_type": "proteinChain", "sequence": "{DESIGN_CHAIN}/{DESIGN_SPAN}"}},
    ],
    "hotspot": ",".join(f"{TARGET_CHAIN}/{{r}}" for r in "{BINDING_SITE}".split(",")),
    "partial_diff": "{DESIGN_CHAIN}/{CDR3_SPAN}",
    "center_method": "hotspot_center",
}}]
spec_path = WORK / "odesign_input.json"
spec_path.write_text(json.dumps(spec, indent=2))
print(spec_path.read_text())
'''),
        md("""
## 4. L'inférence

Les surcharges hydra d'`inference_demo.sh`, sans le script shell. La diffusion
partielle est **activée** : c'est elle qui ne rebruite que le CDR3 et laisse le
reste du complexe en place.
"""),
        code('''
log = OUT / f"odesign_{RUN}.log"
cmd = [py, "./scripts/inference.py",
       f"exp=train_{MODEL}",
       f"data_root_dir={DATA}",
       f"ckpt_root_dir={CKPT}",
       f"exp.infer_model_name={MODEL}",
       "exp.design_modality=protein",
       f"exp.input_json_path={spec_path}",
       f"exp.exp_name={RUN}",
       f"exp.seeds={SEEDS}",
       f"exp.model.sample_diffusion.N_sample={N_SAMPLE}",
       "exp.use_msa=false",
       "exp.num_workers=2",
       "exp.model.inference_noise_schedulers.coordinate.partial_diffusion.enable=true",
       f"exp.model.inference_noise_schedulers.coordinate.partial_diffusion.snr={SNR}"]
step(" ".join(cmd[1:6]) + " ...")
env = dict(os.environ, DATA_ROOT_DIR=str(DATA), CKPT_ROOT_DIR=str(CKPT),
           MPLBACKEND="Agg", HYDRA_FULL_ERROR="1")
with open(log, "w") as fh:
    proc = subprocess.Popen(cmd, cwd=str(WORK / "ODesign"), stdout=fh,
                            stderr=subprocess.STDOUT, env=env)
    while proc.poll() is None:
        time.sleep(120)
        print(f"  [{time.time()-T0:7.0f}s] en cours, log {log.stat().st_size/1024:.0f} ko",
              flush=True)
print(open(log).read()[-4000:])
if proc.returncode != 0:
    raise SystemExit(f"inference en echec (code {proc.returncode}), voir {log}")
step("inference terminee")
'''),
        md("""
## 5. Récupérer les séquences

ODesign écrit des mmCIF, pas de CSV. Le même extracteur que le carnet de
conception BoltzGen : les CA, puis la chaîne **qui varie** d'un fichier à
l'autre — la cible étant identique partout. Prendre « la plus courte » prendrait
la cible (115 résidus contre 127).

Le criblage par nos portes se fait **en local** : ANARCI n'est pas ici.
"""),
        code('''
import glob
AA3 = dict(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E",
           GLY="G", HIS="H", ILE="I", LEU="L", LYS="K", MET="M", PHE="F",
           PRO="P", SER="S", THR="T", TRP="W", TYR="Y", VAL="V")
def chains_of(path):
    out, seen = dict(), set()
    text = pathlib.Path(path).read_text(errors="ignore")
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
    return out

roots = sorted(glob.glob(str(WORK / "ODesign" / "outputs" / "**" / "*.cif"),
                         recursive=True))
print(f"{len(roots)} mmCIF produits")
per_file = [(f, chains_of(f)) for f in roots]
per_file = [(f, c) for f, c in per_file if c]
counts = dict()
for _, ch in per_file:
    for seq in ch.values():
        counts[seq] = counts.get(seq, 0) + 1
constant = set(s for s, n in counts.items() if n > 1 and n >= len(per_file) * 0.5)
rows = []
for f, ch in per_file:
    seq = ch.get("B")
    if seq is None:
        usable = [v for v in ch.values() if v not in constant]
        seq = max(usable, key=len) if usable else None
    if seq and 90 <= len(seq) <= 150 and "X" not in seq:
        rows.append(dict(designed_chain_sequence=seq, id=pathlib.Path(f).stem))
import pandas as pd
res = pd.DataFrame(rows).drop_duplicates("designed_chain_sequence")
dest = OUT / f"designs_{RUN}.csv"
res.to_csv(dest, index=False)
print(f"{len(per_file)} structures, {len(constant)} sequence(s) constante(s) ecartee(s)")
print(f"{len(res)} sequences distinctes -> {dest}")
if len(res):
    s = res.designed_chain_sequence
    print(f"longueurs : {sorted(set(s.str.len()))}")
    print("\\nA FAIRE EN LOCAL, avant toute lecture de scores :")
    print(f"  python -m scripts.verify_designs --file results/designs_{RUN}.csv")
res.head(10)
'''),
    ]
    return {"cells": cells, "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python",
                       "name": "python3"},
        "language_info": {"name": "python"}}, "nbformat": 4, "nbformat_minor": 5}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--nova-dir", type=Path, default=Path("/home/user/nova-upstream"))
    parser.add_argument("--run", default="odesign_cdr3_8aok")
    parser.add_argument("--n-sample", type=int, default=40,
                        help="echantillons par graine")
    parser.add_argument("--seeds", default="42,43",
                        help="graines separees par des virgules")
    parser.add_argument("--snr", type=float, default=0.1,
                        help="rapport signal/bruit de la diffusion partielle")
    parser.add_argument("--out", type=Path,
                        default=REPO_ROOT / "notebooks" / "kaggle_odesign.ipynb")
    args = parser.parse_args()

    structure = args.nova_dir / "data" / "structures" / "8AOK.cif"
    if not structure.exists():
        raise SystemExit(f"structure absente : {structure}. Passer --nova-dir.")
    seeds = [int(s) for s in args.seeds.split(",") if s.strip()]
    nb = build(args.nova_dir, args.run, args.n_sample, seeds, args.snr)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(nb, indent=1))
    size = args.out.stat().st_size / 1024
    print(f"wrote {args.out} ({size:.0f} kB)")
    print(f"  run {args.run}, {args.n_sample} echantillons x {len(seeds)} graines")
    if size > 1024:
        print("  ATTENTION : Kaggle refuse un carnet de plus de 1 Mo a l'import")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
