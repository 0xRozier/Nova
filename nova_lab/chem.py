"""Small-molecule validity checks, mirroring upstream molecule_validity.py."""
import math

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, rdFingerprintGenerator

from . import DATA_DIR
from .config import MAX_SIMILARITY_TO_HISTORICAL, MIN_ENTROPY, NUM_MOLECULES
from .filters import first_failure
from .reactions import name_to_smiles

RDLogger.DisableLog("rdApp.*")

ARCHIVE_URL = "https://huggingface.co/datasets/Metanova/Submission-Archive/resolve/main/{target}_{kind}.csv"
MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)  # novelty check
ATOM_PAIR = rdFingerprintGenerator.GetAtomPairGenerator(fpSize=2048)  # entropy


def heavy_atoms(smiles: str) -> int:
    mol = Chem.MolFromSmiles(smiles)
    return mol.GetNumHeavyAtoms() if mol else 0


def boltz_score(prob_binary: float, pred_value: float, n_heavy: int) -> float:
    """Validator score of one molecule ('heavy_atom_normalization', higher is better)."""
    return (prob_binary - pred_value) / n_heavy if n_heavy else -math.inf


def is_boltz_safe(smiles: str) -> tuple[bool, str | None]:
    """Boltz atom names must fit in 4 chars and RDKit must embed a conformer."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "RDKit failed to parse SMILES"
        mol = Chem.AddHs(mol)
        for atom, rank in zip(mol.GetAtoms(), Chem.CanonicalRankAtoms(mol)):
            if len(atom.GetSymbol().upper() + str(rank + 1)) > 4:
                return False, "atom name would exceed 4 chars"
        params = AllChem.ETKDGv3()
        if AllChem.EmbedMolecule(mol, params) == -1:
            params.useRandomCoords = True
            if AllChem.EmbedMolecule(mol, params) == -1:
                return False, "no conformer"
        return True, None
    except Exception as e:
        return False, str(e)


def fingerprint_entropy(smiles_list: list[str]) -> float | None:
    bits, n = np.zeros(2048), 0
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            bits += np.array(ATOM_PAIR.GetFingerprint(mol))
            n += 1
    if n == 0:
        return None
    p = bits / n
    with np.errstate(divide="ignore", invalid="ignore"):
        h = np.where((p > 0) & (p < 1), -p * np.log2(p) - (1 - p) * np.log2(1 - p), 0.0)
    return float(h.mean())


def load_archive(target: str, kind: str = "molecules", refresh: bool = False) -> pd.DataFrame:
    """Everything already submitted for `target` (the validator's 'historical submissions')."""
    path = DATA_DIR / "archive" / f"{target}_{kind}.csv"
    if refresh or not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        pd.read_csv(ARCHIVE_URL.format(target=target, kind=kind)).to_csv(path, index=False)
    return pd.read_csv(path)


class HistoricalIndex:
    """Archive fingerprints for the uniqueness and novelty checks."""

    def __init__(self, target: str, refresh: bool = False):
        df = load_archive(target, refresh=refresh)
        self.smiles = set(df["SMILES"])
        fps = [MORGAN.GetFingerprint(m) for m in map(Chem.MolFromSmiles, df["SMILES"]) if m is not None]
        self.fps = fps

    def max_similarity(self, fp) -> float:
        return max(DataStructs.BulkTanimotoSimilarity(fp, self.fps), default=0.0)

    def add(self, smiles: str):
        self.smiles.add(smiles)
        self.fps.append(MORGAN.GetFingerprint(Chem.MolFromSmiles(smiles)))


def check_molecule(name: str, allowed_reaction: str | None = None,
                   history: HistoricalIndex | None = None, smiles: str | None = None) -> tuple[bool, str]:
    """Per-molecule checks in the validator's order. Returns (ok, reason)."""
    if "~" in name:
        return False, "null token"
    if allowed_reaction and not name.startswith(allowed_reaction + ":"):
        return False, f"reaction is not {allowed_reaction}"
    smiles = smiles or name_to_smiles(name)
    if not smiles:
        return False, "no SMILES"
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "RDKit cannot parse"
    failure = first_failure(mol)
    if failure:
        return False, f"filter: {failure}"
    ok, reason = is_boltz_safe(smiles)
    if not ok:
        return False, f"boltz-unsafe: {reason}"
    if history is not None:
        if smiles in history.smiles:
            return False, "already submitted"
        sim = history.max_similarity(MORGAN.GetFingerprint(mol))
        if sim >= MAX_SIMILARITY_TO_HISTORICAL:
            return False, f"too similar to archive ({sim:.2f})"
    return True, ""


def check_submission(names: list[str], allowed_reaction: str | None = None,
                     history: HistoricalIndex | None = None) -> list[str]:
    """Whole-submission checks. Returns the problems found (empty list = valid)."""
    problems = []
    if len(names) != NUM_MOLECULES:
        problems.append(f"{len(names)} molecules, {NUM_MOLECULES} required")
    if len(set(names)) != len(names):
        problems.append("duplicate names")
    smiles = []
    for name in names:
        smi = name_to_smiles(name)
        ok, reason = check_molecule(name, allowed_reaction, history, smi)
        if not ok:
            problems.append(f"{name}: {reason}")
        if smi:
            smiles.append(smi)
    inchikeys = [Chem.MolToInchiKey(Chem.MolFromSmiles(s)) for s in smiles]
    if len(set(inchikeys)) != len(inchikeys):
        problems.append("chemically identical molecules (same InChIKey)")
    entropy = fingerprint_entropy(smiles) if len(smiles) > 1 else None
    if entropy is not None and entropy < MIN_ENTROPY:
        problems.append(f"entropy {entropy:.3f} < {MIN_ENTROPY}")
    return problems
