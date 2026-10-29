"""Molecule featurization for the CPU surrogate of the Boltz-2 score."""
import numpy as np
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdFingerprintGenerator, rdMolDescriptors

N_BITS = 2048
_MORGAN_COUNT = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=N_BITS)

DESCRIPTORS = {
    "heavy_atoms": lambda m: m.GetNumHeavyAtoms(),  # the score is divided by it
    "mw": Descriptors.MolWt,
    "logp": Crippen.MolLogP,
    "tpsa": rdMolDescriptors.CalcTPSA,
    "hbd": Lipinski.NumHDonors,
    "hba": Lipinski.NumHAcceptors,
    "rot_bonds": Descriptors.NumRotatableBonds,
    "rings": Descriptors.RingCount,
    "aromatic_rings": Descriptors.NumAromaticRings,
    "fsp3": rdMolDescriptors.CalcFractionCSP3,
    "charge": Chem.GetFormalCharge,
    "n_N": lambda m: sum(a.GetSymbol() == "N" for a in m.GetAtoms()),
    "n_O": lambda m: sum(a.GetSymbol() == "O" for a in m.GetAtoms()),
    "n_halogen": lambda m: sum(a.GetSymbol() in ("F", "Cl", "Br", "I") for a in m.GetAtoms()),
}
FEATURE_NAMES = [f"ecfp4_{i}" for i in range(N_BITS)] + list(DESCRIPTORS)


def featurize(smiles: str) -> np.ndarray | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = _MORGAN_COUNT.GetCountFingerprintAsNumPy(mol).astype(np.float32)
    desc = np.array([f(mol) for f in DESCRIPTORS.values()], dtype=np.float32)
    return np.concatenate([fp, desc])


def featurize_many(smiles_list) -> tuple[np.ndarray, np.ndarray]:
    """Returns (X, ok_mask) with one row per parseable SMILES."""
    rows, ok = [], []
    for smi in smiles_list:
        x = featurize(smi) if isinstance(smi, str) else None
        ok.append(x is not None)
        if x is not None:
            rows.append(x)
    return (np.vstack(rows) if rows else np.zeros((0, len(FEATURE_NAMES)), np.float32)), np.array(ok)
