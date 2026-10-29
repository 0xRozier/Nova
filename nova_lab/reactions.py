"""Resolve combinatorial-DB molecule names (rxn:R:A:B[:C][:kK]) to SMILES.

Reimplements upstream combinatorial_db/reactions.py; product order and the
k-th product selection must stay identical or names would map to other molecules.
"""
import sqlite3
from functools import lru_cache

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from . import DB_PATH

RDLogger.DisableLog("rdApp.*")

AMIDE_SMARTS = "[C:1](=O)[OH].[N:2]>>[C:1](=O)[N:2]"
SUZUKI_BR_SMARTS = "[#6:1][Br].[#6:2][B]([OH])[OH]>>[#6:1][#6:2]"
SUZUKI_CL_SMARTS = "[#6:1][Cl].[#6:2][B]([OH])[OH]>>[#6:1][#6:2]"


def as_mask(value) -> int:
    return int.from_bytes(bytes(value), "big") if value is not None else 0


def has_role(role_mask: int, role: int) -> bool:
    """A molecule fills a role only when it carries every bit the role requires."""
    return role != 0 and (role_mask & role) == role


@lru_cache(maxsize=1)
def _conn():
    return sqlite3.connect(f"file:{DB_PATH}?mode=ro", uri=True, check_same_thread=False)


@lru_cache(maxsize=None)
def reaction_info(rxn_id: int) -> tuple | None:
    """(name, smarts, roleA, roleB, roleC)."""
    row = _conn().execute(
        "SELECT name, smarts, roleA, roleB, roleC FROM reactions WHERE rxn_id = ?", (str(rxn_id),)
    ).fetchone()
    if row is None:
        return None
    name, smarts, a, b, c = row
    return name, smarts, as_mask(a), as_mask(b), as_mask(c)


def reaction_ids() -> list[int]:
    return sorted(int(r[0]) for r in _conn().execute("SELECT rxn_id FROM reactions"))


@lru_cache(maxsize=200_000)
def molecule(mol_id: int) -> tuple | None:
    """(smiles, role_mask)."""
    row = _conn().execute("SELECT smiles, role_mask FROM molecules WHERE mol_id = ?", (mol_id,)).fetchone()
    return (row[0], as_mask(row[1])) if row else None


def molecules_for_role(role: int) -> list[tuple[int, str]]:
    """[(mol_id, smiles)] of every reagent able to fill `role`."""
    return [(mid, smi) for mid, smi, mask in _conn().execute("SELECT mol_id, smiles, role_mask FROM molecules")
            if has_role(as_mask(mask), role)]


@lru_cache(maxsize=64)
def _rxn(smarts: str):
    return AllChem.ReactionFromSmarts(smarts)


def enumerate_smarts_products(smiles1: str, smiles2: str, smarts: str) -> list[str]:
    """Distinct products in RunReactants order (not sorted, so k stays stable)."""
    rxn = _rxn(smarts)
    mol1, mol2 = Chem.MolFromSmiles(smiles1), Chem.MolFromSmiles(smiles2)
    if rxn is None or mol1 is None or mol2 is None:
        return []
    found = []
    for product_set in rxn.RunReactants((mol1, mol2)):
        if not product_set:
            continue
        candidate = product_set[0]
        try:
            Chem.SanitizeMol(candidate)
        except Exception:
            continue
        smiles = Chem.MolToSmiles(candidate)
        if smiles not in found:
            found.append(smiles)
    return found


def smarts_product(smiles1: str, smiles2: str, smarts: str, k: int = 0) -> str | None:
    products = enumerate_smarts_products(smiles1, smiles2, smarts)
    return products[k] if k < len(products) else None


def combine_triazole_synthons(azide_smiles: str, alkyne_smiles: str) -> str | None:
    """Join the [1*] azide synthon and the [2*] alkyne synthon."""
    try:
        m1 = Chem.RWMol(Chem.MolFromSmiles(azide_smiles))
        m2 = Chem.RWMol(Chem.MolFromSmiles(alkyne_smiles))
        a1 = next((i for i, at in enumerate(m1.GetAtoms()) if at.GetSymbol() == "*" and at.GetIsotope() == 1), None)
        a2 = next((i for i, at in enumerate(m2.GetAtoms()) if at.GetSymbol() == "*" and at.GetIsotope() == 2), None)
        if a1 is None or a2 is None:
            return None
        n1 = m1.GetAtomWithIdx(a1).GetNeighbors()[0].GetIdx()
        n2 = m2.GetAtomWithIdx(a2).GetNeighbors()[0].GetIdx()

        combined = Chem.RWMol(m1)
        mapping = {}
        for i, atom in enumerate(m2.GetAtoms()):
            if i != a2:
                mapping[i] = combined.AddAtom(atom)
        for bond in m2.GetBonds():
            b, e = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if a2 not in (b, e):
                combined.AddBond(mapping[b], mapping[e], bond.GetBondType())

        combined.RemoveAtom(a1)
        n1_adj = n1 - (1 if n1 > a1 else 0)
        n2_adj = mapping[n2] - (1 if mapping[n2] > a1 else 0)
        combined.AddBond(n1_adj, n2_adj, Chem.BondType.SINGLE)
        Chem.SanitizeMol(combined)
        return Chem.MolToSmiles(combined)
    except Exception:
        return None


def order_reactants(m1: tuple, m2: tuple, role_a: int, role_b: int) -> tuple | None:
    (s1, r1), (s2, r2) = m1, m2
    if has_role(r1, role_a) and has_role(r2, role_b):
        return s1, s2
    if has_role(r1, role_b) and has_role(r2, role_a):
        return s2, s1
    return None


def react(rxn_id: int, mol_ids: list[int], k: int = 0) -> str | None:
    info = reaction_info(rxn_id)
    mols = [molecule(m) for m in mol_ids]
    if info is None or not all(mols):
        return None
    _, smarts, role_a, role_b, role_c = info
    ordered = order_reactants(mols[0], mols[1], role_a, role_b)
    if ordered is None:
        return None
    r1, r2 = ordered

    if len(mol_ids) == 2:
        if rxn_id == 1:  # triazole
            return combine_triazole_synthons(r1, r2) if k == 0 else None
        return smarts_product(r1, r2, smarts, k)

    r3, mask3 = mols[2]
    if not has_role(mask3, role_c):
        return None
    if rxn_id == 3:  # click_amide_cascade
        triazole_acid = combine_triazole_synthons(r1, r2)
        return smarts_product(triazole_acid, r3, AMIDE_SMARTS, k) if triazole_acid else None
    if rxn_id == 5:  # suzuki_bromide_then_chloride
        intermediate = smarts_product(r1, r2, SUZUKI_BR_SMARTS)
        return smarts_product(intermediate, r3, SUZUKI_CL_SMARTS, k) if intermediate else None
    return None


def parse_name(name: str) -> tuple | None:
    """rxn:<rxn>:<m1>:<m2>[:<m3>][:k<k>] -> (rxn_id, [mol_ids], k)."""
    parts = name.replace("'", "").replace('"', "").split(":")
    if len(parts) < 4 or parts[0] != "rxn":
        return None
    k = 0
    if parts[-1][:1] == "k":
        if not parts[-1][1:].isdigit():
            return None
        k = int(parts[-1][1:])
        parts = parts[:-1]
    if len(parts) not in (4, 5):
        return None
    try:
        return int(parts[1]), [int(p) for p in parts[2:]], k
    except ValueError:
        return None


def name_to_smiles(name: str) -> str | None:
    parsed = parse_name(name)
    if parsed is None:
        return None
    rxn_id, mol_ids, k = parsed
    try:
        return react(rxn_id, mol_ids, k)
    except Exception:
        return None


def make_name(rxn_id: int, mol_ids: list[int], k: int = 0) -> str:
    name = f"rxn:{rxn_id}:" + ":".join(str(m) for m in mol_ids)
    return f"{name}:k{k}" if k else name
