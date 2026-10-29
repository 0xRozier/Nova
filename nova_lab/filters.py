"""Molecular plausibility filters, reimplemented from upstream utils/filters.py.

Evaluated in the upstream order so the first failure reported is the same one
the validator logs.
"""
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, GraphDescriptors, Lipinski, rdMolDescriptors

from .config import FILTERS as THRESHOLDS

METALS = frozenset("Fe Zn Cu Mg Ca Mn Mo Ni Co Pt Au Hg Pd Ru W V Cr Cd Pb As Sb Bi Sn Al Li Na K".split())
HALOGENS = frozenset("F Cl Br I".split())
PATTERNS = {k: Chem.MolFromSmarts(v) for k, v in {
    "nitrile": "[CX2]#[NX1]",
    "azide": "[NX1]~[NX2]~[NX2,NX3]",
    "isocyanide": "[CX1-]#[NX2+]",
    "allene": "[CX3]=[CX2]=[CX3]",
    "hydroxyl": "[OX2H]",
}.items()}


def _counts(mol):
    counts, aromatic, aromatic_n, aromatic_c = {}, 0, 0, 0
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        counts[sym] = counts.get(sym, 0) + 1
        if atom.GetIsAromatic():
            aromatic += 1
            aromatic_n += sym == "N"
            aromatic_c += sym == "C"
    return counts, aromatic, aromatic_n, aromatic_c


def _group(mol, key):
    return len(mol.GetSubstructMatches(PATTERNS[key]))


def _checks(mol):
    """Yield (name, value, ok) lazily, in upstream order."""
    t = THRESHOLDS
    counts, aromatic, aromatic_n, aromatic_c = _counts(mol)
    el = lambda s: counts.get(s, 0)  # noqa: E731
    heavy = mol.GetNumHeavyAtoms()
    hetero_frac = (heavy - el("C")) / max(heavy, 1)
    largest_ring = max((len(r) for r in mol.GetRingInfo().AtomRings()), default=0)
    triple = sum(1 for b in mol.GetBonds() if b.GetBondType() == Chem.BondType.TRIPLE)
    mw = Descriptors.MolWt(mol)
    rot = Descriptors.NumRotatableBonds(mol)

    yield "kier kappa1", GraphDescriptors.Kappa1(mol), lambda v: v <= t["max_kappa1"]
    yield "molecular weight upper", mw, lambda v: v <= t["max_molecular_weight"]
    yield "heavy atom count", heavy, lambda v: v <= t["max_heavy_atoms"]
    yield "minimum heavy atoms", heavy, lambda v: v >= t["min_heavy_atoms"]
    yield "molecular weight lower", mw, lambda v: v >= t["min_molecular_weight"]
    yield "aromatic ring count", Descriptors.NumAromaticRings(mol), lambda v: v <= t["max_aromatic_rings"]
    yield "aliphatic ring count", Descriptors.NumAliphaticRings(mol), lambda v: v <= t["max_aliphatic_rings"]
    yield "ring count minimum", Descriptors.RingCount(mol), lambda v: v >= t["min_rings"]
    yield "fsp3", rdMolDescriptors.CalcFractionCSP3(mol), lambda v: v <= t["max_fsp3"]
    yield "spiro atoms", rdMolDescriptors.CalcNumSpiroAtoms(mol), lambda v: v <= t["max_spiro_atoms"]
    yield "macrocycle", largest_ring, lambda v: v <= t["max_ring_size"]
    yield "aromatic nitrogen", aromatic_n, lambda v: v <= t["max_aromatic_nitrogen"]
    yield "aromatic atom fraction", aromatic / max(heavy, 1), lambda v: v <= t["max_aromatic_fraction"]
    yield "aromatic carbon", aromatic_c, lambda v: v <= t["max_aromatic_carbon"]
    yield "total halogens", sum(n for s, n in counts.items() if s in HALOGENS), lambda v: v <= t["max_halogens"]
    yield "nitrogen count", el("N"), lambda v: v <= t["max_nitrogen"]
    yield "oxygen count", el("O"), lambda v: v <= t["max_oxygen"]
    yield "sulfur count", el("S"), lambda v: v <= t["max_sulfur"]
    yield "phosphorus count", el("P"), lambda v: v <= t["max_phosphorus"]
    yield "iodine count", el("I"), lambda v: v <= t["max_iodine"]
    yield "bromine count", el("Br"), lambda v: v <= t["max_bromine"]
    yield "fluorine count", el("F"), lambda v: v <= t["max_fluorine"]
    yield "metal atoms", sum(n for s, n in counts.items() if s in METALS), lambda v: v == 0
    yield "boron", el("B"), lambda v: v == 0
    yield "silicon", el("Si"), lambda v: v == 0
    yield "selenium", el("Se"), lambda v: v == 0
    yield "heteroatom fraction", hetero_frac, \
        lambda v: t["min_heteroatom_fraction"] <= v <= t["max_heteroatom_fraction"]
    yield "h-bond donors", Lipinski.NumHDonors(mol), lambda v: v <= t["max_hbd"]
    yield "amide bond count", rdMolDescriptors.CalcNumAmideBonds(mol), lambda v: v <= t["max_amide_bonds"]
    yield "nitrile count", _group(mol, "nitrile"), lambda v: v <= t["max_nitrile"]
    yield "hydroxyl count", _group(mol, "hydroxyl"), lambda v: v <= t["max_hydroxyl"]
    yield "azide", _group(mol, "azide"), lambda v: v == 0
    yield "isocyanide", _group(mol, "isocyanide"), lambda v: v == 0
    yield "allene", _group(mol, "allene"), lambda v: v == 0
    yield "triple bond count", triple, lambda v: v <= t["max_triple_bonds"]
    yield "rotatable bonds lower", rot, lambda v: v >= t["min_rotatable_bonds"]
    yield "rotatable bonds upper", rot, lambda v: v <= t["max_rotatable_bonds"]
    yield "bertzCT complexity", GraphDescriptors.BertzCT(mol), lambda v: v <= t["max_bertz_ct"]
    yield "h-bond acceptors", Lipinski.NumHAcceptors(mol), lambda v: v <= t["max_hba"]
    yield "kier kappa2", GraphDescriptors.Kappa2(mol), lambda v: v <= t["max_kappa2"]
    yield "kier kappa3", GraphDescriptors.Kappa3(mol), lambda v: v <= t["max_kappa3"]
    yield "stereocenter count", len(Chem.FindMolChiralCenters(mol, includeUnassigned=True)), \
        lambda v: v <= t["max_stereocenters"]
    yield "logP (crippen)", Crippen.MolLogP(mol), lambda v: v <= t["max_logp"]


def first_failure(mol) -> str | None:
    for name, value, ok in _checks(mol):
        if not ok(value):
            return f"{name} = {value:.4g}"
    return None
