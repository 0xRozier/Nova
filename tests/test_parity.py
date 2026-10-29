"""Parity checks against what validators actually computed.

The public archive stores, for each submitted molecule, the SMILES the
validator derived from its name; every archived molecule also passed the
validator's filters (with the thresholds active when it was submitted).

    python -m pytest tests/ -q
"""
import random

import pytest
from rdkit import Chem

from nova_lab.chem import check_molecule, load_archive
from nova_lab.reactions import name_to_smiles

SAMPLE = 400


@pytest.fixture(scope="module")
def archive():
    df = load_archive("O15427")
    rows = list(df[df["Molecule_ID"].str.startswith("rxn:")].itertuples(index=False))
    random.Random(0).shuffle(rows)
    return rows[:SAMPLE]


def canon(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol else None


def test_names_resolve_to_validator_smiles(archive):
    mismatches = [(r.Molecule_ID, r.SMILES, name_to_smiles(r.Molecule_ID)) for r in archive
                  if canon(name_to_smiles(r.Molecule_ID) or "") != canon(r.SMILES)]
    assert not mismatches, f"{len(mismatches)}/{len(archive)} mismatches, e.g. {mismatches[:3]}"


def test_archived_molecules_pass_filters(archive):
    failures = [(r.Molecule_ID, reason) for r in archive
                for ok, reason in [check_molecule(r.Molecule_ID, smiles=r.SMILES)] if not ok]
    # Thresholds changed over time, so a few old molecules may fail today's filters.
    assert len(failures) <= 0.02 * len(archive), failures[:10]
