"""Tests for the CDR3 permutation batch.

Tests the batch CSV structure, arm properties, CDR3 permutation invariants,
and the analyse/verdict logic. Run after generating the batch:

    python -m scripts.build_scramble_batch --out inputs/kaggle_scramble.csv
    python -m pytest tests/test_scramble_batch.py -v
"""
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from nova_lab import REPO_ROOT
from nova_lab.nanobodies import regions
from nova_lab.provenance import (assert_current_target, runs_on_current_target,
                                 target_of)
from scripts.build_scramble_batch import (
    score_path,
    JACCARD_CEILING,
    N_ANCHOR,
    N_HIGH,
    N_LOW,
    N_PERM_HIGH,
    N_PERM_LOW,
    OPERATIONAL_THRESHOLD,
    PRESERVATION_PASS,
    RHO_THRESHOLD,
    SIGMA_2,
    jaccard_5mer,
    kmers,
    permute_cdr3,
    verdict,
)

BATCH_PATH = REPO_ROOT / "inputs" / "kaggle_scramble.csv"
SCORE_GLOB = "scores_nanobody_*.csv"


@pytest.fixture(scope="module")
def batch():
    if not BATCH_PATH.exists():
        pytest.skip("batch not generated yet: run scripts.build_scramble_batch first")
    return pd.read_csv(BATCH_PATH)


# ---------------------------------------------------------------------------
# structural tests
# ---------------------------------------------------------------------------
def test_expected_columns(batch):
    required = {"sequence", "arm", "label", "group", "pair_id", "donor_iiptm",
                "donor_ptm", "cdr3_original", "cdr3_permuted", "cdr3_length",
                "jaccard_5mer", "layer1"}
    assert required <= set(batch.columns), f"missing: {required - set(batch.columns)}"


def test_total_size(batch):
    expected = N_HIGH + N_HIGH * N_PERM_HIGH + N_LOW + N_LOW * N_PERM_LOW + N_ANCHOR
    assert len(batch) == expected, f"expected {expected}, got {len(batch)}"


def test_arm_sizes(batch):
    counts = batch["arm"].value_counts().to_dict()
    assert counts.get("donor_high", 0) == N_HIGH
    assert counts.get("permuted_high", 0) == N_HIGH * N_PERM_HIGH
    assert counts.get("donor_low", 0) == N_LOW
    assert counts.get("permuted_low", 0) == N_LOW * N_PERM_LOW
    assert counts.get("anchor", 0) == N_ANCHOR


def test_no_duplicate_sequences(batch):
    dups = batch[batch["sequence"].duplicated()]
    assert dups.empty, f"{len(dups)} duplicate sequences"


def test_all_sequences_are_strings(batch):
    # pandas 3 lit une colonne de texte en StringDtype, pandas 2 en object :
    # c'est le contenu qui compte, pas le dtype que la version choisit.
    assert all(isinstance(s, str) and len(s) >= 80 for s in batch["sequence"])


def test_interleave_order(batch):
    """The batch should interleave arms so a partial run still compares."""
    first_arms = set(batch["arm"].iloc[:5])
    assert len(first_arms) >= 3, f"first 5 rows only cover {first_arms}"


# ---------------------------------------------------------------------------
# CDR3 permutation invariants
# ---------------------------------------------------------------------------
def test_permutations_share_cdr3_composition(batch):
    """Each permuted CDR3 must have the same amino acid composition as its donor's CDR3."""
    perm_rows = batch[batch["arm"].isin(["permuted_high", "permuted_low"])]
    for _, row in perm_rows.iterrows():
        orig = str(row["cdr3_original"])
        perm = str(row["cdr3_permuted"])
        if not orig or not perm or orig == "nan" or perm == "nan":
            continue
        assert sorted(orig) == sorted(perm), (
            f"composition mismatch for pair {row['pair_id']}: "
            f"{sorted(orig)} vs {sorted(perm)}")


def test_cdr3_head_tail_preserved(batch):
    """First 3 and last 3 residues of CDR3 must match between original and permutation."""
    perm_rows = batch[batch["arm"].isin(["permuted_high", "permuted_low"])]
    for _, row in perm_rows.iterrows():
        orig = str(row["cdr3_original"])
        perm = str(row["cdr3_permuted"])
        if not orig or not perm or orig == "nan" or perm == "nan":
            continue
        assert orig[:3] == perm[:3], (
            f"head mismatch for pair {row['pair_id']}: {orig[:3]} vs {perm[:3]}")
        assert orig[-3:] == perm[-3:], (
            f"tail mismatch for pair {row['pair_id']}: {orig[-3:]} vs {perm[-3:]}")


def test_permutations_differ_from_original(batch):
    """Each permuted CDR3 must differ from the original."""
    perm_rows = batch[batch["arm"].isin(["permuted_high", "permuted_low"])]
    for _, row in perm_rows.iterrows():
        orig = str(row["cdr3_original"])
        perm = str(row["cdr3_permuted"])
        if not orig or not perm or orig == "nan" or perm == "nan":
            continue
        assert orig != perm, f"permutation == original for pair {row['pair_id']}: {orig}"


def test_high_donors_have_3_permutations(batch):
    """Each high donor must have exactly N_PERM_HIGH permutations."""
    high_donors = batch[batch["arm"] == "donor_high"]
    high_perms = batch[batch["arm"] == "permuted_high"]
    for _, donor in high_donors.iterrows():
        pid = donor["pair_id"]
        n = len(high_perms[high_perms["pair_id"] == pid])
        assert n == N_PERM_HIGH, f"donor pair_id {pid} has {n} permutations, expected {N_PERM_HIGH}"


def test_low_donors_have_2_permutations(batch):
    """Each low donor must have exactly N_PERM_LOW permutations."""
    low_donors = batch[batch["arm"] == "donor_low"]
    low_perms = batch[batch["arm"] == "permuted_low"]
    for _, donor in low_donors.iterrows():
        pid = donor["pair_id"]
        n = len(low_perms[low_perms["pair_id"] == pid])
        assert n == N_PERM_LOW, f"donor pair_id {pid} has {n} permutations, expected {N_PERM_LOW}"


# ---------------------------------------------------------------------------
# score-range tests
# ---------------------------------------------------------------------------
def test_high_donors_above_070(batch):
    """High donors must have donor_iiptm > 0.70."""
    high = batch[batch["arm"] == "donor_high"]
    below = high[high["donor_iiptm"] <= 0.70]
    assert below.empty, f"{len(below)} high donors with iiptm ≤ 0.70: {below['donor_iiptm'].tolist()}"


def test_low_donors_below_045(batch):
    """Low donors must have donor_iiptm < 0.45."""
    low = batch[batch["arm"] == "donor_low"]
    above = low[low["donor_iiptm"] >= 0.45]
    assert above.empty, f"{len(above)} low donors with iiptm ≥ 0.45: {above['donor_iiptm'].tolist()}"


# ---------------------------------------------------------------------------
# diversity tests
# ---------------------------------------------------------------------------
def test_jaccard_below_ceiling(batch):
    """All permutation Jaccard 5-mer values must be below the ceiling."""
    perm_rows = batch[batch["arm"].isin(["permuted_high", "permuted_low"])]
    jaccards = perm_rows["jaccard_5mer"].dropna()
    above = jaccards[jaccards >= JACCARD_CEILING]
    assert above.empty, f"{len(above)} permutations with Jaccard ≥ {JACCARD_CEILING}: {above.tolist()}"


# ---------------------------------------------------------------------------
# anchors
# ---------------------------------------------------------------------------
def _score_files():
    """Tous les CSV de scores, dans results/ comme a la racine."""
    seen: dict[str, Path] = {}
    for d in (REPO_ROOT / "results", REPO_ROOT):
        for p in sorted(d.glob(SCORE_GLOB)):
            seen.setdefault(p.name, p)
    return list(seen.values())


def test_anchors_are_real_repeats(batch):
    """Every anchor must have been scored in a previous run."""
    anchor_seqs = set(batch[batch["arm"] == "anchor"]["sequence"])
    scored_before: set[str] = set()
    for path in _score_files():
        try:
            df = pd.read_csv(path, usecols=["sequence"])
        except (ValueError, KeyError):
            continue
        scored_before.update(df["sequence"])
    missing = anchor_seqs - scored_before
    assert not missing, f"{len(missing)} anchors not in any previous score CSV"


def test_anchor_sources_are_all_on_the_live_target():
    """Le vivier d'ancres ne doit contenir que des runs de la cible vivante.

    Les quatre ancres du lot `scramble` venaient de l'ancienne cible : elles
    n'ont mesure aucun bruit entre runs, parce que select_anchors lisait une
    liste de fichiers ecrite en dur (FINDINGS 90).
    """
    names = runs_on_current_target()
    assert names, "aucun run de la cible vivante : la table de provenance est vide"
    assert_current_target([score_path(n) for n in names], "le vivier d'ancres")


def test_every_score_file_has_an_established_target():
    """Aucun fichier de scores ne doit avoir de cible inconnue.

    C'est la condition pour qu'un constructeur de lot puisse refuser de puiser
    dans la mauvaise : un fichier ajoute sans cible etablie doit casser ici,
    pas apres 110 min de GPU.
    """
    unknown = [p.name for p in _score_files() if target_of(p) is None]
    assert not unknown, (f"cible non etablie pour {unknown} — "
                         f"l'inscrire dans nova_lab/provenance.RUN_TARGET")


# ---------------------------------------------------------------------------
# unit tests: kmers / jaccard
# ---------------------------------------------------------------------------
def test_kmers_basic():
    assert kmers("ABCDE", 3) == {"ABC", "BCD", "CDE"}
    assert kmers("AB", 5) == set()
    assert len(kmers("ABCDEFGH", 5)) == 4


def test_jaccard_identical():
    assert jaccard_5mer("ABCDEFGH", "ABCDEFGH") == 1.0


def test_jaccard_different():
    j = jaccard_5mer("ABCDEFGHIJ", "KLMNOPQRST")
    assert j == 0.0


# ---------------------------------------------------------------------------
# unit tests: permute_cdr3
# ---------------------------------------------------------------------------
def test_permute_cdr3_basic():
    """permute_cdr3 on a known sequence produces valid permutations."""
    # Use a real-ish nanobody-like sequence; we just need regions() to work.
    # We'll use a donor from the batch if available, else skip.
    if not BATCH_PATH.exists():
        pytest.skip("batch not generated")
    batch = pd.read_csv(BATCH_PATH)
    donor = batch[batch["arm"] == "donor_high"].iloc[0]
    seq = donor["sequence"]
    perms = permute_cdr3(seq, 2, seed=99)
    assert len(perms) == 2
    assert all(p != seq for p in perms)
    assert len(set(perms)) == 2  # all distinct


# ---------------------------------------------------------------------------
# verdict tests with synthetic data
# ---------------------------------------------------------------------------
def test_verdict_strategy_works():
    """When reading 1 and 3 pass, verdict says strategy works."""
    result = {"reading1_pass": True, "reading3_pass": True,
              "reading2_rho": 0.1, "reading4_pass": True}
    v = verdict(result)
    assert "STRATEGY WORKS" in v


def test_verdict_strategy_dead():
    """When reading 1 fails, verdict says strategy dead."""
    result = {"reading1_pass": False, "reading3_pass": True,
              "reading2_rho": -0.1, "reading4_pass": True}
    v = verdict(result)
    assert "STRATEGY DEAD" in v


def test_verdict_strategy_marginal():
    """When reading 1 passes but 3 fails, verdict says marginal."""
    result = {"reading1_pass": True, "reading3_pass": False,
              "reading2_rho": 0.1, "reading4_pass": True}
    v = verdict(result)
    assert "MARGINAL" in v


def test_verdict_inconclusive():
    """When reading 1 is None, verdict says inconclusive."""
    result = {"reading1_pass": None, "reading3_pass": None,
              "reading2_rho": None, "reading4_pass": None}
    v = verdict(result)
    assert "INCONCLUSIVE" in v


def test_verdict_displacement_warning():
    """When reading 2 rho > threshold, verdict includes warning."""
    result = {"reading1_pass": True, "reading3_pass": True,
              "reading2_rho": 0.5, "reading4_pass": True}
    v = verdict(result)
    assert "WARNING" in v
    assert "grows with score" in v


def test_verdict_diversity_fail():
    """When reading 4 fails, verdict includes diversity failure."""
    result = {"reading1_pass": True, "reading3_pass": True,
              "reading2_rho": 0.1, "reading4_pass": False}
    v = verdict(result)
    assert "DIVERSITY FAILED" in v


# ---------------------------------------------------------------------------
# analyse round-trip with synthetic data
# ---------------------------------------------------------------------------
def test_analyse_round_trip(tmp_path):
    """analyse() recovers a planted signal from synthetic data."""
    from scripts.build_scramble_batch import analyse

    # Build a synthetic batch
    batch_rows = []
    # 2 high donors with 3 permutations each
    for i in range(2):
        batch_rows.append({"sequence": f"DONOR_HIGH_{i}", "arm": "donor_high",
                           "pair_id": i, "donor_iiptm": 0.75 + i * 0.02,
                           "jaccard_5mer": 1.0, "label": f"dh{i}", "group": "c"})
        for j in range(3):
            batch_rows.append({"sequence": f"PERM_HIGH_{i}_{j}", "arm": "permuted_high",
                               "pair_id": i, "donor_iiptm": 0.75 + i * 0.02,
                               "jaccard_5mer": 0.85, "label": f"ph{i}{j}", "group": "c"})

    # Scores: permutations score close to their donor (displacement < 0.027)
    score_rows = []
    for i in range(2):
        score_rows.append({"sequence": f"DONOR_HIGH_{i}",
                           "design_iiptm": 0.75 + i * 0.02})
        for j in range(3):
            score_rows.append({"sequence": f"PERM_HIGH_{i}_{j}",
                               "design_iiptm": 0.75 + i * 0.02 + (j - 1) * 0.005})

    batch_df = pd.DataFrame(batch_rows)
    scores_df = pd.DataFrame(score_rows)

    batch_path = tmp_path / "batch.csv"
    scores_path = tmp_path / "scores.csv"
    batch_df.to_csv(batch_path, index=False)
    scores_df.to_csv(scores_path, index=False)

    result = analyse(str(scores_path), str(batch_path), log=lambda *a: None)

    # Both donors should have displacement < 0.027
    assert result["reading1_n_pass"] == 2
    # Median of high permutations should be > 0.70
    assert result["reading3_pass"] == True
