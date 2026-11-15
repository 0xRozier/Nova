"""Ce que `scripts.repair_binder` promet, et que le reste du depot ne couvre pas.

Trois choses sont verifiables sans GPU, et ce sont les trois que le script
revendique :

1. **l'arithmetique** des deux murs. Le nombre de substitutions dans le CDR3
   n'est pas un reglage : il tombe de `cdr_similarity`, qui vaut
   1 - mismatches / longueur a longueur egale. Le test l'epingle contre la
   fonction du validateur elle-meme plutot que contre une formule recopiee ;
2. **ou** les substitutions sont autorisees a tomber. C'est tout l'ecart avec
   `scripts.humanise`, qui deplacait 12 residus de charpente en moyenne (§50) :
   ici l'apex du CDR3 (IMGT 107-114), CDR1, CDR2, FR1, FR2 et FR3 ne sont
   jamais touches, et chaque substitution va vers un residu germinal humain ;
3. **l'honnetete de la coupure**. L'archive ne porte pas de date, donc un
   donneur non date n'est admissible que si `nova_lab.archive_ledger` l'avait
   deja vu. C'est le trou que `tests/test_humanise_cutoff.py` epingle ouvert
   dans `scripts.humanise` ; ici il doit etre ferme.

Le refus est teste autant que l'acceptation : un donneur qui passe deja le
plancher J n'a plus rien de germinal a deplacer dans sa queue, et le script
doit le refuser au lieu d'inventer une mutation d'apex.

    python -m pytest tests/test_repair_binder.py -q
"""
import pandas as pd
import pytest

import scripts.repair_binder as mod
from nova_lab.nanobodies import cdr_similarity, germline_identity, regions

# un bon lieur de l'archive, longueur de CDR3 13, V 82,7 et J 78,6 : le J est a
# reparer, et c'est cette reparation qui finance la nouveaute
DONOR = ("DVMLAESGGGLVQAGGSLRLSCAASGFTFSSYAMSWVRQHPGKSLEWVSDVNTGGDSTGYADLVKGRFTI"
         "SRDNAKNTLYLQMNSLRPEDTAVYYCSLGPKYIWTHLGVHGQGTTVTVSS")
# un lieur deja conforme (J 92,9) : sa fenetre J est deja germinale
COMPLIANT = ("QQQLVESGGGLIQPGGSLRLSCSASGFTFSSYSMSWVRQPPGQQRDWVSEISSDGEPTDYADSVKGRFT"
             "ISRDSAKNTLYLQMNSLRAEDTAVYYCARTYYRSWHYRDSSNFDVWGQGTTVTVSS")


def test_the_cdr3_quota_is_the_validators_own_arithmetic():
    """`cdr3_mutations_needed` doit etre le plus petit m qui passe la porte."""
    for length in range(6, 25):
        for threshold in (0.85, 0.80):
            m = mod.cdr3_mutations_needed(length, threshold)
            motif = "A" * length
            mutated = "C" * m + "A" * (length - m)
            assert cdr_similarity({"CDR3": motif}, {"CDR3": mutated})["CDR3"] < threshold
            if m:
                almost = "C" * (m - 1) + "A" * (length - m + 1)
                assert cdr_similarity({"CDR3": motif},
                                      {"CDR3": almost})["CDR3"] >= threshold


def test_a_120mer_needs_seven_substitutions_for_the_whole_sequence_wall():
    """La borne de 0,95 de `is_duplicate`, lue sur une vraie sequence.

    Le donneur est lui-meme archive a 1,0 : la nouveaute commence par le faire
    descendre sous 0,95 contre lui-meme.
    """
    from nova_lab.nanobodies import whole_identity
    assert len(DONOR) == 120
    for k in range(4, 7):
        mutant = "".join("W" if i < k else a for i, a in enumerate(DONOR))
        assert whole_identity(DONOR, mutant) >= 0.95, k
    mutant = "".join("W" if i < 7 else a for i, a in enumerate(DONOR))
    assert whole_identity(DONOR, mutant) < 0.95


def test_every_substitution_is_germline_directed_and_outside_the_paratope():
    pos_map = mod.imgt_map(DONOR)
    assert pos_map is not None
    seq, info = mod.repair(DONOR, "IGHV3-23", cdr3_need=3, budget=8, pos_map=pos_map)
    assert seq is not None, info
    assert len(seq) == len(DONOR)
    allowed = set(mod.STALK) | set(mod.TAIL) | set(range(118, 129))
    moved = [p for p in pos_map if seq[pos_map[p]] != DONOR[pos_map[p]]]
    assert moved, "nothing moved"
    assert set(moved) <= allowed, sorted(set(moved) - allowed)
    assert info["n_apex"] == 0 and info["n_framework"] == 0
    assert info["n_mut"] >= 8 and info["n_cdr3"] >= 3
    # germline-directed: the new residue is the one the cited allele carries
    jref = mod._positions(mod.JREF[info["j_allele"]])
    vref = mod._positions(mod.VREF[info["v_allele"]])
    for p in moved:
        want = jref.get(p) if p >= 115 else vref.get(p)
        assert seq[pos_map[p]] == want, (p, seq[pos_map[p]], want)
    # the CDRs that bind are untouched
    before, after = regions(DONOR), regions(seq)
    assert (before["CDR1"], before["CDR2"]) == (after["CDR1"], after["CDR2"])
    assert before["CDR3"][2:-3] == after["CDR3"][2:-3]


def test_the_repair_lifts_the_germline_floors_it_was_built_for():
    from nova_lab.config import MIN_HUMAN_J_IDENTITY, MIN_HUMAN_V_IDENTITY
    seq, info = mod.repair(DONOR, "IGHV3-23", cdr3_need=3, budget=8)
    assert seq is not None, info
    before, after = germline_identity(DONOR), germline_identity(seq)
    assert before["J"]["identity_pct"] < MIN_HUMAN_J_IDENTITY
    assert after["J"]["identity_pct"] >= MIN_HUMAN_J_IDENTITY
    assert after["V"]["identity_pct"] >= MIN_HUMAN_V_IDENTITY
    assert after["V"]["identity_pct"] >= before["V"]["identity_pct"]


def test_an_already_compliant_donor_is_refused_not_mutated_in_the_apex():
    """Sa fenetre J est germinale : il ne reste rien a deplacer hors de l'apex."""
    seq, why = mod.repair(COMPLIANT, "IGHV3-74", cdr3_need=3, budget=8)
    assert seq is None
    assert "CDR3" in why or "stalk" in why


@pytest.fixture
def pool(tmp_path, monkeypatch):
    """Un vivier dont un seul donneur non date est situable dans le temps."""
    monkeypatch.setattr(mod, "DATA_DIR", tmp_path)
    import scripts.propose_nanobody as propose_mod
    monkeypatch.setattr(propose_mod, "DATA_DIR", tmp_path)
    import nova_lab.archive_ledger as ledger_mod
    monkeypatch.setattr(ledger_mod, "DATA_DIR", tmp_path)
    (tmp_path / "dataset").mkdir(parents=True)
    (tmp_path / "prebuilt").mkdir(parents=True)
    pd.DataFrame({"sequence": ["OLD_ARCHIVE"]}).to_csv(
        tmp_path / "prebuilt" / "T_regions.csv", index=False)
    pd.DataFrame({
        "sequence": ["OLD_ARCHIVE", "NEW_ARCHIVE", "SUBMITTED_BEFORE", "SUBMITTED_AFTER"],
        "epoch": [None, None, 26200.0, 26400.0],
        "design_iiptm": [0.73, 0.78, 0.74, 0.80],
        "cdr3": ["A" * 13] * 4,
        "v_genes": ["IGHV3-23"] * 4,
    }).to_parquet(tmp_path / "dataset" / "T_nanobodies.parquet")
    return tmp_path


def test_the_cutoff_refuses_an_undated_sequence_it_cannot_place_in_time(pool):
    kept = set(mod.donors("T", 13, 0.72, "IGHV3-23", cutoff=26300)["sequence"])
    assert kept == {"OLD_ARCHIVE", "SUBMITTED_BEFORE"}
    assert "NEW_ARCHIVE" not in kept     # le trou de scripts.humanise
    assert "SUBMITTED_AFTER" not in kept


def test_without_a_cutoff_the_whole_pool_is_visible(pool):
    kept = set(mod.donors("T", 13, 0.72, "IGHV3-23")["sequence"])
    assert len(kept) == 4
