"""Subnet parameters, copied from upstream config/config.yaml (UPSTREAM_COMMIT)."""

EPOCH_LENGTH = 361  # Tempo (360) + 1
NO_SUBMISSION_BLOCKS = 10

SMALL_MOLECULE_TARGET = ["O15427"]
# La cible nanobody a change le 20 novembre 2026 : Q9NZQ7 (PD-L1, structure
# 8AOK) laisse la place a P05231 (interleukine-6, structure 4O9H). Tout ce qui
# depend de la cible doit lire ces constantes et non des litteraux : la premiere
# bascule a montre que `8AOK` etait en dur dans six scripts.
NANOBODY_TARGET = ["P05231"]

# L'intervalle de la sequence UniProt que le validateur soumet au repliement.
NANOBODY_CLIP_INTERVAL = [[27, 212]]

# La structure de reference, sa chaine cible, l'intervalle de residus et le site
# de liaison -- copies de la configuration amont.
#
# Les numeros du site de liaison sont ceux du `label_seq_id` du mmCIF, c'est-a-dire
# l'index canonique sequentiel que BoltzGen lit, et NON la numerotation auteur du
# cristal, qui est inferieure de 2 sur 4O9H. Verifie en calculant les contacts :
# au decalage -2, les 8 residus du site sont tous a moins de 5 A du Fab lie dans
# la structure, contre 2 sur 8 si on les prend pour de la numerotation auteur.
# Ce controle est a refaire a chaque changement de cible (FINDINGS 66 et 67).
NANOBODY_STRUCTURE = "4O9H"
NANOBODY_STRUCTURE_CHAIN = "A"
NANOBODY_STRUCTURE_RES_INDEX = "21..186"
NANOBODY_BINDING_SITE = "24,77,80,82,131,184..186"

# La cible precedente, gardee parce que tout l'historique de mesure de ce depot
# porte sur elle : les sections 1 a 66 de FINDINGS, les lots de results/, le
# vivier de donneurs et le substitut.
PREVIOUS_NANOBODY = {
    "target": "Q9NZQ7",
    "structure": "8AOK",
    "chain": "A",
    "res_index": "1..115",
    "binding_site": "37,39,41..44,46,49,51,56,59,94..96,98,100,104..108",
}
NANOBODY_WEIGHT = 0.4

NUM_MOLECULES = 20
MIN_ENTROPY = 0.25
MAX_SIMILARITY_TO_HISTORICAL = 0.5
BOLTZ_METRICS = ("affinity_probability_binary", "affinity_pred_value")

FILTERS = {
    "max_heavy_atoms": 50,
    "min_heavy_atoms": 10,
    "max_molecular_weight": 700,
    "min_molecular_weight": 100,
    "max_hbd": 5,
    "max_hba": 10,
    "max_logp": 7.0,
    "min_rotatable_bonds": 1,
    "max_rotatable_bonds": 10,
    "min_rings": 1,
    "max_ring_size": 12,
    "max_spiro_atoms": 2,
    "max_aromatic_rings": 5,
    "max_aliphatic_rings": 5,
    "max_fsp3": 0.8,
    "max_triple_bonds": 1,
    "max_amide_bonds": 5,
    "max_stereocenters": 6,
    "max_halogens": 6,
    "max_iodine": 1,
    "max_bromine": 2,
    "max_fluorine": 9,
    "max_phosphorus": 2,
    "max_nitrogen": 9,
    "max_oxygen": 9,
    "max_sulfur": 2,
    "max_aromatic_nitrogen": 6,
    "max_aromatic_carbon": 28,
    "max_aromatic_fraction": 0.92,
    "min_heteroatom_fraction": 0.05,
    "max_heteroatom_fraction": 0.65,
    "max_bertz_ct": 1950,
    "max_kappa1": 35,
    "max_kappa2": 17,
    "max_kappa3": 12,
    "max_nitrile": 2,
    "max_hydroxyl": 5,
}

V_NEIGHBORHOOD_CATEGORIES = {
    "V1": ["IGHV3-23"],
    "V2": ["IGHV3-66"],
    "V3": ["IGHV3-74"],
    "V4": ["IGHV3-7", "IGHV3-11", "IGHV3-43", "IGHV3-48", "IGHV3-53", "IGHV3-64"],
}
CDR3_LENGTH_RANGE = range(5, 23 + 1)

# ---------------------------------------------------------------------------
# What changes at epoch 26331 (upstream PR 102, branch rxn-group-filtering)
# ---------------------------------------------------------------------------
# Announced on the subnet's Discord on 2026-11-13: PR 102 merges at epoch 26331
# and the branch's config values take effect with it. The announcement warns
# nanobody miners that "this change will also affect how valid germlines and
# CDR3 lengths are chosen each epoch due to RNG state consumption".
#
# READ THE CODE AND IT DOES NOT. Both the old and the new reaction selection are
# pure arithmetic on `seed` and never touch the `rng` object:
#
#   old:  allowed_option = seed % total_reactions
#   new:  for group in rxn_groups: allowed_rxn_index = seed % len(group)
#
# `rng` is only consumed afterwards, by the allowed_v and allowed_cdr3_length
# draws, in the same order as before. Checked on 200 000 block hashes: the
# (V gene, CDR3 length) pair comes out identical under both functions, 0
# disagreements (tests/test_challenge_params.py).
#
# So either the announcement is cautious rather than exact, or what gets merged
# differs from the branch head we read. We do not have to choose: a nanobody
# submission is only accepted if its V gene and CDR3 length match the epoch's,
# so the field itself tells us whether our derivation is right, every epoch,
# without trusting either source. That is `scripts.verify_challenge`, and the
# miner should not submit when it disagrees.
RXN_GROUP_FILTERING_FROM_EPOCH = 26331
# config.yaml on the branch; together these cover reactions 1-44.
RXN_GROUPS = (
    (2, 3, 4, 5, 6, 11, 13, 15, 20, 23, 27, 29, 33, 44),
    (1, 7, 8, 9, 12, 14, 19, 21, 22, 28, 41, 42),
    (10, 16, 17, 18, 24, 25, 26, 30, 31, 32, 34, 35, 36, 37, 38, 39, 40, 43),
)
# The molecule filters tighten sharply at the same time. We do not mine
# molecules, so FILTERS above is left on the values in force until 26331 and
# these are kept beside it rather than swapped in: tests/test_parity.py checks
# that archived molecules pass today's filters, and the archive was built under
# the looser ones, so switching would fail that test for the right reason.
FILTERS_FROM_26331 = {
    **FILTERS,
    "max_logp": 5.5,
    "max_ring_size": 7,
    "max_stereocenters": 5,
    "max_iodine": 0,
    "max_bromine": 0,
    "max_fluorine": 3,
    "max_phosphorus": 0,
    "max_sulfur": 1,
    "max_aromatic_nitrogen": 5,
    "max_aromatic_carbon": 24,
    "max_nitrile": 1,
}

# nanobody_requirements + germline_filtering in upstream config/config.yaml
ALLOWED_AAS = set("ACDEFGHIKLMNPQRSTVWY")
NUM_SEQUENCES = 1
MIN_SEQUENCE_LENGTH = 90
MAX_SEQUENCE_LENGTH = 150
MAX_HOMOPOLYMER_RUN = 6
MAX_DI_REPEAT_PAIRS = 4
SP_WINDOW = 12
SP_HYDRO_MIN_IN_WINDOW = 8
SP_SCAN_PREFIX = 30
# config.yaml: a rejection far stricter than TNP's own red threshold of 155.47,
# added upstream on 2026-10-20 and missing from this lab until production hit it.
# It IS reproducible: on the install the validator runs, the 36 submissions it
# accepted in epoch 26174 all read at or under 89 (median 83.5). See
# nova_lab.psh_calibration for the probability model and nova_lab.toolchain for
# what "the install it runs" means.
PSH_STRICT_REJECTION_THRESHOLD = 89
MIN_HUMAN_V_IDENTITY = 82
MIN_HUMAN_J_IDENTITY = 86
# NOVA-nanobody-filter metanano/config.py developability red flags
TOTAL_CDR_LENGTH_RANGE = (20, 39)
