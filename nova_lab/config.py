"""Subnet parameters, copied from upstream config/config.yaml (UPSTREAM_COMMIT)."""

EPOCH_LENGTH = 361  # Tempo (360) + 1
NO_SUBMISSION_BLOCKS = 10

SMALL_MOLECULE_TARGET = ["O15427"]
NANOBODY_TARGET = ["Q9NZQ7"]
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
# Our reproduction of the input value is NOT calibrated: 62% of real submissions
# the validators accepted exceed 89 on the PSH we compute, so treat a breach as
# a warning rather than a verdict.
PSH_STRICT_REJECTION_THRESHOLD = 89
MIN_HUMAN_V_IDENTITY = 82
MIN_HUMAN_J_IDENTITY = 86
PSH_STRICT_REJECTION_THRESHOLD = 89
# NOVA-nanobody-filter metanano/config.py developability red flags
TOTAL_CDR_LENGTH_RANGE = (20, 39)
