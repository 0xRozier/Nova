"""Sequence features for predicting a nanobody's BoltzGen metrics.

Region-aware: the framework is nearly constant across submissions, so most of
the signal has to come from the CDRs (and CDR3 above all). Everything here is
computed from the sequence and its IMGT regions, on CPU.
"""
import numpy as np

from .nanobodies import germline_identity, regions

AAS = "ACDEFGHIKLMNPQRSTVWY"
# Kyte-Doolittle hydropathy, net charge at pH 7, and volume class
HYDROPATHY = dict(zip(AAS, [1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8,
                            1.9, -3.5, -1.6, -3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3]))
CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1, "H": 0.1}
AROMATIC = set("FWY")
POLAR = set("STNQCY")
REGION_KEYS = ("CDR1", "CDR2", "CDR3", "FR2", "FR3")


def _region_features(seq: str) -> tuple[list[float], list[str]]:
    if not seq:
        return [0.0] * 8, []
    hydro = [HYDROPATHY[a] for a in seq]
    values = [
        len(seq),
        float(np.mean(hydro)),
        float(np.max(hydro)),
        sum(CHARGE.get(a, 0) for a in seq),
        sum(a in AROMATIC for a in seq) / len(seq),
        sum(a in POLAR for a in seq) / len(seq),
        sum(a == "G" for a in seq) / len(seq),
        sum(a == "P" for a in seq) / len(seq),
    ]
    return values, ["len", "hydro_mean", "hydro_max", "charge", "aromatic_frac",
                    "polar_frac", "gly_frac", "pro_frac"]


def feature_names() -> list[str]:
    names = []
    for key in REGION_KEYS:
        _, suffixes = _region_features("A")
        names += [f"{key}_{s}" for s in suffixes]
        names += [f"{key}_aa_{a}" for a in AAS]
    names += ["total_len", "total_charge", "total_hydro_mean", "cdr_total_len",
              "v_identity", "j_identity", "n_cys"]
    return names


def featurize(seq: str, region_dict: dict | None = None, germline: dict | None = None) -> np.ndarray | None:
    """Feature vector, or None when the sequence cannot be IMGT-numbered."""
    r = region_dict if region_dict is not None else regions(seq)
    if r is None:
        return None
    values = []
    for key in REGION_KEYS:
        part = r.get(key, "")
        region_values, _ = _region_features(part)
        values += region_values
        values += [part.count(a) / max(len(part), 1) for a in AAS]

    g = germline if germline is not None else germline_identity(seq)
    values += [
        len(seq),
        sum(CHARGE.get(a, 0) for a in seq),
        float(np.mean([HYDROPATHY[a] for a in seq])),
        sum(len(r.get(k, "")) for k in ("CDR1", "CDR2", "CDR3")),
        g["V"]["identity_pct"] if g else np.nan,
        g["J"]["identity_pct"] if g else np.nan,
        seq.count("C"),
    ]
    return np.array(values, dtype=np.float32)
