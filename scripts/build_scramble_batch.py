"""Build the CDR3 permutation batch — does shuffling the loop interior preserve score?

HYPOTHESIS: permuting the CDR3 interior (keeping first 3 + last 3 residues,
shuffling the middle) preserves iiptm even at high score (0.70–0.79). §77
measured only 0.004 displacement on a low-scoring donor (iiptm 0.21), but
§87 noted this was not measured where it matters.

DESIGN: 56 sequences in five arms:
  donor_high (10)    — high-scoring field sequences (0.70–0.79 iiptm)
  permuted_high (30) — 3 CDR3 interior permutations per high donor
  donor_low (4)      — low-scoring archive sequences (0.20–0.40 iiptm)
  permuted_low (8)   — 2 CDR3 interior permutations per low donor
  anchor (4)         — repeats from previous runs for between-run noise

PRE-REGISTERED READINGS (written before the run):
  1. Score preservation at high iiptm. For each of the 10 high donors, compute
     |median(permutation_iiptm) − original_iiptm|. Pass: at least 7/10 donors
     have displacement < 2σ (0.027).
  2. Score-dependence of displacement. ρ(original_iiptm,
     mean_|permutation_iiptm − original_iiptm|) across all 14 donors. If
     ρ > +0.3, displacement grows with score and the §77 null was specific
     to the low range.
  3. Operational threshold. Median iiptm of the 30 high permutations must
     exceed 0.70.
  4. Diversity confirmation. Whole-sequence Jaccard 5-mer between each
     permutation and its original. Expected: ~0.854. All must be < 0.90.

Usage:
    python -m scripts.build_scramble_batch --out inputs/kaggle_scramble.csv
    python -m scripts.build_scramble_batch --analyse scores.csv --batch inputs/kaggle_scramble.csv
"""
import argparse

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from nova_lab import REPO_ROOT
from nova_lab.nanobodies import _regions_of, check_sequence, regions

SCORE_GLOB = "scores_nanobody_*.csv"


def _check_sequence_safe(seq: str) -> str:
    """Run check_sequence, returning 'ok' or the problems. Tolerates missing hmmscan."""
    try:
        warnings: list[str] = []
        problems = check_sequence(seq, warnings=warnings)
        return "; ".join(problems) if problems else "ok"
    except (FileNotFoundError, OSError):
        return "check unavailable (hmmscan missing)"


def _regions_safe(seq: str) -> dict | None:
    """Extract IMGT regions, using anarcii (Python-only) when hmmscan is unavailable."""
    r = regions(seq)
    if r is not None:
        return r
    # Fallback: use anarcii via Chain.batch(use_anarcii=True)
    from abnumber import Chain
    try:
        chains, _ = Chain.batch({"0": seq}, scheme="imgt", cdr_definition="imgt",
                                use_anarcii=True)
        c = chains.get("0")
        return _regions_of(c) if c is not None else None
    except Exception:
        return None


def _regions_batch_safe(sequences: list[str]) -> list[dict | None]:
    """Batch-extract IMGT regions, using anarcii when hmmscan is unavailable."""
    from abnumber import Chain
    try:
        chains, _ = Chain.batch({str(i): s for i, s in enumerate(sequences)},
                                scheme="imgt", cdr_definition="imgt", use_anarcii=True)
    except Exception:
        return [_regions_safe(s) for s in sequences]
    out = []
    for i in range(len(sequences)):
        c = chains.get(str(i))
        out.append(_regions_of(c) if c is not None else None)
    return out

# Live-target runs with high-scoring field sequences
LIVE_TARGET_RUNS = [
    "scores_nanobody_arms_26299_26306.csv",
    "scores_nanobody_deep_26304_26345_26347.csv",
    "scores_nanobody_design3arms_26345.csv",
    "scores_nanobody_k1sample_26345.csv",
    "scores_nanobody_wide4arms_26345.csv",
]

N_HIGH = 10
N_PERM_HIGH = 3
N_LOW = 4
N_PERM_LOW = 2
N_ANCHOR = 4

# Pre-registered thresholds
SIGMA_2 = 0.027            # 2× between-run σ (0.0133)
PRESERVATION_PASS = 7      # at least 7/10 donors pass displacement < 2σ
RHO_THRESHOLD = 0.3        # score-dependence of displacement
OPERATIONAL_THRESHOLD = 0.70  # median of high permutations must exceed
JACCARD_CEILING = 0.90     # diversity: all must be below


# ---------------------------------------------------------------------------
# k-mer Jaccard
# ---------------------------------------------------------------------------
def kmers(seq: str, k: int = 5) -> set[str]:
    """Return the set of k-mers in a sequence."""
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def jaccard_5mer(a: str, b: str) -> float:
    """Jaccard similarity on the 5-mer sets of two sequences."""
    sa, sb = kmers(a, 5), kmers(b, 5)
    if not sa and not sb:
        return 1.0
    return len(sa & sb) / len(sa | sb)


# ---------------------------------------------------------------------------
# load live-target scores
# ---------------------------------------------------------------------------
def load_live_target_scores(log=print) -> pd.DataFrame:
    """Load all scores from live-target runs, deduplicate by sequence keeping best iiptm."""
    frames = []
    for name in LIVE_TARGET_RUNS:
        path = REPO_ROOT / name
        if not path.exists():
            log(f"  WARNING: {name} not found, skipping")
            continue
        df = pd.read_csv(path)
        df["source_file"] = name
        frames.append(df)
    if not frames:
        raise RuntimeError("no live-target score files found")
    all_df = pd.concat(frames, ignore_index=True)
    # Deduplicate: keep the row with the best iiptm per sequence
    dedup = (all_df.sort_values("design_iiptm", ascending=False)
             .drop_duplicates(subset=["sequence"], keep="first")
             .reset_index(drop=True))
    log(f"loaded {len(all_df)} rows from {len(frames)} files, "
        f"{len(dedup)} unique sequences")
    return dedup


# ---------------------------------------------------------------------------
# CDR3 permutation
# ---------------------------------------------------------------------------
def permute_cdr3(sequence: str, n: int, seed: int, log=print) -> list[str]:
    """Generate n distinct CDR3 interior permutations of a sequence.

    Keeps cdr3[:3] and cdr3[-3:] fixed, shuffles the interior.
    Returns full sequences with the permuted CDR3 swapped in.
    """
    r = _regions_safe(sequence)
    if r is None:
        raise ValueError(f"cannot parse regions for sequence: {sequence[:30]}...")
    cdr3 = r["CDR3"]
    if len(cdr3) < 9:
        raise ValueError(f"CDR3 too short for permutation: {cdr3} (len {len(cdr3)})")

    head = cdr3[:3]
    middle = list(cdr3[3:-3])
    tail = cdr3[-3:]

    rng = np.random.default_rng(seed)
    seen = {cdr3}  # reject the original
    permutations = []

    for _ in range(n * 200):
        if len(permutations) >= n:
            break
        rng.shuffle(middle)
        new_cdr3 = head + "".join(middle) + tail
        if new_cdr3 in seen:
            continue
        seen.add(new_cdr3)
        # Swap CDR3 in the full sequence
        new_seq = sequence.replace(cdr3, new_cdr3, 1)
        if new_seq == sequence:
            raise ValueError(f"CDR3 {cdr3} not found in sequence")
        permutations.append(new_seq)

    if len(permutations) < n:
        raise ValueError(f"only generated {len(permutations)}/{n} distinct permutations "
                         f"for CDR3 {cdr3} (interior len {len(middle)})")
    return permutations


# ---------------------------------------------------------------------------
# donor selection
# ---------------------------------------------------------------------------
def select_donors(scores: pd.DataFrame, n_high: int = N_HIGH, n_low: int = N_LOW,
                  log=print) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Select diverse high and low donors.

    High donors: from sequences with iiptm > 0.70, pick n_high with diverse
    CDR3s spanning different lengths and the full 0.70–0.79 score range.
    CDR3 must be ≥ 9 residues (need ≥ 3 interior positions).

    Low donors: from scores_nanobody_ptm.csv, pick n_low with iiptm 0.20–0.40.
    """
    # --- High donors ---
    high_pool = scores[scores["design_iiptm"] > 0.70].copy()
    log(f"high pool: {len(high_pool)} sequences with iiptm > 0.70")

    # Compute CDR3 for each (batch for speed)
    region_list = _regions_batch_safe(list(high_pool["sequence"]))
    cdr3s = [r["CDR3"] if r else "" for r in region_list]
    high_pool["cdr3_computed"] = cdr3s
    high_pool["cdr3_len"] = high_pool["cdr3_computed"].str.len()

    # Filter: CDR3 ≥ 9 (need ≥ 3 interior positions)
    high_pool = high_pool[high_pool["cdr3_len"] >= 9].copy()
    log(f"  after CDR3 ≥ 9 filter: {len(high_pool)}")

    # Deduplicate by CDR3 (keep best iiptm per unique CDR3)
    high_pool = (high_pool.sort_values("design_iiptm", ascending=False)
                 .drop_duplicates(subset=["cdr3_computed"], keep="first"))
    log(f"  after CDR3 dedup: {len(high_pool)} unique CDR3s")

    if len(high_pool) < n_high:
        raise RuntimeError(f"only {len(high_pool)} diverse high donors, need {n_high}")

    # Select: cover the score range and CDR3 length diversity
    # Sort by score, pick evenly spaced + ensure CDR3 length diversity
    high_pool = high_pool.sort_values("design_iiptm", ascending=False).reset_index(drop=True)
    if len(high_pool) == n_high:
        high_donors = high_pool
    else:
        # Greedy: pick to maximize CDR3 length diversity across score range
        picked = []
        lengths_seen: dict[int, int] = {}
        # First pass: stratified by score (top, bottom, then fill)
        indices = np.linspace(0, len(high_pool) - 1, n_high * 2, dtype=int)
        for idx in indices:
            if len(picked) >= n_high:
                break
            row = high_pool.iloc[idx]
            cdr3_len = row["cdr3_len"]
            # Prefer new CDR3 lengths, but allow repeats if needed
            if cdr3_len not in lengths_seen or lengths_seen[cdr3_len] < 2:
                if idx not in picked:
                    picked.append(idx)
                    lengths_seen[cdr3_len] = lengths_seen.get(cdr3_len, 0) + 1
        # Fill remaining if needed
        for idx in range(len(high_pool)):
            if len(picked) >= n_high:
                break
            if idx not in picked:
                picked.append(idx)
        high_donors = high_pool.iloc[sorted(picked[:n_high])]

    log(f"  selected {len(high_donors)} high donors: "
        f"iiptm {high_donors['design_iiptm'].min():.3f}-{high_donors['design_iiptm'].max():.3f}, "
        f"CDR3 lengths {sorted(high_donors['cdr3_len'].unique())}")

    # --- Low donors ---
    ptm_path = REPO_ROOT / "scores_nanobody_ptm.csv"
    if not ptm_path.exists():
        raise RuntimeError(f"ptm scores not found: {ptm_path}")
    ptm_scores = pd.read_csv(ptm_path)
    low_pool = ptm_scores[(ptm_scores["design_iiptm"] >= 0.20) &
                          (ptm_scores["design_iiptm"] <= 0.40)].copy()
    log(f"low pool: {len(low_pool)} sequences with iiptm 0.20-0.40")

    # Compute CDR3 (batch for speed)
    region_list_low = _regions_batch_safe(list(low_pool["sequence"]))
    cdr3s = [r["CDR3"] if r else "" for r in region_list_low]
    low_pool["cdr3_computed"] = cdr3s
    low_pool["cdr3_len"] = low_pool["cdr3_computed"].str.len()

    # Filter: CDR3 ≥ 9
    low_pool = low_pool[low_pool["cdr3_len"] >= 9].copy()
    log(f"  after CDR3 ≥ 9 filter: {len(low_pool)}")

    # Pick n_low with diverse CDR3s
    low_pool = (low_pool.sort_values("design_iiptm", ascending=False)
                .drop_duplicates(subset=["cdr3_computed"], keep="first")
                .head(n_low))
    if len(low_pool) < n_low:
        raise RuntimeError(f"only {len(low_pool)} diverse low donors, need {n_low}")
    log(f"  selected {len(low_pool)} low donors: "
        f"iiptm {low_pool['design_iiptm'].min():.3f}-{low_pool['design_iiptm'].max():.3f}")

    return high_donors, low_pool


# ---------------------------------------------------------------------------
# anchor selection
# ---------------------------------------------------------------------------
def select_anchors(batch_sequences: set[str], n: int = N_ANCHOR,
                   log=print) -> pd.DataFrame:
    """Pick sequences scored in 2+ previous runs, not already in the batch."""
    seq_runs: dict[str, list[str]] = {}
    seq_iiptm: dict[str, list[float]] = {}
    for path in sorted(REPO_ROOT.glob(SCORE_GLOB)):
        run = path.stem.replace("scores_nanobody_", "")
        try:
            df = pd.read_csv(path, usecols=["sequence", "design_iiptm"])
        except (ValueError, KeyError):
            continue
        for _, row in df.iterrows():
            if pd.notna(row["design_iiptm"]):
                seq_runs.setdefault(row["sequence"], []).append(run)
                seq_iiptm.setdefault(row["sequence"], []).append(row["design_iiptm"])

    repeats = {seq: runs for seq, runs in seq_runs.items()
               if len(runs) >= 2 and seq not in batch_sequences}
    if len(repeats) < n:
        raise RuntimeError(f"only {len(repeats)} repeat sequences not in batch, need {n}")

    candidates = pd.DataFrame([
        {"sequence": seq, "n_runs": len(runs),
         "mean_iiptm": np.mean(seq_iiptm[seq]),
         "previous_runs": ";".join(runs)}
        for seq, runs in repeats.items()
    ]).sort_values("mean_iiptm")

    # Stratified pick: evenly spaced across the score range
    step = max(1, len(candidates) // n)
    picked = candidates.iloc[::step].head(n)
    if len(picked) < n:
        picked = candidates.head(n)
    log(f"anchors: {len(picked)} sequences from {len(repeats)} available, "
        f"iiptm range {picked['mean_iiptm'].min():.3f}-{picked['mean_iiptm'].max():.3f}")
    return picked


# ---------------------------------------------------------------------------
# interleave
# ---------------------------------------------------------------------------
def interleave(batch: pd.DataFrame) -> pd.DataFrame:
    """Interleave by [depth, arm] so a partial run still compares arms."""
    arm_order = {"donor_high": 0, "permuted_high": 1, "donor_low": 2,
                 "permuted_low": 3, "anchor": 4}
    batch = batch.copy()
    batch["_arm_order"] = batch["arm"].map(arm_order)
    batch["_depth"] = batch.groupby("arm").cumcount()
    batch = batch.sort_values(["_depth", "_arm_order"]).reset_index(drop=True)
    return batch.drop(columns=["_arm_order", "_depth"])


# ---------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------
def build(out: str, log=print) -> pd.DataFrame:
    """Build the CDR3 permutation batch CSV."""
    # 1. Load live-target scores
    scores = load_live_target_scores(log=log)

    # 2. Select donors
    high_donors, low_donors = select_donors(scores, log=log)

    # 3. Build rows
    rows = []

    # --- donor_high + permuted_high ---
    # Pre-compute CDR3 from the donor selection (already done in select_donors)
    for i, (_, donor) in enumerate(high_donors.iterrows()):
        pair_id = i
        seq = donor["sequence"]
        iiptm = donor["design_iiptm"]
        cdr3 = donor["cdr3_computed"]
        ptm = donor.get("design_ptm", np.nan)

        layer1 = _check_sequence_safe(seq)

        # Donor row
        rows.append({
            "sequence": seq, "arm": "donor_high", "label": f"donor_high#{i:02d}",
            "group": "candidate", "pair_id": pair_id,
            "donor_iiptm": iiptm, "donor_ptm": float(ptm) if pd.notna(ptm) else np.nan,
            "cdr3_original": cdr3, "cdr3_permuted": "",
            "cdr3_length": len(cdr3), "jaccard_5mer": 1.0,
            "layer1": layer1,
        })

        # Permutations
        perm_seqs = permute_cdr3(seq, N_PERM_HIGH, seed=42 + i, log=log)
        perm_regions = _regions_batch_safe(perm_seqs)
        for j, (perm_seq, r_perm) in enumerate(zip(perm_seqs, perm_regions)):
            perm_cdr3 = r_perm["CDR3"] if r_perm else ""
            jac = jaccard_5mer(seq, perm_seq)

            rows.append({
                "sequence": perm_seq, "arm": "permuted_high",
                "label": f"permuted_high#{i:02d}.{j}",
                "group": "candidate", "pair_id": pair_id,
                "donor_iiptm": iiptm,
                "donor_ptm": float(ptm) if pd.notna(ptm) else np.nan,
                "cdr3_original": cdr3, "cdr3_permuted": perm_cdr3,
                "cdr3_length": len(cdr3), "jaccard_5mer": jac,
                "layer1": _check_sequence_safe(perm_seq),
            })

    # --- donor_low + permuted_low ---
    for i, (_, donor) in enumerate(low_donors.iterrows()):
        pair_id = N_HIGH + i
        seq = donor["sequence"]
        iiptm = donor["design_iiptm"]
        cdr3 = donor["cdr3_computed"]
        ptm = donor.get("design_ptm", np.nan)

        rows.append({
            "sequence": seq, "arm": "donor_low", "label": f"donor_low#{i:02d}",
            "group": "candidate", "pair_id": pair_id,
            "donor_iiptm": iiptm, "donor_ptm": float(ptm) if pd.notna(ptm) else np.nan,
            "cdr3_original": cdr3, "cdr3_permuted": "",
            "cdr3_length": len(cdr3), "jaccard_5mer": 1.0,
            "layer1": _check_sequence_safe(seq),
        })

        perm_seqs = permute_cdr3(seq, N_PERM_LOW, seed=1000 + i, log=log)
        perm_regions = _regions_batch_safe(perm_seqs)
        for j, (perm_seq, r_perm) in enumerate(zip(perm_seqs, perm_regions)):
            perm_cdr3 = r_perm["CDR3"] if r_perm else ""
            jac = jaccard_5mer(seq, perm_seq)

            rows.append({
                "sequence": perm_seq, "arm": "permuted_low",
                "label": f"permuted_low#{i:02d}.{j}",
                "group": "candidate", "pair_id": pair_id,
                "donor_iiptm": iiptm,
                "donor_ptm": float(ptm) if pd.notna(ptm) else np.nan,
                "cdr3_original": cdr3, "cdr3_permuted": perm_cdr3,
                "cdr3_length": len(cdr3), "jaccard_5mer": jac,
                "layer1": _check_sequence_safe(perm_seq),
            })

    # --- anchors ---
    batch_seqs = {r["sequence"] for r in rows}
    anchor_df = select_anchors(batch_seqs, log=log)
    anchor_regions = _regions_batch_safe(list(anchor_df["sequence"]))
    for (_, anchor), r in zip(anchor_df.iterrows(), anchor_regions):
        seq = anchor["sequence"]
        cdr3 = r["CDR3"] if r else ""

        rows.append({
            "sequence": seq, "arm": "anchor", "label": "anchor",
            "group": "candidate", "pair_id": np.nan,
            "donor_iiptm": np.nan, "donor_ptm": np.nan,
            "cdr3_original": cdr3, "cdr3_permuted": "",
            "cdr3_length": len(cdr3), "jaccard_5mer": np.nan,
            "layer1": _check_sequence_safe(seq),
        })

    # 4. Assemble and interleave
    batch = pd.DataFrame(rows)
    assert len(batch) == N_HIGH + N_HIGH * N_PERM_HIGH + N_LOW + N_LOW * N_PERM_LOW + N_ANCHOR, \
        f"expected 56 rows, got {len(batch)}"

    # Check no duplicate sequences
    dups = batch[batch["sequence"].duplicated(keep=False)]
    if not dups.empty:
        log(f"WARNING: {len(dups)} duplicate sequences; dropping later copies")
        batch = batch.drop_duplicates(subset=["sequence"], keep="first")

    batch = interleave(batch)

    # 5. Write
    out_cols = ["sequence", "arm", "label", "group", "pair_id", "donor_iiptm",
                "donor_ptm", "cdr3_original", "cdr3_permuted", "cdr3_length",
                "jaccard_5mer", "layer1"]
    path = REPO_ROOT / out
    path.parent.mkdir(parents=True, exist_ok=True)
    batch[out_cols].to_csv(path, index=False)
    log(f"\nwrote {len(batch)} sequences to {path}")
    for arm in ["donor_high", "permuted_high", "donor_low", "permuted_low", "anchor"]:
        log(f"  {arm:16s}: {(batch['arm'] == arm).sum()}")
    return batch


# ---------------------------------------------------------------------------
# analyse: pre-registered readings after the GPU run
# ---------------------------------------------------------------------------
def analyse(scores_path: str, batch_path: str, log=print) -> dict:
    """Compute the 4 pre-registered readings."""
    scores = pd.read_csv(scores_path)
    batch = pd.read_csv(batch_path)

    # Merge scores onto batch metadata
    score_cols = ["sequence"]
    for c in ("design_iiptm", "design_ptm"):
        if c in scores.columns:
            score_cols.append(c)
    df = batch.merge(scores[score_cols], on="sequence", how="left")
    scored = df.dropna(subset=["design_iiptm"])
    log(f"{len(scored)} of {len(df)} sequences scored")

    result: dict = {"n_scored": len(scored)}

    # --- Reading 1: Score preservation at high iiptm ---
    high_donors = scored[scored["arm"] == "donor_high"]
    high_perms = scored[scored["arm"] == "permuted_high"]
    result["n_high_donors"] = len(high_donors)
    result["n_high_perms"] = len(high_perms)

    if len(high_donors) >= 1 and len(high_perms) >= 1:
        n_pass = 0
        displacements_high = []
        for _, donor in high_donors.iterrows():
            pid = donor["pair_id"]
            perms = high_perms[high_perms["pair_id"] == pid]
            if perms.empty:
                continue
            median_perm = perms["design_iiptm"].median()
            displacement = abs(median_perm - donor["design_iiptm"])
            displacements_high.append({
                "pair_id": pid, "donor_iiptm": donor["design_iiptm"],
                "median_perm_iiptm": median_perm, "displacement": displacement,
                "pass": displacement < SIGMA_2,
            })
            if displacement < SIGMA_2:
                n_pass += 1

        n_tested = len(displacements_high)
        result["reading1_displacements"] = displacements_high
        result["reading1_n_pass"] = n_pass
        result["reading1_n_tested"] = n_tested
        result["reading1_pass"] = n_pass >= PRESERVATION_PASS
        log(f"\n1. SCORE PRESERVATION: {n_pass}/{n_tested} donors have displacement < {SIGMA_2} "
            f"({'PASS' if result['reading1_pass'] else 'FAIL'}: need {PRESERVATION_PASS}/{N_HIGH})")
        for d in displacements_high:
            log(f"   pair {d['pair_id']:2d}: donor {d['donor_iiptm']:.3f}, "
                f"median perm {d['median_perm_iiptm']:.3f}, "
                f"displacement {d['displacement']:.4f} {'✓' if d['pass'] else '✗'}")
    else:
        result["reading1_pass"] = None
        log("\n1. SCORE PRESERVATION: insufficient data")

    # --- Reading 2: Score-dependence of displacement ---
    low_donors = scored[scored["arm"] == "donor_low"]
    low_perms = scored[scored["arm"] == "permuted_low"]

    all_displacements = []
    for arm_donors, arm_perms in [(high_donors, high_perms), (low_donors, low_perms)]:
        for _, donor in arm_donors.iterrows():
            pid = donor["pair_id"]
            perms = arm_perms[arm_perms["pair_id"] == pid]
            if perms.empty:
                continue
            mean_disp = abs(perms["design_iiptm"] - donor["design_iiptm"]).mean()
            all_displacements.append({
                "original_iiptm": donor["design_iiptm"],
                "mean_displacement": mean_disp,
            })

    if len(all_displacements) >= 4:
        orig = [d["original_iiptm"] for d in all_displacements]
        disp = [d["mean_displacement"] for d in all_displacements]
        rho = spearmanr(orig, disp).statistic
        result["reading2_rho"] = float(rho)
        result["reading2_growing"] = rho > RHO_THRESHOLD
        log(f"\n2. SCORE-DEPENDENCE: ρ(original_iiptm, mean_displacement) = {rho:+.3f} "
            f"({'displacement grows with score' if rho > RHO_THRESHOLD else 'no growth'})")
    else:
        result["reading2_rho"] = None
        log("\n2. SCORE-DEPENDENCE: insufficient data")

    # --- Reading 3: Operational threshold ---
    if len(high_perms) >= 3:
        median_high_perm = high_perms["design_iiptm"].median()
        result["reading3_median"] = float(median_high_perm)
        result["reading3_pass"] = median_high_perm > OPERATIONAL_THRESHOLD
        log(f"\n3. OPERATIONAL THRESHOLD: median iiptm of {len(high_perms)} high permutations = "
            f"{median_high_perm:.4f} "
            f"({'PASS' if result['reading3_pass'] else 'FAIL'}: need > {OPERATIONAL_THRESHOLD})")
    else:
        result["reading3_pass"] = None
        log("\n3. OPERATIONAL THRESHOLD: insufficient data")

    # --- Reading 4: Diversity confirmation ---
    perm_rows = scored[scored["arm"].isin(["permuted_high", "permuted_low"])]
    if "jaccard_5mer" in perm_rows.columns:
        jaccards = perm_rows["jaccard_5mer"].dropna()
        if len(jaccards) > 0:
            result["reading4_max_jaccard"] = float(jaccards.max())
            result["reading4_median_jaccard"] = float(jaccards.median())
            result["reading4_pass"] = bool(jaccards.max() < JACCARD_CEILING)
            log(f"\n4. DIVERSITY: Jaccard 5-mer max={jaccards.max():.3f}, "
                f"median={jaccards.median():.3f} "
                f"({'PASS' if result['reading4_pass'] else 'FAIL'}: need all < {JACCARD_CEILING})")
        else:
            result["reading4_pass"] = None
            log("\n4. DIVERSITY: no Jaccard values in batch")
    else:
        result["reading4_pass"] = None
        log("\n4. DIVERSITY: no jaccard_5mer column")

    # --- Verdict ---
    log(f"\n{'='*60}")
    v = verdict(result)
    log(v)
    return result


def verdict(result: dict) -> str:
    """Return the pre-registered conclusion from analyse() readings."""
    r1 = result.get("reading1_pass")
    r2_rho = result.get("reading2_rho")
    r3 = result.get("reading3_pass")
    r4 = result.get("reading4_pass")

    parts = []

    # Reading 1 + 3 → strategy works
    if r1 is True and r3 is True:
        parts.append("STRATEGY WORKS: permuting the CDR3 interior of a high-scoring "
                      "field sequence produces a new sequence scoring in the field's range.")
    elif r1 is False:
        parts.append("STRATEGY DEAD: CDR3 order matters at high score "
                      "(most displacements > 2σ).")
    elif r1 is True and r3 is False:
        parts.append("STRATEGY MARGINAL: score preserved per donor, but median of "
                      "high permutations falls below 0.70.")
    else:
        parts.append("INCONCLUSIVE: insufficient data to decide.")

    # Reading 2 → degradation warning
    if r2_rho is not None and r2_rho > RHO_THRESHOLD:
        parts.append(f"WARNING: displacement grows with score (ρ = {r2_rho:+.3f}). "
                      "Strategy degrades at higher targets.")

    # Reading 4 → diversity
    if r4 is False:
        parts.append("DIVERSITY FAILED: some permutations too similar to their donor "
                      "(Jaccard 5-mer ≥ 0.90).")
    elif r4 is True:
        parts.append("Diversity confirmed: all permutations below Jaccard 0.90 threshold.")

    return "\n".join(parts)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", default="inputs/kaggle_scramble.csv",
                        help="output batch CSV (default: inputs/kaggle_scramble.csv)")
    parser.add_argument("--analyse", metavar="SCORES_CSV",
                        help="analyse a scored batch instead of building one")
    parser.add_argument("--batch", metavar="BATCH_CSV",
                        help="the batch CSV to compare against (for --analyse)")
    args = parser.parse_args()

    if args.analyse:
        if not args.batch:
            parser.error("--analyse requires --batch")
        analyse(args.analyse, args.batch)
    else:
        build(args.out)


if __name__ == "__main__":
    main()
