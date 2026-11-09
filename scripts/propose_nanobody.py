"""Propose a nanobody sequence for an epoch, on CPU.

The epoch fixes the V gene and the CDR3 length, so a seed rarely matches both.
The generator therefore:

  1. picks framework donors from the archive whose V gene call includes the
     epoch's gene, best rank_sum first;
  2. grafts CDR3s of exactly the required length, taken from the best-scoring
     archived sequences (any V gene), and mutates CDR1/CDR2/CDR3 — the
     framework is left alone so the V gene call survives;
  3. keeps only candidates the offline checker accepts (nova_lab.nanobodies);
  4. ranks them with the per-metric surrogate against a reference pool of
     recent archived sequences, reproducing the validator's rank sum.

The surrogate ranks within a lineage at Spearman ~0.7, so this orders
candidates usefully but does not replace BoltzGen: re-score the shortlist on a
GPU before trusting it. Structure-based developability stays unchecked.

Usage:
    python -m scripts.propose_nanobody                       # current epoch
    python -m scripts.propose_nanobody --epoch 26025 --candidates 4000
"""
import argparse
import random

import numpy as np
import pandas as pd

from nova_lab import DATA_DIR
from nova_lab.cdr3_generators import GENERATORS, J_TAIL, make_generator
from nova_lab.config import MIN_HUMAN_J_IDENTITY, MIN_HUMAN_V_IDENTITY, NANOBODY_TARGET
from nova_lab.nanobodies import (NanobodyArchive, cdr_similarity, check_sequence,
                                 germline_identity, regions)
from nova_lab.nanobody_features import featurize
from nova_lab.surrogate import Surrogate
from scripts.check_nanobody import resolve_epoch
from scripts.train_nanobody_surrogate import METRIC_MODES, rank_sum

AAS = "ACDEFGHIKLMNPQRSTVWY"
# Most of the archive predates the germline floors (6227 of 8578 undated
# sequences are below them) and those old sequences score better, so donors and
# the reference pool must both come from sequences that satisfy today's rules.
GERMLINE_MARGIN = 2.0  # headroom donors need, since mutating CDRs lowers V identity
REFERENCE_POOL = 200
# The validator calls a sequence a duplicate at CDR3 identity >= 0.85, so a
# grafted CDR3 has to be mutated away from every archived one of that length.
MAX_CDR3_IDENTITY = 0.85


def assemble(r: dict) -> str:
    return r["FR1"] + r["CDR1"] + r["FR2"] + r["CDR2"] + r["FR3"] + r["CDR3"] + r["FR4"]


def mutate_region(seq: str, rng: random.Random, n: int, keep_tail: int = 0) -> str:
    """n point substitutions, leaving the last `keep_tail` residues alone."""
    out = list(seq)
    positions = range(max(len(out) - keep_tail, 0))
    if not positions:
        return seq
    for pos in rng.sample(list(positions), min(n, len(positions))):
        out[pos] = rng.choice(AAS)
    return "".join(out)


def generate(donors: pd.DataFrame, cdr3_pool: list[str], n: int, rng: random.Random,
             generator=None) -> list[str]:
    """Framework from a donor, CDR3 of the required length, CDRs mutated lightly.

    CDR1 and CDR2 are scored against the V reference, so changing them eats into
    the germline identity margin: CDR3 carries most of the variation.

    When *generator* is given, CDR3 creation is delegated to it instead of the
    inline uniform-random mutate_region.
    """
    seen, out = set(), []
    attempts = 0
    while len(out) < n and attempts < 40 * n:
        attempts += 1
        donor = donors.iloc[rng.randrange(len(donors))]
        r = dict(donor["regions"])
        donor_cdr3 = rng.choice(cdr3_pool)
        cdr3_length = len(donor_cdr3)
        if generator is not None:
            r["CDR3"] = generator.generate_cdr3(donor_cdr3, cdr3_length, J_TAIL, rng)
        else:
            r["CDR3"] = mutate_region(donor_cdr3, rng, rng.randint(2, 4), keep_tail=J_TAIL)
        if rng.random() < 0.15:
            r["CDR1"] = mutate_region(r["CDR1"], rng, 1)
        if rng.random() < 0.1:
            r["CDR2"] = mutate_region(r["CDR2"], rng, 1)
        seq = assemble(r)
        if seq not in seen:
            seen.add(seq)
            out.append(seq)
    return out


def score(sequences: list[str], models: dict, reference: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    """Predicted metrics plus rank sums, with the reference pool's own rank sums
    computed in the same combined pool so the two are comparable."""
    rows, feats = [], []
    for seq in sequences:
        r = regions(seq)
        g = germline_identity(seq)
        x = featurize(seq, r, g) if r else None
        if x is None or not np.isfinite(x).all():
            continue
        rows.append({"sequence": seq, "cdr3": r["CDR3"],
                     "v_genes": ",".join(g["V"]["genes"]) if g else "",
                     "v_identity": g["V"]["identity_pct"] if g else np.nan})
        feats.append(x)
    if not rows:
        return pd.DataFrame(), np.array([])
    X = np.vstack(feats)
    pred = pd.DataFrame({m: model.predict(X) for m, model in models.items()})
    # Rank candidates and reference sequences together, the way a validator
    # ranks an epoch's submissions against each other.
    pool = pd.concat([reference[list(METRIC_MODES)], pred], ignore_index=True)
    all_ranks = rank_sum(pool).to_numpy()
    out = pd.DataFrame(rows)
    out["predicted_rank_sum"] = all_ranks[len(reference):]
    for m in METRIC_MODES:
        out[m] = pred[m].to_numpy()
    return out.sort_values("predicted_rank_sum"), all_ranks[:len(reference)]


def propose(params: dict, target: str = NANOBODY_TARGET[0], candidates: int = 3000,
            n_donors: int = 40, n_cdr3: int = 60, shortlist: int = 10, seed: int = 0,
            archive_cutoff: int | None = None, generator_name: str = "chimeric",
            skip_surrogate: bool = False, log=print) -> dict:
    """Build and rank candidates for one epoch. Pure CPU, no chain write.

    `params` is what scripts.check_nanobody.resolve_epoch returns: it must carry
    `allowed_v` and `allowed_cdr3_length`. Returns a dict with the shortlist and
    the diagnostics the observation log needs; `shortlist` is empty whenever the
    archive cannot supply the epoch's gene or CDR3 length.

    `archive_cutoff` is mandatory for an honest replay: without it the generator
    grafts from sequences submitted at the replayed epoch or later. Measured on
    the first 21-epoch sweep, that let 23 of 40 candidates carry a CDR3 at
    0.70-0.84 identity to a sequence that did not exist yet — just under the
    0.85 duplicate threshold — and it accounted for 4 of 5 apparent wins.
    """
    allowed_v, cdr3_len = params["allowed_v"], params["allowed_cdr3_length"]
    out: dict = {"epoch": params.get("epoch"), "allowed_v": allowed_v,
                 "allowed_cdr3_length": cdr3_len, "shortlist": pd.DataFrame(),
                 "generated": 0, "valid": 0, "reasons": {}, "status": "ok",
                 "reference_size": 0, "place": None, "best_rank_sum": None}
    log(f"epoch {params.get('epoch')}: V gene {allowed_v}, CDR3 length {cdr3_len}")

    ds = DATA_DIR / "dataset"
    df = pd.read_parquet(ds / f"{target}_nanobodies.parquet")
    if archive_cutoff is not None:
        before = len(df)
        df = df[df["epoch"].isna() | (df["epoch"] < archive_cutoff)].copy()
        log(f"archive restricted to before epoch {archive_cutoff}: {len(df)} of {before} sequences")
    df["true_rank_sum"] = rank_sum(df)
    models = Surrogate.load(f"{target}_nanobody").main

    compliant = df[(df["v_identity"] >= MIN_HUMAN_V_IDENTITY + GERMLINE_MARGIN)
                   & (df["j_identity"] >= MIN_HUMAN_J_IDENTITY)]
    log(f"{len(compliant)} of {len(df)} archived sequences clear today's germline floors with margin")

    # Framework donors: the right V gene, best archived rank sum first.
    donors = compliant[compliant["v_genes"].str.split(",").apply(lambda g: allowed_v in g)]
    donors = donors.nsmallest(n_donors, "true_rank_sum").copy()
    if donors.empty:
        log(f"no archived sequence calls {allowed_v}; cannot build a framework")
        out["status"] = "no_framework_donor"
        return out
    donors["regions"] = [regions(s) for s in donors["sequence"]]
    donors = donors[donors["regions"].notna()]

    # CDR3 donors: exactly the required length, from the best sequences.
    right_length = compliant[compliant["cdr3"].str.len() == cdr3_len]
    if right_length.empty:
        right_length = df[df["cdr3"].str.len() == cdr3_len]
        log(f"no compliant CDR3 of length {cdr3_len}; falling back to the whole archive")
    cdr3_donors = list(right_length.nsmallest(n_cdr3, "true_rank_sum")["cdr3"])
    if not cdr3_donors:
        log(f"no archived CDR3 of length {cdr3_len}; cannot graft")
        out["status"] = "no_cdr3_donor"
        return out
    out["framework_donors"] = len(donors)
    out["cdr3_donors"] = len(cdr3_donors)
    log(f"{len(donors)} framework donors ({allowed_v}), {len(cdr3_donors)} CDR3 donors of length {cdr3_len}")

    # Build the CDR3 generator. PSSM uses all compliant CDR3s (not just the
    # right-length donors) so it has enough data for every column.
    generator = make_generator(generator_name,
                               cdr3_sequences=list(compliant["cdr3"].dropna()))
    out["generator"] = generator_name
    log(f"CDR3 generator: {generator_name}")

    rng = random.Random(seed)
    generated = generate(donors, cdr3_donors, candidates, rng, generator=generator)
    out["generated"] = len(generated)
    log(f"generated {len(generated)} distinct candidates; checking validity...")

    archived_cdr3 = [c for c in df["cdr3"] if len(c) == cdr3_len]

    def cdr3_is_novel(seq: str) -> bool:
        """Below the validator's CDR3 duplicate threshold against the archive."""
        mine = regions(seq)
        if mine is None:
            return False
        return all(cdr_similarity(mine, {"CDR3": other})["CDR3"] < MAX_CDR3_IDENTITY
                   for other in archived_cdr3)

    valid, reasons = [], {}
    for seq in generated:
        problems = check_sequence(seq, allowed_v, cdr3_len)  # archive checked later
        if not problems and not cdr3_is_novel(seq):
            problems = [f"CDR3 >= {MAX_CDR3_IDENTITY:.0%} identical to an archived CDR3"]
        if problems:
            key = problems[0].split(" (")[0][:40]
            reasons[key] = reasons.get(key, 0) + 1
        else:
            valid.append(seq)
    out["valid"], out["reasons"] = len(valid), reasons
    log(f"{len(valid)} pass the offline checks. Main rejection causes:")
    for reason, n in sorted(reasons.items(), key=lambda kv: -kv[1])[:6]:
        log(f"  {n:5d}  {reason}")
    if not valid:
        out["status"] = "no_valid_candidate"
        return out

    archive = NanobodyArchive(target)

    if skip_surrogate:
        # Mining path: PSH + BoltzGen do the real selection, so skip the
        # surrogate entirely and just pass valid, de-duplicated candidates.
        log(f"skip_surrogate=True: passing {len(valid)} valid candidates directly (no surrogate ranking)")
        rows = []
        for seq in valid:
            r = regions(seq)
            g = germline_identity(seq)
            warnings = []
            if check_sequence(seq, allowed_v, cdr3_len, archive, warnings=warnings):
                continue
            rows.append({"sequence": seq, "cdr3": r["CDR3"] if r else "",
                         "predicted_rank_sum": np.nan, "v_identity": g["V"]["identity_pct"] if g else np.nan,
                         "design_iiptm": np.nan, "interaction_pae": np.nan,
                         "archive_warning": "; ".join(warnings)})
            if len(rows) >= shortlist:
                break
        if not rows:
            out["status"] = "all_duplicates"
            return out
        out["shortlist"] = pd.DataFrame(rows)
        return out

    # The real competition: the most recent dated submissions, which all comply.
    dated = df[df["epoch"].notna()]
    reference = dated.nlargest(REFERENCE_POOL, "epoch") if len(dated) else compliant.head(REFERENCE_POOL)
    ranked, reference_ranks = score(valid, models, reference)
    out["reference_size"] = len(reference)
    log(f"\nranking {len(ranked)} candidates against {len(reference)} recent archived sequences")

    rows = []
    for row in ranked.itertuples():
        warnings = []
        if check_sequence(row.sequence, allowed_v, cdr3_len, archive, warnings=warnings):
            continue
        rows.append({"sequence": row.sequence, "cdr3": row.cdr3,
                     "predicted_rank_sum": row.predicted_rank_sum,
                     "v_identity": row.v_identity,
                     "design_iiptm": row.design_iiptm, "interaction_pae": row.interaction_pae,
                     "archive_warning": "; ".join(warnings)})
        if len(rows) >= shortlist:
            break
    if not rows:
        out["status"] = "all_duplicates"
        return out
    table = pd.DataFrame(rows)
    out["shortlist"] = table
    best = table["predicted_rank_sum"].iat[0]
    out["best_rank_sum"] = float(best)
    out["place"] = int((reference_ranks < best).sum()) + 1
    out["reference_best"] = float(reference_ranks.min())
    out["reference_median"] = float(np.median(reference_ranks))
    log(f"\nranked in one pool with {len(reference)} recent real submissions "
        f"(which score {reference_ranks.min():.0f} best, {np.median(reference_ranks):.0f} median):")
    log(f"  our best candidate scores {best:.0f} -> place {out['place']} of {len(reference) + 1}")
    log("  predicted only: confirm with BoltzGen before trusting it")
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--epoch", type=int)
    parser.add_argument("--target", default=NANOBODY_TARGET[0])
    parser.add_argument("--candidates", type=int, default=3000)
    parser.add_argument("--donors", type=int, default=40, help="framework donors to draw from")
    parser.add_argument("--cdr3-pool", type=int, default=60, help="CDR3 donors of the right length")
    parser.add_argument("--shortlist", type=int, default=10)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--archive-cutoff", type=int,
                        help="ignore archived sequences from this epoch onward (replays)")
    parser.add_argument("--generator", choices=sorted(GENERATORS), default="pssm",
                        help="CDR3 generation strategy (default: pssm)")
    args = parser.parse_args()

    p = resolve_epoch(args.epoch)
    result = propose(p, args.target, args.candidates, args.donors, args.cdr3_pool,
                     args.shortlist, args.seed, args.archive_cutoff,
                     generator_name=args.generator)
    out = result["shortlist"]
    if out.empty:
        print(f"nothing to propose ({result['status']})")
        return
    with pd.option_context("display.width", 200, "display.max_colwidth", 30):
        print(out[["cdr3", "predicted_rank_sum", "v_identity", "design_iiptm",
                   "interaction_pae", "archive_warning"]].to_string(index=False))

    out_dir = DATA_DIR / "proposals"
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{args.target}_epoch{p['epoch']}_nanobodies.csv"
    out.to_csv(path, index=False)
    (out_dir / f"{args.target}_epoch{p['epoch']}_submission.txt").write_text(f"~|{out['sequence'].iat[0]}\n")
    print(f"\nshortlist saved to {path}")
    print("Re-score it with BoltzGen on a GPU before submitting: the surrogate only orders candidates.")


if __name__ == "__main__":
    main()
