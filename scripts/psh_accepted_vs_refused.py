"""Does our PSH separate the submissions the validators refused from those they kept?

scripts.psh_by_length established that at CDR3 >= 17 a reading of 92 is what an
ACCEPTED sequence reads: refusing on a flat 89 there discards three quarters of
the population the validator admits. It did not establish why, and two
mechanisms fit the same signature:

  a. our pipeline over-reads as the loop grows (bias b > 0), the validator's own
     measurement being near-deterministic, so accepted and refused split
     sharply in the true value;
  b. the validator's measurement is as noisy as ours (sigma = 4.01, measured),
     so its gate keeps the lucky tail of a population whose true PSH is above
     89, and our independent draw lands back on the true value.

A constant bias shifts both groups equally, so the gap between their medians
does not tell them apart. What does is how much the two groups OVERLAP: under
(a) the split is sharp and our reading separates them well; under (b) the two
groups come from the same population cut by a noisy threshold, and the overlap
is wide.

So: measure our PSH on both groups at the same lengths, and fit

    mu ~ N(m, s)                  true PSH of what miners submit
    validator reads mu + e,       e ~ N(0, sigma_v), accepts iff <= 89
    we read mu + b + f,           f ~ N(0, 4.01)

sigma_v and b are what the two stories disagree about: (a) is sigma_v ~ 0 with
b > 0, (b) is sigma_v = 4.01 with b = 0. The fit is reported with its
uncertainty, not as a verdict.

Refused sequences are refused for any gate, not only PSH, so those our own
offline checker already rejects are dropped first: what remains is the set
refused despite passing everything we can verify.

Usage:
    python -m scripts.psh_accepted_vs_refused                # 180 sequences
    python -m scripts.psh_accepted_vs_refused --parse <log>  # re-read a run
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from nova_lab import DATA_DIR
from nova_lab.config import PSH_STRICT_REJECTION_THRESHOLD as PSH_STRICT
from nova_lab.nanobodies import check_sequence, regions_batch, seq_hash
from scripts.tnp_developability import parse, run

PSH_SIGMA = 4.01  # measured on 10 sequences profiled twice
LONG, SHORT = (17, 25), (5, 8)


def paired_field() -> pd.DataFrame:
    """Every decrypted submission, with whether the validators scored it.

    A submission with metrics cleared every gate. One without, in an epoch
    where others were scored, was refused.
    """
    scores = pd.read_json(DATA_DIR / "scores" / "nanobodies.jsonl", lines=True)
    sequence_of = {}
    with (DATA_DIR / "submissions" / "submissions.jsonl").open() as handle:
        for line in handle:
            row = json.loads(line)
            for seq in row.get("sequences") or []:
                sequence_of[(seq_hash(seq), row["epoch"])] = seq
    scores["sequence"] = [sequence_of.get((h, e))
                          for h, e in zip(scores["hash"], scores["epoch"])]
    scores["accepted"] = scores["metrics"].notna()
    df = scores[scores["sequence"].notna()].drop_duplicates(["sequence", "epoch"]).copy()

    # An epoch nobody was scored in was not processed, so its blanks say nothing.
    scored_epochs = set(df.loc[df["accepted"], "epoch"])
    df = df[df["epoch"].isin(scored_epochs)]

    regs = regions_batch(df["sequence"].tolist())
    df["cdr3"] = [r["CDR3"] if r else "" for r in regs]
    imposed = (df[df["accepted"]].assign(n=df["cdr3"].str.len())
               .groupby("epoch")["n"].agg(lambda s: s.mode().iat[0]))
    df["epoch_length"] = df["epoch"].map(imposed)
    return df[df["epoch_length"].notna() & (df["cdr3"].str.len() == df["epoch_length"])]


def drop_other_gates(df: pd.DataFrame, log=print) -> pd.DataFrame:
    """Keep refusals our own checker cannot already explain."""
    keep, dropped = [], 0
    for row in df.itertuples():
        if row.accepted:
            keep.append(True)
            continue
        # No archive: the uniqueness rule depends on the archive as it stood
        # that epoch, which we cannot reconstruct, so we do not count it.
        problems = check_sequence(row.sequence, allowed_v=None,
                                  allowed_cdr3_length=int(row.epoch_length))
        keep.append(not problems)
        dropped += bool(problems)
    log(f"  {dropped} refusals our own gates already explain, dropped")
    return df[pd.Series(keep, index=df.index)]


def predicted(m: float, s: float, sigma_v: float, bias: float,
              n: int = 40000, seed: int = 0) -> dict:
    """Summary statistics the model implies, by simulation."""
    rng = np.random.default_rng(seed)
    mu = rng.normal(m, s, n)
    accepted = mu + rng.normal(0, sigma_v, n) <= PSH_STRICT
    ours = mu + bias + rng.normal(0, PSH_SIGMA, n)
    if accepted.sum() < 50 or (~accepted).sum() < 50:
        return {}
    return {"median_accepted": float(np.median(ours[accepted])),
            "median_refused": float(np.median(ours[~accepted])),
            "sd_accepted": float(np.std(ours[accepted])),
            "sd_refused": float(np.std(ours[~accepted])),
            "acceptance": float(accepted.mean()),
            "auc": float(auc(ours[~accepted], ours[accepted]))}


def auc(refused: np.ndarray, accepted: np.ndarray) -> float:
    """P(a refused sequence reads higher than an accepted one), ties at half.

    Mann-Whitney U from the rank sum: the pairwise form is a |refused| x
    |accepted| matrix, which at simulation sizes is tens of millions of
    comparisons per grid point and turned a parameter sweep into half an hour.
    """
    both = np.concatenate([refused, accepted])
    ranks = pd.Series(both).rank().to_numpy()  # averages ties, as U requires
    u = ranks[:len(refused)].sum() - len(refused) * (len(refused) + 1) / 2
    return float(u / (len(refused) * len(accepted)))


def observed(group: pd.DataFrame) -> dict:
    acc = group.loc[group["accepted"], "PSH"].to_numpy()
    ref = group.loc[~group["accepted"], "PSH"].to_numpy()
    return {"median_accepted": float(np.median(acc)), "median_refused": float(np.median(ref)),
            "sd_accepted": float(np.std(acc)), "sd_refused": float(np.std(ref)),
            "auc": float(auc(ref, acc)), "n_accepted": len(acc), "n_refused": len(ref)}


# The real base rates at each band, from the whole paired field rather than
# from the sample, so the curve below reads as a probability a miner can use.
BASE_RATES = {"long 17+": (741, 692), "court 5-8": (651, 218)}


def wilson(successes: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Binomial interval that stays inside [0, 1] at small n and extreme p."""
    p = successes / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return centre - half, centre + half


def acceptance_curve(group: pd.DataFrame, band: str,
                     edges=(88, 92, 96, 100)) -> pd.DataFrame:
    """P(the validators accept it | our PSH reads X), at this band's base rate.

    This replaces a four-parameter fit of the measurement model that the data
    cannot identify: see FINDINGS 3ter. A quarter of the refusals at 17+ read
    BELOW the accepted median, so they were refused for something our PSH does
    not see (Aggrescan3D, archive similarity). The fit had no term for that and
    absorbed it as validator noise, which is why it ran to whatever ceiling the
    grid allowed. The question a miner actually needs answered does not need
    the model: given a reading, how often was such a sequence accepted?
    """
    n_accepted, n_refused = BASE_RATES[band]
    sampled_accepted = int(group["accepted"].sum())
    sampled_refused = len(group) - sampled_accepted
    if not (sampled_accepted and sampled_refused):
        return pd.DataFrame()
    w_accepted = n_accepted / sampled_accepted
    w_refused = n_refused / sampled_refused

    def to_base_rate(share: float) -> float:
        """A sample share, reweighted to the band's real acceptance rate."""
        odds = share * w_accepted
        return odds / (odds + (1 - share) * w_refused)

    bounds = [-np.inf, *edges, np.inf]
    rows = []
    for low, high in zip(bounds, bounds[1:]):
        band_rows = group[(group["PSH"] > low) & (group["PSH"] <= high)]
        if band_rows.empty:
            continue
        k, n = int(band_rows["accepted"].sum()), len(band_rows)
        low_ci, high_ci = wilson(k, n)
        label = (f"<= {high:.0f}" if low == -np.inf else
                 f"> {low:.0f}" if high == np.inf else f"{low:.0f}-{high:.0f}")
        rows.append({"notre PSH": label, "n": n, "acceptees": k,
                     "P(accepte)": f"{to_base_rate(k / n):.0%}",
                     "IC95": f"{to_base_rate(low_ci):.0%}-{to_base_rate(high_ci):.0%}"})
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--budget", type=int, default=180)
    parser.add_argument("--ncores", type=int, default=4)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--work", type=Path, default=DATA_DIR / "psh_accepted_vs_refused")
    parser.add_argument("--parse", type=Path)
    parser.add_argument("--out", type=Path, default=Path("results/psh_accepted_vs_refused.csv"))
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    meta_path = args.work / "meta.csv"

    if args.parse and meta_path.exists():
        chosen, text = pd.read_csv(meta_path), args.parse.read_text()
    else:
        field = paired_field()
        print(f"{len(field)} soumissions appariees, "
              f"{int(field['accepted'].sum())} acceptees, "
              f"{int((~field['accepted']).sum())} refusees")
        field = drop_other_gates(field)
        parts = []
        for (lo, hi), share in ((LONG, 0.7), (SHORT, 0.3)):
            band = field[field["epoch_length"].between(lo, hi)]
            per_group = max(1, round(args.budget * share / 2))
            for accepted in (True, False):
                sub = band[band["accepted"] == accepted]
                n = min(len(sub), per_group)
                parts.append(sub.sample(n, random_state=args.seed) if n < len(sub) else sub)
                print(f"  CDR3 {lo}-{hi}, {'acceptees' if accepted else 'refusees'}: "
                      f"{n} de {len(sub)} disponibles")
        chosen = pd.concat(parts).reset_index(drop=True)
        chosen["band"] = np.where(chosen["epoch_length"] >= LONG[0], "long 17+", "court 5-8")
        chosen["name"] = [f"e{r.epoch}_L{int(r.epoch_length)}_"
                          f"{'acc' if r.accepted else 'ref'}_{i}"
                          for i, r in enumerate(chosen.itertuples())]
        fasta = args.work / "sequences.fasta"
        fasta.write_text("".join(f">{r.name}\n{r.sequence}\n" for r in chosen.itertuples()))
        chosen.to_csv(meta_path, index=False)
        print(f"\n{len(chosen)} sequences -> {fasta}")
        text = run(fasta, args.ncores, args.work)

    merged = chosen.merge(parse(text), on="name", how="left")
    merged = merged[merged["PSH"].notna()]
    print(f"\n{len(merged)} sequences profilees")

    for band, group in merged.groupby("band"):
        obs = observed(group)
        print(f"\n=== {band} ===")
        print(f"  acceptees n={obs['n_accepted']}: mediane {obs['median_accepted']:.1f}, "
              f"ecart-type {obs['sd_accepted']:.1f}")
        print(f"  refusees  n={obs['n_refused']}: mediane {obs['median_refused']:.1f}, "
              f"ecart-type {obs['sd_refused']:.1f}")
        print(f"  ecart des medianes {obs['median_refused'] - obs['median_accepted']:+.1f}")
        print(f"  AUC {obs['auc']:.2f}  (0,50 = notre PSH ne separe pas du tout ; "
              f"1,00 = separation parfaite)")
        below = group[~group["accepted"] & (group["PSH"] < obs["median_accepted"])]
        print(f"  {len(below)} des {obs['n_refused']} refusees lisent SOUS la mediane des "
              f"acceptees : refusees pour ce que notre PSH ne voit pas")
        curve = acceptance_curve(group, band)
        if not curve.empty:
            print(f"  P(accepte | notre PSH), au taux de base reel du champ "
                  f"({BASE_RATES[band][0]} acceptees / {BASE_RATES[band][1]} refusees) :")
            print("    " + curve.to_string(index=False).replace("\n", "\n    "))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    cols = [c for c in ("name", "epoch", "uid", "band", "epoch_length", "accepted",
                        "cdr3", "PSH", "PSH_flag", "reds", "red_on", "sequence")
            if c in merged.columns]
    merged[cols].to_csv(args.out, index=False)
    print(f"\nsauve dans {args.out}")


if __name__ == "__main__":
    main()
