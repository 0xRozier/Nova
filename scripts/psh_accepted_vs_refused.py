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


def fit(obs: dict, acceptance: float) -> pd.DataFrame:
    """Grid over (m, s, sigma_v, bias); score on the statistics we measured."""
    rows = []
    for m in np.arange(84, 109, 1.0):
        for s in (2.0, 3.0, 4.0, 6.0, 8.0, 10.0):
            for sigma_v in (0.0, 1.0, 2.0, 4.01, 6.0, 8.0, 10.0, 12.0, 15.0):
                for bias in np.arange(-4.0, 15.0, 1.0):
                    pred = predicted(m, s, sigma_v, bias)
                    if not pred:
                        continue
                    # Each statistic in units of its own measurement scatter.
                    cost = (((pred["median_accepted"] - obs["median_accepted"]) / 1.0) ** 2
                            + ((pred["median_refused"] - obs["median_refused"]) / 1.0) ** 2
                            + ((pred["sd_accepted"] - obs["sd_accepted"]) / 1.5) ** 2
                            + ((pred["sd_refused"] - obs["sd_refused"]) / 1.5) ** 2
                            + ((pred["auc"] - obs["auc"]) / 0.04) ** 2
                            + ((pred["acceptance"] - acceptance) / 0.08) ** 2)
                    rows.append({"m": m, "s": s, "sigma_v": sigma_v, "bias": bias,
                                 "cost": cost, **pred})
    return pd.DataFrame(rows).sort_values("cost").reset_index(drop=True)


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
        acceptance = float(group["accepted"].mean())
        top = fit(obs, acceptance)
        if top.empty:
            continue
        best = top.iloc[0]
        print(f"  ajustement : PSH vraie ~ N({best['m']:.0f}, {best['s']:.0f}), "
              f"bruit du validateur sigma_v = {best['sigma_v']:.2f}, "
              f"notre biais = {best['bias']:+.0f}")
        near = top[top["cost"] <= top["cost"].iloc[0] + 4]
        print(f"  parmi les {len(near)} ajustements comparables : "
              f"sigma_v de {near['sigma_v'].min():.2f} a {near['sigma_v'].max():.2f}, "
              f"biais de {near['bias'].min():+.0f} a {near['bias'].max():+.0f}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    cols = [c for c in ("name", "epoch", "uid", "band", "epoch_length", "accepted",
                        "cdr3", "PSH", "PSH_flag", "reds", "red_on", "sequence")
            if c in merged.columns]
    merged[cols].to_csv(args.out, index=False)
    print(f"\nsauve dans {args.out}")


if __name__ == "__main__":
    main()
