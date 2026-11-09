"""Does the CDR3 strategy actually lower PSH? Measured, not proxied.

scripts.compare_generators reports hydropathy and AA-frequency divergence, but
those are proxies: on 205 sequences with measured PSH, the best sequence-level
hydrophobicity proxy reaches Spearman 0.30 once CDR length is regressed out, and
it changes sign with CDR3 length. A strategy can therefore look much better on
hydropathy and move PSH not at all.

PSH is what the validator rejects on (> 89), so this measures PSH itself, with
TNP, on candidates each strategy actually produces.

Usage:
    python -m scripts.psh_by_generator --epoch 26174 --per-generator 30
"""
import argparse
import random
from pathlib import Path

import pandas as pd

from nova_lab import DATA_DIR
from nova_lab.config import MIN_HUMAN_J_IDENTITY, MIN_HUMAN_V_IDENTITY, NANOBODY_TARGET
from nova_lab.config import PSH_STRICT_REJECTION_THRESHOLD as PSH_STRICT
from nova_lab.cdr3_generators import make_generator
from nova_lab.nanobodies import check_sequence, regions
from scripts.propose_nanobody import GERMLINE_MARGIN, generate
from scripts.train_nanobody_surrogate import rank_sum


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--epoch", type=int, required=True)
    parser.add_argument("--allowed-v", help="skip the chain lookup")
    parser.add_argument("--cdr3", type=int, help="skip the chain lookup")
    parser.add_argument("--target", default=NANOBODY_TARGET[0])
    parser.add_argument("--generators", nargs="+",
                        default=["chimeric", "pssm", "pssm_denovo"])
    parser.add_argument("--per-generator", type=int, default=30)
    parser.add_argument("--donors", type=int, default=40)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    if args.allowed_v and args.cdr3:
        allowed_v, cdr3_len = args.allowed_v, args.cdr3
    else:
        from scripts.check_nanobody import resolve_epoch
        p = resolve_epoch(args.epoch)
        allowed_v, cdr3_len = p["allowed_v"], p["allowed_cdr3_length"]
    print(f"epoch {args.epoch}: V {allowed_v}, CDR3 {cdr3_len}")

    df = pd.read_parquet(DATA_DIR / "dataset" / f"{args.target}_nanobodies.parquet")
    df["true_rank_sum"] = rank_sum(df)
    df = df[df["epoch"].isna() | (df["epoch"] < args.epoch)]
    compliant = df[(df["v_identity"] >= MIN_HUMAN_V_IDENTITY + GERMLINE_MARGIN)
                   & (df["j_identity"] >= MIN_HUMAN_J_IDENTITY)]
    donors = compliant[compliant["v_genes"].str.split(",").apply(lambda g: allowed_v in g)]
    donors = donors.nsmallest(args.donors, "true_rank_sum").copy()
    donors["regions"] = [regions(s) for s in donors["sequence"]]
    donors = donors[donors["regions"].notna()]
    if donors.empty:
        raise SystemExit(f"no donor calls {allowed_v}")
    pool = list(compliant[compliant["cdr3"].str.len() == cdr3_len]
                .nsmallest(60, "true_rank_sum")["cdr3"])
    if not pool:
        raise SystemExit(f"no archived CDR3 of length {cdr3_len}")
    print(f"{len(donors)} donors, {len(pool)} CDR3 donors")

    rows = []
    for name in args.generators:
        gen = make_generator(name, list(compliant["cdr3"].dropna()))
        rng = random.Random(args.seed)
        made, tries = [], 0
        while len(made) < args.per_generator and tries < 60 * args.per_generator:
            tries += 1
            for seq in generate(donors, pool, 1, rng, generator=gen):
                if not check_sequence(seq, allowed_v, cdr3_len):
                    made.append(seq)
        print(f"  {name}: {len(made)} valid candidates in {tries} tries")
        rows += [{"name": f"{name}_{i}", "sequence": s, "group": name,
                  "epoch": args.epoch} for i, s in enumerate(made)]

    work = DATA_DIR / "psh_by_generator"
    work.mkdir(parents=True, exist_ok=True)
    meta = pd.DataFrame(rows).drop_duplicates("sequence")
    meta.to_csv(work / "meta.csv", index=False)
    (work / "sequences.fasta").write_text(
        "".join(f">{r.name}\n{r.sequence}\n" for r in meta.itertuples()))

    from scripts.tnp_developability import parse, run
    df_psh = parse(run(work / "sequences.fasta", 4, work))
    merged = meta.merge(df_psh[["name", "PSH"]], on="name", how="inner")

    print(f"\nPSH measured on {len(merged)} candidates (validator refuses > {PSH_STRICT}):")
    print(f"{'strategy':<14}{'n':>4}{'median':>9}{'min':>8}{'<=89':>8}{'<=84.4':>9}")
    for name, sub in merged.groupby("group"):
        print(f"{name:<14}{len(sub):>4}{sub['PSH'].median():>9.1f}{sub['PSH'].min():>8.1f}"
              f"{(sub['PSH'] <= PSH_STRICT).mean():>7.0%}{(sub['PSH'] <= 84.4).mean():>9.0%}")
    out = work / "psh_by_generator.csv"
    merged.to_csv(out, index=False)
    print(f"\nsaved to {out}")
    print("84.4 is the usable threshold with a 3-fold mean: 89 minus two "
          "deviations of 2.32 (one measurement deviates by 4.01)")


if __name__ == "__main__":
    main()
