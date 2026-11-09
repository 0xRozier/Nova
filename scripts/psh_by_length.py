"""Is our PSH reproduction biased by CDR3 length?

Every submission the validators scored passed their PSH gate, so for each of
them the true PSH is <= 89 by construction. That makes the scored field a
ground-truth set we did not have to build: run our own PSH on it, bucket by
the epoch's imposed CDR3 length, and the share reading over 89 is the rate at
which our reproduction would have refused a sequence the validators accepted,
length by length.

Why it matters. Replaying 21 epochs offline, our best candidate ranked 1st in
six of the nine epochs with a CDR3 of 18 or more -- and in all nine, every
candidate read over 89 on our PSH (95 to 117, against a limit of 89). Two
readings fit that:

  a. long CDR3s really are unable to pass the gate, and a third of the
     calendar is closed to everyone;
  b. our PSH over-reads as the loop grows, and we have been refusing our own
     best epochs on a number we compute badly.

The accepted field separates them. Under (a) the accepted sequences at length
18+ read at or under 89 like ours do not; under (b) they read over it too.

A third fact makes the answer expensive to get wrong: the field thins as the
loop grows. Over 71 scored epochs the median field is 26 submissions at
length 5-8 and 9 at length 17-23 (Spearman rho = -0.63, p = 3e-09). The
epochs we rank best in are also the ones with the fewest rivals.

Usage:
    python -m scripts.psh_by_length                 # ~120 accepted sequences
    python -m scripts.psh_by_length --budget 240    # more, at 3.6 s each
    python -m scripts.psh_by_length --parse <log>   # re-read a finished run
"""
import argparse
from pathlib import Path

import pandas as pd

from nova_lab import DATA_DIR
from nova_lab.config import PSH_STRICT_REJECTION_THRESHOLD as PSH_STRICT
from nova_lab.nanobodies import regions_batch
from scripts.tnp_developability import parse, run

FIELD = DATA_DIR / "submissions" / "Q9NZQ7_nanobody_results.csv"
BUCKETS = [(5, 8, "court 5-8"), (9, 12, "moyen 9-12"),
           (13, 16, "long 13-16"), (17, 25, "tres long 17+")]


def bucket_of(length: int) -> str:
    for lo, hi, label in BUCKETS:
        if lo <= length <= hi:
            return label
    return "hors bornes"


def accepted_field(path: Path = FIELD) -> pd.DataFrame:
    """The scored submissions, with the CDR3 length their epoch imposed.

    A row in this file is a submission the validators scored, which means it
    cleared every gate, PSH included.
    """
    df = pd.read_csv(path)
    regs = regions_batch(df["sequence"].tolist())
    df["cdr3"] = [r["CDR3"] if r else "" for r in regs]
    df["cdr3_length"] = df["cdr3"].str.len()
    df = df[df["cdr3_length"] > 0].drop_duplicates("sequence")
    # The epoch imposes one length; a row that disagrees is an ANARCI
    # disagreement on the boundary, not a different challenge.
    imposed = df.groupby("epoch")["cdr3_length"].agg(lambda s: s.mode().iat[0])
    df["epoch_length"] = df["epoch"].map(imposed)
    df = df[df["cdr3_length"] == df["epoch_length"]]
    df["bucket"] = df["epoch_length"].map(bucket_of)
    return df


def sample(df: pd.DataFrame, budget: int, seed: int = 0) -> pd.DataFrame:
    """Spend the budget on the long buckets first: those are the question."""
    want = {"tres long 17+": 0.45, "long 13-16": 0.25,
            "moyen 9-12": 0.15, "court 5-8": 0.15}
    parts = []
    for label, share in want.items():
        sub = df[df["bucket"] == label]
        n = min(len(sub), max(1, round(budget * share)))
        parts.append(sub.sample(n, random_state=seed) if n < len(sub) else sub)
    return pd.concat(parts).reset_index(drop=True)


def report(tnp: pd.DataFrame, meta: pd.DataFrame, out: Path) -> None:
    m = meta.merge(tnp, on="name", how="left")
    scored = m[m["PSH"].notna()]
    print(f"\n{len(scored)} of {len(m)} sequences profiled by TNP")
    if scored.empty:
        return

    print(f"\n=== PSH de sequences que les validateurs ont ACCEPTEES ===")
    print(f"leur PSH reelle est <= {PSH_STRICT} par construction, donc toute lecture")
    print(f"au-dessus de {PSH_STRICT} est un refus a tort de notre reproduction\n")
    rows = []
    for _, label in sorted({(lo, lab) for lo, _, lab in BUCKETS}):
        sub = scored[scored["bucket"] == label]
        if sub.empty:
            continue
        over = sub["PSH"] > PSH_STRICT
        rows.append({"bucket": label, "n": len(sub),
                     "epoques": sub["epoch"].nunique(),
                     "PSH mediane": round(sub["PSH"].median(), 1),
                     "PSH 90e pct": round(sub["PSH"].quantile(0.9), 1),
                     f"> {PSH_STRICT}": f"{over.mean():.0%}",
                     "refus a tort": int(over.sum())})
    table = pd.DataFrame(rows)
    print(table.to_string(index=False))

    long_label = BUCKETS[-1][2]
    long_rows = scored[scored["bucket"] == long_label]
    if not long_rows.empty:
        share = (long_rows["PSH"] > PSH_STRICT).mean()
        median = long_rows["PSH"].median()
        print(f"\nlecture du resultat, sur {len(long_rows)} sequences a CDR3 >= 17 :")
        print(f"  {share:.0%} des sequences ACCEPTEES lisent au-dessus de {PSH_STRICT} chez"
              f" nous (mediane {median:.1f}).")
        if share > 0.5:
            print(f"  Ce qui est etabli : refuser un candidat parce que notre lecture")
            print(f"  depasse {PSH_STRICT} est faux a cette longueur, puisque c'est la lecture")
            print("  TYPIQUE d'une sequence acceptee. Le seuil plat doit devenir un")
            print("  quantile par longueur.")
            print("  Ce qui n'est PAS etabli : que notre implementation sur-lise. Deux")
            print("  mecanismes produisent la meme signature, et le second suffit")
            print("  quantitativement (voir FINDINGS 3ter) :")
            print("    a. notre pipeline lit trop haut quand la boucle s'allonge ;")
            print("    b. la mesure du validateur est bruitee elle aussi (sigma ~ 4) et sa")
            print("       porte selectionne la queue chanceuse d'une population dont la PSH")
            print("       vraie est vers 98-100 ; notre tirage independant retombe dessus.")
            print("  Les deux impliquent le meme correctif, mais pas le meme espoir : sous")
            print("  (a) un candidat long lu a 92 passe le plus souvent ; sous (b) il passe")
            print("  une fois sur cinq environ. Contre 0 % aujourd'hui, ou nous ne")
            print("  soumettons rien.")
        else:
            print("  notre reproduction tient a cette longueur. Le mur est reel, et nos")
            print("  candidats longs (PSH 95 a 117) en sont tres loin.")

    out.parent.mkdir(parents=True, exist_ok=True)
    cols = [c for c in ("name", "epoch", "uid", "bucket", "epoch_length", "cdr3",
                        "PSH", "PSH_flag", "C", "PPC", "PNC", "reds", "red_on",
                        "sequence") if c in m.columns]
    m[cols].to_csv(out, index=False)
    print(f"\nsauve dans {out}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--budget", type=int, default=120,
                        help="sequences to profile; about 3.6 s each on 4 cores")
    parser.add_argument("--ncores", type=int, default=4)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--work", type=Path, default=DATA_DIR / "psh_by_length")
    parser.add_argument("--parse", type=Path, help="parse an existing TNP log")
    parser.add_argument("--out", type=Path, default=Path("results/psh_by_length.csv"))
    args = parser.parse_args()

    field = accepted_field()
    print(f"{len(field)} soumissions notees, {field['epoch'].nunique()} epoques")
    print(field.groupby("bucket").agg(n=("sequence", "size"),
                                      epoques=("epoch", "nunique")).to_string())

    args.work.mkdir(parents=True, exist_ok=True)
    meta_path = args.work / "meta.csv"
    if args.parse and meta_path.exists():
        meta = pd.read_csv(meta_path)
        text = args.parse.read_text()
    else:
        chosen = sample(field, args.budget, args.seed)
        chosen["name"] = [f"e{r.epoch}_u{r.uid}_L{r.epoch_length}_{i}"
                          for i, r in enumerate(chosen.itertuples())]
        fasta = args.work / "sequences.fasta"
        fasta.write_text("".join(f">{r.name}\n{r.sequence}\n"
                                 for r in chosen.itertuples()))
        chosen.to_csv(meta_path, index=False)
        print(f"\n{len(chosen)} sequences retenues -> {fasta}")
        print(chosen.groupby("bucket").size().to_string())
        text = args.parse.read_text() if args.parse else run(fasta, args.ncores, args.work)
        meta = chosen
    report(parse(text), meta, args.out)


if __name__ == "__main__":
    main()
