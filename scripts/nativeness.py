"""Run the validator's OWN nativeness filter, instead of reconstructing it.

This repo has now reconstructed three validator rules and got two of them
wrong. FINDINGS 54: the novelty gate copied an identity rule from the filter's
search *engine*, which the diversity validator never calls -- the real rule is a
k-mer Jaccard under 0.9 against the top 50, and ours threw away half the
generator's output. FINDINGS 67: `MIN_HUMAN_V_IDENTITY = 82` and
`MIN_HUMAN_J_IDENTITY = 86` are ours, not the filter's; it hard-rejects on V
identity under 0.70, V coverage under 0.80 and a CDR3 outside 6-35, and has no J
floor at all.

So this does not reproduce anything. It calls
`metanano.utils.igblast_nativeness.run` -- the filter's own code, on the filter's
own pre-built camelid and human databases -- and reports what it says.

What it took to make that runnable here, since the filter's installer assumes
conda:

  - IgBLAST 1.22.0 binaries from NCBI, unpacked to /home/user/igblast/bin;
  - the pre-built `camelid_V` and `human_V` protein databases, which ship inside
    NOVA-nanobody-filter at custom_data/igblast_database/database;
  - `internal_data/camelid` from the same repo, AND `internal_data/human` from
    the NCBI tarball. Without the second, every call died on "Germline
    annotation database .../internal_data/human/human_V could not be found" --
    and `NativenessFilter.analyze` swallows that exception and returns None, so
    the failure reads exactly like a refusal. The four known-accepted
    submissions it first scored as refused were this, not a verdict.

Usage:
    python -m scripts.nativeness --file candidates.csv --column sequence
    python -m scripts.nativeness --sequence EVQLVESGGG...
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

FILTER_ROOT = Path(os.environ.get("NOVA_NANOFILTER", "/home/user/nanofilter"))
IGBLAST_DIR = Path(os.environ.get("IGDATA", "/home/user/igblast"))

# What the filter's own config declares. Measured against sequences whose verdict
# we know, on 23 November:
#
#   hard rejects      0 of 60 accepted submissions, 0 of 60 archive binders
#                     -> not what binds, for either target
#   nativeness >= 0.80  0 of 60 ACCEPTED submissions reach it (median 0.507,
#                     max 0.596) -> cannot be the deployed gate, or nothing
#                     would ever be accepted
#   humanness >= 0.75  40 of 40 accepted submissions pass; only 21 % of the
#                     P05231 archive's binders do -> consistent with being
#                     enforced, and it is the one threshold that discriminates
NATIVENESS_THRESHOLD = 0.80
HUMANNESS_THRESHOLD = 0.75
ACCEPTED_NATIVENESS_RANGE = (0.386, 0.596)   # observed on 60 accepted submissions


def available() -> tuple[bool, str]:
    """(usable, why). Checked before scoring so a missing piece is not a refusal."""
    if not (FILTER_ROOT / "metanano").is_dir():
        return False, f"NOVA-nanobody-filter absent de {FILTER_ROOT}"
    if not (IGBLAST_DIR / "bin" / "igblastp").exists():
        return False, f"igblastp absent de {IGBLAST_DIR / 'bin'}"
    for db in ("camelid_V.phr", "human_V.phr"):
        if not (IGBLAST_DIR / "database" / db).exists():
            return False, f"base {db} absente de {IGBLAST_DIR / 'database'}"
    for org in ("camelid", "human"):
        if not (IGBLAST_DIR / "internal_data" / org).is_dir():
            return False, (f"internal_data/{org} absent : igblastp echoue et le "
                           f"filtre rend None, ce qui se lit comme un refus")
    return True, "ok"


def score(sequences: list[str]) -> pd.DataFrame:
    """One row per sequence, straight from the filter's own scorer."""
    usable, why = available()
    if not usable:
        raise SystemExit(f"le filtre n'est pas executable ici : {why}")
    sys.path.insert(0, str(FILTER_ROOT))
    os.environ.setdefault("IGDATA", str(IGBLAST_DIR))
    from metanano.utils import igblast_nativeness as ign

    rows = []
    for seq in sequences:
        try:
            result = ign.run(seq)
            first = result[0] if isinstance(result, list) and result else {}
            rows.append({
                "sequence": seq,
                "nativeness": first.get("vhh_nativeness"),
                "humanness": first.get("human_framework"),
                "hard_reject": bool(first.get("hard_reject")),
                "reject_reason": first.get("reject_reason"),
            })
        except Exception as error:  # noqa: BLE001 - one bad sequence must not stop the batch
            rows.append({"sequence": seq, "nativeness": None, "humanness": None,
                         "hard_reject": True, "reject_reason": f"{type(error).__name__}: {error}"})
    return pd.DataFrame(rows)


def verdict(frame: pd.DataFrame) -> pd.DataFrame:
    """Adds the two readings: the hard rejects, and the humanness gate.

    The nativeness threshold is deliberately NOT applied: no accepted submission
    we have reaches 0.80, so applying it would refuse everything the validators
    took. It is reported, not enforced.
    """
    frame = frame.copy()
    frame["passes_hard"] = ~frame["hard_reject"]
    frame["passes_humanness"] = frame["humanness"].fillna(0) >= HUMANNESS_THRESHOLD
    frame["in_accepted_nativeness_range"] = (
        frame["nativeness"].fillna(-1) >= ACCEPTED_NATIVENESS_RANGE[0])
    return frame


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=Path, help="CSV de sequences")
    parser.add_argument("--column", default="sequence")
    parser.add_argument("--sequence", action="append", default=[])
    parser.add_argument("--sample", type=int, help="n'en noter qu'un echantillon")
    parser.add_argument("--out", type=Path, help="ecrire le detail en CSV")
    args = parser.parse_args()

    usable, why = available()
    print(f"filtre : {why}")
    if not usable:
        return 1

    seqs = list(args.sequence)
    if args.file:
        frame = pd.read_csv(args.file)
        if args.column not in frame.columns:
            raise SystemExit(f"pas de colonne {args.column!r} dans {args.file}")
        col = frame[args.column].dropna().astype(str)
        col = col[col.str.isalpha()].drop_duplicates()
        if args.sample:
            col = col.sample(min(args.sample, len(col)), random_state=0)
        seqs += list(col)
    if not seqs:
        raise SystemExit("aucune sequence : passer --file ou --sequence")

    out = verdict(score(seqs))
    n = len(out)
    print(f"\n{n} sequences notees par le filtre lui-meme\n")
    print(f"  rejets durs                  : {(~out.passes_hard).sum()}/{n}")
    if (~out.passes_hard).any():
        print(f"    motifs : {out[~out.passes_hard].reject_reason.value_counts().head(3).to_dict()}")
    ok = out[out.passes_hard]
    if len(ok):
        print(f"  nativite   : mediane {ok.nativeness.median():.3f} "
              f"(min {ok.nativeness.min():.3f}, max {ok.nativeness.max():.3f})")
        print(f"    seuil annonce {NATIVENESS_THRESHOLD} : {(ok.nativeness>=NATIVENESS_THRESHOLD).sum()}/{len(ok)}"
              f"  -- NON applique, aucune acceptee connue ne l'atteint")
        print(f"    dans la plage des acceptees (>= {ACCEPTED_NATIVENESS_RANGE[0]}) : "
              f"{ok.in_accepted_nativeness_range.sum()}/{len(ok)}")
        print(f"  humanite   : mediane {ok.humanness.median():.3f}")
        print(f"    seuil {HUMANNESS_THRESHOLD} : {ok.passes_humanness.sum()}/{len(ok)} = "
              f"{ok.passes_humanness.mean():.0%}  -- le seul qui discrimine")
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(args.out, index=False)
        print(f"\ndetail -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
