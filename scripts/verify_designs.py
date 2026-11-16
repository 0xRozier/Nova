"""Are the designs nanobodies at all? The check that was missing since §11.

Two runs of BoltzGen design were read as evidence before anyone counted the
cysteines. All 210 designs had **none** (FINDINGS 55): the vendored BoltzGen
lacks the inverse-folding constraint mask, so `filter_cysteine=true` stripped
the canonical Cys23-Cys104 disulfide from positions the spec had frozen. Those
sequences could not fold like a nanobody and `min_cysteines` would have refused
them outright -- yet their low design_iiptm was published as a refutation of
structure-conditioned design.

So no design run is read again until it passes this, and a failure here is a
broken pipeline, never a result about the method.

    python -m scripts.verify_designs --file designs.csv --frozen <scaffold.fasta>
"""
import argparse
from pathlib import Path

import pandas as pd

# The guide's universally conserved positions, as motifs rather than indices:
# IMGT numbering has gaps and the designs are not numbered here.
CANONICAL_CYS_PAIR = 2          # Cys23 + Cys104
MIN_CYSTEINES = 1               # the validator's Layer 1 floor, typically >= 1


def cysteines(sequence: str) -> int:
    return sequence.count("C")


def frozen_intact(design: str, scaffold: str, head: int = 26, tail: int = 11) -> bool:
    """FR1 and FR4 came from the scaffold and the spec froze them."""
    return design[:head] == scaffold[:head] and design[-tail:] == scaffold[-tail:]


def check(designs: pd.DataFrame, column: str, scaffolds: list[str] | None = None) -> dict:
    """Returns a report; `ok` is False as soon as anything is off."""
    sequences = designs[column].dropna().astype(str)
    cys = sequences.map(cysteines)
    report = {
        "n": int(len(sequences)),
        "no_cysteine": int((cys == 0).sum()),
        "below_canonical_pair": int((cys < CANONICAL_CYS_PAIR).sum()),
        "below_validator_floor": int((cys < MIN_CYSTEINES).sum()),
        "median_cysteines": float(cys.median()) if len(cys) else float("nan"),
        "distinct_fr1": int(sequences.str[:26].nunique()),
        "distinct_fr4": int(sequences.str[-11:].nunique()),
    }
    if scaffolds:
        intact = sum(any(frozen_intact(s, sc) for sc in scaffolds) for s in sequences)
        report["frozen_intact"] = int(intact)
        report["frozen_broken"] = int(len(sequences) - intact)
    report["ok"] = (report["below_validator_floor"] == 0
                    and report["below_canonical_pair"] == 0
                    and report.get("frozen_broken", 0) == 0)
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=Path, required=True, help="the design CSV")
    parser.add_argument("--column", default="designed_chain_sequence")
    parser.add_argument("--scaffold", type=Path, action="append", default=[],
                        help="a scaffold the spec froze; repeatable")
    args = parser.parse_args()

    frame = pd.read_csv(args.file)
    if args.column not in frame.columns:
        raise SystemExit(f"no column {args.column!r} in {args.file} "
                         f"(have: {list(frame.columns)[:8]})")
    scaffolds = [p.read_text().split("\n", 1)[-1].replace("\n", "").strip()
                 for p in args.scaffold] or None
    report = check(frame, args.column, scaffolds)

    print(f"{report['n']} designs lus dans {args.file}")
    print(f"  cysteines : mediane {report['median_cysteines']:.0f}")
    print(f"  a zero cysteine                    : {report['no_cysteine']}")
    print(f"  sous le plancher du validateur ({MIN_CYSTEINES}) : {report['below_validator_floor']}")
    print(f"  sans la paire canonique ({CANONICAL_CYS_PAIR})       : {report['below_canonical_pair']}")
    print(f"  variantes distinctes de FR1        : {report['distinct_fr1']} (la spec en impose une par charpente)")
    print(f"  variantes distinctes de FR4        : {report['distinct_fr4']}")
    if "frozen_broken" in report:
        print(f"  regions gelees alterees            : {report['frozen_broken']}")
    if report["ok"]:
        print("\nOK : ces designs sont des nanobodies, leurs scores veulent dire quelque chose.")
        return 0
    print("\nECHEC : ne pas lire les scores de ce run.")
    print("Un design sans son pont disulfure canonique ne se replie pas comme un")
    print("nanobody, et le validateur le refuse sur min_cysteines. Un iiptm faible")
    print("ici ne dit rien sur la methode : il dit que la chaine est cassee.")
    print("Cause connue (FINDINGS 55) : BoltzGen embarquee au lieu de l originale.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
