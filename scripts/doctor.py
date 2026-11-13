"""Is this machine measuring what the validator measures? Run it before anything.

Every failure this project has had was a tool that answered instead of
refusing. The list is short and it is all the same shape:

  - no `hmmscan`, so ANARCI returned None for every sequence, so the training
    set was built on 0 of 9 466 rows and said "saved" (fixed 2026-11-13);
  - DSSP absent or 4.x, so TNP printed one error line and reported a PSH taken
    on another code path -- 89.50 where the right path gives 91.21 -- and
    eleven production submissions were refused on a gate we thought we were
    screening for;
  - ImmuneBuilder's set/dict typo, so 3 sequences in 37 came back with no
    model and therefore no PSH, which reads like "nothing to refuse";
  - `data/scores/` missing, so step 2 of the deployment died on an exception.

So this does not ask whether the binaries exist. It runs them on a known
nanobody and checks the numbers come out, which is the only question that
matters: a tool that is present and wrong is worse than one that is absent.

Usage:
    python -m scripts.doctor             # the install
    python -m scripts.doctor --quick     # skip the end-to-end fold (no TNP run)

Exits 1 if anything the miner depends on is wrong, so it can gate a run.
"""
import argparse
import shutil
import sys
import tempfile
from pathlib import Path

from nova_lab import DATA_DIR, toolchain

# A clinical-stage nanobody (caplacizumab VHH), so the numbers below are not
# ours to choose: any faithful install reproduces them.
PROBE = ("EVQLVESGGGLVQPGGSLRLSCAASGRTFSYNPMGWFRQAPGKGRELVAAISRTGGSTYYPDSVEGRFT"
         "ISRDNAKRMVYLQMNSLRAEDTAVYYCAAAGVRAEDGRVRTLPSEYTFWGQGTQVTVSS")


def line(ok: bool | None, what: str, detail: str) -> None:
    mark = "ok  " if ok else ("-- " if ok is None else "FAIL")
    print(f"  {mark} {what:<26} {detail}")


def check_python_side() -> list[bool]:
    """The CPU-only half: ANARCI, the archive, the dataset, the model."""
    results = []
    print("\nnumbering and germlines (ANARCI + hmmscan)")
    hmmscan = shutil.which("hmmscan")
    line(bool(hmmscan), "hmmscan", hmmscan or
         "absent; ANARCI then returns None for every sequence WITHOUT raising")
    results.append(bool(hmmscan))

    try:
        from nova_lab.nanobodies import germline_identity, regions
        r, g = regions(PROBE), germline_identity(PROBE)
        ok = bool(r) and bool(g) and g["V"]["identity_pct"] > 50
        line(ok, "IMGT numbering",
             f"CDR3 {r['CDR3']}, V {g['V']['identity_pct']:.1f}% {g['V']['genes']}"
             if ok else "returned nothing on a known nanobody")
    except Exception as error:  # noqa: BLE001 - the report is the point
        ok = False
        line(False, "IMGT numbering", f"{type(error).__name__}: {error}")
    results.append(ok)

    print("\nthe donor pool")
    from scripts.bootstrap import state
    from nova_lab.config import NANOBODY_TARGET
    for name, path, _ in state(NANOBODY_TARGET[0]):
        here = path.exists()
        if here and "dataset" in name:
            matrix = DATA_DIR / "dataset" / f"{NANOBODY_TARGET[0]}_nanobodies_X.npy"
            here = matrix.exists()
        line(here, name, f"{path.stat().st_size / 1e6:.1f} MB" if path.exists()
             else "missing; python -m scripts.bootstrap --build")
        results.append(here)
    return results


def check_structural(quick: bool) -> list[bool]:
    """The gate the validator applies before it scores anything."""
    print("\nthe developability toolchain")
    rows = toolchain.check()
    for what, ok, detail in rows:
        line(ok, what, detail)
    results = [ok for _, ok, _ in rows]
    if quick or not all(results):
        if not quick:
            print("     (skipping the end-to-end fold: fix the above first)")
        return results

    print("\nend to end, on a clinical-stage nanobody")
    from scripts.aggrescan_gate import screen as a3d_screen
    from scripts.tnp_developability import parse, run
    work = Path(tempfile.mkdtemp(prefix="nova-doctor-"))
    try:
        fasta = work / "probe.fasta"
        fasta.write_text(f">probe\n{PROBE}\n")
        df = parse(run(fasta, 4, work), work)
        got = df[df["name"] == "probe"]
        if got.empty:
            line(False, "TNP", "no row came back for the probe")
            return results + [False]
        row = got.iloc[0]
        dssp = bool(row.get("dssp", False))
        line(dssp, "DSSP annotation",
             "present, so PSH is on the validator's path" if dssp else
             "ABSENT: this PSH is not the validator's measurement")
        line(True, "TNP metrics",
             f"PSH {row['PSH']:.2f}, C {row['C']:.3f}, L {row['L']:.0f}, "
             f"L3 {row['L3']:.0f}, {int(row['reds'])} red")
        models = work / "TNP_output" / "Final_Models"
        verdict, values = a3d_screen({"probe": PROBE}, models, work / "a3d",
                                     with_values=True)
        a3d_ok = bool(values)
        if a3d_ok:
            avg, peak = values[PROBE]
            line(True, "Aggrescan3D", f"average {avg:.4f}, peak {peak:.4f}")
        else:
            line(False, "Aggrescan3D", "produced nothing on a model TNP built")
        return results + [dssp, a3d_ok]
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true",
                        help="check the tools are right without folding a structure")
    args = parser.parse_args()

    print("NOVA miner: is this machine measuring what the validator measures?")
    results = check_python_side() + check_structural(args.quick)
    bad = results.count(False)
    print(f"\n{len(results) - bad}/{len(results)} checks pass")
    if bad:
        print("Do not start the miner: a tool that is present and wrong reports "
              "numbers that are not the validator's, which is how the first "
              "production run lost twelve epochs.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
