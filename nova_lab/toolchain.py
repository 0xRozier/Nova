"""Where the structural tools live, and whether they are the ones the validator runs.

The developability gate is the one place where *our* number has to be the
validator's number, because PSH > 89 is a refusal and a refusal costs the
epoch. That makes the install itself part of the parity surface, and it had
been treated as a detail.

What the validator runs, read from its own install (upstream `install_deps.sh`,
`NOVA-nanobody-filter/install/install.sh` and `requirements-minimal.txt`):

  - TNP at commit 29dcac72f1380e8538e8870f45a699d3c6156162,
  - with `NOVA-nanobody-filter/install/patches/TNP/process_pdb.py` copied over
    `site-packages/scripts/process_pdb.py`,
  - ImmuneBuilder at 0df4e2ad82a1aa60f37ea9dae335d1198159ef78 (NanoBodyBuilder2),
  - DSSP built from cmbi/dssp, which is the autotools 3.1.4 tree,
  - aggrescan3d-rs and the FreeSASA 2.0.1 pinned beside it, from the filter's
    `install/install-aggrescan3d.sh`.

DSSP is the one that bites. TNP asks Biopython for a DSSP annotation of the
NanoBodyBuilder2 model; mkdssp 4.x refuses that model (no PDB HEADER line, so
it assumes mmCIF and fails) while 3.1.4 accepts it once the filter's patch pads
the lines to 80 characters. TNP does not stop on that failure: it prints

    ERROR: Failed to process output PDB. DSSP failed to produce an output

and goes on to report a PSH computed without the annotation. Measured on one
candidate, the same structure reads **89.50 without DSSP and 91.21 with it**,
and its positive-charge flag moves from AMBER to GREEN. A reading taken on the
wrong path is not a noisy version of the validator's number, it is a different
number, so this module makes that case loud instead of silent.

Paths default to where this project installs them and can be moved with
NOVA_TNP, NOVA_DSSP and NOVA_NANOFILTER, which is what a rented box needs.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

TNP_COMMIT = "29dcac72f1380e8538e8870f45a699d3c6156162"
DSSP_VERSION = "3.1.4"  # cmbi/dssp, the tree the filter builds

TNP = Path(os.environ.get("NOVA_TNP", "/home/user/venv-tnp/bin/TNP"))
NANOFILTER = Path(os.environ.get("NOVA_NANOFILTER", "/home/user/nanofilter"))
A3D_BIN_DIR = NANOFILTER / ".tools" / "bin"


def dssp_bin() -> Path | None:
    """The mkdssp TNP will find, in the order it will find it."""
    if explicit := os.environ.get("NOVA_DSSP"):
        path = Path(explicit)
        return path if path.exists() else None
    for candidate in (TNP.parent / "mkdssp", TNP.parent / "dssp"):
        if candidate.exists():
            return candidate
    for name in ("mkdssp", "dssp"):
        if found := shutil.which(name):
            return Path(found)
    return None


def dssp_version(path: Path | None = None) -> str | None:
    """The version string mkdssp reports, or None if it cannot be asked."""
    path = path or dssp_bin()
    if path is None:
        return None
    for flag in ("--version", "-h"):
        try:
            proc = subprocess.run([str(path), flag], capture_output=True,
                                  text=True, timeout=30)
        except (OSError, subprocess.SubprocessError):
            return None
        text = f"{proc.stdout}\n{proc.stderr}"
        for token in text.split():
            if token[:1].isdigit() and token.count(".") >= 2:
                return token
    return None


def tnp_env() -> dict:
    """The environment TNP needs: NanoBodyBuilder2 and mkdssp on PATH.

    TNP shells out to a bare `NanoBodyBuilder2`, and Biopython shells out to a
    bare `mkdssp`. Both are resolved from PATH, so a run that forgets either
    one still produces numbers -- the wrong ones, in the DSSP case.
    """
    parts = [str(TNP.parent)]
    if (dssp := dssp_bin()) is not None and str(dssp.parent) not in parts:
        parts.append(str(dssp.parent))
    parts.append(os.environ.get("PATH", ""))
    env = {**os.environ, "PATH": ":".join(p for p in parts if p)}

    # TNP keeps its own modules in a top-level package called `scripts`, and so
    # do we. Any PYTHONPATH entry pointing at this repository therefore shadows
    # TNP's package, and TNP dies on
    #     ModuleNotFoundError: No module named 'scripts.region_definitions'
    # having profiled nothing -- 100 sequences in, 0 rows out. It is an easy
    # PYTHONPATH to set (any script importing nova_lab from outside the repo
    # needs it), so the collision is removed here rather than warned about.
    repo = Path(__file__).resolve().parent.parent
    kept = [p for p in env.get("PYTHONPATH", "").split(os.pathsep)
            if p and not _is_repo(p, repo)]
    if kept:
        env["PYTHONPATH"] = os.pathsep.join(kept)
    else:
        env.pop("PYTHONPATH", None)
    return env


def _is_repo(entry: str, repo: Path) -> bool:
    try:
        return Path(entry).resolve() == repo
    except OSError:
        return False


def annotated(out_dir: Path, name: str) -> bool:
    """Did TNP get its DSSP annotation for this sequence?

    `<name>_annotations.csv` is written from the DSSP result and is absent
    exactly when the annotation step failed, which makes it a per-sequence
    verdict on whether PSH was measured the validator's way. The stdout error
    line says the same thing but cannot be attributed to a sequence once TNP
    runs several in parallel.
    """
    raw = out_dir / "Raw_Model_Outputs" / name
    return (raw / f"{name}_annotations.csv").exists()


def check() -> list[tuple[str, bool, str]]:
    """(what, is it right, what it says) for every tool the gate needs."""
    rows = []

    rows.append(("TNP", TNP.exists(),
                 str(TNP) if TNP.exists() else f"absent from {TNP}"))
    builder = TNP.parent / "NanoBodyBuilder2"
    rows.append(("NanoBodyBuilder2", builder.exists(),
                 str(builder) if builder.exists()
                 else f"absent from {builder}; TNP calls it by bare name"))

    site = next((p for p in TNP.parent.parent.glob("lib/python3.*/site-packages")
                 if (p / "scripts" / "process_pdb.py").exists()), None)
    patched = False
    detail = "TNP is not installed, so its process_pdb.py cannot be read"
    if TNP.exists():
        if site is None:
            detail = "scripts/process_pdb.py not found beside TNP"
        else:
            text = (site / "scripts" / "process_pdb.py").read_text()
            patched = "fix_pdb_for_dssp" in text
            detail = ("the filter's patch is in place" if patched else
                      "UNPATCHED: copy NOVA-nanobody-filter/install/patches/"
                      "TNP/process_pdb.py over it, or DSSP will fail")
    rows.append(("TNP process_pdb patch", patched, detail))

    # ImmuneBuilder 1.2 passes a set where OpenMM wants a dict, so the
    # strained-sidechain path dies with "Wrong number or type of arguments for
    # overloaded function 'new_Context'". It is not every sequence: 3 of the 37
    # calibration sequences hit it, and each one came back with no model and so
    # no PSH at all.
    refine = site / "ImmuneBuilder" / "refine.py" if site else None
    if refine is not None and refine.exists():
        bad = "{'Threads', str(n_threads)}" in refine.read_text()
        rows.append(("ImmuneBuilder refine.py", not bad,
                     "the set/dict typo is fixed" if not bad else
                     f"line 346 of {refine} passes a SET to OpenMM; replace "
                     f"{{'Threads', str(n_threads)}} with {{'Threads': str(n_threads)}} "
                     f"or 8% of sequences silently get no model"))
    else:
        rows.append(("ImmuneBuilder refine.py", False,
                     "ImmuneBuilder not found beside TNP"))

    dssp = dssp_bin()
    version = dssp_version(dssp)
    major = (version or "").split(".")[0]
    rows.append(("mkdssp", dssp is not None and major == "3",
                 f"{dssp} reports {version}" if dssp is not None else
                 "absent from PATH and from TNP's directory"
                 if version is None else
                 f"{dssp} reports {version}; 4.x refuses NanoBodyBuilder2 models, "
                 f"build {DSSP_VERSION} from cmbi/dssp"))

    for name in ("aggrescan3d-rs", "freesasa"):
        binary = A3D_BIN_DIR / name
        rows.append((name, binary.exists(),
                     str(binary) if binary.exists()
                     else f"absent from {binary}; run the filter's "
                          f"install/install-aggrescan3d.sh"))
    return rows
