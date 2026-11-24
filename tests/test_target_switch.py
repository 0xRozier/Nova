"""The target can change, and on 20 November 2026 it did.

Q9NZQ7 (PD-L1, structure 8AOK) became P05231 (interleukin-6, structure 4O9H).
The switch found `8AOK` hard-coded in six scripts, so the structure, its chain,
its residue range and its binding site now live in nova_lab.config and the
scripts read them. These tests keep it that way, and pin the one thing that
silently decides whether we design against the right face of the protein.

    python -m pytest tests/test_target_switch.py -q
"""
import re
from pathlib import Path

from nova_lab import config

REPO = Path(__file__).resolve().parent.parent


def test_the_configured_target_and_structure_agree():
    assert config.NANOBODY_TARGET == ["P05231"]
    assert config.NANOBODY_STRUCTURE == "4O9H"
    assert config.NANOBODY_STRUCTURE_CHAIN == "A"
    assert config.NANOBODY_STRUCTURE_RES_INDEX == "21..186"


def test_the_previous_target_is_kept_because_the_measurements_are_its():
    """Sections 1 to 66 of FINDINGS, every batch in results/, the donor pool and
    the surrogate all describe Q9NZQ7. Dropping it would orphan them."""
    assert config.PREVIOUS_NANOBODY["target"] == "Q9NZQ7"
    assert config.PREVIOUS_NANOBODY["structure"] == "8AOK"


def test_the_binding_site_is_inside_the_residue_range():
    lo, hi = (int(x) for x in config.NANOBODY_STRUCTURE_RES_INDEX.split(".."))
    for piece in config.NANOBODY_BINDING_SITE.split(","):
        bounds = [int(x) for x in piece.split("..")] if ".." in piece else [int(piece)]
        for r in bounds:
            assert lo <= r <= hi, f"residu {r} hors de {lo}..{hi}"


def test_the_scoring_and_design_scripts_read_the_config():
    """Not literals. The switch showed 8AOK hard-coded in six scripts, and a
    stale structure constant is silent: it scores against the wrong protein."""
    for name in ("scripts/boltzgen_score.py", "scripts/design_nanobody.py"):
        text = (REPO / name).read_text()
        assert "NANOBODY_STRUCTURE" in text, f"{name} ne lit pas la config"
        literals = re.findall(r'"(8AOK|4O9H)"', text)
        assert not literals, f"{name} porte encore {literals} en dur"


def test_no_script_still_hard_codes_the_old_structure():
    offenders = []
    for path in sorted((REPO / "scripts").glob("*.py")):
        text = path.read_text()
        # a mention in a docstring or comment is history, not configuration
        code = "\n".join(l.split("#")[0] for l in text.splitlines())
        if re.search(r'=\s*"8AOK"', code):
            offenders.append(path.name)
    assert not offenders, f"structure en dur dans {offenders}"


# --- FINDINGS 76: the hardcoded TARGET that survived the switch -------------
# It cost 6 GPU minutes of a 2 h Kaggle run, and nothing disagreed locally
# because the name is only ever used to build a path.

def test_the_scorer_reads_its_target_from_the_config():
    from scripts.boltzgen_score import TARGET
    assert TARGET == config.NANOBODY_TARGET[0]


def test_no_payload_file_pins_the_previous_target():
    """Only the PREVIOUS_NANOBODY record and comments may still name it."""
    from pathlib import Path

    from nova_lab import REPO_ROOT
    from scripts.build_boltzgen_notebook import OUR_FILES

    previous = config.PREVIOUS_NANOBODY["target"]
    for rel in OUR_FILES:
        for number, line in enumerate(Path(REPO_ROOT / rel).read_text().splitlines(), 1):
            code = line.split("#")[0]
            if previous in code and "PREVIOUS_NANOBODY" not in code and '"target"' not in code:
                raise AssertionError(f"{rel}:{number} fixe encore {previous} : {line.strip()}")


def test_the_msa_the_notebook_fetches_is_the_one_the_scorer_opens():
    from scripts.boltzgen_score import TARGET
    from scripts.build_boltzgen_notebook import MSA_URLS
    for url in MSA_URLS:
        assert url.format(commit="x", target=TARGET).endswith(f"{TARGET}.a3m")
