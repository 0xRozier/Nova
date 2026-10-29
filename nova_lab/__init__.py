"""Offline research tooling for NOVA (Bittensor SN68).

The validator rules are reimplemented here from the upstream source
(metanova-labs/nova, see UPSTREAM_COMMIT) so they can run without the GPU stack.
Re-check them against upstream whenever the subnet updates.
"""
from pathlib import Path

UPSTREAM_COMMIT = "db9bf190a6907338a05871b9658f42a089229675"

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
DB_PATH = DATA_DIR / "molecules.sqlite"
DB_URL = "https://huggingface.co/datasets/Metanova/Mol-Rxn-DB/resolve/main/molecules.sqlite"
