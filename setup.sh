#!/usr/bin/env bash
# Offline research environment: CPU only, no GPU needed.
set -euo pipefail
cd "$(dirname "$0")"

python3 -m venv .venv
.venv/bin/pip install -q --upgrade pip
.venv/bin/pip install -q -r requirements.txt

mkdir -p data
if [ ! -f data/molecules.sqlite ]; then
  echo "Downloading the combinatorial database (~120 MB)..."
  curl -sSL -o data/molecules.sqlite \
    "https://huggingface.co/datasets/Metanova/Mol-Rxn-DB/resolve/main/molecules.sqlite"
fi

echo "Done. Activate with: source .venv/bin/activate"
