# SN68 Miner — Déploiement sur Lium (2×RTX 4090)

Guide de déploiement du mineur nanobody sur une machine GPU louée via Lium.

## Pré-requis

- Compte Lium avec solde (20$ minimum pour ~22h)
- CLI Lium configuré (`lium init`)
- Wallet Bittensor avec un hotkey enregistré sur le subnet 68
- Repo de soumissions GitHub (privé, initialisé avec au moins un commit)
- GitHub Personal Access Token (fine-grained) avec accès Contents R/W sur le repo de soumissions

## 1. Provisionner la machine

```bash
# Lister les 2×RTX4090 disponibles, triées par prix
lium ls --gpu RTX4090 --count 2 --sort price_gpu

# Lancer le pod (22h de TTL = budget ~20$)
lium up <index> --name sn68-miner --ttl 22h -y

# Vérifier
lium ps
lium exec sn68-miner "nvidia-smi --query-gpu=name,memory.total --format=csv"
```

## 2. Cloner les repos

```bash
lium exec sn68-miner "cd /root && git clone --recurse-submodules https://github.com/metanova-labs/nova.git nova-upstream"
lium exec sn68-miner "cd /root && git clone https://github.com/0xRozier/Nova.git nova-lab"
```

## 3. Installer les dépendances (upstream)

```bash
lium exec sn68-miner "cd /root/nova-upstream && bash install_deps.sh"
```

Cela installe : uv, Rust, CUDA PyTorch, Boltz-2, BoltzGen, timelock (WASM),
RDKit, bittensor, les poids NanoBodyBuilder2 et Boltz-2, Aggrescan3D, HMMER,
MMseqs2, IgBLAST, et la base combinatoire (~120 MB).

Vérification :
```bash
lium exec sn68-miner "cd /root/nova-upstream && source .venv/bin/activate && python3 -c 'import torch; print(torch.cuda.device_count(), torch.cuda.get_device_name(0))'"
# Doit afficher: 2 NVIDIA GeForce RTX 4090
```

## 4. Installer les dépendances (lab)

```bash
lium exec sn68-miner "source /root/nova-upstream/.venv/bin/activate && cd /root/nova-lab && uv pip install -r requirements.txt"
```

## 5. Transférer le hotkey (pas le coldkey)

```bash
lium exec sn68-miner "mkdir -p /root/.bittensor/wallets/Main/hotkeys"
# Copier le contenu du fichier hotkey local :
cat ~/.bittensor/wallets/Main/hotkeys/sub68
# Puis sur le pod :
lium exec sn68-miner "cat > /root/.bittensor/wallets/Main/hotkeys/sub68 << 'EOF'
{...contenu du hotkey...}
EOF"
```

## 6. Configurer le .env

```bash
lium exec sn68-miner "cat > /root/nova-upstream/.env << 'ENV_EOF'
SUBTENSOR_NETWORK=\"wss://entrypoint-finney.opentensor.ai:443\"
BT_NO_PARSE_CLI_ARGS=false
NOVA_SKIP_HISTORICAL_CHECKS=0

# GPU settings pour 2×RTX4090
NOVA_BOLTZ_SHARDS=4
NOVA_BOLTZ_GPUS=0,1,0,1
NOVA_BOLTZGEN_GPUS=0,1
NOVA_BOLTZ_THREADS=1
NOVA_MAX_CONCURRENT_TNP=8

# GitHub
GITHUB_TOKEN=\"<ton_token>\"
GITHUB_REPO_NAME=\"sn68-submissions\"
GITHUB_REPO_BRANCH=\"main\"
GITHUB_REPO_OWNER=\"0xRozier\"
GITHUB_REPO_PATH=\"\"

# Wallet
BT_WALLET_NAME=\"Main\"
BT_WALLET_HOTKEY=\"sub68\"
AUTO_UPDATE=\"0\"
ENV_EOF"
```

## 7. Liens et données

```bash
# Lier la base combinatoire pour le lab
lium exec sn68-miner "mkdir -p /root/nova-lab/data && ln -sf /root/nova-upstream/combinatorial_db/molecules.sqlite /root/nova-lab/data/molecules.sqlite"
```

## 8. Construire le dataset et le modèle

```bash
# Dataset nanobody (~5 min)
lium exec sn68-miner "cd /root/nova-lab && source /root/nova-upstream/.venv/bin/activate && export \$(grep -v '^#' /root/nova-upstream/.env | xargs) && python3 -m scripts.build_nanobody_dataset"

# Surrogate model (~45 min sur 2×RTX4090 — CPU-bound)
lium exec sn68-miner "cd /root/nova-lab && source /root/nova-upstream/.venv/bin/activate && export \$(grep -v '^#' /root/nova-upstream/.env | xargs) && python3 -m scripts.train_nanobody_surrogate"
```

Note : si le dataset n'a pas d'epochs associées (pas de `fetch_submissions`
préalable), le script ne fait que le random-split et sauvegarde le modèle.
Le surrogate sert uniquement de pré-filtre (Spearman ~0.2) ; BoltzGen fait
le vrai classement.

## 9. Créer le script de lancement

```bash
lium exec sn68-miner "cat > /root/run_mine.sh << 'SCRIPT'
#!/bin/bash
cd /root/nova-lab
source /root/nova-upstream/.venv/bin/activate
set -a; source /root/nova-upstream/.env; set +a
export PYTHONUNBUFFERED=1
python3 -m scripts.mine --wallet Main --hotkey sub68 --nova-dir /root/nova-upstream --scorer surrogate \"\$@\" >> /root/mine.log 2>&1
SCRIPT
chmod +x /root/run_mine.sh"
```

## 10. Dry run

```bash
lium exec sn68-miner "bash /root/run_mine.sh --once"
# Vérifier le log :
lium exec sn68-miner "cat /root/mine.log"
# Doit se terminer par "dry run: nothing committed on chain"
```

## 11. Lancement live

```bash
lium exec sn68-miner "bash /root/run_mine.sh --live"
# Surveiller :
lium exec sn68-miner "tail -20 /root/mine.log"
lium exec sn68-miner "cat /root/nova-lab/data/mining/journal.jsonl"
```

## Coûts

| Configuration | $/h | $/jour | Budget 20$ |
|---|---|---|---|
| 2×RTX 4090 | ~0.90 | ~21.6 | ~22h |

Seuil de rentabilité : 0.3% de taux de victoire (~1 victoire tous les 3 jours).

## Dépannage

- **`FileNotFoundError: ...nanobodies.parquet`** : lancer `scripts.build_nanobody_dataset`
- **`FileNotFoundError: ...surrogate.joblib`** : lancer `scripts.train_nanobody_surrogate`
- **`bt.wallet` / `bt.subtensor` AttributeError** : bittensor v10 utilise `bt.Wallet` / `bt.Subtensor` (capitalisé)
- **lium exec & background** : `lium exec` ne supporte pas `&` ; utiliser un wrapper script avec `nohup` ou `tmux`
- **Config drift** : `python3 -m scripts.check_config` avant chaque session GPU
