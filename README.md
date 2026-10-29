# NOVA (SN68) — laboratoire hors ligne

Outils de recherche **sans GPU** pour préparer une participation au subnet 68
de Bittensor ([metanova-labs/nova](https://github.com/metanova-labs/nova)) avant
d'engager du matériel ou des TAO.

## Idée

Le concours de petites molécules (60 % des récompenses, un seul gagnant par
époque de 361 blocs) note 20 molécules avec Boltz-2 :
`score = (affinity_probability_binary − affinity_pred_value) / nb_atomes_lourds`.

Les validateurs publient tous leurs scores Boltz-2 sur une API publique
(score sharing). On s'en sert pour entraîner un **modèle de substitution sur
CPU**, qui sert à explorer la base combinatoire et à proposer des candidats.
Boltz-2 sur GPU ne sert ensuite qu'à confirmer une courte liste.

## Installation

```bash
./setup.sh
source .venv/bin/activate
```

## Utilisation

| Étape | Commande |
|---|---|
| Paramètres de l'époque en cours (réaction imposée, gène V, longueur CDR3) | `python -m scripts.epoch_info` |
| Télécharger tous les scores publiés (~280k molécules, ~1 h) | `python -m scripts.fetch_scores` |
| Statistiques : niveau des gagnants, réactions favorables | `python -m scripts.analyze_scores` |
| Déchiffrer les soumissions passées de chaque mineur (nœud archive + drand) | `python -m scripts.fetch_submissions --last 50` |
| Scores exacts des gagnants par époque, qui gagne | `python -m scripts.analyze_winners` |
| Construire le jeu d'entraînement | `python -m scripts.build_dataset` |
| Entraîner et backtester le modèle | `python -m scripts.train_surrogate` |
| Proposer 20 molécules pour une réaction | `python -m scripts.propose_molecules --reaction rxn:41` |
| Noter un lot avec Boltz-2 comme les validateurs (GPU ; CPU trop lent : >1 h 30 par molécule) | `python -m scripts.boltz_score --file inputs/kaggle_batch.csv --accelerator gpu` |
| Idem sur GPU Kaggle gratuit (recherche uniquement) | importer `notebooks/kaggle_boltz.ipynb` dans Kaggle |
| Vérifier la parité avec le validateur | `python -m pytest tests -q` |

## Fidélité au validateur

`nova_lab/` réimplémente les règles du validateur (résolution des noms
`rxn:R:A:B`, filtres RDKit, compatibilité Boltz, entropie, nouveauté par
rapport à l'archive) à partir du code amont au commit indiqué dans
`nova_lab/__init__.py`. `tests/test_parity.py` vérifie que les SMILES calculés
sont identiques à ceux de l'archive publique des validateurs. **À revérifier à
chaque mise à jour du subnet** (`nova_lab/config.py` reprend `config.yaml`).

## Limites

- Le modèle de substitution n'est pas Boltz-2. Il trie, il ne remplace pas.
- Le backtest ne choisit que parmi des molécules déjà soumises (et déjà
  plutôt bonnes) : c'est une indication, pas une preuve.
- Le concours nanobodies (40 %) n'est pas encore couvert.
