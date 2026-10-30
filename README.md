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
| Idem pour les nanobodies (classement par somme de rangs BoltzGen) | `python -m scripts.analyze_nanobodies` |
| Construire le jeu d'entraînement | `python -m scripts.build_dataset` |
| Entraîner et backtester le modèle | `python -m scripts.train_surrogate` |
| Proposer 20 molécules pour une réaction | `python -m scripts.propose_molecules --reaction rxn:41` |
| Noter un lot avec Boltz-2 comme les validateurs (GPU ; CPU trop lent : >1 h 30 par molécule) | `python -m scripts.boltz_score --file inputs/kaggle_batch.csv --accelerator gpu` |
| Idem sur GPU Kaggle gratuit (recherche uniquement) : notebook autonome à importer dans Kaggle | `python -m scripts.build_kaggle_notebook` puis importer `notebooks/kaggle_boltz.ipynb` |
| **Nanobodies** : vérifier des séquences candidates avant de soumettre (CPU) | `python -m scripts.check_nanobody --file candidats.txt` |
| **Nanobodies** : jeu d'entraînement depuis l'archive (~10 min) | `python -m scripts.build_nanobody_dataset` |
| **Nanobodies** : entraîner le modèle et le backtester | `python -m scripts.train_nanobody_surrogate` |
| **Nanobodies** : proposer des séquences pour l'époque | `python -m scripts.propose_nanobody` |
| **Nanobodies** : noter un lot avec la BoltzGen des validateurs (GPU) | `python -m scripts.boltzgen_score --file inputs/kaggle_nanobody_batch.csv` |
| Idem sur GPU Kaggle : notebook autonome (nécessite le clone amont) | `python -m scripts.build_boltzgen_notebook` puis importer `notebooks/kaggle_boltzgen.ipynb` |
| Vérifier la parité avec le validateur | `python -m pytest tests -q` |

## Fidélité au validateur

`nova_lab/` réimplémente les règles du validateur (résolution des noms
`rxn:R:A:B`, filtres RDKit, compatibilité Boltz, entropie, nouveauté par
rapport à l'archive) à partir du code amont au commit indiqué dans
`nova_lab/__init__.py`. `tests/test_parity.py` vérifie que les SMILES calculés
sont identiques à ceux de l'archive publique des validateurs. **À revérifier à
chaque mise à jour du subnet** (`nova_lab/config.py` reprend `config.yaml`).

## Résultats Boltz-2 (Kaggle, `results/`)

`kaggle_batch1` (T4, 32 bits, 4,2 min/molécule) :
- **Calibration** : corrélation 0,976 (Spearman 0,988), écart moyen 0,012 avec les
  scores publiés → notre Boltz reproduit fidèlement celui des validateurs.
- **Proposition rxn:33 du générateur** : total réel **1,13** pour 2,93 prédit ;
  moyenne 0,05 par molécule, sous la moyenne des molécules soumises (0,10).
  Spearman prédiction/réel ≈ 0,08 sur les 40 candidats : hors de la zone des
  données d'entraînement, le modèle de substitution ne vaut rien. Or la règle de
  nouveauté (Tanimoto < 0,5 avec l'archive) oblige justement à sortir de cette zone.

## Stratégies observées (époques 25956–26025)

**Molécules** (`scripts.analyze_winners`) : `zeromaj/nova` gagne 41 époques sur 51.
Score moyen par molécule 0,132 contre 0,074 pour les autres, à taille égale
(~24 atomes lourds) et à distance de l'archive égale (Tanimoto max ~0,35).
Ni astuce de taille ni stratégie au bord de la limite de nouveauté : une bien
meilleure sélection, vraisemblablement beaucoup d'évaluations Boltz par époque.

**Nanobodies** (`scripts.analyze_nanobodies`) : `frankji-groundcontrol` gagne
30 époques sur 56 (71 % des époques où il soumet), avec un seul hotkey et
100 % de soumissions valides, alors que **52 % de toutes les soumissions sont
invalides**. Il fait évoluer une même lignée (identité ~0,87 d'une victoire à la
suivante, CDR2 conservé, CDR3 réadapté à la longueur imposée). Les gagnants se
distinguent surtout sur les métriques de confiance (design_iiptm 0,48 contre 0,26,
PAE d'interaction 14 contre 20). La règle de nouveauté est bien plus souple que pour
les molécules (rejet à ≥ 95 % d'identité globale ou ≥ 85 % sur le CDR3) : on peut
optimiser localement autour des meilleures séquences connues, là où un modèle
appris sur les 9 000 séquences notées de l'archive reste dans sa zone de validité.

## Vérificateur de validité nanobody (`nova_lab/nanobodies.py`)

52 % des soumissions nanobodies sont rejetées par les validateurs. Le
vérificateur reproduit hors ligne, sur CPU, tous les contrôles qui ne demandent
pas de prédiction de structure : jeton nul, longueur, acides aminés, homopolymères,
di-répétitions, peptide signal, identité germinale humaine V/J (ANARCI),
gène V imposé par l'époque, classification VHH avec les résidus marqueurs IMGT
(H49 ∈ {E,Q}, H50 = R), cystéines canoniques H23/H104, longueur de CDR3 imposée,
et les signaux rouges de développabilité calculables sur la séquence
(longueur totale des CDR, longueur du CDR3).

**Non couvert** : les contrôles de développabilité structuraux (TNP : surface
hydrophobe, plaques de charge, compacité du CDR3 ; Aggrescan3D), qui exigent
une prédiction de structure. Une séquence acceptée ici est donc « valide
sauf rejet de développabilité ».

**Nuance sur la similarité à l'archive** : nous comparons à toute l'archive,
alors que le validateur ne compare qu'aux candidats retenus par son pré-filtre
par k-mers. Nous sommes donc plus stricts que lui ; cette similarité est
rapportée comme avertissement, pas comme rejet.

Dépendances supplémentaires : `abnumber`, `anarci` et le binaire `hmmscan`
(paquet HMMER). `tests/test_nanobody_parity.py` confronte le vérificateur aux
verdicts réels des validateurs sur les soumissions déchiffrées.

## Générateur de nanobodies

**Le modèle a du signal, contrairement au cas des molécules.** Un modèle
LightGBM par métrique, sur des descripteurs par région IMGT (composition,
hydropathie, charge, identité germinale) :

| Évaluation | Spearman sur la somme de rangs |
|---|---|
| Aléatoire (optimiste) | 0,94 |
| **Au sein d'une lignée** (régime du générateur) | **0,70** |
| Temporelle stricte, dans l'époque | 0,55 |

En test temporel, choisir par prédiction aurait désigné le gagnant réel de
**8 époques sur 14** (place médiane 1 sur 22 ; au hasard : place 10,5).
**Nuance capitale** : le modèle choisit la séquence de frankji dans 10 époques
sur 10 où il soumet. Il reconnaît donc d'abord une bonne lignée. Les 0,70
au sein d'une lignée disent qu'il discrimine aussi à l'intérieur, ce dont le
générateur a besoin, mais c'est plus faible.

**Le générateur** (`scripts.propose_nanobody`) part de charpentes d'archive
appelant le gène V de l'époque, y greffe des CDR3 de la longueur imposée issus
des meilleures séquences, mute les CDR, puis filtre et classe. Deux
enseignements tirés des mesures :

- **L'archive est trompeuse.** 6 227 des 8 578 séquences non datées sont sous
  les seuils germinaux actuels (antérieures à la règle) et scorent mieux que les
  conformes. Les charpentes doivent venir des ~1 000 séquences conformes.
- **La fin du CDR3 est codée par le gène J.** La muter faisait échouer 285
  candidats sur 800 ; en la protégeant, le rendement passe de 48 % à 79 %.

Les quasi-doublons (≥ 95 % d'identité globale, ou CDR3 ≥ 85 %) sont écartés
comme rejets certains ; les règles plus subtiles restent des avertissements.

**Non validé** : aucune séquence proposée n'a encore été notée par BoltzGen.
Le classement affiché est une prédiction. `inputs/kaggle_nanobody_batch.csv`
contient un lot prêt à noter (candidats + séquences de calibration dont on
connaît le classement réel), et `notebooks/kaggle_boltzgen.ipynb` le notebook
pour le faire tourner.

### Notebook BoltzGen (`scripts.build_boltzgen_notebook`)

Le scoring utilise **la BoltzGen embarquée par les validateurs**, pas celle de
PyPI, avec leur structure cible (8AOK, chaîne A, résidus 1..115), leur site de
liaison, leur protocole `nanobody-anything` et leurs échelles (step 2.0,
noise 0.88). Le notebook embarque cette source (~300 Ko compressés), la
structure et le code ; la MSA vient de jsDelivr vérifiée par empreinte, les
poids (~6 Go) de Hugging Face. Les noyaux cuEquivariance se désactivent d'eux-
mêmes sur T4 (capacité 7.5 < 8). Les 103 dépendances épinglées se résolvent.

**Jamais exécuté de bout en bout** : la durée par séquence et le comportement
réel sur T4 restent à mesurer (cet environnement n'a ni GPU ni les ~6 Go
d'espace nécessaires).

## Limites

- Le modèle de substitution n'est pas Boltz-2 et ne généralise pas aux molécules nouvelles (voir ci-dessus).
- Le backtest ne choisit que parmi des molécules déjà soumises (et déjà
  plutôt bonnes) : c'est une indication, pas une preuve.
- Nanobodies : le générateur n'est pas encore confronté à BoltzGen ; sa sortie est une prédiction, pas un score.
- La développabilité structurale reste hors de portée hors ligne (prédiction de structure requise).
