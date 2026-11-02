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
| **Mineur en mode observation** : journaliser à chaque époque ce qu'on *aurait* soumis | `python -m scripts.observe --watch` |
| Rejouer des époques passées en observation | `python -m scripts.observe --epochs 25996-26006` |
| État du journal d'observation | `python -m scripts.observe --report` |
| Lot GPU unique couvrant toutes les époques journalisées | `python -m scripts.observe --combined-batch` |
| Verdict réel d'un lot d'observation noté par BoltzGen | `python -m scripts.review_observation data/observation/batches/...` |
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

**Validé par BoltzGen** (`results/boltzgen_batch2_epoch26025.csv`, 11 séquences
notées sur 2×T4 en 21 min, soit **1,9 min/séquence** — le coût par séquence
baisse avec la taille du lot) :

| | |
|---|---|
| Fidélité : notre classement contre celui des validateurs (6 séquences de calibration) | **Spearman 1,00** |
| Nos 5 candidats, dans un panier de 11 | rangs **36 à 64** |
| Les 6 soumissions réelles de l'époque 26025 | rangs **27 à 80** |
| Notre meilleur candidat | **3ᵉ sur 11**, 9 points derrière le gagnant réel |
| Prédiction du modèle contre rang BoltzGen réel, sur ces 5 candidats | **Spearman 0,31** |

Trois conclusions nettes :

1. **L'outil de mesure est bon.** Notre BoltzGen ordonne les 6 séquences de
   calibration exactement comme les validateurs. Ce qu'il dit d'une séquence est
   fiable.
2. **Les candidats sont compétitifs, pas gagnants.** Ils battent 3 des
   6 soumissions réelles, mais pas les deux meilleures (`design_iiptm` 0,58 au
   mieux contre 0,72). Nous aurions perdu cette époque de peu.
3. **Le modèle de substitution ne sait pas classer ses propres candidats**
   (0,31 sur 5 points, donc peu concluant mais loin du 0,70 intra-lignée). Il
   amène dans la bonne zone ; il ne désigne pas le meilleur. C'est le même mode
   d'échec que pour les molécules, en moins brutal.

**Expérience du volume** (`results/boltzgen_batch3_epoch26025.csv`, 40 candidats
échantillonnés sur 5 strates de prédiction + les 6 mêmes ancres,
46 séquences en 74 min soit **1,6 min/séquence**) :

| | |
|---|---|
| Meilleur candidat sur 40 | **rang 33** |
| Meilleure soumission réelle de l'époque | rang 45 |
| Candidats battant la meilleure soumission réelle | **1 sur 40** |
| Candidats battant la soumission réelle médiane | 13 sur 40 |
| `design_iiptm` : notre meilleur / le leur | 0,683 / 0,707 |
| `interaction_pae` de notre meilleur | **10,02, le plus bas du lot** |

**Le meilleur d'un lot de 40 atteint le niveau des gagnants** : il bat les
6 soumissions réelles, dont celle qui a gagné l'époque 26025. Sa lignée est
une variante évoluée de celle de QTTECH (CDR3 `AKFMHRIKDFDY` contre
`AKFPHRYKDFDY`), ce qui valide la stratégie de mutation locale.

**Trois réserves, à ne pas perdre de vue :**

1. **Un seul tirage.** La fidélité de calibration tombe à 0,94 ici, et les deux
   meilleures séquences de QTTECH se sont inversées (rang réel 31 → notre 52,
   34 → notre 45). Notre avance de 12 points est du même ordre que ce bruit :
   « aurait gagné » est plausible, pas établi. Il faut des répliques.
2. **Un sur 40, c'est mince.** La médiane des candidats est à 119, loin derrière.
   Le pipeline ne produit pas des gagnants, il produit un gagnant *possible*
   par lot de 40 — à rejouer sur plusieurs époques pour savoir si ça tient.
3. **Le modèle de substitution ne sert pas à choisir** (Spearman 0,26 sur
   40 points ; le meilleur candidat venait de la strate médiane). Les 40 points
   mesurés sont en revanche de quoi le réentraîner sur son angle mort.

**Réplique indépendante** (`results/boltzgen_batch3_rep2_epoch26025.csv`,
même lot de 46, nouveau tirage, 67,9 min soit 1,5 min/séquence) :

| | Tirage 1 | Tirage 2 |
|---|---|---|
| Fidélité au classement des validateurs (6 ancres) | 0,94 | **1,00** |
| Meilleure soumission réelle | 45 | 45 |
| **Notre meilleur candidat** | **33** | **34** |
| Candidats devant *toutes* les soumissions réelles | 1 | **2** |
| `design_iiptm` de notre meilleur | 0,683 | 0,719 |

**Le résultat tient.** Notre meilleur candidat (CDR3 `AKFMHRIKDFDY`) reste
premier du lot aux deux tirages, devant les six soumissions réelles dont celle
qui a gagné l'époque. Au second tirage son `design_iiptm` monte à 0,719, soit
le niveau du gagnant (0,714), et un deuxième candidat passe aussi devant tout
le monde.

**L'inversion du tirage 1 était bien du bruit** : au tirage 2, notre classement
reproduit exactement celui des validateurs sur les ancres (Spearman 1,00).

**Bruit de BoltzGen mesuré** sur les 46 séquences : Spearman 0,95 entre les deux
tirages, écart moyen de 14,5 places sur la somme de rangs, 0,028 sur
`design_iiptm` (jusqu'à 0,25 dans les cas extrêmes). C'est pourquoi les
validateurs moyennent plusieurs runs ; pour choisir entre candidats proches, il
faudra répliquer.

**Le modèle de substitution reste inutile pour trier** : Spearman 0,26 aux deux
tirages, à l'identique. Ce n'est donc pas du bruit, c'est une limite réelle.

**Conséquence opérationnelle** : 40 évaluations tiennent dans une époque de
72 min sur 2 GPU modestes (1,6 min/séquence, et le coût baisse encore avec la
taille du lot). Le concours nanobody est donc techniquement jouable avec un
matériel accessible — ce qui n'était pas le cas des molécules.

### Notebook BoltzGen (`scripts.build_boltzgen_notebook`)

Le scoring utilise **la BoltzGen embarquée par les validateurs**, pas celle de
PyPI, avec leur structure cible (8AOK, chaîne A, résidus 1..115), leur site de
liaison, leur protocole `nanobody-anything` et leurs échelles (step 2.0,
noise 0.88). Le notebook embarque cette source (~300 Ko compressés), la
structure et le code ; la MSA vient de jsDelivr vérifiée par empreinte, les
poids (~6 Go) de Hugging Face. Les noyaux cuEquivariance se désactivent d'eux-
mêmes sur T4 (capacité 7.5 < 8). Les 103 dépendances épinglées se résolvent.

**Exécuté sur Kaggle (2×T4) : 3,7 min/séquence.** Trois obstacles
d'environnement levés au passage : `MPLBACKEND` imposé par Kaggle, l'import de
`bittensor` par la BoltzGen de NOVA (substitut dans `kaggle/bittensor.py`), et
le `TOP` qui ne prenait que les premières lignes du lot.

**Fidélité mesurée** (`results/boltzgen_batch1_calibration.csv`, 2 séquences
dont on connaît les métriques publiées) :

| Métriques | Écart à la valeur publiée |
|---|---|
| `liability_score`, `liability_num_violations` | exact |
| confiance (`design_ptm`, `design_to_target_iptm`, `interaction_pae`, PAE min) | 0,004 à 0,18 |
| `design_iiptm` | ±0,02 à 0,05 |
| `plip_hbonds_refolded` | −4 sur les deux séquences |
| `plip_saltbridge_refolded`, `delta_sasa_refolded` | 0 à −3 ; −65 à +13 |

Les métriques d'interaction se calculent sur la **structure prédite**, et
l'étape de conception de BoltzGen est stochastique (aucune graine dans la
config du subnet). Notre run est donc un tirage de la même distribution,
là où les valeurs publiées sont une moyenne sur plusieurs validateurs — ce qui
explique que les métriques robustes (confiance, développabilité) collent et que
les métriques de comptage (PLIP, SASA) s'écartent. Ce sont précisément les
trois que le modèle de substitution prédit le moins bien. Pour classer des
candidats sérieusement, il faudra répliquer ou accepter ce bruit.

## Économie du subnet (`results/economics.txt`, mesuré le 1ᵉʳ nov. 2026)

Chiffres on-chain : 1 alpha émis par bloc, répartition 41 % mineurs / 41 %
validateurs / 18 % propriétaire (confirmée par le métagraphe, mineurs et
validateurs exactement à égalité), alpha à 0,0166 TAO, enregistrement à
0,026 TAO.

| | Par époque (72 min) |
|---|---|
| Pool des mineurs | 148 alpha |
| **Gagnant nanobody (40 %)** | **59,2 alpha = 0,985 TAO** |
| Gagnant molécules (60 %) | 88,8 alpha = 1,477 TAO |

Revenu du concours nanobody selon le taux de victoire, et coût GPU (40 évaluations
BoltzGen par époque demandent ~64 min de 2×T4, donc 2 GPU en continu) :

| Taux de victoire | TAO/jour | Coût GPU 2×RTX4090 | Coût GPU 2×A100 |
|---|---|---|---|
| 10 % | 2,0 | 0,06 TAO/jour | 0,19 TAO/jour |
| 30 % | 5,9 | idem | idem |
| 71 % (niveau de frankji) | 13,9 | idem | idem |

**Le seuil de rentabilité est à 0,3 % de victoires** (2×RTX4090) ou **1 %**
(2×A100) : une victoire tous les trois jours suffit à payer le calcul. Le coût
GPU est négligeable devant la récompense, quel que soit le prix du TAO.

**Trois réserves sérieuses :**

1. **La valeur est au prix de marché, pas garantie.** Le subnet émet
   7 184 alpha/jour mais ne reçoit que 16,5 TAO/jour. Le prix tient parce que
   les détenteurs conservent leur alpha, pas grâce à un flux entrant. Vendre une
   victoire (59 alpha) ne déplace le prix que de 0,0025 %, mais l'ensemble des
   mineurs vendant chaque jour représenterait 0,3 % du pool par jour.
2. **Le prix du TAO est incertain** : mes recherches donnent 250–260 $ fin
   octobre 2026, sans cotation fiable du jour. À vérifier avant de décider.
3. **Les versements passent par une API** (`payout.enabled`, UID 42) : les
   validateurs mettent tout le poids sur un UID relais, et Metanova reverse
   ensuite l'alpha à la coldkey du gagnant. Cela ajoute une dépendance à un
   tiers que nous n'avons pas auditée.

## Quelles tâches le subnet tire-t-il ? (`scripts.epoch_tasks`)

Le couple (gène V, longueur CDR3) sort du hash de bloc via `random.Random`, donc
sa loi découle de l'algorithme : une catégorie parmi quatre, puis un gène dans
la catégorie, puis une longueur uniforme dans 5..23. Vérifié sur 200 000 tirages
simulés.

| Gène V imposé | Fréquence | Charpentes conformes en archive |
|---|---|---|
| IGHV3-23 | 25 % | 392 |
| IGHV3-66 | 25 % | 298 |
| IGHV3-74 | 25 % | 177 |
| IGHV3-53 | 4,2 % | 77 |
| IGHV3-7 | 4,2 % | 50 |
| IGHV3-11 | 4,2 % | 52 |
| IGHV3-48 | 4,2 % | 44 |
| IGHV3-64 | 4,2 % | 23 |
| **IGHV3-43** | **4,2 %** | **4** |

Trois gènes occupent chacun une catégorie entière et sortent donc un quart du
temps ; les six de la catégorie V4 se partagent le dernier quart. La longueur de
CDR3 est uniforme (5,26 % par valeur).

**Nous sommes équipés pour 95,8 % des époques** (au moins 10 charpentes et
5 CDR3 de la bonne longueur). Le seul trou est **IGHV3-43**, tiré une époque sur
24, pour lequel l'archive ne contient que 4 séquences conformes. La fréquence
réelle observée sur mes 56 époques colle à la théorie (IGHV3-23 18 fois,
IGHV3-43 une seule).

Table complète par combinaison : `data/submissions/Q9NZQ7_epoch_coverage.csv`.

### Test sur le pire cas (époque 25989, IGHV3-43, CDR3 de 15)

4 charpentes seulement, combinaison marquée « non équipée ». Résultat :
**1 658 candidats valides sur 4 000** (41 %, contre 81 % sur une époque bien
couverte). **Le nombre de charpentes n'est donc pas la contrainte** : la
variation vient du CDR3, et 4 charpentes suffisent à en produire des milliers.
Mon seuil de 10 charpentes était mal choisi.

Deux différences réelles apparaissent en revanche :

- **Plus de collisions de CDR3** : 1 110 rejets pour CDR3 ≥ 85 % identique à
  l'archive, parce que seuls 79 CDR3 de longueur 15 existent et que les muter
  les ramène souvent sur un voisin connu.
- **Aucune marge germinale** : tous les candidats sortent à 82,8–83,8 %
  d'identité V, pour un seuil à 82. Les charpentes IGHV3-43 de l'archive sont
  elles-mêmes à la limite, donc la moindre mutation fait basculer.

**Verdict BoltzGen** (`results/boltzgen_batch4_epoch25989.csv`, 24 séquences
en 40,6 min) :

| | |
|---|---|
| Gagnant réel (frankji, rang vrai 15) | rang 27, `design_iiptm` 0,561 |
| **Notre meilleur candidat** | **rang 25**, `design_iiptm` 0,588 |
| Candidats devant *toutes* les soumissions réelles | 1 sur 20 |
| Candidats devant la médiane réelle | 14 sur 20 |

**C'est une égalité, pas une victoire.** Notre candidat passe devant, mais de
**2 places**, alors que le bruit de BoltzGen mesuré sur deux tirages est de
14,5 places en moyenne. Sur l'époque 26025 l'avance était de 11 à 12 places et
se reproduisait sur deux tirages ; ici nous n'avons qu'un tirage et un écart
dans le bruit. La conclusion honnête : **à parité avec le gagnant sur une époque
mal couverte**, ce qui reste un bon résultat mais n'est pas une victoire
démontrée.

Trois observations utiles :

- **La tâche était dure pour tout le monde** : frankji gagne cette époque avec
  0,561 de `design_iiptm`, contre 0,714 pour le gagnant de l'époque 26025. Le
  niveau requis baisse avec la difficulté, donc s'abstenir serait une erreur.
- **Diversité de charpente très faible** : 18 de nos 20 candidats partagent la
  même charpente (3 distinctes sur 4 disponibles). Si cette charpente avait été
  mauvaise, tout le lot l'aurait été. C'est le vrai risque de ces époques, et il
  ne se voit pas dans le rendement de validité.
- **Aucune marge germinale** : 82,8–83,8 % d'identité V pour un seuil à 82.
- Le modèle de substitution tombe à Spearman 0,13 ici (contre 0,26), confirmant
  qu'il ne sert pas à trier.

## Mineur en mode observation (`scripts.observe`)

Avant d'engager du TAO, le mineur tourne à vide. `scripts.observe` ne touche
**ni portefeuille, ni chaîne, ni dépôt GitHub** : il lit les paramètres de
l'époque, fait tourner le générateur sur CPU, et écrit dans
`data/observation/journal.jsonl` ce qu'il aurait soumis (séquence, empreinte,
gène V imposé, rendement de validité, causes de rejet). Le champ `submitted`
vaut toujours `false` ; le chiffrement, la publication GitHub et l'ancrage sur
la chaîne ne sont volontairement pas dans ce fichier.

Chaque époque produit aussi un lot GPU (`data/observation/batches/`) qui
**entrelace nos candidats et les soumissions réelles de cette époque**. C'est le
point clé : BoltzGen est stochastique, et l'écart mesuré entre deux tirages des
mêmes séquences atteint une quinzaine de places sur 46 — des métriques venant de
deux runs différents ne sont pas comparables. Il faut donc noter les deux camps
dans le même run. Pour une époque en cours, les concurrents sont encore
chiffrés : le lot part avec nos seuls candidats, et
`scripts.observe --rebuild-batches` y réinjecte les rivaux une fois
`scripts.fetch_submissions` et `scripts.analyze_nanobodies` passés.

La cible ne change jamais (Q9NZQ7 / 8AOK), donc 21 runs Kaggle séparés seraient
du gâchis : `--combined-batch` réunit toutes les époques journalisées en un seul
lot, avec une colonne `epoch` que le scorer conserve. Les lignes y sont rangées
par profondeur et non par époque — premier tour : le meilleur candidat et le
meilleur rival de chaque époque ; deuxième tour : les deuxièmes. Un run
interrompu donne donc un verdict peu profond sur **toutes** les époques plutôt
qu'un verdict complet sur les premières. À profondeur 2, les 21 époques tiennent
en 84 séquences, soit ~5 h sur 2×T4 : une seule session Kaggle.

`scripts.review_observation` prend les scores BoltzGen d'un lot, applique le
classement du validateur (rang dense par métrique, somme des rangs) sur le pool
réellement noté, et dit à quelle place notre candidat serait arrivé. Il classe
**chaque époque séparément** quand le lot en couvre plusieurs, comme le fait un
validateur. `--summary` agrège les époques déjà jugées.

### Marche à suivre sur Kaggle

1. `python -m scripts.observe --combined-batch` puis
   `python -m scripts.build_boltzgen_notebook --batch data/observation/batches/Q9NZQ7_combined.csv --run observation_combined`
2. Importer `notebooks/kaggle_boltzgen.ipynb` dans Kaggle, régler *Accelerator*
   sur **GPU T4 x2** et *Internet* sur **On**.
3. Test rapide dans l'éditeur : mettre `TOP = 4` et **Run All** (~15 min, vérifie
   que les poids se téléchargent et que BoltzGen tourne).
4. Lot complet : remettre `TOP = None`, puis **Save Version → Save & Run All
   (Commit)**. Compter ~5 h ; le quota Kaggle est de 30 h de GPU par semaine et
   12 h par session, donc le lot passe en une fois.
5. Récupérer `scores_nanobody_observation_combined.csv` dans les *Output* du run,
   le placer dans `results/`, puis
   `python -m scripts.review_observation results/scores_nanobody_observation_combined.csv`
   et `python -m scripts.review_observation --summary`.

Pour approfondir ensuite : `--depth 4` donne 168 séquences (~10 h, deux
sessions) et un `RUN` différent est **obligatoire** à chaque lot, BoltzGen
reprenant sinon les résultats du répertoire existant au lieu de recalculer.

Ce que cette boucle ne peut pas dire :

- la présélection dans le lot vient du modèle de substitution, qui ne trie pas
  (Spearman ~0,2) : notre meilleur candidat n'est probablement pas dans les
  4 envoyés au GPU ;
- seuls les rivaux portés dans le lot sont dans le pool, et seulement ceux que
  le validateur a lui-même acceptés ;
- la développabilité structurale reste invisible hors ligne, donc un candidat
  « valide » peut encore être refusé.

## Limites

- Le modèle de substitution n'est pas Boltz-2 et ne généralise pas aux molécules nouvelles (voir ci-dessus).
- Le backtest ne choisit que parmi des molécules déjà soumises (et déjà
  plutôt bonnes) : c'est une indication, pas une preuve.
- Nanobodies : hors des lots déjà notés, la sortie du générateur est une
  prédiction, pas un score. Le mode observation existe pour transformer cela en
  taux de victoire mesuré, mais il faut d'abord noter les lots sur GPU.
- La développabilité structurale reste hors de portée hors ligne (prédiction de structure requise).
