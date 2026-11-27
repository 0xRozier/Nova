# Où en est le projet

Passation pour le prochain agent. Mis à jour le 2026-11-27, branche
`claude/eloquent-feynman-ipdup2`.

**Sur l'état des tests, et ne pas s'y laisser prendre.** Dans un conteneur neuf
(dépôt fraîchement cloné, `requirements.txt` installé, rien d'autre) :
218 passent, 20 échouent, 15 en erreur. **Aucun de ces rouges n'est une
régression du code** : ils viennent tous de ce que `.gitignore` exclut `data/*`,
donc le cache d'archive (`data/archive/<cible>_regions.csv`) et les pièces de
nativité (`hmmscan`, `internal_data/human`) sont absents après chaque
reconstruction. Les modules touchés : `test_ptm_batch`, `test_calibration_batch`,
`test_learning_batch`, `test_regions_cache`, `test_analyze_target`,
`test_nativeness`, `test_repair_binder`, `test_parity`, `test_esm_surrogate`.
Le « 206 tests verts » des passations précédentes supposait une machine où ces
données avaient été régénérées. **Avant de conclure qu'un test est cassé, lire
son erreur : « No such file or directory: data/archive/… » n'est pas un bug.**
Et `scripts.nativeness.available()` d'abord, toujours (§6).

`docs/FINDINGS.md` est le journal complet, 91 sections, dans l'ordre
chronologique.

**À lire avant tout le reste : le §91.** L'instrument n'a jamais été étalonné
depuis la bascule de cible du 20 novembre, et l'archive publique de P05231
porte 46 118 séquences notées par les validateurs à 0,654 de médiane quand nos
327 mesures plafonnent à 0,446 avec 0,181 de médiane. Tant que la lecture 0 du
lot `inputs/kaggle_instrument.csv` n'a pas tranché, **toute comparaison de nos
valeurs absolues à celles du champ est suspendue**, et le « résultat central »
du §2 ci-dessous avec elle. Ce fichier-ci n'en est pas un résumé neutre : c'est ce qu'il faut
savoir pour ne pas refaire un travail déjà fermé par la mesure.

---

## 1. Le contexte en dix lignes

Laboratoire hors ligne (CPU) pour miner le sous-réseau Bittensor 68
(`metanova-labs/nova`), concours nanobody. Une époque dure 361 blocs (~72 min),
le gagnant prend tout, les soumissions ferment 10 blocs avant la fin.
`num_sequences: 1` — **une seule séquence par soumission**, donc rien à
diversifier à l'intérieur d'une soumission.

Une soumission franchit deux barrières puis est notée :

1. **Couche 1**, lue dans `config/config.yaml` amont : longueur 90-150, acides
   aminés valides, homopolymères ≤ 6, di-répétitions ≤ 4, pas de peptide signal,
   `enforce_vhh_hallmarks: true` (E/Q en IMGT 49, R en 50),
   `min_human_v_identity: 82`, `min_human_j_identity: 86`, gène V et longueur de
   CDR3 imposés par l'époque, unicité globale.
2. **Couche 2**, le filtre NOVA-nanobody-filter : diversité (Jaccard 5-mer < 0,9
   contre le top 50), nativité, développabilité.
3. **Notation** BoltzGen : rang-somme sur dix métriques.

**Règle du dépôt à respecter : ne pas engager de TAO pendant la recherche.** Le
dépôt a 1 époque journalisée en dry-run et **0 soumission**. Un jeton GitHub ne
doit jamais passer par la conversation (`docs/DEPLOY.md`).

---

## 2. Le résultat central, et il change tout — mais il est suspendu depuis le §91

**Le score des validateurs ne mesure pas la liaison biologique** (§77).

> **Suspendu le 27 novembre (§91).** Ce résultat repose sur l'échelle de
> calibration du §77, mesurée sur la cible vivante avec un instrument qui n'a
> jamais été étalonné contre l'archive depuis la bascule. Si cet instrument
> déprime tout à 0,18 — et 61H7, le liant cristallographique de la structure que
> la config nomme, y note 0,178, sous 92,9 % des soumissions ordinaires de
> l'archive — alors l'échelle n'a pas été mesurée, elle a été aplatie, et un
> classement dans du bruit explique le ρ = −0,232 aussi bien que « le score
> ignore la liaison ». 14 requêtes tranchent. **Ne rien bâtir sur ce qui suit
> avant de les avoir dépensées.** Ce qui reste vrai quoi qu'il arrive : la
> reproduction du score était exacte sur l'ANCIENNE cible (ρ +0,984, décalage
> +0,0012 sur 84 séquences), et les lectures internes à un run survivent.

Nous reproduisons le score exactement : leur build de BoltzGen vendorisée, leur
`design.yaml`, leur protocole `nanobody-anything`, leurs dix métriques, leur
structure, leur MSA, leur site de liaison. Les deux couches de barrières sont
reproduites aussi, la seconde en exécutant le code du filtre lui-même. **Il n'y
a aucun critère du validateur que nous ne sachions calculer.**

Mais une échelle de calibration de six liants cristallographiques d'IL-6, de
8/8 à 0/8 de recouvrement avec l'épitope de la config, donne
`Spearman(recouvrement, iiptm) = −0,232` — le mauvais signe. Le liant prouvé de
l'épitope exact (61H7, 4O9H) note 0,178 et 9 des 15 non-liants le dépassent. Un
VHH qui lie **une autre protéine** (PD-L1) est la mieux notée des ancres.

Deux conséquences à garder en tête en permanence :

- **chercher un liant réel est hors sujet.** Quatre impasses — §50 (la liaison
  ne voyage pas avec le motif), §68 (la greffe cristallographique est le pire
  bras), §72 (la conformation ne prédit pas), §73/§77 (le liant prouvé note sous
  les non-liants) — ne sont pas quatre échecs mais un seul : elles posaient
  toutes une question que l'instrument n'entend pas ;
- **le rang-somme n'est pas dix métriques mais une mesure en quatre copies.**
  `design_to_target_iptm`, `design_iiptm`, `interaction_pae` et
  `min_design_to_target_pae` corrèlent de 0,90 à 0,94 avec le rang-somme et de
  0,87 à 0,96 entre elles ; les six autres sont sous 0,30. On ne peut donc pas
  moyenner le bruit sur dix métriques.

Le problème est **une optimisation d'un objectif coûteux à évaluer** (2,2 min
la requête GPU), pas de la biologie structurale. Attention au vocabulaire : les
dix métriques sont **documentées, nommées et calculables**, et le §87 les
décompose une par une. Ce qui est opaque, c'est la carte séquence → valeur
**sans lancer BoltzGen** — et pour `design_ptm` nous savons en partie la faire.
J'ai écrit « boîte noire » dans un sens trop large et cela a retardé la
décomposition qui a trouvé le levier du §87.

---

## 3. Les chiffres à ne pas réestimer

| grandeur | valeur | où |
|---|---|---|
| bruit d'une mesure **entre runs** | 0,0121 (ancienne cible, 248 répétitions) ; 0,0133 (cible vivante, 8 répétitions) | §78, §80 |
| biais entre runs | −0,0008 : les runs sont bruités, **pas** biaisés | §80 |
| bruit **à l'intérieur** d'un run | ≤ 0,006 | §77 bis |
| écart typique entre rivaux dans une époque | 0,014 | §49 |
| effet du choix du CDR3 à charpente fixe | écart-type 0,0218 | §77 |
| effet d'une **permutation** de l'intérieur du CDR3 | écart-type 0,0600, moyenne signée −0,004 (38 permutations) | §90 |
| meilleur iiptm jamais mesuré sur la cible vivante | 0,4457 sur 311 séquences uniques, **mais non soumettable** (V 79,6 < 82) | §90 |
| meilleur iiptm **soumettable** sur la cible vivante | **0,3869** | §90 |
| part du top 24 de la cible qui passe la couche 1 | **12 sur 24** | §90 |
| transfert entre cibles, métrique d'**interface** | **ρ = +0,003 ± 0,15** | §80 |
| transfert entre cibles, `design_ptm` (repliement du monomère) | **ρ = +0,487** | §87 |
| coût d'une requête | 2,0 à 2,2 min sur 2× T4 | §77, §80 |

**Tous les ρ de validation croisée de ce dépôt sont des majorants, pas des
estimations** (§89). Deux fois le prospectif a démenti la validation croisée, et
deux fois il a donné **le même 0,145** : §83 (croisée 0,260) et §89 (croisée
0,661). Le regroupement par Jaccard empêche les quasi-jumelles de traverser un
pli, mais pas la distribution d'entraînement de différer de celle du
déploiement — notre jeu est dominé par nos propres familles de candidats. **Ne
jamais agir sur un ρ croisé sans l'avoir confirmé en prospectif.**

**Comparer entre runs est beaucoup moins fiable qu'à l'intérieur d'un run.**
C'est pourquoi `scripts/build_arms_batch.py` impose que tout ce qu'on veut
comparer tienne dans un seul run. Ne pas contourner cette règle.

---

## 4. Le substitut : état et plan

`scripts/train_score_surrogate.py`. Il prédit `design_iiptm`, d'où le rang-somme
découle à ρ −0,92. Validation croisée **par grappes** (Jaccard 5-mer ≥ 0,80,
union-find) : une charpente entière est toujours d'un seul côté du pli.

**Le chiffre opératoire est ρ ≈ 0,145 en prospectif** — pas 0,25, et pas 0,58. Le §78 annonçait 0,580 sur
l'ancienne cible ; le §83 a montré que c'était de l'**amplitude**, pas de la
prédiction. Cette cible mélangeait des conceptions de novo (mauvaises) et des
donneurs d'archive (bons), d'où un rapport signal/bruit de 15,6 ; restreinte à
la bande réaliste elle donne ρ 0,243, et la cible vivante donne 0,260. Les deux
concordent.

| | requêtes | ρ (bande réaliste) | part de l'écart hasard→oracle captée sur 30 candidats |
|---|---|---|---|
| ancienne cible | 1 490 | 0,243 | — |
| cible vivante | 273 | **0,260** | **18 %**, soit +0,55 écart-type |

Et il **ne grandit presque pas** avec les requêtes : 0,235 à 132, 0,260 à 202.
La limite n'est pas la quantité de données, c'est que nos descripteurs ne
captent pas ce que fait l'oracle à l'intérieur d'une bande étroite. La « courbe
d'apprentissage » du §80 bis (200 requêtes → 0,42) est **invalidée** par le même
artefact d'amplitude.

Un **test prospectif** a été fait (§83) et il est négatif : 25 séquences choisies
par le substitut contre 45 choisies pour la couverture, p = 0,209, et
ρ(prédit, obtenu) = +0,145 sur 70 séquences réellement nouvelles. Les deux
meilleures du lot venaient du bras couverture. **Une validation croisée à 0,25
ne se traduit pas en choix fiables.**

Structure à retenir : le substitut prédit le couple (séquence, cible). La cible
change chaque **semaine** (§80 : transfert ρ +0,003 pour l'interface ; mais
§87 : ρ +0,487 pour `design_ptm`, qui n'implique pas la cible et se mutualise),
les contraintes d'époque
chaque **époque**. Il se rachète donc par semaine.

### Le levier dominant est le GPU, pas le substitut

Comparaison sur la distribution de la cible vivante, où l'écart typique entre
rivaux (0,014) vaut 0,35 écart-type :

| | gain |
|---|---|
| passer de 30 à 100 candidats interrogés | **+1,19 écart-type** (3,4 écarts-rivaux) |
| le substitut, 1 choisi parmi 30 | +0,55 écart-type (1,6 écart-rival) |

Le GPU vaut **environ deux fois** le substitut. Le §78 affirmait l'inverse (« un
facteur 100 ») : c'est faux, et c'était fondé sur le ρ gonflé.

## 5. La cible vivante

`nova_lab/config.py` : `P05231` (interleukine-6), structure `4O9H`, chaîne `A`,
`res_index 21..186`, site `24,77,80,82,131,184..186`. Vérifié par
`scripts/analyze_target.py` (§73, §75) :

- `label_seq_id = auth + 2`, décalage unique. **La config est en numérotation
  mmcif canonique**, qui est celle que BoltzGen attend pour `res_index` ;
- la fenêtre couvre exactement les 152 résidus observés, deux boucles
  désordonnées absentes (label 51-62 et 107-108), hors épitope ;
- **8/8 du site sont au contact** du Fab co-cristallisé, 2,67 à 4,00 Å ;
- le site est l'empreinte de la **chaîne lourde** : 6 des 8 touchés par H seule,
  et les 7 résidus que seule la chaîne légère touche sont exclus. La géométrie
  demandée est donc atteignable par un domaine unique.

`python -m scripts.analyze_target` refait toute cette analyse sur n'importe
quelle cible, plus le relevé PDB des liants connus et les mutations hors CDR qui
rendraient chacun soumettable. **À lancer en premier à chaque changement de
cible** — c'est le §5 du guide amont, et c'est gratuit.

---

## 6. Pièges qui ont déjà coûté

- **Une constante survit à un changement de cible, et ça a mordu deux fois.**
  La deuxième (§90) a coûté une heure de GPU : `build_scramble_batch.py` tirait
  ses donneurs d'une liste de fichiers nommée `LIVE_TARGET_RUNS`, écrite avant
  la bascule du 20 novembre et devenue fausse ce jour-là sans que rien ne le
  dise. Dix séquences du champ choisies pour leur iiptm de 0,715-0,795 **sur
  PD-L1** notaient 0,145-0,224 sur IL-6, et 40 des 56 séquences du lot ont été
  payées sur une prémisse vide. La cause racine : **les CSV de `results/`
  n'enregistraient pas la cible**. Désormais `boltzgen_score` écrit une colonne
  `target`, `nova_lab/provenance.py` porte la table des quinze fichiers
  antérieurs (établie sur deux signatures, pas sur la date de commit, qui ne
  prouve rien dans le sens utile), et `assert_current_target()` refuse tout
  fichier dont la cible n'est pas établie. **À appeler dans tout constructeur
  de lot qui puise une séquence, un score ou une ancre dans l'historique.**
- **Une constante survit à un changement de cible.** `TARGET = "Q9NZQ7"` était
  en dur dans `boltzgen_score.py` alors que quatre autres valeurs venaient de la
  config : 6 min de GPU perdues au milieu d'un run de 2 h, parce que le nom ne
  sert qu'à construire un chemin et ne casse rien en local.
  `verify_boltzgen_notebook` compare désormais le `TARGET` de la copie
  **embarquée** au `MSA_NAME` du carnet, et `tests/test_target_switch.py` scanne
  les fichiers de la charge utile.
- **BoltzGen nomme ses entrées par l'empreinte SHA-256 de la séquence.** Une
  séquence répétée dans un lot est **une** mesure payée deux fois. Mon premier
  lot de calibration avait 3 répliques par ancre et aurait gâché un tiers de
  l'heure GPU. `verify_boltzgen_notebook` le détecte.
- **`hmmscan` absent fait passer la couche 1 pour « indisponible » sans que
  rien ne s'arrête.** Le lot du §90 a été construit sans lui : ses 56 lignes
  `layer1` disent « ANARCI did not return exactly one heavy-chain domain », et
  personne n'a su avant le 27 novembre que son record (0,4457) avait une
  identité V de 79,6 pour un plancher à 82 — donc n'était pas soumettable, ni
  lui ni aucune permutation de sa charpente. Sur le top 24 de la cible,
  **12 séquences sur 24 ne passent pas la couche 1**. Un vivier trié sur
  l'iiptm seul est à moitié inutilisable. `apt-get install -y hmmer` après
  `apt-get update` remet `regions()` et `germline_identity()` en marche ; sans
  lui, **filtrer la couche 1 avant l'iiptm est impossible, et tout choix de
  donneur est suspect**.
- **Le filtre de nativité ne lève pas quand une pièce manque** : il répond
  `hard_reject=True` pour **tout**, y compris pour des soumissions acceptées.
  Arrivé deux fois (`internal_data/human` au §69, puis `abnumber`/`anarci` au
  §73). `scripts.nativeness.available()` vérifie chaque pièce par son nom.
  **Toujours l'appeler avant de conclure quoi que ce soit d'un refus.**
- **Le conteneur est éphémère** et perd tous les paquets pip.
  `pip install --ignore-installed wheel -r requirements.txt` (le drapeau est
  nécessaire : `wheel` vient de Debian sans fichier RECORD et pip abandonne
  sans lui). Détails dans `docs/DEPLOY.md`.
- **Jaccard 5-mer ≥ 0,90 est la règle du validateur, pas un seuil
  anti-fuite.** Deux permutations d'un même CDR3 sont à 0,854 : distinctes pour
  le subnet, jumelles pour une validation croisée. Un ρ de 0,627 est tombé à
  0,580 en passant le regroupement à 0,80.
- **La recherche RCSB rend une limitation de débit sous forme de résultat
  vide**, pas d'erreur — ce qui se lit « cette protéine n'a aucun liant connu ».
  `scripts/analyze_target.py` réessaie et rapporte un vide persistant comme un
  échec réseau.
- **Un échantillon pris dans une seule bande de score atténue par restriction
  d'étendue.** Le lot de transfert stratifie sur huit bandes pour cette raison,
  et un test le vérifie.

---

## 7. Discipline de mesure en vigueur

Elle a attrapé plusieurs de mes propres erreurs ; la garder.

- **Écrire les seuils de lecture avant le run**, dans le docstring du module qui
  construit le lot, et une fonction `analyse` / `verdict` testée sur un
  aller-retour synthétique. Voir `scripts/build_transfer_batch.py` : le verdict
  ne pouvait pas être choisi après coup. Fait de même aux §72, §74, §79.
- **Un témoin dans chaque lot.** Le bras `scramble` du §74 (même charpente, même
  composition de CDR3, ordre mélangé) a donné le bruit intra-run qu'on n'avait
  jamais mesuré. Le bras `repeat` du §79 a donné le bruit de la cible vivante.
- **Vérifier la provenance des séquences avant de construire le lot.** Tout
  score historique qui a servi à *choisir* une séquence doit porter la cible
  vivante : `nova_lab.provenance.assert_current_target(paths)`. C'est gratuit,
  et son absence a coûté une heure de GPU au §90. `verify_boltzgen_notebook` ne
  le fait pas : il vérifie la charge utile, pas la question posée. Un lot peut
  être techniquement parfait et poser une question qui n'existe plus.
- **Vérifier le carnet en local avant l'upload** :
  `python -m scripts.verify_boltzgen_notebook <carnet> --batch <lot> --run <nom>`.
  Il a attrapé les deux pannes les plus chères du projet.
- **Consigner les réfutations aussi soigneusement que les résultats.** Plusieurs
  sections de FINDINGS sont des retraits (§52, §57, §58, §61, §77 bis, §78).

---

## 8. Ce qui est fermé par la mesure — ne pas rouvrir sans raison neuve

| voie | verdict | où |
|---|---|---|
| conception BoltzGen de novo | borne 0,698, sous le vivier de donneurs | §64 |
| ODesign en de novo | fermé | §65 |
| greffe d'un motif cristallographique sur notre charpente | 0,338, le pire bras | §68 |
| crible conformationnel sur CPU | tous les \|ρ\| ≤ 0,041 | §72 |
| garder la charpente du liant, redessiner le CDR3 | non soutenu, p = 0,84 | §77 |
| descripteur unique prédisant le score | tous \|ρ\| ≤ 0,19 | §77 |
| réutiliser les requêtes d'une cible sur une autre, pour l'interface | ρ = +0,003 | §80, §82 |
| une représentation de séquence plus riche (ESM-2) | égale les 3-mers partout | §86 |
| `design_ptm` comme cible de substitution | prospectif ρ +0,145, iiptm en baisse | §89 |
| « permuter le CDR3 donne une jumelle diverse **au même score** » | faux : déplacement 10× le bruit intra-run, 29/38 permutations au-delà | §90 |
| la permutation du CDR3 à 0,70-0,79 | **impossible à poser** : 0 séquence sur 311 au-dessus de 0,45 sur cette cible | §90 |

---

## 9. La suite, par ordre de rapport

1. **Fait, et négatif.** Les 70 requêtes du §81 ont été dépensées (§83) : le
   palier n'existe pas, le test prospectif est négatif, et le ρ 0,58 était de
   l'amplitude. Ne pas racheter de requêtes en espérant un meilleur substitut
   sans changer de représentation.
2. **Fermé par le §89.** Le levier `design_ptm` du §87 a été testé (lot de 60,
   dose-réponse sur quatre bandes) et il tombe au premier critère : le
   prédicteur ne généralise pas, ρ(prédit, mesuré) = **+0,145** là où il
   fallait 0,40. Le `ptm` réel est **plat** à 0,778 dans les quatre bandes
   alors que la prédiction s'étalait sur 0,058, et l'iiptm **baisse** quand le
   `ptm` prédit monte. La restriction d'étendue est écartée : le `ptm` mesuré a
   un écart-type de 0,0150, soit 4,4 fois le bruit.

3. **Le seul axe ouvert, et il est sur le générateur.** Le §90 a trouvé que la
   permutation de l'intérieur du CDR3 est une source de variance **gratuite**
   (aucune requête pour produire le candidat), réelle (écart-type 0,060 contre
   un bruit de mesure de 0,006), sans direction (moyenne signée −0,004), qui
   préserve la couche 1 et laisse le Jaccard à 0,80. Le record de la cible
   (0,4457 sur 311 séquences) est une permutation d'un donneur à 0,2551. Sur un
   objectif aussi tassé, une variance gratuite est ce que consomme une
   recherche par meilleur-de-k. **Ce qui n'est pas fait, et qui décide :
   comparer à budget GPU égal**, dans un seul run, k permutations d'un donneur
   contre k candidats neufs du générateur habituel, verdict écrit avant. Sans
   cette comparaison, on ne sait pas si la permutation vaut mieux que
   simplement tirer k candidats de plus — et c'est la seule question qui
   compte. Ne pas rebâtir le lot du §90 : son bras haut est impossible à poser
   sur cette cible.
4. **La question stratégique, qui n'est pas technique.** Le dépôt n'a jamais
   soumis. Les validateurs notent gratuitement et exactement ce qu'on leur
   envoie, là où Kaggle ne donne que des approximations payantes — mais soumettre
   engage du TAO, ce que la règle du dépôt interdit pendant la recherche.
   **C'est une décision pour le propriétaire du dépôt, pas pour un agent.**
5. **Le GPU, si la question revient.** Une époque dure 72 min, donc 2× T4
   n'évaluent qu'une trentaine de candidats dedans : avec Kaggle on ne peut pas
   *concourir*, seulement rejouer des époques passées, puisqu'on ne commande pas
   quand la session démarre. Passer de 30 à 300 candidats interrogés vaut
   +0,017 d'iiptm, soit 1,2 fois l'écart typique entre rivaux ; de 30 à 3 000,
   +0,030, soit 2,2 fois. Le rendement est logarithmique. Et il faut en plus
   ~7 h par semaine pour le substitut de la nouvelle cible.

---

## 10. Repères de code

| fichier | rôle |
|---|---|
| `docs/FINDINGS.md` | le journal, 87 sections, chronologique |
| `docs/DEPLOY.md` | remonter l'environnement, secrets |
| `nova_lab/config.py` | cible, structure, site, seuils ; `PREVIOUS_NANOBODY` garde l'ancienne |
| `nova_lab/provenance.py` | **contre quelle cible chaque fichier de `results/` a été noté** ; `assert_current_target` avant tout lot |
| `nova_lab/nanobodies.py` | couche 1 (`check_sequence`), régions, germinal, `is_vhh` |
| `nova_lab/structures.py` | lecture mmCIF, deux numérotations, contacts |
| `scripts/nativeness.py` | exécute le **vrai** filtre du validateur ; `available()` d'abord |
| `scripts/analyze_target.py` | l'analyse de cible complète, à lancer à chaque changement |
| `scripts/boltzgen_score.py` | la notation, telle que les validateurs la font |
| `scripts/build_boltzgen_notebook.py` | carnet Kaggle autonome (< 1 Mo) |
| `scripts/verify_boltzgen_notebook.py` | **à passer avant tout upload** |
| `scripts/build_arms_batch.py` | plusieurs bras dans **un seul** run |
| `scripts/train_score_surrogate.py` | le substitut, validation par grappes |
| `scripts/build_transfer_batch.py` | le lot du §80, et le motif « lecture écrite avant le run » |
| `scripts/mine.py` | le mineur : timelock, GitHub, on-chain |
| `results/` | tous les scores GPU payés, un fichier par run |

Développer et pousser **uniquement** sur `claude/eloquent-feynman-ipdup2`. Ne
pas créer de pull request sans demande explicite.
