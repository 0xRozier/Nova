# Où en est le projet

Passation pour le prochain agent. Écrit le 2026-11-25, branche
`claude/eloquent-feynman-ipdup2`, commit `e15e1de`, 206 tests verts (hors
`tests/test_nanobody_parity.py`, qui dure 21 min).

`docs/FINDINGS.md` est le journal complet, 87 sections, dans l'ordre
chronologique. Ce fichier-ci n'en est pas un résumé neutre : c'est ce qu'il faut
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

## 2. Le résultat central, et il change tout

**Le score des validateurs ne mesure pas la liaison biologique** (§77).

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

Le problème est **une optimisation en boîte noire d'un objectif coûteux**
(2,2 min la requête GPU), pas de la biologie structurale.

---

## 3. Les chiffres à ne pas réestimer

| grandeur | valeur | où |
|---|---|---|
| bruit d'une mesure **entre runs** | 0,0121 (ancienne cible, 248 répétitions) ; 0,0133 (cible vivante, 8 répétitions) | §78, §80 |
| biais entre runs | −0,0008 : les runs sont bruités, **pas** biaisés | §80 |
| bruit **à l'intérieur** d'un run | ≤ 0,006 | §77 bis |
| écart typique entre rivaux dans une époque | 0,014 | §49 |
| effet du choix du CDR3 à charpente fixe | écart-type 0,0218 | §77 |
| transfert entre cibles | **ρ = +0,003 ± 0,15** | §80 |
| coût d'une requête | 2,0 à 2,2 min sur 2× T4 | §77, §80 |

**Comparer entre runs est beaucoup moins fiable qu'à l'intérieur d'un run.**
C'est pourquoi `scripts/build_arms_batch.py` impose que tout ce qu'on veut
comparer tienne dans un seul run. Ne pas contourner cette règle.

---

## 4. Le substitut : état et plan

`scripts/train_score_surrogate.py`. Il prédit `design_iiptm`, d'où le rang-somme
découle à ρ −0,92. Validation croisée **par grappes** (Jaccard 5-mer ≥ 0,80,
union-find) : une charpente entière est toujours d'un seul côté du pli.

| | requêtes | ρ | part de l'écart hasard→oracle captée sur 30 candidats |
|---|---|---|---|
| ancienne cible `8AOK` | 1 490 | +0,580 | 73 % (et 55 % dans le haut du vivier) |
| **cible vivante `4O9H`** | **141** | **+0,235 ± 0,14** | 20 % |

**Le transfert entre cibles est nul** (§80), donc les 1 490 requêtes de
l'ancienne cible sont de l'histoire, pas un actif. La méthode reste validée : le
ρ 0,580 dit qu'un substitut *peut* atteindre ce niveau à partir d'un millier de
requêtes **sur une même cible**.

Courbe d'apprentissage, mesurée en sous-échantillonnant par grappes (§80 bis) :

| requêtes | heures GPU | ρ |
|---|---|---|
| 100 | 3,7 h | +0,27 ± 0,22 |
| **200** | **7,3 h** | **+0,42 ± 0,10** |
| 700 | 25,7 h | +0,47 ± 0,06 |
| 1 000 | 36,7 h | +0,54 ± 0,01 |

**200 requêtes donnent les trois quarts de ce que 1 157 atteignent.** Nous en
avons 141 sur la cible vivante ; environ 70 de plus (2,5 h) atteignent le palier.

Structure à retenir : le substitut prédit le couple (séquence, cible). La cible
change chaque **semaine**, les contraintes d'époque (gène V, longueur de CDR3)
chaque **époque**. Le substitut se rachète donc **par semaine** et sert à toutes
les époques de la semaine. Le guide amont prévoit ~24 h de préavis sur la
nouvelle cible : c'est la fenêtre pour le construire.

---

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
| réutiliser les requêtes d'une cible sur une autre | ρ = +0,003 | §80 |

---

## 9. La suite, par ordre de rapport

1. **Compléter le substitut de la cible vivante — ~70 requêtes, 2,5 h de GPU.**
   Les choisir pour l'information apportée (étaler l'espace que le modèle
   connaît mal), pas au hasard. C'est le palier des 200 et ρ ≈ 0,42.
2. **Brancher la génération sur ce substitut — gratuit.** `propose_nanobody`
   propose, un substitut obsolète classe. Le remplacer par celui du §78,
   réentraîné sur la cible vivante, et générer beaucoup pour n'interroger que
   le haut du classement.
3. **La question stratégique, qui n'est pas technique.** Le dépôt n'a jamais
   soumis. Les validateurs notent gratuitement et exactement ce qu'on leur
   envoie, là où Kaggle ne donne que des approximations payantes — mais soumettre
   engage du TAO, ce que la règle du dépôt interdit pendant la recherche.
   **C'est une décision pour le propriétaire du dépôt, pas pour un agent.**
4. **Le GPU, si la question revient.** Une époque dure 72 min, donc 2× T4
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
