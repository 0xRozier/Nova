# Audit du rejeu : toutes les fuites de la même classe que la coupure d'archive

Audit conduit le 15 novembre 2026, après le correctif `7a90b69` (« la coupure
ignorait les lignes non datées »). La question posée n'est pas « le correctif
est-il bon » mais « **où ailleurs** le rejeu voit-il le futur, et où mesure-t-il
autre chose que ce que le mineur fait ».

Les époques rejouées par les §48, §49 et §51 sont **26299-26306**. Le seul
instantané daté de l'archive publique est `data/prebuilt/Q9NZQ7_regions.csv`,
figé à 9466 séquences, soit tout ce qui précède **l'époque 26025**. Dans tout ce
document, « du futur » veut dire « hors de cet instantané et non daté, donc
impossible à situer avant l'époque rejouée » — et l'erreur va dans un seul sens :
une ligne hors instantané est arrivée *après* 26025, jamais avant.

Aucun code de production n'a été modifié. Les mesures ci-dessous ont été faites
avec des scripts jetables (hors dépôt) et les données présentes sur disque au
15 novembre 14 h UTC.

## Tableau récapitulatif, du plus gonflant au moins

| # | Fuite / écart | Où | Verdict | Taille mesurée |
|---|---|---|---|---|
| 1 | `scripts/humanise.py` garde les lignes non datées : **le trou exact du correctif, non corrigé**. Tout le bras `humanised` du §47/§49 est bâti sur des donneurs qui n'existaient pas. | `humanise.py:120` | **fuite confirmée** | **117 donneurs → 0** sur les 8 époques. Le pic du bras (0,694 à 26301, 92 % du meilleur rival) vient de donneurs post-instantané : honnêtement, le bras n'a **aucun** donneur à ces époques. |
| 2 | Le substitut `design_iiptm` qui choisit les candidats du mineur est entraîné **à 100 % sur des étiquettes postérieures** à toute époque rejouée. | `train_iiptm_surrogate.py:155-160`, `mine.py:244` | **fuite confirmée** | 438/438 lignes d'entraînement : 348 du champ (époques 26299-26348), 208 de nos propres candidats notés (époques 25996-26025). Zéro ligne antérieure à 26299. |
| 3 | La « place » journalisée est calculée contre un vivier de repli qui n'a rien d'un champ : `dated` est vide sous toute coupure ≤ 26306, et `compliant.head(200)` prend silencieusement le relais. | `propose_nanobody.py:429-430` | **fuite confirmée** | `surrogate_place_vs_reference` = **1 à 24 observations sur 24**, contre **0 victoire sur 8** au GPU. Médiane `design_iiptm` du vivier de repli : **0,205** contre 0,710 pour le champ. Même les vraies soumissions du champ y arrivent 1ʳᵉ : la métrique ne discrimine rien. |
| 4 | Le substitut 10-métriques qui pré-filtre le vivier du rejeu est entraîné sur **tout le jeu de données**, à dessein. | `train_nanobody_surrogate.py:93`, `propose_nanobody.py:428` | **fuite confirmée** | ≥ 1070 des 10536 lignes d'entraînement sont postérieures à 26025 et 152 sont datées 26306-26348 ; aucune des 9466 autres n'est datable plus finement. Le modèle honnête n'existe pas sur disque. |
| 5 | Le rejeu et le mineur ne donnent **pas le même vivier** au même sélecteur : `observe` appelle `propose` sans `skip_surrogate`, `mine` avec. | `observe.py:227`, `mine.py:718-724` | **écart confirmé** | Le rejeu remet au sélecteur les `k` **meilleurs** du substitut (fuite n° 4), le mineur les `k` **premiers** dans l'ordre de génération. Recouvrement attendu `k²/V` : **0,1 à 5,5 séquences sur 13 à 66** — les deux viviers sont quasi disjoints. |
| 6 | `scripts/psh_by_generator.py` porte le même trou `isna()`. | `psh_by_generator.py:53` | **fuite confirmée** | Même 197 lignes non datées admises ; l'A/B des générateurs (`results/psh_generator_ab_cdr3_18.csv`) a été mesuré sur un vivier qui voyait le futur. |
| 7 | `PSH_BIAS = 2.0` : une de ses trois estimations est un seuil de Youden ajusté sur **les 296 verdicts des époques rejouées**. | `developability.py:52-53` | **fuite confirmée, faible** | Sans cette estimation le biais vaut 1,6 à 3,1 au lieu de 2,0. Sur les 8 PSH retenues du journal, P(refus) moyenne passe de 15,9 % à 14,2 % (biais 1,6) ou 18,0 % (biais 2,5) : **≤ 4 points**, et l'ordre des candidats est inchangé à qualité égale. |
| 8 | Paramètres par défaut divergents entre rejeu et mineur. | `observe.py:300,303` vs `mine.py:866,869` | **écart confirmé, mineur** | `--candidates` 3000 (rejeu) contre 5000 (mineur) ; `--shortlist` 10 contre 8. Le premier affaiblit le rejeu, le second lui donne deux tirages GPU de plus. |
| 9 | **La porte de nouveauté du rejeu est plus permissive que celle du mineur** : `propose` passe `warnings=[]`, `mine.archive_screen` ne passe rien, et `check_sequence` fait de tout coup non-`certain` une rejection quand `warnings is None`. | `nanobodies.py:287`, `propose_nanobody.py:438` vs `mine.py:385` | **écart confirmé, dans le sens qui gonfle** | Un coup à 0,90 ≤ identité < 0,95 sans donnée CDR est **gardé avec un avertissement** par le rejeu et **refusé** par le mineur : **3 des 20 candidats échantillonnés (15 %)** survivent dans le rejeu et seraient jetés par le mineur. |
| 9b | La porte de nouveauté juge contre l'archive **d'aujourd'hui** : `archive_exclude` n'est passé par aucun chemin de production. | `propose_nanobody.py:393,412,438`, `mine.py:383` | **pas de fuite** (anachronisme pessimiste) | 0 des 80 candidats journalisés atteint 0,85 de similarité CDR3 contre une ligne post-instantané (max **0,833**), et 0 des 20 échantillonnés atteint 0,95 d'identité entière (max **0,948**, donc le zéro tient à peu). |
| 10 | Le cache de régions est empoisonné : 6740 des 17244 lignes de l'archive sont marquées définitivement « non numérotables » alors qu'elles le sont. | `nanobodies.py:312-343` | **pas une fuite** — défaut qui rétrécit le vivier, mais il **vide la prémisse du §51** | Les 196 soumissions de 26299-26305 sont absentes du vivier bien qu'elles soient dans l'archive et numérotables. Le « vivier conforme de 256 » du §51 en vaut 75 aujourd'hui, et **0** sont dans l'instantané. |
| 11 | `scripts/design_nanobody.py` n'a aucune option de coupure. | `design_nanobody.py:86-95` | **fuite par construction** | Tout rejeu passant par ce générateur greffe depuis le vivier entier. Non utilisé par `observe`/`mine` : dette, pas dégât mesuré. |
| 12 | `review_observation.py` choisit les rivaux du lot par leur **rang réel** (vision du futur) — mais prend les **4 meilleurs**. | `observe.py:87,99`, `review_observation.py:84-86` | **pas de fuite** | La sélection par hindsight durcit le pool au lieu de l'adoucir : battre les 4 meilleurs d'une époque implique la gagner. Le script avertit déjà qu'il ne voit que les soumissions acceptées. |
| 13 | `PSH_PASS_RATE`, `SIGMA = 2.7`, `ACCEPTED_MEDIAN = 83.5`, `WANT_UNDER_THRESHOLD`, les σ Aggrescan3D. | `mine.py:84`, `psh_calibration.py:95-100`, `preflight.py:66-67` | **pas de fuite** | Calibrés sur l'époque **26174** (36 acceptées + 1 refusée), sur l'époque **26238**, et sur 100 candidats frais par longueur — tous antérieurs à 26299. Aucune verdict d'époque rejouée n'y entre. |
| 14 | La coupure porte sur l'époque de **soumission**, alors que les métriques qui classent les donneurs sont publiées plus tard. | `observe.py:216` | **fuite potentielle, nulle aujourd'hui** | Aucune ligne datée n'existe entre 26025 et 26299 : déplacer la coupure change **0 ligne** sur les huit rejeux publiés. Le délai de publication de l'API de scores n'a pas pu être mesuré hors réseau. |
| 15 | Les paramètres d'époque (gène V, longueur CDR3). | `check_nanobody.py:23-34`, `epoch_info.py:30-53` | **pas de fuite** | Dérivés du hash du **premier bloc de l'époque rejouée**, lu sur le nœud d'archive. `len(reaction_ids())` est bien d'aujourd'hui, mais il n'intervient qu'avant que `rng` soit touché, donc il ne déplace ni le gène ni la longueur (`tests/test_challenge_params.py`). |

---

## 1. Le modèle de substitution

Il y en a **deux**, et les deux fuient.

### 1a. `Q9NZQ7_iiptm_surrogate.joblib` — la moitié « qualité » du choix du mineur

C'est celui que `scripts.mine.surrogate_rank` charge, et c'est lui qui décide
quels candidats sont mesurés puis soumis : `select_by_expected_value` multiplie
sa prédiction par P(non-refus). `scripts.observe` appelle **la même fonction**
(`observe.py:225`), donc le rejeu hérite du modèle tel quel.

Son jeu d'entraînement (`train_iiptm_surrogate.labelled`) est la réunion de :

- `data/submissions/Q9NZQ7_nanobody_results.csv` — **348 lignes**, et leurs
  époques sont 26299, 26300, 26301, 26302, 26303, 26304, 26305, 26306, 26345,
  26346, 26347, 26348. **Toutes** sont donc à l'époque rejouée ou après ;
- `results/scores_nanobody_observation_*.csv` — **208 lignes**, qui sont *nos
  propres candidats* des époques 25996-26025, avec leur `design_iiptm` réel
  mesuré par BoltzGen.

Après featurisation il reste `n_train = 438` (la valeur stockée dans le
bundle). **Aucune de ces 438 lignes n'est antérieure à 26299.** Pour un rejeu de
26299-26306, 100 % des étiquettes sont dans le futur. Et pour un rejeu de
25996-26025 — ce que le §46 a fait — les 208 lignes « nos candidats » sont une
mémorisation littérale : le modèle a vu la note réelle des candidats que le
rejeu régénère.

Le commentaire du code assume le choix (« il n'y a pas de raison de retenir des
données du modèle livré », `train_iiptm_surrogate.py:155`). C'est défendable
pour le mineur en production. C'est indéfendable dans un rejeu, et rien dans
`observe.py` ne le signale.

Vérification que ce n'est pas une mémorisation directe des candidats
d'aujourd'hui : des 80 candidats des lots `data/observation/batches/
Q9NZQ7_epoch2630*.csv`, **0** est dans les 438 lignes d'entraînement. La fuite
est donc distributionnelle, pas une recopie — pour ces huit époques.

Mesure de l'écart, sur les 63 candidats à nous dont BoltzGen a mesuré le
`design_iiptm` à ces mêmes époques (`results/scores_nanobody_arms_26299_26306.csv`
et `..._repaired_...csv`), aucun présent dans l'un ou l'autre jeu
d'entraînement :

| modèle | lignes d'entraînement | ρ de Spearman contre le vrai `design_iiptm` | son top-1 lit | son top-8 lit (médiane) |
|---|---|---|---|---|
| **celui du dépôt** (fuite totale) | 556 (champ 26299-26348 + nos candidats notés) | **0,244** (p = 0,017, n = 95) | **0,179** | 0,543 |
| **honnête** (échantillon de 3000 de l'instantané ≤ 26025) | 3000 | **0,418** (p = 2,5·10⁻⁵, n = 95) | **0,581** | 0,535 |
| le vivier lui-même | — | — | médiane 0,270, max 0,696 | 8 au hasard : 0,312 |

**Le résultat va dans le sens inverse de l'intuition, et il faut le dire tel
quel : la fuite ne rapporte rien.** Le modèle honnête ordonne *mieux* nos
candidats (ρ = 0,418 contre 0,244) et son premier choix lit un vrai
`design_iiptm` de 0,581 là où celui du dépôt lit 0,179 — pire que 8 candidats
tirés au hasard (0,312). L'explication plausible est que le champ de novembre
(médiane 0,71) ne ressemble pas à nos greffes (médiane 0,27), donc apprendre sur
le champ transfère mal, tandis que l'archive ancienne couvre la même gamme que
nos candidats.

**Réserve honnête sur cette mesure** : les deux modèles n'ont pas le même nombre
de lignes (3000 contre 556), donc une part de l'écart peut venir de la taille
d'échantillon et non de l'honnêteté. Ce que la mesure établit solidement est le
**signe** : la fuite n'achète aucune performance mesurable sur cette tâche, donc
**aucun chiffre publié n'est gonflé par elle**. Ce qui reste cassé est autre
chose : le rejeu fait tourner le sélecteur du mineur avec un modèle qui a vu le
futur, donc il ne mesure pas le sélecteur que le mineur aura à l'époque suivante.
Les 95 candidats du test incluent les 31 `humanised`, dont la liaison est très
basse, ce qui élargit la gamme et facilite les deux modèles également.

### 1b. `Q9NZQ7_nanobody_surrogate.joblib` — le substitut 10 métriques

`propose_nanobody.propose` le charge (`Surrogate.load`) pour ordonner les
candidats et calculer la `place`. `train_nanobody_surrogate.py:93` l'entraîne
sur `np.arange(len(df))` — **toutes** les lignes du parquet — et le commentaire
juste au-dessus explique pourquoi c'est volontaire. Le modèle du split temporel
honnête est calculé juste après, affiché, puis **jeté** : il n'est jamais
sauvegardé, et la dernière ligne du script le dit (« le modèle du backtest n'est
pas celui sur le disque »).

Datation du jeu d'entraînement (parquet du 15 novembre, 10536 lignes) :

| tranche | lignes |
|---|---|
| dans l'instantané daté (≤ 26025) | 9 466 |
| hors instantané, donc postérieures à 26025 | **1 070** |
| — dont datées 26306, 26345-26348 | 152 |
| — dont non datables | 918 |

Soit **au moins 1070 lignes sur 10536 (10,2 %) à l'époque rejouée ou après**, et
le reste n'est pas datable plus finement que « ≤ 26025 ». Les 152 datées sont
toutes à 26306 ou plus tard : pour un rejeu de 26299 à 26305 ce sont purement et
simplement des soumissions futures, et elles portent les meilleures métriques du
jeu.

L'effet n'est pas seulement sur la `place` (voir §4) : comme `observe` appelle
`propose` **sans** `skip_surrogate`, ce modèle choisit aussi *quels* candidats
arrivent au sélecteur du mineur (voir §5).

**Ce qu'il faudrait pour trancher proprement** : un modèle entraîné par époque
rejouée sur les seules lignes antérieures. Le script sait le faire
(`train_mask = (df["epoch"] < split) | df["epoch"].isna()` — notons au passage
que ce masque porte lui aussi le `isna()` du correctif, donc même le « split
strict » du §1 entraîne sur des lignes non datées du futur). Il ne l'enregistre
pas, et `Surrogate.load` ne prend pas de coupure.

### Verdict §1

**Fuite confirmée, sur les deux modèles.** Pour le substitut `design_iiptm` qui
décide des soumissions, 438 étiquettes sur 438 sont postérieures à toute époque
rejouée : la moitié « qualité » du sélecteur du mineur n'a, dans un rejeu,
**aucune** information d'avant l'époque. Pour le substitut 10-métriques,
≥ 1070/10536 lignes sont postérieures à l'instantané et 152 sont datées après
l'époque rejouée. Le `train_mask` du split « strict » porte lui-même le bug
`isna()` qui vient d'être corrigé ailleurs.

---

## 2. L'écran PSH / développabilité

### Ce qui a été vérifié, constante par constante

| constante | source déclarée | époque de la source | postérieure à 26299 ? |
|---|---|---|---|
| `psh_calibration.THRESHOLD = 89` | code du validateur (`nanobody_validity.py`) | — | non, ce n'est pas une mesure |
| `psh_calibration.SIGMA = 2.7` | 37 séquences × 3 plis, `results/psh_calibration_epoch26174.csv` + nos 20 candidats v8/v9 | **26174** et nos lots | **non** |
| `psh_calibration.ACCEPTED_MEDIAN = 83.5`, `ACCEPTED_MAX = 88.9` | 36 soumissions acceptées, 3 plis | **26174** | **non** |
| `developability.BOUNDS`, `BANDS = 2.0` | `metanano/config.py` en amont | — | non, ce n'est pas une mesure |
| `developability.SIGMA` (C, PPC, PNC) | 37 séquences × 3 plis | **26174** | **non** |
| `preflight.A3D_SIGMA_AVERAGE/PEAK` | 37 séquences × 3 plis | **26174** | **non** |
| `preflight` jeu de calibration | `results/tnp_epoch26174.csv` | **26174** | **non** |
| `mine.PSH_PASS_RATE = {8: .75, 11: .50, 18: .15}` | 100 **candidats frais à nous** par longueur, 1 pli, `results/psh_fresh_candidates_by_length.csv` | nos candidats, pas le champ | **non** |
| `mine.WANT_UNDER_THRESHOLD = 10` | « la développabilité a refusé 1 sur 8 » | **26238** | **non** |
| **`developability.PSH_BIAS = 2.0`** | trois estimations, dont « le seuil qui sépare le mieux leurs 296 verdicts sur notre échelle (**époques 26299-26306**, Youden) » | **26299-26306** | **OUI** |

`PSH_PASS_RATE` mérite d'être dit explicitement parce que c'est le piège
attendu : il n'est pas calibré sur des verdicts de validateur du tout, mais sur
le taux de nos propres candidats sous 89. Il ne peut donc pas encoder le
résultat d'une époque rejouée. Il ne sert d'ailleurs qu'à **dimensionner** le
vivier PSH, jamais à accepter ou refuser.

### La seule constante fautive, et sa taille

`PSH_BIAS` combine trois lectures : −1,6 (la seule séquence dont nous avons le
log du validateur), ≤ −3,1 (une acceptée que nous lisons à 70,3 alors que leur
borne basse est 73,4), et **−2,5 (Youden sur les 296 verdicts de 26299-26306)**.
La troisième est de l'information de l'époque rejouée, utilisée pendant son
rejeu.

Taille de l'effet. En retirant la troisième estimation, l'intervalle restant est
1,6 à 3,1 ; le choix « conservateur » documenté donnerait 1,6 (ou ~2,35 en
milieu d'intervalle). Sur les 8 valeurs `psh_best` effectivement retenues par le
rejeu (77,04 ; 77,89 ; 78,34 ; 79,44 ; 80,27 ; 83,09 ; 86,97 ; 88,12) et 3 plis :

| biais | P(refus) moyenne sur les 8 retenus | P(refus) à PSH 83 / 86 / 88,6 |
|---|---|---|
| 1,6 (sans l'estimation fautive) | 14,2 % | 8 % / 33 % / 65 % |
| **2,0 (en vigueur)** | **15,9 %** | **10 % / 37 % / 70 %** |
| 2,5 | 18,0 % | 13 % / 44 % / 75 % |
| 3,0 | — | 17 % / 50 % / 80 % |

`refusal_probability` est strictement croissante en PSH pour un biais fixé, donc
à qualité égale l'**ordre** des candidats ne change pas du tout quand le biais
bouge : seul le niveau absolu du risque bouge, de ≤ 4 points. Le biais n'entre
pas dans une porte binaire (il n'y a plus de plafond PSH dans ce dépôt), il
entre dans un prix.

### Verdict §2

**Une fuite confirmée, de taille négligeable** : `PSH_BIAS = 2.0` appuie une de
ses trois estimations sur les verdicts des époques 26299-26306. L'effet mesuré
est de ≤ 4 points de P(refus) et **aucun changement d'ordre** dans la sélection.
Toutes les autres constantes PSH et de développabilité — `SIGMA`,
`ACCEPTED_MEDIAN`, `PSH_PASS_RATE`, les bornes, les bandes, les σ Aggrescan3D,
le jeu de calibration de `preflight` — sont mesurées sur les époques 26174 et
26238 ou sur nos propres candidats, **toutes antérieures** aux époques rejouées.
Pas de fuite là.

---

## 3. La porte de nouveauté / doublons pendant un rejeu

### Ce que fait le code

`archive_exclude` existe dans `check_sequence` et `NanobodyArchive.duplicate_of`
pour exactement ce besoin (son docstring parle des `exclude_ids` du validateur).
Il n'est passé **par aucun appelant de production** :

| appel | coupure appliquée ? |
|---|---|
| `propose_nanobody.py:377` `check_sequence(seq, allowed_v, cdr3_len)` | sans archive du tout |
| `propose_nanobody.py:367-373` `cdr3_is_novel` | **oui** : `archived_cdr3` vient du `df` déjà coupé |
| `propose_nanobody.py:412` (chemin mineur, `skip_archive=False`) | **non** : `NanobodyArchive(target)` = archive d'aujourd'hui |
| `propose_nanobody.py:438` (chemin rejeu) | **non** |
| `mine.py:383` `archive_screen` | **non** (mais en production c'est correct) |
| `check_nanobody.py:72`, `preflight.py:318` | **non** |
| `tests/test_nanobody_parity.py` | **oui**, c'est le seul endroit du dépôt qui construit l'exclusion |

Le chemin « avertissements » a été vérifié : dans `check_sequence`, un coup
`hit` n'est une **rejection** que si `certain`, sinon il va dans `warnings` et le
candidat est conservé avec `archive_warning` rempli. Ce comportement est le même
avec et sans coupure, donc il n'introduit pas d'asymétrie entre rejeu et
production. `scripts/observe.py` et `scripts/review_observation.py` ne refont
aucun test de nouveauté : ils consomment la liste que `propose` a déjà filtrée.

### La direction de l'erreur

Juger la nouveauté contre l'archive d'aujourd'hui est **plus strict** que ce que
le validateur faisait à l'époque : un candidat peut être refusé pour sa
ressemblance à une séquence arrivée après. Cela **ne peut pas** gonfler le
résultat, cela ne peut que perdre des candidats. C'est d'ailleurs ce que
`tests/test_nanobody_parity.py::novelty_is_replayable` dit depuis longtemps :
après 26025 la porte n'est pas rejouable, et le test s'abstient d'affirmer sur
elle.

### La taille de ce pessimisme, mesurée

Sur les 80 candidats des lots journalisés de 26299-26306, similarité CDR3
maximale (règle `is_duplicate` à 0,85) contre l'archive numérotée, séparée en
instantané (9466) et post-instantané (1038 lignes numérotées) :

| | contre l'instantané | contre le post-instantané |
|---|---|---|
| médiane | 0,750 | 0,437 |
| maximum | **0,833** | **0,833** |
| candidats ≥ 0,85 | **0** | **0** |

Et sur l'identité entière (règle à 0,95, et le repli à 0,90 pour les lignes
sans donnée CDR), sur un échantillon de 20 des 80 candidats, contre les 9466 de
l'instantané et les 7778 post-instantané :

| identité entière maximale | contre l'instantané | contre le post-instantané |
|---|---|---|
| médiane | 0,864 | 0,915 |
| maximum | 0,941 | **0,948** |
| candidats ≥ 0,95 (rejection dure) | **0** | **0** |
| candidats ≥ 0,90 **seulement** contre le post-instantané | — | **12 / 20** |

Donc aucune rejection dure anachronique — mais 12 candidats sur 20 franchissent
0,90 contre une ligne que l'archive a gagnée depuis, et le plus proche est à
0,948, à deux millièmes de la rejection. Ce zéro tient à peu.

### L'asymétrie entre le rejeu et le mineur sur cette même porte

Le repli à 0,90 ne se déclenche que pour les lignes d'archive **sans donnée
CDR** — c'est-à-dire, à cause du cache empoisonné (§6c), pour 6740 des 7778
lignes post-instantané. Et il ne produit pas le même verdict des deux côtés :

```python
# nova_lab/nanobodies.py:287
hit, reason, certain = archive.duplicate_of(seq, archive_exclude)
if hit and (certain or warnings is None):
    problems.append(...)          # rejection
elif hit:
    warnings.append(...)          # avertissement, candidat conserve
```

`propose` (le rejeu) passe `warnings=[]`, donc un coup non-`certain` est un
**avertissement** et le candidat survit. `mine.archive_screen` (`mine.py:385`)
appelle `check_sequence(seq, allowed_v, cdr3_len, archive)` **sans** liste
d'avertissements, donc `warnings is None` et le même coup est une **rejection**.
Vérifié directement : `is_duplicate(0.93, None)` rend
`(True, "identity 93% without CDR data", False)`, qui devient « avertissement,
gardé » avec `warnings=[]` et « rejection » avec `warnings=None`.

Taille, sur les mêmes 20 candidats, identité maximale contre les **6740 lignes
non numérotées** (toutes post-instantané) :

| | résultat |
|---|---|
| ≥ 0,90 — gardé par le rejeu, **refusé par le mineur** | **3 / 20 (15 %)** |
| ≥ 0,95 — refusé des deux côtés | 0 / 20 |

**Environ un candidat sur sept survit dans le rejeu et serait jeté par le
mineur.** C'est un écart dans le sens qui gonfle, et il est indépendant de la
coupure : il existe aussi en direct, entre les deux appels.

### Verdict §3

**Deux verdicts distincts.**

*Sur l'anachronisme* (juger la nouveauté contre l'archive d'aujourd'hui) :
**pas de fuite**, l'erreur va dans le sens strict, et elle ne produit **aucune
rejection dure** à ces époques — 0 des 80 candidats atteint 0,85 de similarité
CDR3 contre une ligne post-instantané (max 0,833) et 0 des 20 échantillonnés
atteint 0,95 d'identité entière (max 0,948). Mais `archive_exclude` n'est
appliqué **nulle part** en production, donc ce zéro est une chance de ces huit
époques, pas une propriété du code. Réserve : seules 1038 des 7778 lignes
post-instantané sont numérotées dans le cache (§6c), donc la règle CDR3 n'a pu
être évaluée que sur 13 % d'entre elles.

*Sur l'asymétrie rejeu / mineur* : **écart confirmé, dans le sens qui gonfle, et
déclenché.** `propose` passe `warnings=[]` et `mine.archive_screen` ne passe
rien, donc un coup à 0,90 ≤ identité < 0,95 sans donnée CDR est un avertissement
d'un côté et une rejection de l'autre. **3 des 20 candidats échantillonnés
(15 %) sont dans ce cas** : le rejeu les garde, le mineur les jette. Le rejeu
mesure donc un mineur dont la porte de nouveauté est plus lâche que la vraie.

---

## 4. Le vivier de référence et la « place » rapportée

### Le code

```python
# propose_nanobody.py:427-430
# The real competition: the most recent dated submissions, which all comply.
models = Surrogate.load(f"{target}_nanobody").main
dated = df[df["epoch"].notna()]
reference = dated.nlargest(REFERENCE_POOL, "epoch") if len(dated) else compliant.head(REFERENCE_POOL)
```

`df` est déjà coupé à ce stade, donc `nlargest(200, "epoch")` ne peut pas
remonter dans le futur : **sur la question posée, il n'y a pas de fuite**.

Le défaut est ailleurs, et il est pire. Les seules lignes datées du parquet sont
**26306 et 26345-26348**. Donc pour toute coupure ≤ 26306 — c'est-à-dire pour
les huit époques rejouées — `dated` est **vide**, le ternaire bascule sans un mot
sur `compliant.head(REFERENCE_POOL)`, et le « vrai champ » annoncé par le
commentaire devient *les 200 premières lignes conformes de l'archive dans
l'ordre du fichier*.

Mesure, coupure par coupure :

| coupure | `len(dated)` | vivier retenu | médiane `design_iiptm` du vivier |
|---|---|---|---|
| 26299 | 0 | `compliant.head(200)` (**repli**) | **0,205** |
| 26306 | 0 | `compliant.head(200)` (**repli**) | **0,205** |
| aucune (mineur en direct) | 152 | `dated.nlargest(200, epoch)` | 0,698 |

0,205 contre 0,710 pour la médiane du champ : le vivier de repli n'est pas une
compétition, c'est un épouvantail.

### Ce que la « place » vaut alors

Le journal (`data/observation/journal.jsonl`) donne
`surrogate_place_vs_reference = 1` pour **24 observations sur 24**, aux huit
époques, sur les trois passes. `scripts.observe --report` en tire « place
médiane 1 », et c'est ce chiffre que les docs citent. La vérité GPU, au même
moment, est **0 victoire sur 8 et dernière place huit fois** (§46, §49).

La cause n'est même pas le substitut. En remplaçant les prédictions par les
**métriques réelles** et en classant dans le même vivier de repli :

| ce qu'on classe (métriques réelles) | meilleur `rank_sum` | place / 201 | place de la médiane |
|---|---|---|---|
| nos 64 candidats notés | 68 | **1** | 6 |
| les 32 vraies soumissions des époques rejouées | 65 | **1** | 1 |
| le vivier de repli lui-même | 137 (meilleur), 916 (médiane) | — | — |

Les deux côtés arrivent premiers. La métrique ne distingue pas un candidat qui
perd 8 fois sur 8 d'un gagnant d'époque : **elle n'a aucun pouvoir de
discrimination dans un rejeu**, et le `1` qu'elle affiche n'est pas une
performance, c'est la faiblesse du vivier.

Deux aggravations à noter :

1. `score()` met dans le même `rank_sum` les métriques **prédites** de nos
   candidats et les métriques **réelles** du vivier (`propose_nanobody.py:209`).
   Comparer une prédiction à une mesure dans un classement commun n'est pas
   comparable, et le vivier étant dans le jeu d'entraînement du substitut
   (§1b), le modèle est calibré exactement sur la distribution contre laquelle
   on le fait concourir.
2. Sans coupure (`--no-cutoff`, le chemin « biais assumé »), le vivier devient
   26306 + 26345-26348, c'est-à-dire du futur pur pour un rejeu de 26299. Le
   chemin est documenté comme biaisé, mais c'est la deuxième fuite qu'il porte,
   pas seulement celle du greffage.

### Verdict §4

**Le vivier ne contient pas les soumissions de l'époque rejouée ni de plus
tard** (il est construit après la coupure) : sur la fuite soupçonnée, pas de
fuite. Mais la `place` qu'il sert est **gonflée au maximum possible** : le repli
silencieux vers `compliant.head(200)` donne un vivier à 0,205 de médiane, et
`surrogate_place_vs_reference` vaut 1 à 24 observations sur 24 — y compris si on
lui soumet les métriques réelles de candidats qui ont perdu 8 fois sur 8. Le
commentaire du code (« the real competition: the most recent dated submissions »)
décrit une branche qui ne s'exécute jamais dans un rejeu.

---

## 5. `scripts/observe.py` et `scripts/review_observation.py` de bout en bout

### 5a. Le lot écrit contient-il du futur ?

Oui, et c'est voulu : `write_batch` ajoute `CALIBRATION = 4` vraies soumissions
de l'époque (`real_submissions(epoch, target)`). Ces séquences ne sont jamais
données au générateur ni au substitut — elles ne servent qu'à être notées dans
le même run BoltzGen, ce qui est la seule façon de comparer (deux tirages
diffèrent de ~14,5 places). `rebuild_batches` fait la même chose plus tard pour
les époques dont les rivaux étaient encore chiffrés. **Pas de fuite** : rien ne
remonte du lot vers la génération.

### 5b. La revue classe-t-elle contre un vivier assemblé avec du recul ?

Oui, mais **dans le sens dur**. `real_submissions` trie par `rank_sum` réel et
`write_batch` garde `head(4)` : les rivaux portés dans le lot sont les **4
meilleurs** de l'époque, choisis en sachant leur résultat. `review_observation`
calcule ensuite `place_among_rivals` contre ces 4-là seulement. Battre les 4
meilleurs d'une époque implique la gagner, donc `place == 1` reste une victoire
réelle ; l'inverse n'est pas vrai, mais c'est de la prudence, pas de la fuite.
Le script imprime déjà les deux réserves utiles : seuls les rivaux portés dans le
lot sont dans le pool, et seuls ceux que le validateur a **acceptés** ont un
`rank_sum`, donc les refus restent invisibles.

Une conséquence à garder en tête : `would_have_won` est un booléen sur un pool de
5, et `summary()` en fait un « taux de victoire » ; le script avertit sous 20
époques, mais l'avertissement porte sur la taille de l'échantillon, pas sur le
fait que le pool n'est pas le champ.

### 5c. Le rejeu choisit-il comme le mineur ? Non.

C'est l'écart le plus structurel du fichier, et il survit au correctif
`d99ccd5` qui prétendait l'avoir réglé.

```python
# mine.py:718
result = propose(..., shortlist=max(psh_pool, args.shortlist), seed=args.seed + epoch,
                 generator_name=args.generator, skip_surrogate=True, skip_archive=True, log=log)

# observe.py:227
result = propose(params, target, candidates=args.candidates, n_donors=args.donors,
                 n_cdr3=args.cdr3_pool,
                 shortlist=max(psh_pool, args.shortlist) if args.psh else args.shortlist,
                 seed=args.seed + epoch, archive_cutoff=cutoff,
                 generator_name=args.generator)
```

`observe` ne passe **ni `skip_surrogate=True` ni `skip_archive=True`**. Les deux
branches de `propose` sont alors différentes :

| | mineur (`skip_surrogate=True`) | rejeu (défaut `False`) |
|---|---|---|
| vivier remis au sélecteur | les `psh_pool` **premiers** candidats valides, dans l'ordre de génération (`propose_nanobody.py:407-420`) | les `psh_pool` **meilleurs** selon le substitut 10-métriques (`propose_nanobody.py:436-446`) |
| règle d'archive à ce stade | différée (`skip_archive=True`), appliquée ensuite aux survivants de la PSH | appliquée à tout le vivier |
| substitut 10-métriques | pas chargé | chargé, et c'est celui du §1b |

Le commentaire d'`observe.py:220-224` affirme le contraire (« le vivier sort donc
large et `scripts.mine.select_by_expected_value` le coupe, exactement comme en
production »). Il est faux : en production le vivier est **non trié**, dans le
rejeu il est **trié par un modèle entraîné sur le futur**. Le rejeu mesure donc
un mineur qui dispose d'un pré-filtre que le mineur n'a pas, et ce pré-filtre est
alimenté par la fuite n° 4.

Taille du recouvrement entre les deux viviers. Elle se calcule exactement, sans
relancer `propose` : le vivier du mineur est les `k` **premiers** valides et
celui du rejeu les `k` **mieux notés** par le substitut, donc si l'ordre de
génération est indépendant de la note — il l'est, la génération tire au hasard
dans `donors` et `cdr3_donors` — le recouvrement attendu est `k² / V`. À
l'époque 26303, le journal donne `V = 1769` candidats valides et
`psh_pool = 66`. Sur la dernière observation journalisée de chaque époque :

| époque | CDR3 | valides `V` | `psh_pool` `k` | recouvrement attendu `k²/V` |
|---|---|---|---|---|
| 26299 | 14 | 1260 | 28 | **0,6 sur 28** |
| 26300 | 8 | 1520 | 13 | **0,11 sur 13** |
| 26301 | 20 | 793 | 66 | **5,5 sur 66** |
| 26302 | 11 | 2785 | 20 | **0,14 sur 20** |
| 26303 | 18 | 1786 | 66 | **2,4 sur 66** |
| 26304 | 14 | 1825 | 28 | **0,4 sur 28** |
| 26305 | 10 | 2863 | 17 | **0,10 sur 17** |
| 26306 | 17 | 1737 | 50 | **1,4 sur 50** |

Autrement dit : **les deux viviers sont pratiquement disjoints**, de 0,1 à 5,5
séquences communes sur 13 à 66. Le rejeu ne mesure pas le sélecteur du mineur sur
le vivier du mineur ; il le mesure sur un vivier choisi par un modèle entraîné
sur le futur. (Un run de confirmation de `propose` dans les deux modes a été
lancé puis abandonné : la machine était saturée par un autre agent, et le calcul
ci-dessus ne dépend d'aucune mesure.)

### 5d. Autres écarts de paramètres

| paramètre | `mine.py` | `observe.py` | effet sur le rejeu |
|---|---|---|---|
| `--candidates` | 5000 | **3000** | affaiblit le rejeu |
| `--shortlist` | 8 | **10** | deux tirages GPU de plus → gonfle légèrement |
| `--generator` | `pssm` | `pssm` | aligné (corrigé le 14/11) |
| `--donors`, `--cdr3-pool`, `--psh-replicates` | 40 / 60 / 3 | 40 / 60 / 3 | alignés |
| rafraîchissement du vivier en cours d'époque (`donor_pool_lag`, `refresh_donor_pool`) | **oui** | **non** | le rejeu mesure un mineur à vivier figé |

### Verdict §5

**Le lot ne fuit pas** (les rivaux portés n'alimentent ni le générateur ni le
substitut) et **la revue ne fuit pas** (le recul sert à durcir le pool, pas à
l'adoucir). Mais le rejeu **ne choisit pas comme le mineur** : `observe` laisse
`skip_surrogate` et `skip_archive` à `False`, donc le vivier remis au sélecteur
est pré-trié par le substitut 10-métriques entraîné sur tout, là où le mineur
prend les candidats dans l'ordre de génération. Plus trois écarts de paramètres
par défaut (3000 vs 5000 candidats, 10 vs 8 de liste courte, pas de
rafraîchissement de vivier).

---

## 6. Le reste : ce que le ratissage a trouvé

Le `grep` demandé (`nlargest`, `epoch`, `max()`, `isna()`, `notna()` autour des
chargements de jeu de données) a sorti **deux autres occurrences de la forme
exacte du bug corrigé** — un `isna()` conservé à côté d'une comparaison
d'époque — et aucune des deux n'est corrigée.

### 6a. `scripts/humanise.py:120` — la plus grosse fuite du dépôt

```python
df = pd.read_parquet(DATA_DIR / "dataset" / f"{target}_nanobodies.parquet")
if archive_cutoff is not None:
    df = df[df["epoch"].isna() | (df["epoch"] < archive_cutoff)]
good = df[df["design_iiptm"] >= floor].copy()
```

C'est mot pour mot la ligne que `7a90b69` a remplacée dans
`propose_nanobody.py`, et elle est ici sur le chemin dont le dépôt a fait sa
**piste 2** : `humanise` prend « les bons lieurs de l'archive » (`floor = 0.72`
par défaut) et les humanise. Or la population que le correctif a mesurée comme
entièrement future est précisément *les donneurs non datés qui lient à 0,72 ou
mieux*.

Mesure, avec le filtre buggé et avec le filtre corrigé (instantané daté), à la
longueur de CDR3 que chaque époque impose :

| époque | CDR3 | donneurs ≥ 0,72 (filtre buggé) | meilleur | donneurs (filtre honnête) | meilleur |
|---|---|---|---|---|---|
| 26299 | 14 | 0 | — | 0 | — |
| 26300 | 8 | 1 | 0,721 | **0** | — |
| 26301 | 20 | 23 | 0,767 | **0** | — |
| 26302 | 11 | 41 | 0,785 | **0** | — |
| 26303 | 18 | 52 | 0,781 | **0** | — |
| 26304 | 14 | 0 | — | 0 | — |
| 26305 | 10 | 0 | — | 0 | — |
| 26306 | 17 | 0 | — | 0 | — |
| **total** | | **117** | **0,785** | **0** | **—** |

**Les 117 donneurs du bras `humanised` sont tous du futur. Honnêtement, ce bras
n'a aucun donneur à ces huit époques.**

Ce que ça fait au résultat publié (`results/scores_nanobody_arms_26299_26306.csv`) :

| bras | n | `design_iiptm` médian | max | époques gagnées |
|---|---|---|---|---|
| `real` | 32 | 0,710 | 0,778 | — |
| `grafted` | 32 | 0,499 | 0,696 | 0/8 |
| `humanised` | 31 | **0,185** | **0,694** | 0/8 |

Le bras était déjà réfuté, donc **la fuite n'a pas produit de fausse victoire** :
elle a produit une mesure vide. Le seul point où `humanised` approche le champ
(0,694 à l'époque 26301, soit 92 % du meilleur rival 0,755) sort des donneurs à
0,746-0,766 de la longueur 20 — les 23 donneurs post-instantané du tableau
ci-dessus. Sans eux il n'y avait rien à mesurer à cette époque.

Et le §51 est atteint au même endroit. Son chiffre d'affiche — « le vivier
conforme n'est plus de 20 séquences mais de **256** », avec un plafond à 0,788
au-dessus du maximum du champ — se décompose aujourd'hui ainsi :

| | lignes |
|---|---|
| conformes et `design_iiptm` ≥ 0,72, total | 75 |
| — dans l'instantané daté (≤ 26025) | **0** |
| — hors instantané (futur pour 26299-26306) | **75** |
| plafond `design_iiptm` conforme, instantané | **0,704** |
| plafond `design_iiptm` conforme, hors instantané | 0,781 |

**Zéro** des donneurs conformes au niveau du champ existait aux époques rejouées,
et le tableau par époque du §51 (15 / 1 / 23 / 41 / 52 / 15 / 0 / 20 donneurs)
devient 0 partout. Le correctif `7a90b69` annonce que les §48, §49 et §51 sont
viciés ; ce qui n'était pas dit, c'est que le **code** qui a produit la piste 2
porte encore le bug.

**Verdict 6a : fuite confirmée, la plus grande du dépôt — 117 donneurs du futur
contre 0 honnêtes. Inoffensive quant à la conclusion (le bras perdait déjà 8 fois
sur 8), mais elle vide la mesure.**

### 6b. `scripts/psh_by_generator.py:53` — même trou

```python
df = df[df["epoch"].isna() | (df["epoch"] < args.epoch)]
```

Même ligne, sur l'A/B des trois générateurs (`chimeric`, `pssm`, `pssm_denovo`,
`results/psh_generator_ab_cdr3_18.csv`). Le script admet donc les mêmes 197
lignes non datées du futur dans ses viviers de donneurs et de CDR3. L'A/B compare
trois générateurs **sur le même vivier**, donc la comparaison relative résiste
probablement ; les niveaux absolus de PSH qu'il rapporte, non.
**Verdict 6b : fuite confirmée, effet sur le classement relatif non mesuré
(il faudrait un run PSH, ~30 min de CPU par génerateur).**

### 6c. Le cache de régions est empoisonné — pas une fuite, mais il vide le §51

`NanobodyArchive._cached_regions` mémorise « non numérotable » comme un
`CDR3` vide, et ne le réessaie **jamais** (`missing` ne contient que les
séquences absentes du cache, pas celles marquées vides). État actuel de
`data/archive/Q9NZQ7_regions.csv` :

| | lignes |
|---|---|
| archive téléchargée | 17 244 |
| numérotées dans le cache | 10 504 |
| marquées « non numérotable » | **6 740** |
| dont réellement non numérotables (échantillon de 70 : 40 des 196 soumissions de 26299-26305 + 30 tirées au hasard des 6740) | **0** |

Les lignes nulles forment un bloc contigu : index 9466 à ~16 800, c'est-à-dire
**tout ce que l'archive a gagné depuis l'instantané**, sauf les ~450 dernières.
Conséquences mesurées :

- les **196 soumissions des époques 26299-26305** sont dans l'archive, avec
  leurs 10 métriques complètes, et **aucune** n'est dans le parquet
  (`data/dataset/Q9NZQ7_nanobodies.parquet`, 10 536 lignes) ;
- `duplicate_of` passe donc pour 6740 lignes d'archive par la branche dégradée
  `elif identity >= 0.90` au lieu de la règle CDR3 ;
- le « vivier conforme de 256 » du §51 en vaut **75** sur le parquet
  d'aujourd'hui.

La direction est pessimiste : le vivier est trop petit, pour le rejeu **et** pour
le mineur en direct, qui lisent le même fichier. Ce n'est donc pas une fuite.
Mais c'est un défaut de mesure de premier ordre, il est silencieux, et il rend
le chiffre central du §51 non reproductible aujourd'hui.

**Verdict 6c : pas une fuite ; défaut de mesure confirmé, 6740 lignes sur 17 244
perdues à tort, dont les 196 soumissions des époques rejouées.**

### 6d. `scripts/design_nanobody.py` n'a pas de coupure du tout

`framework_donors` lit le parquet entier et n'accepte aucun `archive_cutoff`
(`design_nanobody.py:86-95`). Tout rejeu passant par ce générateur greffe depuis
le vivier d'aujourd'hui. Il n'est pas sur le chemin de `observe`/`mine`, donc
aucun résultat publié n'en dépend à ma connaissance.
**Verdict 6d : fuite par construction, non exploitée. Dette.**

### 6e. Le `train_mask` du split « strict » porte lui aussi le bug

`train_nanobody_surrogate.py:107` :

```python
train_mask = (df["epoch"] < split) | df["epoch"].isna()
```

avec le commentaire « les séquences sans époque viennent de tours antérieurs du
même cible, donc elles sont sûres pour l'entraînement ». Depuis le
re-téléchargement de l'archive, c'est faux : une ligne sans époque peut être
d'hier. Le modèle concerné n'est pas sauvegardé — c'est celui du backtest
affiché — donc la fuite ne touche aucun artefact, mais elle touche **le chiffre
que le backtest publie** (« notre pick gagne n/k époques »).
**Verdict 6e : fuite confirmée dans le backtest du substitut, taille non
mesurée (il faudrait réentraîner avec l'instantané ; ~10 min de CPU).**

### 6g. La coupure porte sur l'époque de **soumission**, pas sur celle de **publication**

Une subtilité que la coupure ne traite pas, et qui est de la même classe.
`observe.py:216` pose `cutoff = epoch`, donc une soumission de l'époque `E-1`
est admissible pour un rejeu de `E`. Le déverrouillage timelock le permet : le
round qui ouvre les soumissions de `E-1` tombe dix blocs avant la fin de `E-1`
(`timelock.reveal_round`), donc la **séquence** était bien lisible avant que `E`
commence.

Mais le générateur ne classe pas ses donneurs sur la séquence : il les classe sur
`true_rank_sum`, donc sur les **dix métriques BoltzGen**, qui viennent de l'API
de partage des validateurs (`analyze_nanobodies.API`) et n'apparaissent qu'un
certain temps après la fin de l'époque. Si ce délai dépasse une époque — ce qui
est probable pour un pipeline GPU, mais je n'ai pas de mesure du délai de
publication sans réseau — alors la coupure `epoch` est trop laxiste d'autant.

**Ce que ça coûte aujourd'hui : rien, et c'est mesurable.** Les seules lignes
datées de tout le dépôt sont les époques 26299-26306 et 26345-26348, dans
`data/scores/nanobodies.jsonl` comme dans le parquet. **Aucune ligne datée
n'existe entre 26025 et 26299.** Donc pour les huit rejeux publiés, `dated` est
vide quelle que soit la coupure, et déplacer la coupure de `E` à `E-k` ne change
**0 ligne**. Le problème est devant nous, pas derrière : il mordra dès que le
registre `nova_lab/archive_ledger.py` et le cache de scores se rempliront.

**Verdict 6g : fuite potentielle confirmée dans le code, de taille nulle
aujourd'hui (0 ligne concernée). Le délai de publication des métriques n'a pas
pu être mesuré hors réseau — c'est le « je ne peux pas dire » de cet audit.**

### 6f. Ce qui a été vérifié et ne fuit pas

- **Paramètres d'époque** : `resolve_epoch` lit le hash du premier bloc de
  l'époque **rejouée** sur `wss://archive.chain.opentensor.ai` et dérive gène V
  et longueur CDR3 de ce hash seul. `len(reaction_ids())` est d'aujourd'hui, mais
  il n'est consommé qu'avant que `rng` soit touché, donc il ne déplace ni l'un ni
  l'autre — `tests/test_challenge_params.py` l'épingle, et
  `data/submissions/epoch_params_cache.json` ne met en cache qu'un hash de bloc,
  qui est immuable.
- **`archived_cdr3` / `cdr3_is_novel`** (`propose_nanobody.py:365-373`) : bien
  construits sur le `df` coupé.
- **Vivier de donneurs et de CDR3** (`nsmallest` sur `true_rank_sum`,
  `propose_nanobody.py:302,333`) : construits sur le `df` coupé, donc corrects
  depuis `7a90b69`.
- **`observe.py` `real_submissions` / `write_batch` / `rebuild_batches`** : ne
  réalimentent jamais la génération (§5a).
- **`scripts/analyze_field_method.py`** (non suivi par git, écrit pendant cet
  audit) lit les 348 soumissions de 26299-26348 : c'est un script d'analyse
  rétrospective, pas un chemin de rejeu. Hors périmètre.

---

## Ce qu'il faudrait pour que les rejeux redeviennent utiles

Par ordre de rapport qualité/prix, et sans rien promettre sur les résultats :

1. **`scripts/humanise.py:120`** : appliquer `dated_archive_snapshot`, comme
   `propose_nanobody.py` le fait maintenant. Deux lignes. Après quoi la piste 2
   n'a plus de donneur aux époques 26299-26306 et il faut rejouer ailleurs ou
   attendre que le registre (`nova_lab/archive_ledger.py`) se remplisse.
2. **`scripts/psh_by_generator.py:53`** : idem.
3. **Un substitut par coupure.** `Surrogate.load` et `surrogate_rank` doivent
   prendre une époque et charger (ou entraîner) un modèle n'ayant vu que des
   lignes antérieures. Sans cela, le chiffre « notre sélecteur choisit mieux que
   le hasard » n'est pas rejouable : aujourd'hui 438/438 de ses étiquettes sont
   dans le futur.
4. **Le repli du vivier de référence** doit échouer bruyamment au lieu de prendre
   `compliant.head(200)`. En l'état `surrogate_place_vs_reference` devrait être
   retiré du journal et des rapports : il vaut 1 quoi qu'on lui donne.
5. **`observe` doit appeler `propose` avec `skip_surrogate=True, skip_archive=True`**,
   comme `mine`, et aligner `--candidates` sur 5000 et `--shortlist` sur 8.
6. **Le cache de régions** doit réessayer les lignes marquées vides, ou être
   invalidé. 6740 lignes sur 17 244 sont perdues à tort, dont les 196 soumissions
   des époques rejouées.
7. **`PSH_BIAS`** : documenter que la troisième estimation vient des époques
   rejouées, et la retirer des rejeux de 26299-26306. L'effet est de ≤ 4 points
   de P(refus), donc c'est de l'hygiène, pas une urgence.
