# Ce que nous avons mesuré sur SN68 — synthèse

## Les données brutes derrière chaque chiffre

| Fichier | Contenu |
|---|---|
| `results/psh_by_generator.csv` | 540 candidats, 3 stratégies × 3 longueurs de CDR3, PSH mesurée (§3) |
| `results/place_vs_psh.csv` | 21 époques : place dans le champ complet, PSH, admissibilité (§3bis) |
| `results/psh_by_length.csv` | 120 séquences **acceptées** par les validateurs, PSH par longueur de CDR3 (§3ter) |
| `results/tnp_epoch26174.csv` | les 36 acceptées de l'époque 26174 + notre rejetée, 6 métriques TNP (§2) |
| `results/aggrescan_epoch26174.csv` | Aggrescan3D sur ce même champ (§2) |
| `results/tnp_developability.csv` | TNP sur nos 84 candidats et 84 soumissions réelles (§4) |
| `results/aggrescan_gate.csv` | Aggrescan3D sur les mêmes 168 (§4) |
| `results/scores_nanobody_observation_depth4.csv` | les 168 scores BoltzGen du lot à profondeur 4 (§7) |
| `results/economics.txt` | seuil de rentabilité, dimensionnement GPU |
| `journal_final.jsonl`, `mine_live_final.log` | le run de production (11 époques) |

`data/` est ignoré par git : tout ce qui soutient une conclusion est copié dans
`results/`.

État au 2026-11-09. Ce document rassemble les faits mesurés, les erreurs
commises et ce qui reste ouvert. Il est organisé par sujet, pas
chronologiquement. `docs/OBSERVATIONS.md` et `docs/POST_MORTEM.md` gardent le
déroulé.

---

## 1. Le fait central : l'admission, pas le classement

Sur 71 époques, **le taux de rejet médian est de 50 %** et seulement 2 époques
sur 61 passent sous 10 %. Il varie énormément selon l'opérateur, et c'est lui
qui prédit les victoires :

| Mineur | Soumissions | Notées | Victoires |
|---|---|---|---|
| frankji-groundcontrol | 54 | **83 %** | **30** |
| thomasvangurp | 384 | 78 % | 5 |
| QTTECH2023 | 156 | 73 % | 13 |
| bitty-labs | 205 | 68 % | 2 |
| saon-miner | 338 | 59 % | 5 |
| cczulaika | 97 | 40 % | 2 |
| uchiha0818 | 599 | 26 % | 1 |
| aichorx | 91 | 12 % | 0 |
| **nous (UID 22)** | **11** | **0 %** | 0 |

uchiha0818 a soumis 599 fois pour une victoire ; frankji 54 fois pour trente.

**Toutes mes estimations hors ligne mesuraient le classement parmi les
soumissions déjà admises.** Il fallait les multiplier par une probabilité
d'admission que je n'ai jamais mesurée. C'est l'erreur de fond de toute la
phase de recherche, et elle était présente depuis le premier jour : mon jeu de
données ne contenait que les soumissions ayant des métriques, donc les 50 %
rejetées étaient invisibles — un biais de sélection intégré aux données.

---

## 2. Pourquoi nos 11 soumissions ont été rejetées

Le log W&B du validateur le dit directement :

```
UID 22: contains sequences that are rejected by PSH strict rejection: [(..., 89.5651)]
```

`config.yaml` fixe `psh_strict_rejection_threshold: 89`, ajouté le 2026-10-20.
C'est **1,7 fois plus strict** que le seuil rouge de TNP (155,47), et je ne
l'avais pas implémenté.

### La PSH n'est pas mal calibrée, elle est bruitée

Ce n'est pas une différence d'échelle entre notre TNP et le leur. La PSH se
calcule sur une **structure prédite**, et NanoBodyBuilder2 n'est pas
déterministe.

| Mesure du bruit | Valeur |
|---|---|
| 3 répliques d'une séquence | étendue 3,46 (sous-estimation) |
| **10 séquences mesurées deux fois** | **écart-type d'une mesure 4,01** |
| écart médian entre deux lectures | 3,3 points |
| accord de classement entre deux lectures | Spearman 0,89 |

**Pour un seuil à 89, une mesure unique dévie de 4.** Notre 85 et leur 89,57
sont deux tirages de la même distribution. Cela explique aussi pourquoi 44 %
des soumissions *acceptées* mesurent plus de 89 chez nous.

**Conséquence : aucune calibration n'est possible.** On ne calibre pas contre
du bruit. Le correctif proposé dans le post-mortem — abaisser le seuil local à
75 — rejetterait **100 %** des séquences qui passent réellement chez le
validateur (elles s'étalent de 75,3 à 99,3 sur notre échelle).

### Ce qui marche à la place

Moyenner k plis divise l'écart-type par √k :

| Plis | σ de la moyenne | Marge 2σ | Seuil utilisable | Nos candidats retenus |
|---|---|---|---|---|
| 1 | 4,01 | 8,0 | 81,0 | 12 % |
| **3** | **2,32** | **4,6** | **84,4** | **17 %** |
| 5 | 1,79 | 3,6 | 85,4 | 17 % |

k=3 est le bon compromis. Le criblage se fait en deux étages : une passe sur
tout le vivier (3,6 s/candidat sur 4 cœurs), puis deux passes sur la tête de
classement. Pour 300 candidats dont on garde 40 : ~23 min, contre 54 si l'on
repliait tout trois fois.

---

## 3. La PSH ne se prédit pas depuis la séquence

Testé sur 205 PSH mesurées :

| Indicateur | Spearman |
|---|---|
| longueur des CDR | **+0,59** |
| nombre de résidus hydrophobes en CDR | +0,41 |
| somme Kyte-Doolittle positive | +0,38 |
| hydrophobicité moyenne | −0,09 |

Après correction de la longueur, le signal de composition tombe à **+0,30** et
**change de signe** selon la longueur du CDR3 : −0,22 à longueur 15, +0,58 à
longueur 23.

**Donc une pénalité d'hydrophobicité à la génération agirait au hasard dans la
moitié des époques.** Et la longueur n'est pas actionnable : CDR1+CDR2 est figé
à 15–16 résidus pour nous comme pour tous les concurrents.

### Les stratégies de génération ne changent rien à la PSH

540 candidats, 60 par stratégie et par longueur, PSH mesurée avec TNP :

| CDR3 | chimeric | pssm | pssm_denovo |
|---|---|---|---|
| 8 | 88,5 / **53 %** | 89,9 / 43 % | 88,2 / 57 % |
| 11 | 93,2 / **28 %** | 92,7 / 23 % | 93,1 / 33 % |
| 18 | 100,5 / **10 %** | 96,8 / 13 % | 100,2 / 3 % |

*(PSH médiane / part sous 89)*

Aucune différence significative : à la longueur 8, chimeric contre pssm donne
p = 0,36, contre pssm_denovo p = 0,85. **La PSSM ne corrige pas la PSH.** Elle
change bien la distribution des acides aminés — c'est réel et c'est ce que
mesure `scripts.compare_generators` — mais cela ne se traduit pas en PSH.

> Un premier essai à n=30 montrait pssm à 30 % et pssm_denovo à 57 %, et je
> l'ai rapporté comme un effet. En doublant l'effectif l'écart a disparu. Je
> l'avais pourtant prédit : le signal de composition ne vaut que +0,30 après
> correction de la longueur.

### Ce qui pilote la PSH, c'est la longueur imposée

| CDR3 | PSH médiane | ≤89 | ≤84,4 |
|---|---|---|---|
| 8 | 88,5 | 53 % | 20 % |
| 11 | 93,2 | 28 % | 8 % |
| 18 | 100,5 | 10 % | 2 % |

Corrélation longueur/PSH : **+0,81**. Et nous ne choisissons pas la longueur.

> **Leçon de méthode.** `scripts.compare_generators` mesure l'hydropathie et la
> divergence de Jensen-Shannon — des proxys non validés. Sur ces proxys, la
> PSSM paraît excellente. Sur la PSH mesurée, elle est pire que la référence.
> Un proxy non validé peut inverser une décision.

---

## 3bis. La tension structurelle : nous gagnons là où nous ne sommes pas admis

C'est le résultat le plus important de tout le projet.

Sur les 21 époques du jeu noté à profondeur 4, en croisant la place dans le
champ complet et l'admissibilité PSH de nos candidats :

| Bande de CDR3 | Époques | Gagnées (classement) | Au moins un candidat sous 89 | Avec marge (≤84,4) |
|---|---|---|---|---|
| court ≤10 | 4 | 1 | 4 | 3 |
| moyen 11–16 | 8 | 2 | 6 | 3 |
| **long ≥17** | **9** | **6** | **0** | **0** |

- longueur CDR3 vs notre place : Spearman **−0,44** (on classe **mieux** sur les longs)
- longueur CDR3 vs PSH : Spearman **+0,81** (la PSH **monte** avec la longueur)

**Les deux effets s'opposent exactement.** Sur les 9 époques que nous
gagnerions au classement, **aucune n'a de candidat sous 89** si elle est longue ;
toutes bandes confondues, seules 3 ont un candidat admissible et **2** en ont un
avec marge.

| | Époques | Taux |
|---|---|---|
| gagnées sur le classement seul | 9/21 | 43 % |
| admissibles avec marge seule | 6/21 | 29 % |
| **les deux à la fois** | **2/21** | **10 %** |

**Le taux de victoire réellement atteignable est donc ~10 %**, IC 95 % [1 %,
30 %] — pas les 29 % que j'annonçais, qui ne mesuraient que le classement.

Cela explique aussi pourquoi aucun réglage du générateur n'aide : **nous sommes
bons là où nous ne passons pas, et nous passons là où nous sommes mauvais.**
Changer la stratégie de CDR3 ne touche ni l'un ni l'autre.

Face au seuil de rentabilité de 0,3–1 %, 10 % reste au-dessus. Mais la marge
d'erreur a fondu, et la borne basse de l'intervalle (1 %) touche le seuil.

### Et le champ se vide quand la boucle s'allonge

Mesuré sur les 71 époques notées du champ réel (1252 soumissions,
`data/submissions/Q9NZQ7_nanobody_results.csv`) :

| CDR3 imposé | Époques | Champ médian | Champ moyen | Minimum |
|---|---|---|---|---|
| court 5–8 | 19 | 26 | 24,2 | 11 |
| moyen 9–12 | 11 | 20 | 20,5 | 10 |
| long 13–16 | 18 | 22 | 18,9 | 6 |
| **très long 17–23** | **23** | **9** | 9,9 | 1 |

Spearman longueur vs taille du champ : **ρ = −0,63, p = 3×10⁻⁹**.

Cela double l'enjeu du tableau précédent : les époques où nous classons le
mieux sont aussi **celles où il y a le moins de rivaux** — 9 concurrents au
lieu de 26, sur un tiers du calendrier (23 des 71 époques). Ce n'est pas une
coïncidence : si presque personne ne soumet quand la boucle est longue, c'est
vraisemblablement que presque personne n'y passe les portes non plus.

D'où la question que la mesure suivante tranche : la PSH des séquences que les
validateurs ont **acceptées** à CDR3 ≥ 17. Leur PSH réelle est ≤ 89 par
construction. Si elles lisent elles aussi au-dessus de 89 chez nous, le mur
est en grande partie le nôtre (`scripts/psh_by_length.py`).

---

## 3ter. Le seuil PSH plat de 89 est faux tel que nous l'appliquons

Mesure directe, `results/psh_by_length.csv` : 120 séquences que les validateurs
ont **acceptées**, profilées par notre propre TNP, regroupées par la longueur de
CDR3 que leur époque imposait. Chacune a franchi la porte `psh > 89 → rejet`,
donc toute lecture au-dessus de 89 chez nous est un refus que nous aurions
prononcé à tort.

| Bande de CDR3 | n | Époques | PSH médiane | 90ᵉ pct | > 89 | Refusées à tort |
|---|---|---|---|---|---|---|
| court 5–8 | 18 | 13 | 85,1 | 92,6 | 28 % | 5 |
| moyen 9–12 | 18 | 9 | 86,8 | 95,1 | 33 % | 6 |
| long 13–16 | 30 | 13 | 90,6 | 97,0 | 60 % | 18 |
| **très long 17+** | **54** | **20** | **92,0** | **100,5** | **76 %** | **41** |

Pente : **+0,63 point de PSH par résidu de CDR3** (p = 2×10⁻⁸, R² = 0,23),
Spearman ρ = 0,46 (p = 2×10⁻⁷).

**Ce qui est établi.** À CDR3 ≥ 17, lire 92 est la lecture *typique* d'une
séquence acceptée. Refuser un candidat parce que notre chiffre dépasse 89 y
écarte 76 % de la population que le validateur admet. Le seuil plat doit
devenir un quantile par longueur. C'est vrai quel qu'en soit le mécanisme.

**Ce qui n'est pas établi — et que j'ai d'abord écrit à tort.** Que notre
implémentation sur-lise. Deux mécanismes produisent la même signature :

- **(a) biais** : notre chaîne lit trop haut à mesure que la boucle s'allonge
  (NanoBodyBuilder2 sur CPU chez nous, sur GPU chez le validateur ; la
  minimisation OpenMM reste sur CPU des deux côtés, donc l'écart viendrait du
  réseau, pas du raffinement) ;
- **(b) sélection** : la mesure du validateur est bruitée elle aussi
  (σ = 4,01 mesuré) et sa porte retient la **queue chanceuse** d'une population
  dont la PSH vraie est vers 98–100 ; notre tirage indépendant retombe sur la
  vraie valeur, au-dessus de 89.

J'ai simulé (b) : population vraie de moyenne 98–100 et de dispersion ~4,
bruit σ = 4,01 des deux côtés → ensemble accepté de médiane ~92, 90ᵉ
percentile ~99, part au-dessus de 89 ≈ 76 %. **(b) suffit quantitativement**,
sans aucun biais de notre part. Le même jeu de paramètres implique un taux
d'admission de 2 à 3 % par soumission longue — ce qui explique aussi, tout
seul, pourquoi le champ n'a que 9 soumissions à ces longueurs.

**Même correctif, espoirs différents.** Sous (a), un candidat long lu à 92
passe le plus souvent ; sous (b), il passe environ une fois sur cinq. Les deux
valent mieux que les **0 %** d'aujourd'hui, où nous refusons de soumettre dans
un tiers du calendrier — celui où nous classons 1ᵉʳ contre 9 rivaux.

Ce qui trancherait : notre PSH sur les soumissions que le validateur a
**refusées** aux mêmes longueurs. Sous (a) elle sépare nettement refusées et
acceptées ; sous (b) elle les sépare mal, les deux groupes venant de la même
population coupée par un seuil bruité. Ces séquences sont dans
`data/submissions/submissions.jsonl` (soumise mais non notée = refusée), donc
la mesure ne coûte que du CPU.

---

## 4. Les verrous du validateur, reproduits et vérifiés

Tout ceci tourne hors ligne, sans GPU :

| Verrou | Reproduit ? | Note |
|---|---|---|
| longueur, acides aminés, homopolymères, di-répétitions | oui | parité testée |
| peptide signal | oui | |
| identité germinale V ≥ 82, J ≥ 86 | oui | ANARCI, 128 états |
| gène V imposé par l'époque | oui | |
| classification VHH, cystéines canoniques | oui | |
| longueur de CDR3 imposée | oui | |
| unicité et similarité à l'archive | oui | plus strict que le validateur |
| TNP : 6 drapeaux (L, L3, C, PSH, PPC, PNC) | oui | seuils identiques au bit |
| Aggrescan3D (moyenne ≤ −0,680, pic ≤ 1,580) | oui | scorer Rust + FreeSASA 2.0.1 épinglé |
| **PSH strict ≤ 89** | **valeur bruitée** | voir §2 |

Contrôles positifs : sur les 84 soumissions réelles acceptées, TNP rejette 0 %
et A3D 4 %. Sur nos candidats, TNP 2,4 % et A3D 39 %.

**Correction d'une affirmation que j'ai répétée :** cette porte de
développabilité n'est **pas** invérifiable hors ligne. TNP est installable par
pip, prédit la structure sur CPU en quelques secondes, et les 168 séquences ont
été profilées en dix minutes sur 4 cœurs. Les « 52 % de rejets inconnaissables »
que j'invoquais pour nuancer chaque conclusion étaient mesurables depuis le
début.

---

## 5. Le vivier de donneurs était périmé

En production le générateur greffait depuis un jeu de données arrêté à l'époque
**26025** alors qu'il jouait **26174** : 149 époques de retard, 7,5 jours, 26 %
de l'archive absente.

**Aucune des 55 soumissions notées des époques 26173-26174 n'était dans le
vivier**, les cinq premières de chaque époque comprises.

| | Donneurs contemporains (20 époques précédentes) |
|---|---|
| mes rejeux hors ligne | 143 à 316 |
| **production** | **0** |

Trois défauts de ma part s'y renforçaient : le mineur ne vérifiait pas l'âge,
`bootstrap` affichait « ok » pour un fichier d'une semaine, et j'avais livré ce
fichier figé dans `data/prebuilt/` en vantant le gain de temps. **J'avais
optimisé la vitesse de chargement d'un artefact dont la fraîcheur était la
propriété critique.** Corrigé : le mineur refuse au-delà de `--max-donor-lag`
époques (40 par défaut).

---

## 5ter. Deux interprètes Python dans le conteneur, un seul sait soumettre

Le 9 novembre, en préparant l'essai à blanc, la première lecture de la chaîne
a échoué :

```
AttributeError: 'Subtensor' object has no attribute 'get_current_block'
```

Diagnostic : le conteneur a deux interprètes, et `python` nu désigne le
**mauvais**.

| Interprète | `bittensor` | API de soumission |
|---|---|---|
| `/home/user/venv/bin/python` (le projet, le pin) | 10.4.1 | complète |
| `/usr/bin/python` (système, sur le `PATH` par défaut) | 11.3.0 | **absente** |

Dans la 11.3.0 l'objet `Subtensor` a été réécrit : les trois appels que fait le
mineur y ont tous disparu — `get_current_block`,
`get_uid_for_hotkey_on_subnet`, `set_commitment` — et `bittensor_drand`, donc
le verrou temporel, n'y est pas installé du tout.

**La production n'a jamais été touchée** : le mineur tournait sous le venv.
Ce que cela casse, c'est un lancement fait sans activer le venv — et le
message, un `AttributeError` nu sur la première lecture de la chaîne, ne dit
rien de la cause. J'ai moi-même passé vingt minutes à conclure à une
régression de l'API avant de regarder `which python`.

Deux garde-fous, dans `nova_lab/chain.py` :

- `current_block(st)` lit la tête de chaîne via la méthode (10.x) ou la
  propriété (11.x), donc toutes les lectures survivent aux deux API ;
- `require_api(st)` refuse de démarrer, **avant** l'époque, quand la version
  chargée ne peut pas soumettre, en nommant le pin à réinstaller ; en essai à
  blanc c'est un simple avertissement.

Lancer le mineur avec `/home/user/venv/bin/python -m scripts.mine`, jamais
`python -m scripts.mine`.

---

## 6. Hypothèses que j'ai formulées et qui étaient fausses

| Hypothèse | Comment je l'ai réfutée |
|---|---|
| le seuil PSH n'est pas la cause | faux : le log W&B le nomme. J'avais conclu « pas la cause » parce que ma reproduction de la *valeur* était mauvaise — un mécanisme ne cesse pas d'être la cause parce qu'on le mesure mal |
| dérive du temps de bloc (11,34 s) | j'avais pris les horodatages du journal pour une horloge de blocs, alors que le mineur avait redémarré. La chaîne est à 12,00 s, horloge synchrone, et le verrou temporel a **+1,9 min de marge** |
| engagement trop précoce dans l'époque | 45 % de notation pour les précoces contre 54 % pour les tardifs, p = 0,17 |
| écart-type de la PSH = 1,83 | estimé sur 3 répliques d'une seule séquence ; la vraie valeur est 4,01 |
| les 10 framework donors suffisent | invalidé par le test sur l'époque 25989 |

---

## 7. Mesures qui tiennent toujours

- **BoltzGen est un classeur fidèle** : 86 % d'accord par paires avec les
  validateurs sur 118 paires (p = 1,7×10⁻¹⁶), reproduit sur trois runs
  indépendants (85 %, 85 %, 86 %).
- **Le modèle de substitution ne classe pas** : Spearman −0,20 sur nos
  candidats, et son n°1 est notre meilleur 4 fois sur 21 — moins bien que le
  hasard (25 %). Son backtest temporel est pourtant flatteur (Spearman 0,55,
  place médiane 1 sur 22) : il généralise dans le temps, pas d'un générateur à
  l'autre.
- **`liability_score` et `liability_num_violations` sont exactement
  déterministes** : écart 0,0000 sur 84 séquences notées deux fois, là où les
  8 autres métriques varient. 2 métriques sur 10, optimisables sans GPU.
- **Le départage favorise le plus rapide** : somme de rangs → catégories →
  bloc d'engagement → horodatage GitHub → UID. Noter plus de candidats retarde
  la soumission et perd les égalités.
- **La cible change tous les ~28 jours en médiane** (extrêmes 7 et 53), et les
  validateurs se mettent à jour toutes les heures via l'auto-updater.
- **Le verrou timelock fonctionne** : round visé correct, +1,9 min de marge, et
  drand répond HTTP 425 avant l'échéance.
- **Le subnet est plein** : 64 UIDs, dont 7 validateurs ; enregistrement à
  τ0,04 ; immunité 5000 blocs (~14 époques).

---

## 7bis. Le budget d'époque, enfin mesuré

Époque 26238, départ à 2 blocs de la frontière, **67 min de budget** (fin
d'époque moins les 30 blocs sans soumission moins 12 blocs de sécurité). Premier
passage complet de la chaîne.

| Étage | Temps | Part du budget |
|---|---|---|
| `propose` | **1853 s (30,9 min)** | 46 % |
| `psh_select` (300 candidats, puis 3 replis sur 24) | 1058 s (17,6 min) | 26 % |
| `developability` (TNP + A3D sur 8) | 33 s | 1 % |
| `boltzgen` | indisponible ici | — |
| **total** | **2946 s (49,1 min)** | 73 %, **18 min restants** |

Donc le budget tient — mais BoltzGen n'a pas tourné (binaire absent de ce
conteneur, pas de GPU). À ~1,6 min par candidat sur 2×T4, 7 candidats coûtent
~11 min : 49 + 11 = 60 min sur 67. **7 min de marge, c'est trop peu.**

### 28 des 31 minutes de `propose` servaient à rien

Mesuré directement, par opération et par candidat :

| opération | ms/candidat | pour 300 |
|---|---|---|
| `regions` (ANARCI) | 18 | 0,1 min |
| `germline_identity` (hmmer) | 15 | 0,1 min |
| `check_sequence` sans archive | <1 | 0,0 min |
| **`duplicate_of` (archive)** | **5 595** | **28,0 min** |

La règle d'archive aligne le candidat contre **chacune des 14 424** séquences
archivées avec `difflib`. `propose` l'appliquait au lot entier — 300 candidats —
alors que la PSH en jette ensuite 97 %. Et le coût a grandi de 52 % aujourd'hui
même, quand l'archive est passée de 9 466 à 14 424 séquences : un étage dont le
prix suit la taille de l'archive, dans une boucle qui tourne à chaque époque.

Corrigé : `propose(skip_archive=True)` et la règle tourne dans le mineur sur la
poignée qui survit à la PSH. **44 s pour 8 candidats**, mesuré. La séquence
réellement soumise reste intégralement vérifiée, et un doublon reste un refus
certain — pas de repli « soumettre quand même », contrairement au crible de
développabilité.

### Mesuré après les deux correctifs

Même époque, même graine, mêmes donneurs :

| | avant | après |
|---|---|---|
| `propose` | 1853 s | **145 à 189 s** (×10 à ×12,8) |
| candidats valides sur 5000 | 829 | **2163** (×2,6) |
| refus sur l'identité J | 2919 | **0** |

(Le `propose` corrigé a été remesuré quatre fois : 145, 171, 181 et 189 s. Le
coût de la sonde d'appariement varie avec le nombre de donneurs sondés ; la
dispersion de ±15 % est à retenir plutôt que la première valeur.)

Budget d'époque recomposé : propose 145 s + PSH 1058 s + archive 44 s +
développabilité 33 s + BoltzGen ~660 s ≈ **32 min sur 67**. Au lieu de 7 min
de marge, **35** — et de quoi tripler le pool PSH sur les époques longues.

---

## 7ter. Le greffon cassait un segment germinal : un cas sur deux corrigé

Deux défauts symétriques, trouvés en croisant les donneurs plutôt qu'en
raisonnant. Les deux mécanismes sont réels et mesurés. **Un seul des deux
correctifs améliore le résultat**, et c'est tout l'intérêt de l'histoire.

### Le gène J : corrigé, ×2,6 sur le vivier de candidats

Le segment J couvre les **3 derniers résidus du CDR3 plus les 11 de la FR4** —
quatorze positions. Le chimérisme prend la queue chez un donneur et la FR4 chez
un autre : la séquence porte alors deux gènes IGHJ, lit 12/14 = **85,71 %**
contre un plancher de 86, et le validateur la refuse pour trois dixièmes de
point. Époque 26238 : **2919 candidats sur 5000**, tous sur cette même valeur.

Croisement de neuf queues × quatre variantes de FR4, corps fixe :

| queue | L | M | Q | T |
|---|---|---|---|---|
| FDY | **100** | 92,9 | 92,9 | 92,9 |
| FDV | 92,9 | **100** | 92,9 | 92,9 |
| MDV | 92,9 | 92,9 | 92,9 | **100** |
| RDV | *85,7* | 92,9 | *85,7* | 92,9 |
| DDY | 92,9 | *85,7* | *85,7* | *85,7* |
| MDG | *85,7* | *85,7* | *85,7* | 92,9 |

Validé sur la vérité terrain : **69 % des 2536 soumissions acceptées** utilisent
une des quatre queues qui passent partout (`FDY`, `FDV`, `FDP`, `MDV`), `FDY`
seule en faisant la moitié. Les séquences admises portent des extrémités
germino-compatibles parce qu'elles n'ont pas le choix.

Correctif : le donneur de CDR3 est tiré parmi ceux compatibles avec la FR4 de la
charpente. Mesuré à graine et donneurs identiques : **38 % → 100 %** de passage
du plancher J sur 400 candidats.

### Le gène V : mécanisme réel, correctif inutile — retiré

Le segment V couvre CDR1, CDR2 et le **début** du CDR3 : IGHV3-66 finit par
`…YYCAR`, donc la tête du CDR3 est codée par le V. En assemblant 144 paires
**sans aucune mutation**, 38 % perdent le gène V de l'époque, et la coupure est
franche : tête `AR` → 100 % des charpentes gardent l'appel, tête `AK` → 25 %.

J'ai donc étendu la sonde au V. Résultat :

| variante | valides sur 5000 |
|---|---|
| origine | 829 |
| **J seul** | **2163** |
| J + V, avec repli vers le vivier entier | 1955 |
| J + V, sans repli | 2129 |

**Les refus sur le V restent à ~1077 sur 5000 avec le criblage en place.**
Cribler la paire ne corrige que l'assemblage non muté ; ce qui casse l'appel V,
c'est la suite — CDR1 muté 15 % du temps, CDR2 10 %, et le corps du CDR3
librement. Le criblage réduisait les paires utilisables de 2400 à 654 et
écartait trois charpentes pour ne rien gagner. Retiré.

C'est la même leçon que le σ de la PSH et que « le seuil PSH n'est pas la
cause » : **un mécanisme établi ne dit pas quelle part du résultat il
explique.** Seule la mesure de bout en bout le dit.

### Quatre tentatives sur le gène V, quatre rejets par la mesure

| tentative | valides sur 5000 (CDR3 = 7) |
|---|---|
| **référence : appariement J seul** | **2163** |
| criblage des paires sur le V, avec repli | 1955 |
| criblage des paires sur le V, sans repli | 2129 |
| CDR1/CDR2 jamais mutés | inchangé (28 % de refus V dans les 4 variantes) |
| donneurs de CDR3 filtrés sur le gène V | 1751 |
| tête du CDR3 protégée (`V_HEAD = 2`) | **0** |

Le mécanisme, lui, a fini par être trouvé et il est solide : **le deuxième
résidu du CDR3 porte la signature du gène V** (IGHV3-23 finit par `…YYCAK`,
IGHV3-66 par `…YYCAR`). Sur 1200 candidats, une tête en `R` est refusée sur
l'appel V dans 0 à 2 % des cas, une tête en `K` dans 65 à 77 % ; le *premier*
résidu est indifférent. Dans l'archive, 48 % des séquences appelant IGHV3-23
portent `K` là, contre 3 % de celles appelant IGHV3-66.

Et les deux pièces n'ont de sens qu'**ensemble** : un donneur du bon gène V
(pour que la tête soit juste) *et* une tête que le PSSM ne peut pas toucher
(parce qu'il est ajusté sur tous les gènes V confondus et y remet un `K`).
Chacune seule échoue.

**Mais à CDR3 court il n'y a pas de place pour les appliquer.** Longueur 7 :
2 résidus codés par le V + 3 par le J (+1 de marge dans `J_TAIL`) = **une seule
position libre**, donc chaque candidat est son propre donneur et 3361 sont
refusés comme doublons d'archive. Les ~25 % de refus sur le V sont le prix de
la diversité à ces longueurs. Ce n'est pas un défaut de notre générateur, c'est
la géométrie du problème — et cela explique vraisemblablement pourquoi le champ
est le plus fourni aux CDR3 courts (26 soumissions médianes contre 9 aux
longueurs longues) : personne n'y peut différencier beaucoup sa séquence.

À CDR3 long, où les mêmes contraintes laissent 14 positions libres, la paire
de correctifs fait ce qu'elle promet — mesuré à longueur 20, `IGHV3-23` :

| variante | valides sur 1500 | refus V |
|---|---|---|
| témoin | 1353 (90 %) | 142 |
| tête protégée seule | 1362 (91 %) | 136 |
| donneurs du bon V seuls | 1411 (94 %) | 85 |
| **les deux ensemble** | **1441 (96 %)** | **54** |

**Non livré malgré tout.** +6 points de validité dans un étage qui coûte 145 s
et n'est pas le goulot (la PSH en prend 17,6 min), contre un comportement
conditionnel à la longueur dans deux modules. Le vivier de candidats valides
n'a jamais été ce qui nous limite aux longueurs longues. Le mécanisme est
documenté ici ; si un jour la génération devient le goulot, le correctif est à
reprendre tel quel.

### Et une affirmation du code qui était fausse

Le commentaire de `generate` disait depuis toujours que muter CDR1/CDR2 entame
la marge d'identité germinale. Chiffré pour la première fois : couper
**entièrement** ces mutations laisse les refus V à 28 % et la validité à 71,5 %
contre 69,9 % — un écart dans le bruit (± 1,2 point à n = 1500). L'affirmation
est fausse en pratique.

---

## 8. Ce qui reste ouvert

*Révisé au 2026-11-10, après la nuit du 9 au 10 (§9).*

**Fermé depuis :** le preflight calibré (36/36 sur les acceptées de l'époque
26174), le chronométrage par étage (§7bis), la persistance des scores BoltzGen,
la reprise sur HTTP 429 et la vérification du dépôt avant l'engagement du
pointeur (les trois par le propriétaire du dépôt), et la question du générateur
— `pssm` contre `pssm_denovo` ne se distinguent pas à n = 60 (p = 0,36 et 0,85).

**Ce qui demande une décision, pas une mesure :**

1. **Porter le seuil PSH par longueur** dans `preflight.py` et `mine.py`. La
   mesure (§3ter) dit qu'à CDR3 ≥ 17 un candidat lu ≤ 96 est admis 68 % du
   temps et un candidat lu > 96, 28 %. Le seuil plat de 89 nous fait refuser de
   soumettre dans un tiers du calendrier — celui où nous classons 1ᵉʳ contre 9
   rivaux. C'est un changement de politique de soumission, donc une décision du
   propriétaire, pas une correction de bug.
2. **Un essai à blanc avec BoltzGen**, qui demande un GPU loué. Tous les autres
   étages sont chronométrés ; BoltzGen est estimé à ~11 min pour 7 candidats
   sur 2×T4 et n'a jamais tourné dans une époque réelle.

**Ce que le générateur ne sait pas faire :** nos candidats sont à 98,1 de PSH
médiane contre 90,0 pour les acceptées. Le criblage exploite la queue de
distribution ; il ne rend pas le générateur meilleur. Et §7ter montre pourquoi
il ne peut pas beaucoup mieux faire aux CDR3 courts : 5 des 7 positions y sont
réclamées par les segments germinaux.

**Ce que le générateur ne sait pas faire :** nos candidats sont à 98,1 de PSH
médiane contre 90,0 pour les acceptées. Le criblage exploite la queue de
distribution ; il ne rend pas le générateur meilleur.

---

## 9. Bilan de la nuit du 9 au 10 novembre 2026

Huit heures de travail hors ligne, sans soumission, sans dépense et sans GPU.
Vingt-trois mesures, douze commits. Résumé par nature de résultat.

### Deux questions ouvertes depuis des jours, tranchées

| | résultat |
|---|---|
| **Le budget d'époque tient-il ?** | Oui, mesuré pour la première fois : 49,1 min sur 67 avant correctifs, **~32 min après**. Le mineur avait échoué trois fois sur `too_late` sans qu'on sache pourquoi. |
| **Le seuil PSH est-il à la bonne place ?** | Non : la rupture est à **96, pas 89**. À CDR3 ≥ 17, lire ≤ 96 vaut 68 % d'admission [57-77], lire > 96 vaut 28 % [18-42], Fisher p = 2×10⁻⁵. |

### Quatre défauts corrigés, dont trois auraient coûté des époques

1. **`propose` passait 28 de ses 31 minutes** à cribler l'archive sur 300
   candidats dont la PSH jette 97 % — et ce coût suivait la taille de
   l'archive, donc il empirait tout seul (il a grandi de 52 % dans la même
   nuit). Déplacé après la coupe PSH : 28 min → 44 s.
2. **Le chimérisme coupait le gène J en deux** : 58 % de la production perdue
   pour trois dixièmes de point d'identité. Appariement queue/FR4 : 38 % →
   100 % de passage, vivier valide 829 → 2163.
3. **`--once` ne terminait jamais** (fil websocket non démon) : un superviseur
   aurait attendu indéfiniment un mineur ayant fini son travail en 0,8 s.
4. **Le rebuild du vivier que le mineur conseillait lui-même** ne
   re-téléchargeait pas l'archive, donc ne corrigeait rien. Je suis tombé dans
   ce piège avant de lire le constructeur. Vivier : 63 époques de retard → 3.

### Ce qui a échoué, et pourquoi c'est aussi un résultat

**Quatre tentatives sur les refus du gène V, quatre rejets par la mesure**
(§7ter). Le mécanisme a fini par être trouvé et il est solide — la position 1
du CDR3 porte la signature du gène V. Mais aucune de ses corrections ne gagne :
à CDR3 court il n'y a pas de place (5 positions sur 7 sont germinales), et à
CDR3 long le gain est réel mais marginal (+6 points de validité) dans un étage
qui n'est pas le goulot. **Non livré.** Ce que l'enquête laisse n'est pas un
correctif mais une contrainte structurelle, qui explique pourquoi le champ est
le plus fourni là où personne ne peut différencier sa séquence.

### Trois erreurs de ma part, attrapées avant d'être publiées comme des faits

| ce que j'allais écrire | ce qui était vrai |
|---|---|
| « notre PSH sur-lit avec la longueur » | la sélection au seuil suffit à produire la même signature ; simulé et vérifié |
| un ajustement donnant σ_v et notre biais | modèle mal spécifié : un quart des refus n'est pas dû à la PSH, les paramètres ne sont pas identifiables. Retiré, remplacé par une courbe d'admission sans modèle |
| « muter CDR1/CDR2 entame la marge germinale » (dans le code depuis l'origine) | faux en pratique : couper ces mutations ne change rien (28 % de refus V dans les deux cas) |

Les trois relèvent de la même erreur de forme, déjà commise sur le σ de la PSH
et sur « le seuil PSH n'est pas la cause » : **un mécanisme établi ne dit pas
quelle part du résultat il explique.** Seule la mesure de bout en bout le dit,
et elle a dit non quatre fois cette nuit.

### Méthode qui a marché, à réutiliser

Croiser systématiquement les donneurs deux à deux et mesurer, plutôt que
raisonner sur le mécanisme. C'est ce croisement qui a donné le gène J (après
deux hypothèses fausses de ma part) et la position 1 du CDR3 (après trois).
Et valider chaque contrainte trouvée contre le champ réel accepté : 69 % des
2536 acceptées portent une des quatre queues J compatibles, ce qui a confirmé
que la contrainte était réelle et non un artefact de notre chaîne.

---

## 10. Le champ s'est amélioré pendant que nous réparions le mineur

Mesuré le 2026-11-10 sur le run Kaggle de l'étage BoltzGen (8 candidats réels
de l'époque 26238, 2×T4) confronté aux 970 soumissions notées des 35 dernières
époques.

### D'abord, le chiffre qui manquait au budget

| étape de BoltzGen | temps |
|---|---|
| design | 392,9 s |
| folding | 547,5 s |
| analysis | 36,1 s |
| **total, 8 séquences** | **20,0 min = 150 s/séquence** sur 2×T4 |

**150 s, pas les 96 que j'estimais.** Budget d'époque complet : propose ~3 min
+ PSH 3,6 min + archive 0,7 + développabilité 0,6 + BoltzGen 20 = **~28 min
sur 67**. Le budget tient, avec 39 min de marge.

### Ensuite, le résultat qui compte

Le champ réel a changé de régime entre les époques 26025 et 26170 :

| | octobre (25956-26025) | maintenant (26200+) |
|---|---|---|
| `design_iiptm` médian | **0,239** | **0,664** |
| 99ᵉ percentile | 0,630 | 0,778 |
| maximum | 0,704 | 0,786 |
| part au-dessus de 0,6 | **2 %** | **75 %** |
| meilleur par époque (médiane) | 0,543 | 0,757 |
| taux de refus | 52 % | 33 % |

**Notre meilleur candidat : 0,375.** Soit dans les **12 % les plus faibles** du
champ actuel, et **0 des 35 dernières époques** aurait été gagnée avec cette
valeur.

### Deux explications écartées par la mesure

1. **Changement d'échelle de la notation.** Il aurait déplacé le maximum comme
   la médiane. Le maximum n'a bougé que de 12 % (0,704 → 0,786) quand la
   médiane triplait. De plus le dépôt amont n'a pas changé depuis le 26 octobre
   (db9bf19, `HEAD..origin/main` = 0), donc le code de notation est le même.
2. **Survie** — un durcissement des filtres éliminerait les soumissions faibles
   et relèverait la médiane des acceptées sans toucher le plafond, exactement
   la signature observée. Mais le taux de refus a **baissé**, de 52 % à 33 %.
   Les filtres ne se sont pas durcis.

Il reste donc que **le champ s'est réellement amélioré** : en octobre 2 % des
soumissions dépassaient 0,6, aujourd'hui 75 %. Le plafond était déjà atteignable
en octobre (quelqu'un faisait 0,704) ; maintenant presque tout le monde
l'atteint.

### La bascule datée : une courbe d'adoption, pas un palier

Les 141 époques manquantes (26026-26169) ont été récupérées. Le champ noté
couvre désormais **5833 soumissions sur 260 époques sans trou**, et la
transition est **graduelle** :

| bloc de 20 époques | jours avant 26234 | part > 0,6 | médiane |
|---|---|---|---|
| 25960 | 13,7 | 0,4 % | 0,220 |
| 26000 | 11,7 | 1,8 % | 0,241 |
| 26020 | 10,7 | 9,2 % | 0,301 |
| 26040 | 9,7 | 23 % | 0,468 |
| 26080 | 7,7 | 35 % | 0,527 |
| 26120 | 5,7 | 58 % | 0,623 |
| 26180 | 2,7 | 67 % | 0,656 |
| 26200 | 1,7 | **82 %** | 0,683 |

**C'est décisif contre l'hypothèse d'un changement de notation**, qui serait un
palier net à une date unique et identique pour tous. Une diffusion sur 14 jours
ne peut être qu'une méthode qui se propage.

L'ordre des adoptants le confirme — première époque au-dessus de 0,6 par dépôt :

| dépôt | première > 0,6 | soumissions récentes |
|---|---|---|
| `frankji-groundcontrol` | 25969 | 29 |
| `QTTECH2023` | 26004 | **248** |
| `cczulaika` | 26019 | 3 |
| `saon-miner` | 26028 | **168** |
| `bitty-labs` (9 dépôts) | 26040-26076 | 130 au total |
| `thomasvangurp` | 26040 | 22 |

Les neuf dépôts de `bitty-labs` basculent en huit époques : un seul opérateur
qui déploie. **Quinze dépôts ont acquis la méthode en deux semaines, donc elle
est copiable** — ce n'est pas un avantage propriétaire. C'est cohérent avec
BoltzGen-1, qui est public sur Hugging Face et dont les exemples de conception
de nanobodies sont dans le dépôt des validateurs lui-même (§11).

### Ce que cela fait aux estimations de ce document

**Les taux de victoire de 14-29 % (§7) puis 10 % (§3bis) ont été mesurés contre
le champ d'octobre. Ce champ n'existe plus.** Rejoué contre les époques
récentes, notre meilleur candidat se place 17ᵉ à 33ᵉ sur des champs de 22 à 37,
là où il se plaçait 1ᵉʳ dans 7 des 18 époques d'octobre à CDR3 court.

L'écart n'est pas un écart de réglage : greffer des CDR3 archivés sur des
charpentes archivées produit 0,375, et il faut 0,75 pour gagner. Aucun des
correctifs de cette nuit — appariement J, seuil PSH recalibré, budget d'époque
— ne touche à cela. Ils ont rendu le mineur capable de soumettre ; ils ne
rendent pas la séquence compétitive.

**Conséquence directe : ne pas engager de TAO.** L'espérance de gain mesurée
contre le champ actuel est de zéro victoire sur 35 époques.
