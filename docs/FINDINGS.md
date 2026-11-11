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

## 0. Où nous en sommes, et où est le problème

*Ajouté le 2026-11-10. Les sections 1 à 10 ont été écrites avant ce qui suit et
restent exactes sur leur objet, mais la section 1 n'annonce plus « le fait
central » : l'admission a été comprise et corrigée, et ce qui bloque
aujourd'hui est ailleurs.*

### Ce qu'il faut pour gagner, et ce que nous produisons

| | `design_iiptm` |
|---|---|
| requis pour gagner une époque (CDR3 = 16) | **0,777** |
| médiane du champ à cette longueur | 0,697 |
| **notre meilleur mesuré** | **0,673** |
| notre production en sélection PSH (production réelle) | 0,375 |

### Les cinq mesures qui situent le problème (§15)

| candidat | verdict |
|---|---|
| **le GPU** | **non** — notre distribution est **bornée à 0,674** (0/148 au-dessus de 0,68 ; queue de forme c = −0,82). Noter plus ne trouve rien de plus. |
| **les donneurs** | **non** — ils scorent 0,729-0,794 et viennent tous du régime actuel |
| **la distance à l'archive** | **non** — ρ = −0,07, p = 0,4 |
| **trop peu de surface de conception** | **non** — trois CDR conçus font pire qu'un (0,169 contre 0,214) |
| **la recombinaison** | **oui** — charpente d'une séquence à 0,745 + CDR3 d'une à 0,729 → 0,48. Les pièces ne se transfèrent pas. |

### Ce qui a été réparé, et ce que ça valait

| | avant | après |
|---|---|---|
| sélection des candidats (§12, §13) | PSH, **0,375** | substitut `design_iiptm`, **0,616** |
| budget d'époque (§7bis) | 49 min sur 67 | ~28 min |
| vivier de donneurs (§5ter) | 63 époques de retard | 3 |
| identité J du générateur (§7ter) | 58 % de refus | 0 |

La sélection valait **+0,241** et était gratuite. Elle nous fait passer des 14 %
les plus faibles au **38ᵉ percentile** du champ — et **0 époque gagnée sur 8**,
comme avant.

### La leçon de méthode, qui a tenu à chaque fois

**Un mécanisme établi ne dit pas quelle part du résultat il explique.** Quatre
mécanismes réels et bien mesurés n'ont rien rapporté ; une simulation à moi a
rendu une valeur impossible (1,042 pour une métrique bornée par 1) ; un plateau
prédit à 0,756 s'est révélé sur-promis de 0,137 hors de sa plage de calibration.
Seule la mesure de bout en bout a tranché, à chaque fois.

---

## 1. L'admission : le premier obstacle, compris et corrigé

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

---

## 11. La conception BoltzGen contrainte : testée, et elle ne transfère pas

Lancée sur Kaggle (2×T4), 20 designs, pipeline complet en 5 étapes sans échec :
design 775,7 s + inverse_folding 33,1 s + folding 1099,6 s + analysis 39,4 s +
filtering 15,0 s ≈ **32 min**. Données : `results/boltzgen_design_attempt1.csv`.

### La mécanique a marché

| ce qui était demandé | obtenu |
|---|---|
| CDR3 de longueur 16 | **20/20** exactement 16 |
| tête `AR` figée (gène V) | **20/20** |
| identité V ≥ 82 | **20/20**, médiane 87,8 |

L'arithmétique `exclude` + `design_insertions` était donc juste, et l'indexation
aussi. L'outillage est bon.

### Le score, lui, est pire que le greffage

| | `design_iiptm` |
|---|---|
| **designs BoltzGen contraints** | **médiane 0,214, max 0,311** |
| notre greffage | 0,375 |
| médiane du champ réel | 0,664 |
| meilleur par époque du champ | 0,757 |

Aucun design au-dessus de 0,664. **La méthode telle que je l'ai contrainte ne
transfère pas** — elle fait moins bien que ce qu'elle devait remplacer.

### Le mécanisme de l'échec est identifié

Le filtrage interne de BoltzGen a écarté les 20 designs sur un seul critère :

```
filter_rmsd (seuil 2,5) : 0/20 passent   — mediane 9,23 A
filter_rmsd_design      : 7/20 passent   — mediane 2,66 A
```

L'écart entre les deux est l'explication : le CDR3 conçu se replie **localement**
comme prévu (`bb_rmsd_design` 1,04 à 2,68 Å), mais **le nanobody ne se repose pas
là où la conception l'avait placé** (`bb_rmsd` 3,6 à 9,5 Å). La pose de
l'interface n'est pas reproduite au repliement indépendant — et c'est exactement
ce que `design_iiptm` mesure, d'où 0,21.

### Et une hypothèse de ma part qui était fausse

Je croyais que figer une région dans la spécification figeait sa **séquence**.
C'est faux : l'étape `inverse_folding` réassigne la séquence de toute la chaîne.

- les queues de CDR3 obtenues sont `ADY, FAY, FDY, FEY, FPY, LDY, MDY, VDY` —
  **pas** les `FDP`/`FDH` des donneurs que j'avais déclarés fixes ;
- l'identité V monte à 87,8 contre 84-86 chez les donneurs : le modèle a rendu
  la charpente *plus* germinale qu'elle n'était ;
- conséquence directe : **seuls 6/20 passent le plancher J de 86** (médiane 85,7),
  le défaut même que §7ter avait corrigé côté générateur.

Le masque de conception contraint la **structure**, pas la séquence finale.

### Deux pistes, par ordre de coût

1. **Trop peu de liberté.** Nous n'avons laissé que 11 résidus conçus, en gelant
   CDR1, CDR2 et les deux extrémités du CDR3. L'exemple amont
   `nanobody_scaffolds/7eow.yaml` conçoit **les trois CDR** (26..34, 52..59,
   98..118). Avec une seule boucle courte, le modèle ne peut pas construire
   d'interface. Et §7ter a mesuré que muter CDR1/CDR2 ne coûte rien de
   mesurable sur l'identité germinale — donc les geler coûtait sans rapporter.
2. **Vingt designs, quand le CLI dit « typiquement 10 000 ».** Le filtrage est
   conçu pour tailler dans des milliers ; le meilleur de 20 n'est pas le
   meilleur de 10 000. À 1,6 min/design, 200 designs font ~5 h de Kaggle.

La piste 1 est à tester d'abord : elle ne coûte pas plus cher et elle corrige en
même temps l'erreur sur `inverse_folding` (si la séquence est réassignée de toute
façon, geler CDR1/CDR2 n'avait aucun effet protecteur).

---

## 12. L'erreur n'était pas dans ce que nous fabriquions, mais dans ce que nous choisissions

### Le fait qui retourne le diagnostic

Nos candidats greffés, notés par BoltzGen à l'échelle
(`results/scores_nanobody_observation_*.csv`) :

| lot | n | médiane | **max** |
|---|---|---|---|
| combined | 40 | 0,404 | 0,599 |
| cutoff | 42 | 0,326 | 0,662 |
| depth4 | 84 | 0,343 | **0,692** |

**0,692 sur 84 tirages, au-dessus de la médiane du champ (0,664).** Or les 8
candidats de l'essai à blanc de l'époque 26238 plafonnaient à **0,375**. La
différence n'est pas le générateur : c'est que la coupe de milliers à huit se
faisait sur la **PSH**, dont le §3ter avait déjà établi qu'elle ne prédit pas le
verdict sous 96. **Nous jetions les bons candidats avant de les regarder.**

### Le substitut `design_iiptm`

5956 séquences étiquetées (champ entier + nos candidats notés), 147 traits.
Validé trois fois, par ordre de confiance décroissante :

| validation | ρ |
|---|---|
| 5 plis aléatoires sur tout | 0,913 — **gonflé** |
| 5 plis groupés **par époque** | 0,882 (les quasi-doublons ne valent que 0,03) |
| entraîné sur le champ, **classant nos candidats** | **0,569** (p = 6×10⁻¹², n = 124) |

Le troisième est celui à citer. Un modèle qui distinguerait seulement
« séquence du champ de novembre » (0,66) de « greffe » (0,34) réussirait les
deux premiers et serait inutile sur notre vivier — l'erreur de forme commise
trois fois ailleurs dans ce document.

Effet sur la sélection, sur nos 124 candidats notés :

| 8 choisis | `design_iiptm` réel |
|---|---|
| **par le substitut** | **médiane 0,572, max 0,649** |
| au hasard | 0,298 |
| le vivier entier | 0,335, max 0,667 |

### Le coût, mesuré

| | |
|---|---|
| classer 1751 candidats par le substitut | **27 s** = 15,6 ms/candidat |
| mesurer la PSH des mêmes | **105 min** |

C'est ce rapport qui permet de cribler le vivier **entier** au lieu de 31
candidats. Pipeline recâblé en conséquence : générer → classer tout par le
substitut → **admissibilité** PSH sur les meilleurs → archive →
développabilité → BoltzGen. La PSH redevient un seuil, plus un critère de choix.

### Ce que cela ne règle pas

0,649 reste sous le meilleur par époque du champ (0,757). Ce correctif nous fait
passer des 12 % les plus faibles à environ la médiane. **Nécessaire, pas
suffisant.**

Deux réserves sur le substitut lui-même :

- il est entraîné sur un champ qui vient de changer de régime (82 % au-dessus de
  0,6 contre 2 % en octobre), donc il est **périssable** : il faut relancer
  `fetch_submissions` + `analyze_nanobodies` + l'entraînement régulièrement ;
- ρ = 0,569 veut dire qu'il se trompe souvent entre candidats proches. Ce qu'il
  fait bien, c'est écarter le bas de la distribution — assez pour ne plus
  soumettre un 0,34 quand le vivier contient un 0,67, pas assez pour garantir le
  meilleur des 1751.

### Les trois tentatives écartées avant celle-ci

Toutes cherchaient à améliorer la **fabrication** quand le défaut était dans la
**sélection** :

| tentative | résultat mesuré |
|---|---|
| conception BoltzGen contrainte (§11) | 0,214 — pire que le greffage |
| criblage des donneurs de CDR3 par gène V (§7ter) | 1751 valides contre 2163 |
| protection de la tête du CDR3 (§7ter) | **0** valide |

C'est la leçon de méthode de ce projet : **un mécanisme établi ne dit pas quelle
part du résultat il explique.** Les trois mécanismes étaient réels et mesurés.
Aucun n'était le goulot.

### 12bis. Une simulation à moi, invalide, et ce qui la remplace

Pour chiffrer le gain du criblage large j'ai simulé la malédiction du
vainqueur : vivier de taille n, prédictions tirées d'une gaussienne ajustée sur
nos 208 points, valeur réelle = régression + bruit. Résultat annoncé à
n = 40 000 : **1,042**.

`design_iiptm` est borné par 1,0 et le maximum jamais observé sur 6000
soumissions du champ est 0,794. Le chiffre est impossible, et la cause est la
modélisation : **une gaussienne a des queues non bornées**, donc le maximum de
n tirages croît sans limite avec n. La simulation ne mesurait rien.

Remplacée par une lecture non paramétrique, qui n'extrapole pas au-delà du
mesuré — pour chaque tranche de prédiction, les valeurs réelles observées :

| prédit | n | réel médian | réel max |
|---|---|---|---|
| 0,12–0,41 | 104 | 0,243 | 0,438 |
| 0,41–0,51 | 41 | 0,491 | 0,602 |
| 0,51–0,57 | 31 | 0,553 | 0,667 |
| 0,57–0,62 | 21 | 0,566 | 0,662 |
| **0,62–0,67** | **11** | **0,643** | **0,711** |

Deux choses utiles en sortent :

1. **le substitut est bien calibré** dans la plage mesurée — prédire 0,64 donne
   un réel médian de 0,643, sans biais appréciable ;
2. **sa tranche haute contient du 0,711 réel**, proche du meilleur par époque du
   champ (0,757).

Mais la prédiction la plus haute jamais observée sur nos candidats est 0,667.
**Ce que vaut réellement un candidat prédit au-dessus de ça ne s'extrapole pas :
il faut le noter.** C'est précisément ce que l'étape Kaggle doit mesurer, et
c'est la seule façon de savoir si le criblage large nous met au niveau du champ.

---

## 13. La sélection est réparée ; le plafond du générateur, non

24 candidats prédits 0,726–0,752, notés par le BoltzGen des validateurs sur
Kaggle (40,6 min, 1,7 min/séquence). `results/scores_surrogate_top24.csv`.

| | `design_iiptm` |
|---|---|
| prédit par le substitut | 0,726 – 0,752 |
| **réel** | **médiane 0,616, max 0,673, min 0,296** |
| biais | **le substitut sur-promet de +0,137** |
| ρ prédit/réel **dans** la tranche | **+0,31** |

### Ce qui est gagné, et c'est réel

| sélection | `design_iiptm` obtenu |
|---|---|
| par la PSH (essai de l'époque 26238) | 0,375 |
| **par le substitut** | **médiane 0,616** |

**+0,241.** Le correctif de sélection fonctionne, et il est gratuit.

### Ce qui ne suffit pas, et c'est décisif

Comparé au champ **à CDR3 = 16** (200 soumissions, 8 époques récentes) :

| | |
|---|---|
| médiane du champ | 0,697 |
| 90ᵉ percentile | 0,771 |
| meilleur par époque (médiane) | **0,777** |
| **notre max : 0,673** | **38ᵉ percentile** |
| époques gagnables avec 0,673 | **0 / 8** |
| époques gagnables avec l'ancien 0,375 | 0 / 8 |

Nous passons des 14 % les plus faibles au 38ᵉ percentile. **Toujours sous la
médiane du champ, et zéro époque gagnée.**

### Trois choses que ce lot établit

1. **Le plateau prédit à 0,756 était une illusion.** Le substitut sur-promet de
   0,137 dans sa tranche haute, hors de la plage où il était calibré (≤ 0,667).
   Refuser d'extrapoler était justifié — le §12bis avait raison de retirer la
   simulation.
2. **Dans la tranche haute, le substitut ne sait plus rien** (ρ = +0,31 contre
   0,569 globalement). Affiner la sélection *à l'intérieur* de cette tranche ne
   rapportera rien : il n'y a plus d'information à exploiter.
3. **Le vrai plafond du greffage est ~0,67**, pas 0,756. Il faut 0,777 pour
   gagner une époque à cette longueur. **L'écart restant, ~0,10, n'est plus un
   écart de sélection — c'est un écart de qualité du générateur.**

### Ce que cela ordonne pour la suite

La sélection est au bout de ce qu'elle peut donner. Le seul levier restant est
le **générateur**, et la tentative de conception BoltzGen du §11 a échoué pour
une raison identifiée : nous ne lui avions laissé que **11 résidus** en gelant
CDR1, CDR2 et les deux extrémités du CDR3, quand l'exemple amont conçoit **les
trois CDR**. Ce gel ne protégeait rien — `inverse_folding` réassigne la séquence
de toute façon (§11), et §7ter avait mesuré que muter CDR1/CDR2 ne coûte rien en
identité germinale.

**C'est maintenant la seule expérience qui vaille**, et son coût est une session
Kaggle. Les 24 mesures de ce lot rejoignent l'entraînement du substitut : elles
sont 24 observations neuves dans la bande où il n'en avait que 11, et où il vient
de se tromper de 0,137.

---

## 14. Trois CDR conçus au lieu d'un : réfutation nette

60 designs, pipeline complet sur 2×T4 (design 2076 s + inverse_folding 59 s +
folding 3023 s + analysis 61 s + filtering 18 s ≈ **1 h 25**).
`results/boltzgen_design_3cdr.csv`.

27 résidus conçus au lieu de 11 — CDR1 + CDR2 + milieu du CDR3, la forme de
l'exemple amont.

| | `design_iiptm` |
|---|---|
| conception **1 CDR** (§11) | médiane 0,214, max 0,311 |
| conception **3 CDR** | **médiane 0,169, max 0,367** |
| greffage + sélection par substitut (§13) | **0,616, max 0,673** |
| champ à CDR3 = 16 | 0,697, meilleur par époque 0,777 |

Part des designs au-dessus de notre propre greffage sélectionné (0,616) : **0 %**.

### L'hypothèse causale du §11 était fausse, et c'est le résultat utile

Le §11 attribuait le faible score au défaut de pose : `filter_rmsd` médian à
9,23 Å pour un seuil de 2,5, « le CDR3 se replie comme prévu mais le nanobody
ne se pose pas où la conception l'a mis ». Avec trois CDR la pose s'améliore
réellement — **5/60 passent `filter_rmsd` contre 0/20**, et
`filter_rmsd_design` 24/60 contre 7/20.

**Mais ces 5 designs à pose correcte scorent 0,200 à 0,356.** Corriger la pose
ne remonte pas `design_iiptm` d'un point. Les deux grandeurs sont indépendantes,
et l'explication causale du §11 — comme l'hypothèse « trop peu de surface » qui
en découlait — est réfutée par la mesure.

### Ce que cela dit de la méthode des meneurs

**La conception BoltzGen contrainte à une charpente de nanobody ne produit pas
de `design_iiptm` élevé, quelle que soit la surface qu'on lui laisse.** La
prémisse du §11 — « les meneurs conçoivent avec le modèle qui note » — n'est
donc pas vérifiée dans cette forme. Deux hypothèses survivent, et aucune n'est
testable dans nos moyens actuels :

1. **conception libre à grande échelle** : sans charpente imposée, puis filtrage
   sur la conformité germinale. Cohérent avec le `--num_designs` « typiquement
   10 000 » du CLI — mais 10 000 designs font ~10 jours de Kaggle ;
2. **notation massive** : générer largement et noter des centaines de candidats
   par époque avec de vrais GPU, ce que notre substitut approxime gratuitement
   et imparfaitement.

### État honnête de la compétitivité

| approche | mesurée | gagne une époque ? |
|---|---|---|
| greffage + sélection PSH (production) | 0,375 | 0/8 |
| **greffage + substitut** | **0,616 / max 0,673** | **0/8** |
| conception 1 CDR | 0,214 | non |
| conception 3 CDR | 0,169 | non |
| **ce qu'il faut** | **0,777** | — |

**Nous ne sommes pas compétitifs, et je n'ai pas de chemin mesuré pour le
devenir.** Le correctif de sélection valait +0,241 et était gratuit ; il nous
met au 38ᵉ percentile du champ. Les ~0,10 restants demandent une qualité de
générateur que ni le greffage ni la conception contrainte n'atteignent.

---

## 15. Où est exactement le problème

Question posée après le §14 : le goulot est-il le GPU ? Cinq mesures le situent,
et aucune ne le met là où je l'avais cherché.

### 1. Ce n'est pas le GPU

Nos 148 candidats notés par BoltzGen :

| seuil | nos candidats |
|---|---|
| > 0,650 | 7 / 148 |
| > 0,680 | **0 / 148** |
| > 0,777 (requis) | **0 / 148** |

Queue ajustée par une Pareto généralisée au-delà du 80ᵉ percentile : forme
**c = −0,82**, fortement négative, donc **bornée**, borne implicite **0,674**.

Le GPU serait le goulot si la queue atteignait 0,777 *rarement* : il faudrait
alors du débit pour tirer assez. Mais sept candidats s'agglutinent entre 0,65 et
0,673 puis plus rien — signature d'un plafond, pas d'une queue qui s'amincit.
**Noter dix mille candidats d'une distribution bornée à 0,674 donne 0,674.**

### 2. Ce ne sont pas les donneurs

| | `design_iiptm` |
|---|---|
| nos charpentes donneuses (40) | médiane **0,745**, max 0,765 |
| nos donneurs de CDR3 (60) | médiane **0,729**, max 0,794 |
| **notre production** | **médiane 0,48, max 0,674** |

100 % des donneurs viennent du nouveau régime (époque ≥ 26170). Les choisir par
`design_iiptm` au lieu de `rank_sum` ne change rien (0,745 contre 0,746).

### 3. Ce n'est pas la distance à l'archive

Spearman entre l'identité du CDR3 au plus proche CDR3 archivé et le
`design_iiptm` : **ρ = −0,07, p = 0,4**. La bande la plus proche autorisée
(0,80-0,85) score comme la plus lointaine (0,500 contre 0,490).

### 4. C'est la recombinaison

Nous prenons la charpente d'une séquence à 0,745 et le CDR3 d'une séquence à
0,729, et nous produisons 0,48. **Les pièces ne se transfèrent pas** : la
liaison d'un nanobody dépend de l'appariement *spécifique* entre sa charpente et
son CDR3 — la boucle est présentée par la charpente qui la porte. Mélanger deux
bons parents donne une chimère qui lie moins bien que chacun d'eux.

Cela explique aussi l'échec de la conception (§11, §14) : concevoir un CDR3 sur
une charpente qui ne l'a jamais porté est le même problème, et concevoir les
trois CDR l'aggrave puisque plus rien n'est apparié.

### 5. Et la règle rend le correctif évident illégal

| règle | ce qu'elle impose |
|---|---|
| identité séquence entière ≥ 0,95 refusée | ≥ **7 résidus** différents sur 123 |
| identité CDR3 ≥ 0,85 refusée | ≥ 3 résidus sur 16 |

Les deux s'appliquent ensemble. Placer ces 7 différences dans la charpente
entame l'identité germinale ; les placer dans le CDR3 en réécrit **44 %** — soit
exactement là où la liaison se décide.

**La structure du piège :** de bonnes séquences existent dans l'archive, nous
n'avons pas le droit de soumettre leurs quasi-copies, et la perturbation que la
règle exige tombe sur le site de liaison.

### Ce qui n'a jamais été essayé

Notre générateur n'a **jamais** tenté la perturbation minimale d'**un seul**
parent : il mélange toujours deux parents. Ce sont deux choses différentes —
mélanger casse l'appariement charpente/CDR3 *et* perturbe ; perturber
minimalement conserve l'appariement et ne change que 7 résidus.

C'est la cinquième hypothèse de cette enquête et les quatre précédentes ont été
réfutées. Mais elle est la seule qui reste à coût d'une session, et elle est
dimensionnable : 8 parents à ~0,78, leurs variantes à 7 résidus près, notées.

- variantes à ~0,75 → notre générateur se trompait de stratégie, pas de moyens ;
- variantes à ~0,60 → la règle nous plafonne structurellement, et seule une
  vraie conception monte plus haut — ce que nous avons mesuré ne pas savoir faire.

---

## 16. Perturbation minimale : la PSH chiffre le prix

61 variantes minimales d'un seul parent, criblées hors GPU contre leurs 24
parents **acceptés par les validateurs**, dans la même chaîne.

| | n | PSH médiane | > 96 (rupture) | min |
|---|---|---|---|---|
| **parents acceptés** | 24 | **89,4** | 25 % | 82,3 |
| v1, substitutions uniformes sur 20 aa | 25 | 100,6 | 72 % | 88,4 |
| **v2, substitutions hydrophiles** | 52 | **97,1** | **60 %** | 84,5 |

Les parents lisant 89,4 chez nous, la chaîne est calibrée à leur niveau : les
écarts sont des dégradations réelles, pas la sur-lecture du §3ter.

### Décomposition du coût

| | points de PSH |
|---|---|
| mon tirage uniforme (bug) | **+3,5** |
| perturber, quel que soit le résidu | **+7,7** |
| total v1 | +11,2 |

Le tirage uniforme sur les 20 acides aminés mettait un résidu hydrophobe dans
la boucle deux fois sur cinq, alors que la PSH *est* le score de patch
hydrophobe de la surface des CDR. Corrigé (`DEKNQRSTH`), cela récupère 3,5
points. Les 7,7 autres sont structurels.

### Ce que la stratégie vaut encore

**40 % des v2 passent sous la rupture de 96** (minimum 84,5), soit une
vingtaine de variantes soumissibles sur 52. La stratégie n'est donc pas morte :
elle a un rendement, et ce rendement est mesurable hors GPU gratuitement.

Reste la question qui coûte une session : **le `design_iiptm` survit-il à la
perturbation ?** La PSH dit que la développabilité se dégrade mais reste
franchissable pour une variante sur deux. Elle ne dit rien de la liaison.

### Le point de méthode, pour la troisième fois de la journée

Sans les parents mesurés **comme témoin dans la même chaîne**, une PSH de 100,6
se lisait comme une propriété de la stratégie. J'aurais abandonné la piste pour
un défaut de trois caractères dans mon générateur. C'est la même discipline qui
a évité une session Kaggle sur 24 séquences toutes rejetées (§15, niveau 3 de la
règle de doublon) et qui a fait retirer la simulation du §12bis.

## 17. Le reste de la portée, gratuitement : A3D, V, J — et le lot GPU prêt

Les 61 variantes v2 ont été appariées à leur parent sur **tout ce qui se mesure
sans GPU**. La PSH venait des journaux de preflight ; l'Aggrescan3D a été
*récupérée* du cache de résidus (`data/preflight/run/a3d{1,2,3}/residues/<sha256>/
residue-scores.csv`) au lieu d'être recalculée — la moyenne et le pic sont la
moyenne et le maximum de la colonne de score, vérifiés à sept décimales contre
`summary/aggrescan3d.csv`. Le fichier de résumé n'avait gardé que le dernier
lot ; le cache de résidus, lui, garde les 151 séquences. L'identité V/J est
recalculée par ANARCI, ce qui est gratuit.
(`scripts/minimal_variant_deltas.py`, sortie dans `data/measurements/minimal2_deltas.csv`.)

| groupe | n | soumissibles | ΔPSH méd. | ΔA3D moy. | ΔA3D pic | ΔV | ΔJ |
|---|---|---|---|---|---|---|---|
| IGHV3-74/16 | 13 | 4 | +6,0 | −0,030 | +0,014 | +2,0 | 0,0 |
| IGHV3-74/18 | 15 | 8 | +4,5 | −0,086 | −0,024 | +1,0 | 0,0 |
| IGHV3-66/12 | 17 | 4 | +4,4 | −0,012 | +0,009 | 0,0 | 0,0 |
| IGHV3-64/18 | 16 | 4 | +13,1 | −0,016 | +0,003 | +0,5 | 0,0 |

**La PSH est le seul prix.** ΔPSH apparié : médiane +5,7, 39 variantes sur 52
dégradées, pire +24,8, meilleure −11,1. En face :

- **A3D : neutre ou favorable.** Les médianes de ΔA3D moyenne sont toutes
  négatives, c'est-à-dire *plus loin* du plafond (−0,680). Neuf variantes ont
  quand même été refusées, dont six sur A3D — mais à la marge, pas par
  effondrement. L'A3D n'est pas le goulot de cette stratégie.
- **J : rigoureusement inchangée.** ΔJ médian de 0,0 dans les quatre groupes.
  C'était attendu — la perturbation ne touche que la tête du CDR3 et deux
  positions de charpente, jamais les 3 derniers résidus ni FR4 — mais c'est
  maintenant mesuré, après deux hypothèses fausses sur le segment J (§7ter).
- **V : légèrement meilleure** (+0 à +2 points), toutes au-dessus du plancher
  de 82 (82,5 à 86,7).

**20 variantes sur 61 franchissent toutes les portes** (PSH < 96 et A3D sous les
deux plafonds), issues de 14 parents distincts.

### Le lot GPU

`inputs/minimal2_gpu.csv` : 34 séquences = les 20 variantes **et leurs 14
parents**, rangées parent d'abord. Les parents ne sont pas du remplissage :
chacun porte son `design_iiptm` archivé, et **s'ils ne re-notent pas près de
cette valeur, la session a mesuré la session et pas les séquences**. ~52 min de
GPU à 1,5 min/séquence. Carnet : `kaggle_minimal2_gpu.ipynb`.

Ce que le lot tranche, et c'est pour cela qu'il valait la peine :

- variantes près de **0,78** → notre générateur avait la mauvaise *stratégie*,
  pas les mauvais moyens : partir d'un candidat qui lie et le perturber au
  minimum légal préserve la liaison, et il reste à industrialiser le rendement
  de 33 %.
- variantes près de **0,60** → la règle d'unicité nous plafonne
  structurellement : tout ce qui est assez loin de l'archive pour être accepté
  est trop loin du parent pour lier. Ce serait la réponse la plus utile des
  deux, parce qu'elle fermerait définitivement la piste de la recombinaison
  (§13) comme celle de la perturbation.

### `scripts/verify_boltzgen_notebook.py`

Trois sessions Kaggle ont été perdues sur des défauts qu'un contrôle local
attrapait (§14). Ils sont maintenant des assertions : chaque cellule s'analyse,
`PAYLOAD` est défini avant son usage, la charge utile se déballe vraiment, le
lot n'a pas de doublon, **chaque variante a son parent dans le lot**, les gardes
« pas de GPU » et « run déjà existant » sont présentes, et les enfants sont
forcés à `MPLBACKEND=Agg` sous le backend de Kaggle. Le contrôle a été testé
par mutation : garde GPU retirée, cellule PAYLOAD déplacée, parent retiré du
lot — il refuse les trois. Il ne lance rien et n'a pas besoin de GPU.

## 18. Bilan de la nuit du 10 au 11 novembre

Cinq commits, aucune dépense, aucun GPU loué, aucune soumission. Rien n'a été
lancé sur Kaggle : seul l'utilisateur le fait.

### Ce qui a été mesuré

1. **Le niveau 3 de la règle de doublon lie** (§15). Mon arithmétique n'avait
   compté que deux niveaux sur trois : il faut 4 substitutions dans le CDR3, pas
   3. Les 24 premières variantes étaient toutes irrecevables. Attrapé hors GPU.
2. **La PSH est le seul prix de la perturbation minimale** (§16, §17). +11,2
   points en v1, dont **+3,5 de mon propre défaut** (tirage uniforme sur les 20
   acides aminés, donc un résidu hydrophobe dans la boucle deux fois sur cinq,
   là où la PSH *est* le score de patch hydrophobe des CDR) et **+7,7
   structurels**, le prix de perturber quel que soit le résidu. L'A3D, le V et
   le J ne bougent pas — le ΔJ médian est de 0,0 dans les quatre groupes.
3. **Le rendement : 20 variantes sur 61** franchissent toutes les portes
   (PSH < 96, A3D sous les deux plafonds), issues de 14 parents distincts.
   Mesuré gratuitement, hors GPU.
4. **L'A3D était récupérable du cache** plutôt que recalculable : des heures de
   re-pliage évitées en lisant `residues/<sha256>/residue-scores.csv` (§17).

### Ce qui a échoué, et pourquoi

- **Mon générateur de variantes injectait du hydrophobe au hasard.** Trois
  caractères (`rng.choice(AAS)` → `HYDROPHILIC`). Sans les parents mesurés
  comme témoin dans la même chaîne, une PSH de 100,6 se serait lue comme une
  propriété de la stratégie et j'aurais abandonné la piste.
- **La vérification du CDR3 par deux niveaux au lieu de trois.** Même cause
  racine que le reste : avoir raisonné au lieu de relire le code qui juge.
- **Le premier criblage a été tué à la limite de temps**, la vérification
  d'archive étant dans la boucle de réessai — le défaut exact que le §7bis
  avait corrigé 40 minutes plus tôt. Déplacé hors de la boucle.

### L'état du lot, prêt à lancer

`kaggle_minimal2_gpu.ipynb` + `inputs/minimal2_gpu.csv` : **34 séquences**, les
20 variantes et leurs 14 parents, parent d'abord. **Réglages Kaggle :
*Accelerator* GPU T4 x2, *Internet* On**, puis *Save Version → Save & Run All
(Commit)*. ~52 min de GPU à 1,5 min/séquence. `TOP = 2` note une paire
parent/variante pour un test rapide. Vérifié en local par
`scripts/verify_boltzgen_notebook.py`, lui-même testé par mutation.

### Ce que le lot tranchera

- **Variantes autour de 0,78**, c'est-à-dire au niveau de leurs parents → notre
  générateur se trompait de *stratégie* et non de moyens. Partir d'un candidat
  qui lie et le perturber au minimum légal préserve la liaison ; il reste à
  industrialiser un rendement de 33 %.
- **Variantes autour de 0,60** → la règle d'unicité nous plafonne
  structurellement : ce qui est assez loin de l'archive pour être accepté est
  trop loin du parent pour lier. Il ne resterait alors de chemin que par une
  vraie conception *de novo*, que nous avons mesuré **deux fois** ne pas savoir
  faire (BoltzGen en mode design, §11 ; les trois CDR, §12).

### À quoi s'en tenir

C'est **la cinquième hypothèse** sur l'origine de notre écart aux 0,78 du
plateau. Les quatre précédentes — le GPU, les donneurs, la distance à
l'archive, la surface de conception — ont toutes été réfutées par la mesure
(§0). Celle-ci a pour elle d'être la première à expliquer le chiffre plutôt que
de l'excuser (§13 : une charpente à 0,745 et un CDR3 à 0,729 donnent 0,48), et
elle a contre elle qu'une session de GPU peut la réfuter aussi. Rien n'est
acquis avant ces 34 mesures.

## 19. Le verdict : la règle d'unicité est bien le plafond. Cinquième hypothèse réfutée.

34 séquences, 52 min de GPU, 1,5 min/séquence — le chiffre du budget d'époque
est confirmé une troisième fois (design 1231 s, pliage 1794 s, analyse 53 s).
Mesures dans `data/measurements/minimal2_gpu_scores.csv`.

### Le témoin tient

Les 14 parents re-notent leur valeur archivée : **médiane −0,008**, écart-type
0,022, pire −0,062. La chaîne mesure donc les séquences et non la session, et
toute chute au-delà de ~0,05 est réelle.

### Les 20 variantes : 20 sur 20 en dessous de leur parent

| | |
|---|---|
| Δ médian variante − parent (re-mesuré) | **−0,172** |
| pire | −0,617 |
| meilleure | −0,021 |
| `design_iiptm` des variantes | médiane **0,562**, max **0,718**, min 0,149 |
| `design_iiptm` des parents | médiane 0,756, max 0,788 |
| variantes ≥ 0,74 | **0 / 20** |
| effondrements ≤ 0,25 | **6 / 20** |

C'est le scénario « vers 0,60 » annoncé au §18, et il faut le dire sans
l'arrondir : **la médiane est à 0,562**.

### Le coup de grâce, à longueur de CDR3 égale

Dans le régime récent (époque ≥ 26170, gagnant médian 0,753), en ne comparant
qu'aux époques qui imposent la longueur de la séquence :

| | notre meilleur | époques gagnées |
|---|---|---|
| **parents** (archivés, donc **non soumissibles**) | 0,788 à CDR3 = 18 | **6 / 6** |
| **variantes minimales** (soumissibles) | 0,718 à CDR3 = 18 | **0 / 6** |

Les séquences qui gagneraient sont exactement celles que nous ne pouvons pas
soumettre, et **la distance légale minimale à ces séquences coûte déjà plus que
la marge**. La règle d'unicité n'est pas un obstacle contournable à la marge :
c'est le plafond lui-même.

(Un piège évité : « 4 époques sur 8 » en mélangeant les longueurs. Une séquence
à CDR3 = 18 ne peut être soumise que dans une époque qui impose 18. C'est le
confondant de longueur du §12, et il transformait un 0/6 en 4/8.)

### Rien de mesurable hors GPU ne prédit la chute

| prédicteur | ρ de Spearman vs la chute | p |
|---|---|---|
| ΔPSH | −0,05 | 0,83 |
| PSH absolue | +0,10 | 0,69 |
| identité V | −0,03 | 0,90 |
| A3D moyenne | +0,02 | 0,94 |
| A3D pic | +0,38 | 0,10 |
| score du parent | **−0,47** | **0,04** |

Le seul signal significatif va dans le mauvais sens : **plus le parent est bon,
plus la variante tombe.** Et deux variantes légales du *même* parent à 0,776
atterrissent à 0,711 et 0,182 — une étendue de 0,529. L'issue dépend donc des
résidus tirés, pas du parent choisi, et aucune porte de développabilité ne la
voit venir. Il n'existe pas de criblage hors GPU pour cette stratégie.

### Ce qui reste vrai, et qui est maigre

0,718 est **notre meilleur candidat soumissible jamais mesuré** (contre 0,673
par greffage), soit le 71ᵉ percentile du champ à CDR3 = 18 au lieu du 50ᵉ.
C'est un progrès réel et insuffisant : zéro époque gagnée dans les deux cas.

### Bilan des cinq hypothèses

| # | hypothèse | verdict |
|---|---|---|
| 1 | le GPU (puissance de calcul) | réfutée (§12bis : distribution bornée à 0,674) |
| 2 | les donneurs (pool de CDR3) | réfutée (§10 : donneurs à 0,729–0,794) |
| 3 | la distance à l'archive | réfutée (ρ = −0,07) |
| 4 | la surface de conception (3 CDR) | réfutée (§12 : pire qu'1 CDR) |
| 5 | **la recombinaison / la perturbation minimale** | **réfutée ici : −0,172 médian, 0/6 époques** |

Il ne reste pas de chemin par recombinaison ni par perturbation d'un candidat
existant. Le seul chemin restant est **une vraie conception de novo** — que
nous avons mesuré deux fois ne pas savoir faire (BoltzGen en mode design, §11 ;
les trois CDR, §12) — ou un modèle de liaison à nous, entraîné pour proposer,
pas pour trier. Aucune soumission n'est justifiée en l'état, et aucun TAO ne
doit être engagé.

## 20. Non. Les gagnants ne font pas de conception *de novo* — ils remontent une lignée.

Avant d'engager quoi que ce soit dans le *de novo*, une mesure gratuite :
**comment les gagnants satisfont-ils la règle d'unicité ?** Pour chacun des 63
gagnants d'époque du régime récent, sa distance à l'archive *telle qu'elle
était* — uniquement les soumissions d'époques strictement antérieures.
(`scripts/winners_distance.py`, `scripts/winners_lineage.py`,
`data/measurements/winners_distance.csv`.)

### Le contrôle d'abord, parce que sans lui le chiffre ne vaut rien

| | similarité CDR3 |
|---|---|
| deux soumissions au hasard, **même longueur** | **0,333** |
| gagnant → son plus proche ancêtre archivé | **0,667** |
| écart au hasard | **+0,312**, et 61/63 au-dessus |

Les gagnants ne sortent donc pas de nulle part : ils sont **deux fois plus
proches d'une soumission antérieure que le hasard ne le voudrait**.

### Et ce dont ils sont proches, c'est d'un gagnant

| | |
|---|---|
| gagnants dont le plus proche ancêtre est lui-même un gagnant | **45 / 63** |
| gain sur cet ancêtre | **+0,074** de médiane, 59/63 en hausse |
| gain sur les seules chaînes gagnant → gagnant | +0,049 |

Les lignées se lisent époque par époque :

```
CDR3=18 : 0,759 → 0,762 → 0,762 → 0,772 → 0,785 → 0,772
CDR3=15 : 0,722 → 0,745 → 0,769 → 0,775 → 0,774 → 0,761 → 0,765 → 0,780
```

**C'est une montée de colline, pas une conception.** Chaque époque reprend le
meilleur candidat connu — l'archive est publique, donc celui de n'importe qui —
en change environ un tiers du CDR3, mesure, et soumet ce qui a gagné. Le gain
par pas est de +0,05, et il s'accumule sur des dizaines d'époques.

### Ce que cela corrige de notre §19

Le §19 concluait « la règle d'unicité est le plafond ». **C'est faux, et c'est
ma cinquième erreur de généralisation.** Le bon énoncé est plus étroit :

> perturber **au minimum légal**, **une seule fois**, **sans sélection**,
> détruit le score.

Ce que nous avons fait et ce que font les gagnants diffèrent sur trois axes, et
nous nous sommes trompés sur les trois :

| | nous (§19) | les gagnants |
|---|---|---|
| distance au parent (CDR3) | **0,75–0,80**, le minimum légal | **0,667**, un tiers du mot |
| nombre de tours | **1** | des dizaines, une par époque |
| sélection entre les tours | **aucune** | on ne garde que ce qui monte |

Le plus contre-intuitif est le premier : **perturber *plus* score mieux que
perturber au minimum.** Notre bande 0,75–0,80 est à la fois la plus proche que
la règle autorise et, mesurée, la pire. Hypothèse de mécanisme — à tester, pas
à croire : à 4 substitutions sur 16 la règle force la distance *dans le CDR3*,
donc on touche forcément des positions de contact ; à 6 on a le choix des
positions et on peut garder les contacts en déplaçant la périphérie.

### Le test qui s'ensuit, une seule variable

Mêmes parents, même chaîne, **seul le nombre de substitutions dans le CDR3
change** : 4 (déjà mesuré, §19, qui sert de bras témoin sans coûter un GPU), 6,
8, 10. `inputs/sweep18_s{6,8,10}.csv`, IGHV3-74, CDR3 = 18.

Si le score remonte avec la distance, la bande des gagnants est accessible et
le vrai sujet devient la boucle de sélection — c'est-à-dire un GPU dans les 67
minutes de l'époque, ce qui est un problème d'infrastructure et non de méthode.
S'il descend encore, alors la montée des gagnants vient d'autre chose que de la
distance, et il faudra chercher *où* ils substituent et non *combien*.

**Inversion d'ordre assumée :** cette fois la liaison est mesurée **avant** le
criblage PSH/A3D. La liaison est le discriminant, elle coûte 1,5 min/séquence
sur GPU, et le criblage coûte des heures de CPU pour une information qui ne
sert que si la liaison tient. Ne cribler que ce qui score.

## 21. Le balayage de distance : plus loin est strictement pire. Sixième hypothèse réfutée.

32 séquences, 50,1 min de GPU (1,6 min/séquence). Mesures dans
`data/measurements/sweep18_gpu_scores.csv`.

### Le témoin tient, et cette fois entre deux sessions

Trois parents déjà notés au §19 ont été re-notés ici : +0,014, −0,018, −0,010.
**Écart maximal 0,018.** La chaîne transfère donc d'une session à l'autre, pas
seulement au sein d'une session. Le parent jamais noté sort à 0,781.

### Les quatre bras, mêmes parents, même chaîne, une seule variable

| substitutions | identité CDR3 au parent | n | médiane | max | effondrées (< 0,30) | survivantes : médiane |
|---|---|---|---|---|---|---|
| 4 (minimum légal) | 0,778 | 8 | 0,439 | **0,718** | 4/8 (50 %) | **0,707** |
| 6 (bande des gagnants) | 0,667 | 7 | 0,633 | 0,669 | 3/7 (43 %) | 0,645 |
| 8 | 0,556 | 9 | 0,177 | 0,673 | 5/9 (56 %) | 0,614 |
| 10 | 0,444 | 12 | 0,475 | 0,655 | 5/12 (42 %) | 0,559 |

**L'hypothèse du §20 est réfutée.** Perturber plus ne score pas mieux : parmi
les survivantes, le score descend monotonement avec la distance — 0,707 →
0,645 → 0,614 → 0,559, **ρ = −0,64, p = 0,003**. Et la meilleure variante des
36 perturbations mesurées reste celle du minimum légal, à 0,718.

Donc la similarité de 0,667 des gagnants **n'est pas la cause de leur
réussite**. C'était mon hypothèse de mécanisme au §20 ; elle est fausse, et
c'est la sixième.

### Ce que le balayage montre en plus, et qui est nouveau

La distribution est **bimodale dans tous les bras** : un groupe effondré vers
0,13 et un groupe survivant vers 0,6. Et le taux d'effondrement **ne dépend pas
de la distance** (42 à 56 %, sans tendance) ni du nombre de substitutions
(effondrées 7,3 en moyenne, survivantes 7,5 ; Mann-Whitney p = 0,81).

Il y a donc **deux mécanismes distincts**, et un seul est lié à la distance :

1. un effondrement tout-ou-rien, qui frappe une variante sur deux quelle que
   soit la perturbation — et dont la signature (`design_iiptm` ≈ 0,13,
   `interaction_pae` > 22) est celle du §11 : la pose n'est pas reproduite au
   repliement indépendant ;
2. une dégradation graduelle des survivantes, proportionnelle à la distance.

### Ce que je ne peux PAS conclure, et pourquoi

J'ai calculé par position du CDR3 des écarts de ±0,3 à ±0,43. **Ils ne valent
rien et je ne les reporte pas comme résultat** : avec 36 variantes, 13
positions testables, chaque variante en touchant 4 à 10 (donc des appartenances
fortement corrélées) et une coupure effondrée/survivante à 17/19, toute
position plus souvent touchée chez les effondrées produit mécaniquement un
écart négatif. Le dispositif ne sépare pas un effet de position de « quelle
variante s'est effondrée ».

Pour y répondre il faudrait varier **la position à nombre de substitutions
constant** : même parent, même compte, positions choisies (périphérie du
CDR3 contre apex). C'est un dispositif propre, et il coûte une session.

### Une erreur de ma part à l'étape précédente

J'ai d'abord mesuré « où substituent les gagnants » en comparant chaque gagnant
au **gagnant de l'époque de son plus proche voisin**, et non à ce voisin
lui-même. Cela ne mesure pas son pas, et les chiffres qui en sortaient (47 % du
mot) sont faux. `scripts/winners_distance.py` stocke désormais la séquence du
voisin, et la mesure est relancée. Le résultat principal ci-dessus n'en dépend
pas.

### Où en est la compétitivité, sans arrondir

| | `design_iiptm` | époques gagnées (CDR3 = 18, régime récent) |
|---|---|---|
| meilleur parent archivé, **non soumissible** | 0,788 | 6/6 |
| meilleure variante soumissible, 36 mesurées | **0,718** | **0/6** |

Six hypothèses, six réfutations par la mesure. Aucune soumission justifiée,
aucun TAO engagé.
