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

> **RETIRÉ le 16 novembre (§55).** Les 150 designs de cette section n'ont
> **aucune cystéine** — le pont disulfure canonique Cys23-Cys104 a été détruit
> dans 100 % des cas par l'inverse folding de la version embarquée de BoltzGen,
> qui manque la machinerie de contrainte des résidus. La comparaison porte donc
> sur des protéines cassées, et ces designs n'étaient même pas soumissibles
> (`min_cysteines`). La conception conditionnée par la structure n'est pas
> réfutée : elle n'a pas été testée.

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

### Correction de la mesure des positions (§21)

Relancée avec la séquence du voisin lui-même : **pas médian 32 % du CDR3**, ce
qui recoupe exactement la similarité de 0,667 du §20. Les positions touchées se
répartissent tête 22 % / milieu 43 % / queue 35 %, pour 33 % attendus — et la
sous-représentation de la tête est mécanique, ses deux premiers résidus étant
encodés par le gène V. **Pas de préférence de position décelable chez eux non
plus.**

## 22. Ce qui reste debout : le volume de sélection, et c'est mesurable

Après six réfutations, la liste de ce qui pourrait expliquer l'écart s'est
vidée de tout ce qui touche à *ce que nous fabriquons*. Les gagnants partent
d'aussi bon que nous (nos parents re-notent 0,788, soit le sommet de l'archive),
perturbent à une distance qui nous est mesurée **pire** que la nôtre, et ne
choisissent pas leurs positions autrement que nous.

Ce qui reste est **combien de candidats sont notés avant de choisir**.

C'est la seule asymétrie non réfutée, et elle est quantifiable. Le bras à 4
substitutions est le seul dont les survivantes atteignent 0,707, mais il n'en
compte que **4** : son maximum de 0,718 est le maximum d'un échantillon de 4,
pas le plafond de la bande. Il faut 0,777.

### Le dispositif, et ce qu'il interdit de conclure à l'avance

**37 variantes** toutes à 4 substitutions, mêmes 6 parents, un seul bras
(48 demandées, 37 obtenues : les réessais n'ont pas tous franchi les portes).
Avec ~50 % d'effondrement attendu, cela donne une vingtaine de survivantes au
lieu de 4 — assez pour lire la queue haute de la distribution là où elle compte.

Le lot complet fait **44 séquences, ~1 h 10 de GPU** : 37 variantes, 1 parent en
re-contrôle (la chaîne a déjà été validée deux fois, à 0,022 intra-session au
§19 et 0,018 inter-session au §21) et **6 réplicats** dont il faut parler à part.

### Six séquences pour la question qui changerait le plus de choses

L'effondrement frappe une variante sur deux sans dépendre de la distance ni du
nombre de substitutions (§21). Or l'étape de design de BoltzGen **est
stochastique** — le carnet porte un `REPLICATE` exactement pour cela. Donc :
l'effondrement est-il une propriété de la séquence, ou du tirage ?

Trois survivantes (0,718 / 0,711 / 0,710) et trois effondrées (0,103 / 0,112 /
0,115) déjà mesurées repassent dans ce lot, prises aux extrêmes des deux côtés
de la coupure pour que le verdict ne tienne pas à la zone grise.

- **Une effondrée re-note haut** → l'effondrement est du bruit de tirage. Alors
  le « 20 sur 20 en dessous du parent » du §19 est en partie un artefact, et la
  sélection n'a pas besoin de mieux concevoir : elle a besoin de **re-noter**.
- **Elle s'effondre encore** → c'est la séquence, et la moitié de ce que produit
  la perturbation est structurellement perdue.

Piège qui invaliderait le test : BoltzGen **reprend** un run existant et
rendrait les anciens résultats au lieu de recalculer. Le nom de run est neuf
(`tail18_s4`), vérifié comme tel, et le carnet refuse de démarrer si le
répertoire du run existe déjà.

Et ce qu'on **ne** fera pas avec : ajuster une loi et extrapoler. Le §12bis a
déjà produit un 1,042 sur une grandeur bornée par 1, en traitant des queues
gaussiennes sur une distribution bornée. La lecture sera non paramétrique :
combien de survivantes dépassent 0,72, 0,74, 0,75, et le maximum observé.

Si une vingtaine de tirages à la bonne distance ne dépasse pas 0,74, le volume
de sélection n'est pas l'explication non plus, et la septième hypothèse tombe
comme les six autres. Si la queue monte, alors le sujet devient entièrement
**un GPU dans les 67 minutes de l'époque** — c'est-à-dire ce que l'utilisateur
disait au premier jour, mais mesuré au lieu d'être supposé, et avec un chiffre :
le nombre de candidats à noter par époque pour que le maximum dépasse 0,777.

## 23. Le mécanisme de l'effondrement : la séquence se replie encore, l'interface part d'un coup

Sur les 48 variantes mesurées (19 effondrées, 29 survivantes) et leurs 15
parents, toutes les métriques que BoltzGen rend en même temps que le score.

| métrique | effondrées | survivantes | parents |
|---|---|---|---|
| `design_iiptm` | 0,146 | 0,633 | 0,761 |
| `design_ptm` | **0,735** | **0,797** | 0,824 |
| `design_to_target_iptm` | 0,149 | 0,592 | 0,650 |
| `min_design_to_target_pae` | 18,16 | 3,06 | 2,17 |
| `interaction_pae` | 24,20 | 11,03 | 9,13 |
| `delta_sasa_refolded` | 620 | 887 | 894 |
| `plip_hbonds_refolded` | 3 | 9 | 9 |

### Deux formes de dégradation, et c'est tout le résultat

En cherchant le plus grand trou dans chaque distribution, rapporté à son
étendue :

| métrique | plus grand trou | |
|---|---|---|
| `design_iiptm` | **38 %** de l'étendue, entre 0,243 et 0,475 | **bimodal franc** |
| `interaction_pae` | 34 %, entre 15,9 et 21,3 | bimodal |
| `min_design_to_target_pae` | 28 %, entre 7,5 et 12,8 | bimodal |
| `design_ptm` | **8 %**, entre 0,770 et 0,779 | **continuum** |

**Zéro variante entre 0,30 et 0,45 de `design_iiptm`.** Le trou est réel, pas un
effet de seuil que j'aurais choisi.

Et `design_ptm` — le repliement du nanobody *seul*, sans la cible — ne montre
rien de tel : 48 valeurs étalées régulièrement de 0,711 à 0,827, sans trou.
Donc :

> **la qualité du repliement se dégrade continûment, l'interface s'effondre
> discontinûment.** Et le repliement ne rate jamais : un nanobody qui ne se
> replierait plus donnerait un `ptm` vers 0,3–0,4, alors que le minimum observé
> est 0,711.

C'est la lecture du §11 (« la pose n'est pas retrouvée au repliement
indépendant »), mais obtenue sur des séquences perturbées et non sur des designs
BoltzGen, donc par un chemin indépendant. La séquence se replie ; au-delà d'un
certain point de dégradation du repliement, l'interface basculeen bloc.

### Ce que ces chiffres ne sont PAS

Les deux groupes ne se chevauchent pas en `design_ptm` (effondrées ≤ 0,750,
survivantes ≥ 0,755) et ρ(`design_ptm`, `design_iiptm`) = +0,88 (p = 3·10⁻¹⁶).
**Cela ne fournit aucun criblage.** Trois raisons, qu'il faut garder ensemble :

1. **Ces métriques sortent du même tirage que le score.** On ne connaît le
   `design_ptm` d'une séquence qu'après avoir payé les 1,6 min de GPU qui
   donnent déjà le `design_iiptm`. Un seuil à 0,75 ne fait économiser
   rigoureusement rien.
2. La séparation parfaite tient à **0,004 point** sur un continuum. C'est une
   coïncidence d'ordre sur 48 valeurs, pas une frontière.
3. Une corrélation entre une grandeur continue et une grandeur discontinue
   mesure l'ordre, pas un mécanisme : elle ne dit pas lequel cause l'autre.

Un substitut utilisable devrait prédire l'effondrement **à partir de la séquence
seule**, et c'est exactement ce qu'aucune de nos grandeurs hors GPU ne fait
(§19 : PSH ρ = −0,05, A3D ρ = +0,02, identité V ρ = −0,03).

## 24. Le dispositif de position, à nombre de substitutions constant

Le §21 a mesuré des écarts par position allant jusqu'à 0,43 et a refusé de les
reporter : positions tirées au hasard, 4 à 10 par variante, donc appartenances
corrélées et confondues avec « quelle variante s'est effondrée ».

`scripts/minimal_variants.py` accepte désormais `--cdr3-positions`, qui **impose
l'ensemble** au lieu de le tirer. Quatre façons de se tromper sont refusées
nommément, chacune vérifiée : toucher la tête encodée par le gène V ou la queue
encodée par J (hors de la fenêtre libre 2–14 pour un CDR3 de 18), en demander
moins que le minimum légal — avec le calcul d'identité à l'appui —, et répéter
une position.

Deux bras, mêmes 6 parents, **toujours 4 substitutions** :

- **périphérie** : positions 2, 3, 13, 14 — les bords de la fenêtre libre,
  contre la tête V et la queue J figées ;
- **milieu** : positions 7, 8, 9, 10 — le centre de la boucle.

Si la périphérie préserve et que le milieu détruit, la position est le levier et
le générateur doit cesser de tirer au hasard. Si les deux bras se valent, alors
la position n'explique rien non plus, et il ne restera que le volume de
sélection du §22.

**Générés** : `inputs/pos18_periph.csv` (13 variantes, 6 parents) et
`inputs/pos18_milieu.csv` (9 variantes, 5 parents). Vérifié séquence par
séquence que les positions réellement modifiées sont exactement l'ensemble
imposé, et rien d'autre, dans les deux bras.

Une asymétrie de rendement à ne pas lire comme un résultat : le bras milieu
rend 9 variantes contre 13, et perd un parent entièrement, parce que davantage
de candidats sont refusés par nos propres portes déterministes. Cela *pourrait*
vouloir dire que substituer au centre de la boucle casse plus souvent une porte
(cystéines canoniques, marqueurs VHH) — mais la cause exacte des refus n'a pas
été capturée dans ce lancement, donc la seule chose établie est la différence de
rendement. À vérifier avant d'en tirer quoi que ce soit.

## 25. Bilan de la nuit du 11 au 12 novembre

Quatre commits, aucune dépense, aucun GPU loué, aucune soumission, aucun
lancement Kaggle.

### Ce qui a été mesuré

1. **Le mécanisme de l'effondrement est identifié** (§23). En cherchant le plus
   grand trou de chaque distribution rapporté à son étendue : `design_iiptm`
   38 % (zéro variante entre 0,30 et 0,45), `interaction_pae` 34 %,
   `min_design_to_target_pae` 28 % — mais `design_ptm` seulement 8 %, un
   continuum régulier de 0,711 à 0,827. Donc **le repliement se dégrade
   continûment, l'interface s'effondre discontinûment**, et le repliement ne
   rate jamais (un nanobody déplié donnerait 0,3–0,4). C'est la lecture du §11
   obtenue par un chemin indépendant : sur des séquences perturbées au lieu de
   designs BoltzGen.
2. **Et ce n'est pas un criblage**, ce que la section dit explicitement : ces
   métriques sortent du même tirage que le score, donc un seuil n'économise
   rien ; la séparation parfaite tient à 0,004 point sur un continuum ; et une
   corrélation entre une grandeur continue et une discontinue mesure l'ordre,
   pas un mécanisme.
3. **Le dispositif de position existe enfin** (§24). `--cdr3-positions` impose
   l'ensemble des positions au lieu de le tirer, et refuse nommément quatre
   façons de se tromper — tête encodée par V, queue encodée par J, moins que le
   minimum légal, position répétée — chacune vérifiée.

### Ce qui a échoué, et pourquoi

- **La génération a rendu 37 variantes sur 48 demandées**, et les deux bras de
  position 13 et 9 sur 18 : beaucoup de candidats sont refusés par nos propres
  portes déterministes. Le rendement suffit ici, mais il dit qu'un générateur
  qui tire au hasard travaille contre les portes.
- **Un redémarrage du processus de session** vers 21 h 40 a interrompu une
  suppression de déclencheur en attente d'approbation. Sans conséquence :
  fichiers et commits vérifiés contre `origin`, et le déclencheur était à tir
  unique.
- **Deux travaux faits plus tôt que leur réveil**, donc deux réveils ont trouvé
  l'ouvrage fait. Pas une perte, mais la preuve qu'une tâche programmée doit
  commencer par constater l'état plutôt que par refaire.

### Les lots prêts à lancer

Pour les deux : *Accelerator* **GPU T4 x2**, *Internet* **On**, puis
*Save Version → Save & Run All (Commit)*.

| lot | n | GPU | état |
|---|---|---|---|
| **`kaggle_tail18.ipynb`** | 44 | **~1 h 10** | carnet vérifié, prêt |
| position (périphérie 13 + milieu 9) | 22 + parents | ~40 min | lots générés, carnet à construire |

Le premier passe d'abord : il porte les deux questions qui décident de la suite.

### Ce que le lot de queue haute tranchera

37 variantes à 4 substitutions, la seule bande dont les survivantes atteignent
0,707. Lecture **non paramétrique décidée d'avance** : combien dépassent 0,72 /
0,74 / 0,75, et le maximum observé. Pas d'ajustement de loi — le §12bis a déjà
produit un 1,042 sur une grandeur bornée par 1.

- **Si une vingtaine de tirages ne dépasse pas 0,74** : le volume de sélection
  n'explique pas l'écart, et la septième hypothèse tombe comme les six autres.
- **Si la queue monte vers 0,777** : tout le sujet devient *un GPU dans les 67
  minutes de l'époque*, et il faudra chiffrer combien de candidats noter par
  époque pour que le maximum y arrive.

### Et ce que les six réplicats trancheront, qui compte davantage

Trois survivantes (0,718 / 0,711 / 0,710) et trois effondrées (0,103 / 0,112 /
0,115) repassent sous un nom de run neuf.

- **Une effondrée re-note haut** → l'effondrement est du tirage, pas de la
  séquence. Alors le « 20 sur 20 en dessous du parent » du §19 est en partie un
  artefact, et la sélection n'a pas besoin de mieux concevoir : elle a besoin de
  **re-noter**. Cela changerait la lecture de tout ce qui précède.
- **Elle s'effondre encore** → c'est la séquence, et la moitié de ce que produit
  la perturbation est structurellement perdue.

### À quoi s'en tenir

Sept hypothèses, **six réfutées par la mesure** : le GPU comme puissance, les
donneurs, la distance à l'archive, la surface de conception, la recombinaison,
et la distance de perturbation. La septième — le volume de sélection — est la
seule qui reste, et elle peut tomber comme les autres. Rien n'est acquis avant
ces 44 mesures, et aucun TAO ne doit être engagé.

## 26. La septième hypothèse tient : le volume de sélection est la réponse

44 séquences, 63,9 min de GPU. `data/measurements/tail18_s4_scores.csv`. Témoin
parent : 0,788 → 0,790.

### Les réplicats : l'effondrement est une propriété de la séquence

| | avant | maintenant | |
|---|---|---|---|
| survivante | 0,718 | 0,734 | +0,016 |
| survivante | 0,711 | 0,656 | −0,055 |
| survivante | 0,710 | 0,691 | −0,019 |
| effondrée | 0,112 | 0,113 | +0,001 |
| effondrée | 0,103 | 0,107 | +0,004 |
| effondrée | 0,115 | 0,177 | +0,062 |

**3 effondrées sur 3 restent effondrées** (maximum 0,177), 3 survivantes sur 3
restent survivantes (minimum 0,656). Donc **le « 20 sur 20 en dessous du
parent » du §19 n'est pas du bruit** : l'effondrement est structurel, et la
sélection ne peut pas le contourner en re-notant. Des deux réponses possibles,
c'est celle qui ferme une porte.

À noter quand même : les survivantes bougent de ±0,055 d'un tirage à l'autre,
contre ±0,018 pour les parents. Le bruit de mesure est donc **trois fois plus
grand sur une variante que sur une séquence archivée**, ce qui compte pour toute
sélection par seuil.

### La queue haute : pour la première fois, une variante soumissible gagnerait

37 tirages à 4 substitutions, contre 8 au §21. Lecture non paramétrique décidée
d'avance :

| | |
|---|---|
| maximum | **0,76984** |
| ≥ 0,72 | 4/37 |
| ≥ 0,74 | 1/37 |
| ≥ 0,75 | 1/37 |
| ≥ 0,777 | 0/37 |
| effondrées | 4/37 (11 %) |

Et contre le champ, à longueur de CDR3 égale et dans le régime récent :

| | `design_iiptm` | époques gagnées (CDR3 = 18) |
|---|---|---|
| ancien record soumissible (§19) | 0,718 | **0/6** |
| **ce lot** | **0,770** | **3/6** |
| meilleur parent, non soumissible | 0,790 | 6/6 |

Gagnants récents à CDR3 = 18 : 0,759 / 0,762 / 0,762 / 0,772 / 0,772 / 0,785.
**0,770 en bat trois.** C'est la première fois en sept hypothèses qu'une
séquence que nous pouvons réellement soumettre franchit un seuil de victoire.

### Ce que cela coûte, chiffré — et pourquoi ce n'est pas gagné

0,770 est **un maximum sur 37 tirages**, pas un résultat qu'on refait à volonté.
En fréquence : P(un tirage ≥ 0,75) = 1/37 = 2,7 % [IC95 0,1–14,2 %].

| candidats notés | P(au moins un ≥ 0,75) | GPU |
|---|---|---|
| **31** (ce que l'époque permet) | **57 %** | 0,8 h |
| 50 | 75 % | 1,2 h |
| 84 | 90 % | 2,1 h |

Le budget d'époque laisse ~47 min à BoltzGen, soit **~31 candidats**. On tirerait
donc un candidat ≥ 0,75 environ une époque sur deux — et 0,75 ne suffit pas
toujours, la médiane des gagnants récents à CDR3 = 18 étant 0,767.

**C'est le diagnostic que l'utilisateur donnait au premier jour** — « le goulot
est le GPU » — mais mesuré au lieu d'être supposé, et avec un chiffre : il ne
s'agit pas de puissance brute mais de **débit de sélection dans la fenêtre de
67 minutes**, et il en faudrait 2 à 3 fois plus que ce que l'époque permet pour
viser la victoire de façon fiable.

### Une anomalie que je ne sais pas expliquer

Le taux d'effondrement est de **11 % ici contre 50 %** sur les 8 variantes du
même bras au §21 (Fisher p = 0,024). Même générateur, mêmes parents, seule la
graine diffère (0 contre 42). La différence est réelle et **inexpliquée**. Elle
doit être traitée comme une variance de tirage non comprise, pas comme un
acquis : si le vrai taux est 50 %, le nombre de candidats utiles par époque
tombe de moitié et les probabilités ci-dessus avec lui.

### État des sept hypothèses

| # | hypothèse | verdict |
|---|---|---|
| 1 | le GPU comme puissance de calcul | réfutée (§12bis) |
| 2 | les donneurs | réfutée (§10) |
| 3 | la distance à l'archive | réfutée (ρ = −0,07) |
| 4 | la surface de conception | réfutée (§12) |
| 5 | la recombinaison | réfutée (§19) |
| 6 | la distance de perturbation | réfutée (§21) |
| 7 | **le volume de sélection** | **soutenue (§26)** |

Première hypothèse qui survit à sa mesure. Elle reste à confirmer sur un second
tirage indépendant — un maximum sur 37 n'est pas une capacité — et rien ne
justifie encore d'engager du TAO.

## 27. Le second tirage : prédictions enregistrées AVANT de le lancer

Un maximum sur 37 n'est pas une capacité, et l'anomalie du taux d'effondrement
(§26) n'est pas expliquée. Le second tirage est donc un vrai test, et il ne vaut
que si les prédictions sont écrites d'abord. Elles le sont ici.

### Les valeurs groupées, qui servent de référence

Tout le bras à 4 substitutions, CDR3 = 18, graines 0 et 42 confondues :

| | | IC95 |
|---|---|---|
| effondrées (< 0,30) | 8/45 = **17,8 %** | 8,0–32,1 % |
| ≥ 0,72 | 4/45 = 8,9 % | 2,5–21,2 % |
| ≥ 0,75 | 1/45 = **2,2 %** | 0,1–11,8 % |

L'anomalie du §26 vue de près : 4/8 (graine 0) contre 4/37 (graine 42), Fisher
p = 0,024, **mais sur 8 observations, et c'est un test parmi plusieurs lancés
sans correction**. Les intervalles se chevauchent (15,7–84,3 % contre
3,0–25,4 %). Signal fragile ; pour planifier, la valeur groupée.

Le parent n'explique pas la différence : dans le lot de 37 les effondrements
sont répartis sur trois parents (1, 2 et 1), et dans celui de 8 il n'y avait que
1 à 2 variantes par parent, donc des taux par parent sans contenu.

### Prédictions, sur les nouvelles variantes (graine 1234, mêmes 6 parents)

La génération a rendu **45** variantes, pas 40. Les prédictions ci-dessous
étaient écrites en comptes sur 40 ; elles sont rapportées à 45 **au même taux**,
sans rien changer d'autre — ce sont les taux groupés qui sont la prédiction, pas
les comptes. Sur 45 : effondrées **8 ± 5**, ≥ 0,72 : **4**, ≥ 0,75 : **1**.
Zéro chevauchement de séquences avec le premier tirage, vérifié.

1. **effondrées : 7 ± 5** (17,8 % de 40). Au-delà de 14, le taux de 11 % était la
   chance et les probabilités du §26 tombent de moitié. En dessous de 3, c'est le
   50 % de la graine 0 qui était la chance.
2. **≥ 0,72 : 3 ou 4** (8,9 %). Zéro serait un désaveu net.
3. **≥ 0,75 : 1**, et donc **0 ou 1 est compatible** — c'est la limite de ce
   tirage, qui ne peut pas départager 2,2 % de 0,5 %. Une seconde variante
   ≥ 0,75 confirmerait que la queue est réelle ; aucune ne la réfuterait pas.
4. **maximum attendu : 0,72 à 0,77.** Un maximum au-dessus de 0,78 serait une
   surprise à ne pas croire sans réplicat.

### Les quatre meilleures du §26 repassent dans le lot

C'est l'autre moitié du test, et c'est celle qui décide d'une soumission
éventuelle. Les survivantes bougent de ±0,055 d'un tirage à l'autre (§26), donc
**notre 0,770 pourrait être noté 0,72 par le validateur**, qui ne le note qu'une
fois. Les quatre variantes ≥ 0,72 (0,770 / 0,737 / 0,728 / 0,725) sont donc
re-mesurées ici : deux points chacune.

Si 0,770 re-note au-dessus de 0,75, c'est un candidat. S'il re-note vers 0,72,
alors ce que nous avons n'est pas une séquence à 0,770 mais une séquence dont
l'espérance est plus basse, et le 3/6 époques du §26 est un artefact du
maximum.

**Prédiction : les quatre perdent en moyenne, par régression vers la moyenne —
elles ont été choisies pour être les plus hautes de 37.** Un maximum sélectionné
re-mesuré redescend presque toujours. Si les quatre tiennent leur valeur, c'est
que la variance de ±0,055 du §26 était surestimée sur trois points.

### Le lot

`kaggle_tail18b.ipynb` — **50 séquences, ~1 h 15 de GPU** : 45 nouvelles
variantes, les 4 meilleures du §26 re-mesurées, 1 parent en témoin. Nom de run
neuf (`tail18b_s4`), vérifié contre les trois runs déjà consommés, sans quoi
BoltzGen rendrait les anciens résultats au lieu de recalculer — ce qui
annulerait précisément les 4 re-mesures. *Accelerator* **GPU T4 x2**,
*Internet* **On**.

## 28. Le second tirage : la queue est réelle mais deux fois plus mince, et le 3/6 était un artefact

50 séquences, 76,7 min de GPU. Témoin parent : 0,790 → 0,788.
`data/measurements/tail18b_s4_scores.csv`.

### Les prédictions du §27, vérifiées une par une

| prédiction | observé | verdict |
|---|---|---|
| les 4 meilleures perdent en moyenne | **4/4 en baisse**, moyenne −0,015 | **vérifiée** |
| effondrées 8 ± 5 (3–13) | **14/45 = 31 %** | **démentie**, par le haut |
| ≥ 0,72 : 4 | **1/45** | **démentie**, par le bas |
| ≥ 0,75 : 0 ou 1 | 0/45 | vérifiée (mais non informative, dit d'avance) |
| maximum 0,72–0,77 | **0,723** | au bord inférieur |

Deux prédictions sur quatre démenties, et dans le sens défavorable. C'est le
mérite de les avoir écrites d'abord : sans elles, j'aurais présenté ce lot comme
« cohérent avec le §26 ».

### Les trois tirages du même bras, côte à côte

| graine | n | effondrées | ≥ 0,72 | max |
|---|---|---|---|---|
| 0 | 8 | 50 % | 0 | 0,718 |
| 42 | 37 | **11 %** | 4 | **0,770** |
| 1234 | 45 | 31 % | 1 | 0,723 |
| **groupé** | **90** | **24,4 %** [16–35] | 5,6 % [1,8–12,5] | — |

**C'était le 11 % qui était la chance, pas le 50 %.** Le §26 a été écrit sur le
tirage le plus favorable des trois, et c'est exactement l'erreur contre laquelle
le §27 avait été écrit.

### Le candidat : le 3/6 époques ne survit pas

Notre meilleure séquence a maintenant deux mesures : **0,770 puis 0,756**,
moyenne 0,763.

| | époques gagnées (CDR3 = 18, régime récent) |
|---|---|
| à 0,770 (mesure 1) | 3/6 |
| **à 0,756 (mesure 2)** | **0/6** |

Gagnants : 0,759 / 0,762 / 0,762 / 0,772 / 0,772 / 0,785. La seconde mesure
passe **sous le plus faible d'entre eux**. Le validateur ne note qu'une fois :
soumettre cette séquence, c'est tirer une fois dans une distribution dont le
centre est à 0,763 et l'écart-type par mesure à **0,016** (estimé sur 7 paires
répétées). Elle est donc à la frontière de la plus facile des six époques, et
perd les cinq autres.

**Aucun de nos 90 tirages n'atteint 0,772**, la médiane des gagnants.

### Planification révisée

| candidats notés | P(au moins un ≥ 0,75) | GPU |
|---|---|---|
| **31** (ce que l'époque permet) | **29 %** | 0,8 h |
| 62 | 50 % | 1,6 h |
| 100 | 67 % | 2,5 h |
| 200 | 89 % | 5,0 h |

Le §26 annonçait 57 % pour 31 candidats, puis 50 % après groupement. C'est
**29 %**, et un candidat à 0,75 ne gagne de toute façon aucune des six époques
récentes. Il faudrait viser 0,772, que 0 tirage sur 90 atteint.

### Où en est la septième hypothèse

**Soutenue dans sa forme, démentie dans son ampleur.** Plus de tirages donnent
bien un meilleur maximum — 8 tirages plafonnent à 0,718, 37 atteignent 0,770 —
donc le volume de sélection agit. Mais le taux est 2 à 3 fois plus faible que
l'estimation du §26, et le débit que l'époque autorise reste très loin de ce
qu'il faudrait : il manque un facteur ~6 (200 candidats contre 31) pour
seulement 89 % de chance de tirer un ≥ 0,75, qui ne gagne rien.

Donc : la mesure ne réfute pas l'hypothèse, elle en chiffre l'insuffisance. En
l'état, **aucune soumission n'est justifiée** — notre meilleure séquence perd
cinq époques sur six et joue la sixième à la pièce — et aucun TAO ne doit être
engagé.

### Ce que je retiens sur la méthode

Le §26 a été écrit après un seul tirage, et ses deux chiffres clés (11 %
d'effondrement, 3/6 époques) étaient tous les deux les plus favorables possibles
— l'un un tirage chanceux, l'autre un maximum sélectionné. Les deux sont tombés
au second tirage. **Un maximum et un taux sur un seul lot ne sont pas des
résultats**, et la seule chose qui l'a montré ici est d'avoir écrit les
prédictions avant de regarder.

## 29. Le plan : le pré-calcul lève le goulot, et il faut tester là où c'est le plus fragile

Le §28 a laissé la septième hypothèse debout mais insuffisante : il faudrait
noter ~200 candidats pour 89 % de chance d'un ≥ 0,75, contre ~31 que l'époque
permet. Ce chiffre de 31 vient d'une contrainte précise, et elle est
contournable.

### Ce qui change tout : l'époque impose un couple, et l'espace est énumérable

Vérifié sur les 63 époques récentes : **chacune impose exactement un gène V**
(une intersection non vide des gènes de toutes ses soumissions, dans les 63 cas)
et une longueur de CDR3. Soit 47 couples distincts, dont 35 vus une seule fois —
une queue lourde, mais une tête concentrée.

**Donc rien n'oblige à noter dans la fenêtre de 72 minutes.** On peut pré-noter
des candidats par couple, hors époque, et n'avoir qu'à soumettre le meilleur
déjà mesuré quand le couple revient. Le plafond de 31 candidats était un
artefact de la notation en direct.

| candidats par couple | 31 | 60 | 100 | 150 |
|---|---|---|---|---|
| IGHV3-74 / L=18 (5 époques) | 46 % | 69 % | **86 %** | 95 % |
| IGHV3-23 / L=15 (3 époques) | 47 % | 71 % | 87 % | 95 % |
| IGHV3-23 / L=7 (2 époques) | 88 % | 98 % | 100 % | 100 % |

### Le coût du pré-calcul, sur les 8 couples les plus fréquents

| candidats/couple | GPU total | semaines Kaggle (30 h) | époques couvertes | part des époques gagnées |
|---|---|---|---|---|
| 31 | 6,2 h | 0,2 | 32 % | 13,7 % |
| 60 | 12,0 h | 0,4 | 32 % | 20,3 % |
| **100** | **20,0 h** | **0,7** | **32 %** | **25,4 %** |

C'est à la portée du quota Kaggle, sans louer un seul GPU.

### Les quatre tests, dans cet ordre, et pourquoi cet ordre

**1. Le test de position** (déjà généré, `inputs/pos18_periph.csv` et
`pos18_milieu.csv`, 22 variantes, **~35 min**). Dernier levier non testé sur ce
que nous fabriquons. Si la périphérie préserve et le milieu détruit, tout le
reste du plan devient 2 à 3 fois plus efficace, parce que le générateur cesse de
tirer au hasard.

**2. IGHV3-23 / L=7** (~30 variantes + parents, **~1 h**). C'est le couple qui
affiche 88 %, et c'est aussi **celui dont la justification est la plus faible** :
à L=7 il ne reste que 2 positions libres dans le CDR3, donc **5 des 7
substitutions imposées tombent dans la charpente**. La perturbation y est de
nature différente de tout ce que nous avons mesuré, et la distribution de coût
du §28 (mesurée à L=18, 4 substitutions dans le CDR3) n'a aucune raison de s'y
transporter. Tester d'abord là où le chiffre est le plus beau et le moins fondé.

**3. IGHV3-74 / L=21** (~40 variantes, **~1 h**). L'alternative riche en
liberté : 16 positions libres, 4368 combinaisons, ~10⁸ variantes possibles,
2 substitutions de charpente seulement. Vérifie le transport de la distribution
de coût dans l'autre sens.

**4. La campagne de pré-calcul**, seulement si 2 ou 3 confirment le transport :
100 candidats sur les 8 meilleurs couples, 20 h de GPU étalées.

### Ce qui a été fermé aujourd'hui, sans GPU

- **Les longueurs L=5 et L=6 sont infaisables** : le nombre de substitutions que
  la règle impose dans le CDR3 dépasse le nombre de positions libres.
- **Prédire l'effondrement depuis la séquence seule : trop faible pour servir.**
  Sur 122 variantes à L=18 (36 effondrées), validation croisée groupée par
  parent : AUC 0,610 en régression logistique, 0,579 en forêt aléatoire, contre
  **0,516 pour le témoin à étiquettes mélangées**. Il y a donc un signal
  au-dessus du hasard, mais 0,61 ne permet pas d'écarter un candidat qui coûte
  1,5 min, et cela ne repose que sur 9 parents indépendants. À reprendre vers
  n ≈ 300, pas avant. Les 24 % de budget perdus en effondrements restent perdus.

### Les réserves, nommées

1. **Tous les P(gain) héritent d'une queue estimée sur 82 tirages**, où seuls
   2 ont coûté moins de 0,02. Les intervalles sont larges et le §28 a déjà montré
   qu'un seul lot surestime.
2. **Le transport de la distribution de coût à une autre longueur n'est pas
   vérifié** — c'est précisément l'objet des tests 2 et 3.
3. **L'archive grossit.** Un candidat pré-noté peut devenir doublon si
   quelqu'un soumet près de lui. La vérification d'unicité doit être refaite au
   moment de soumettre ; elle coûte du CPU, pas du GPU.
4. **La queue des couples est lourde** : 8 couples ne couvrent que 32 % des
   époques, et la couverture plafonne vite. Les 68 % restants ne sont pas joués.

## 30. Retrait : « nos portes refusent le milieu deux fois plus souvent » est faux

J'ai rapporté il y a une heure, et commité, que nos portes déterministes
refusent une substitution au centre du CDR3 deux fois plus souvent qu'à ses
bords : 31/42 contre 18/42, Fisher p = 0,0075, sur deux graines indépendantes.

**La mesure directe de ces portes le contredit.** En répétant 60 tentatives par
bras et en comptant les causes de refus, sans la règle d'archive :

| bras | acceptés par les portes déterministes |
|---|---|
| périphérie 2,3,13,14 | **56/60** |
| milieu 7,8,9,10 | **52/60** |

Et les causes sont les mêmes des deux côtés — un gène V hors cible après
substitution de charpente, et deux fois une identité V à 81,6 contre un plancher
de 82. Rien qui distingue le milieu.

Donc **les 74 % contre 43 % ne viennent pas des portes.** Ils viennent d'ailleurs
dans la boucle de génération — la règle d'archive, que ce contrôle exclut, ou
des collisions entre variantes d'un même parent quand les positions sont figées
et l'alphabet réduit à neuf résidus hydrophiles. C'est en cours de mesure.

### Ce que j'ai mal fait

J'ai calculé un Fisher sur deux nombres lus dans des journaux de génération, et
j'ai nommé la cause (« nos propres portes ») sans l'avoir mesurée. Les deux
nombres étaient justes ; l'explication était une supposition présentée comme un
fait, et elle est fausse. Le p = 0,0075 porte sur un rendement de bout en bout,
pas sur les portes.

C'est la même faute que les §7ter, §12 et §21 : un chiffre réel, une cause
devinée. La seule défense qui ait jamais marché ici est de mesurer le mécanisme
séparément, ce que ce contrôle a fait — en une commande, et avant que
l'affirmation ne serve à décider quoi que ce soit.

**Ce qui reste établi** : le rendement de bout en bout du bras milieu est la
moitié de celui de la périphérie, sur deux graines. **Ce qui ne l'est pas** : la
raison. Et rien là-dedans ne dit quoi que ce soit sur la liaison, qui est ce que
la session en cours mesure.

## 31. Le test de position : pas d'effet. Et deux choses que la lecture pré-enregistrée ne couvrait pas.

23 séquences, 36,7 min de GPU. Témoin parent : 0,788 → 0,782.
Analyse par `scripts/analyze_position.py`, **inchangé depuis son commit e4a4381,
écrit avant que les données existent**.

### Le verdict pré-enregistré

| bras | n | médiane | max | effondrées |
|---|---|---|---|---|
| périphérie 2,3,13,14 | 13 | 0,619 | 0,783 | 0/13 |
| milieu 7,8,9,10 | 9 | 0,615 | 0,660 | 0/9 |

**Mann-Whitney p = 0,894, delta de Cliff −0,04.** Les deux médianes sont à
0,004 l'une de l'autre. **Pas d'effet**, et pas seulement « pas d'effet
large » : le delta de Cliff est nul à la précision disponible.

La position n'est donc pas le levier que les écarts du §21 laissaient espérer —
ces écarts étaient bien l'artefact que le §21 avait soupçonné. Le plan du §29 ne
gagne pas le facteur 2 ou 3 que ce test aurait pu lui donner, et le générateur
peut continuer à tirer ses positions au hasard.

### Hors pré-enregistrement 1 : l'effondrement a disparu

| | effondrées |
|---|---|
| positions **imposées** (ce lot) | **0/22** |
| positions **tirées au hasard** dans 2–14, 4 substitutions, L=18 | 18/82 = 22 % |
| | Fisher p = **0,0114** |

Ce n'est **pas** une comparaison pré-enregistrée — elle oppose ce lot à
l'historique, pas un bras à l'autre — donc c'est une hypothèse, pas un résultat.
Mais elle est précise et testable : les deux ensembles imposés couvrent
2,3,7,8,9,10,13,14, et **ne touchent jamais 4, 5, 6, 11, 12**. Si
l'effondrement vient de ces cinq positions, l'éviter rendrait 22 % du budget
GPU — ce que le classifieur du §29 n'a pas su faire (AUC 0,61).

Le test qui trancherait : un bras imposé sur 4,5,6,11 (quatre des cinq
suspectes) contre un bras imposé sur 2,3,13,14, même nombre de substitutions.
Si les suspectes s'effondrent et les autres non, c'est la première règle de
génération utile que nous aurons.

### Hors pré-enregistrement 2 : un candidat à 0,783, et pourquoi je n'y crois pas encore

La meilleure variante du lot sort à **0,78312**, au-dessus de notre précédent
record soumissible (0,770) et au niveau du parent re-noté dans la même session
(0,782). Elle battrait **5 des 6** gagnants récents à CDR3 = 18.

**C'est un maximum, et le §28 a déjà montré ce que valent nos maximums** : le
0,770 du §26 a re-noté 0,756 au tirage suivant, et les quatre meilleures ont
toutes baissé, de 0,015 en médiane. À 0,783 − 0,015 ≈ 0,769, le compte retombe
à 3 époques sur 6.

Donc : **cette séquence doit être re-notée avant qu'on en parle comme d'un
candidat**, et c'est six minutes de GPU dans n'importe quel prochain lot. Rien
ne justifie une soumission sur un maximum non répliqué — c'est exactement
l'erreur du §26.

### Et une note sur le rendement du bras milieu

Le §30 a retiré l'explication par les portes. Les mesures suivantes ont écarté
les deux autres causes faciles : **aucune collision** (240 variantes distinctes
sur 240 tentatives dans les deux bras) et les portes passent même légèrement
mieux au milieu (211/240 contre 206/240). Il ne reste que la règle d'archive,
dont la mesure tourne encore. Le fait brut — rendement de bout en bout deux fois
plus faible au milieu — reste établi, sa cause reste inconnue, et le présent
§31 montre qu'elle n'a aucun rapport avec la liaison.

## 32. Les deux hypothèses du §31, testées gratuitement : l'une morte, l'autre confondue

### 1. Les positions 4,5,6,11,12 ne causent pas l'effondrement — l'inverse, même

Testé sur les 82 variantes à positions aléatoires déjà mesurées, en un seul test
sur le groupement que le §31 avait énoncé d'avance :

| positions suspectes touchées | n | effondrées |
|---|---|---|
| **aucune** | 8 | **3 = 38 %** |
| au moins une | 74 | 15 = 20 % |
| | | Fisher p = 0,364 |

Aucune association, et la tendance va dans le sens contraire de l'hypothèse.
**Piste morte, pour zéro minute de GPU** — là où j'allais proposer une session
entière pour la tester.

### 2. Le 0/22 est confondu avec la session

| session | effondrées (positions aléatoires, 4 substitutions, L=18) |
|---|---|
| tail18 | 4/37 = 11 % |
| tail18b | 14/45 = 31 % |
| pos18 (positions imposées) | **0/22 = 0 %** |
| | **χ² d'homogénéité : p = 0,0029** |

**Le taux d'effondrement n'est pas stable d'une session à l'autre.** Le 0/22 du
§31 ne peut donc pas être attribué aux positions imposées : les deux
explications — mode de génération et session — sont inséparables dans les
données existantes.

Et cela a une conséquence qui dépasse la question : **les taux d'effondrement
qui servent au §29 sont instables**, donc les probabilités de gain par époque
qu'ils alimentent sont plus incertaines que leurs intervalles ne le disaient.
Ces intervalles supposaient un taux unique.

### Ce que le prochain lot doit faire, et c'est la seule façon de trancher

Les deux modes de génération **dans la même session**, une seule variable :

- bras **aléatoire** : positions tirées dans 2–14, 4 substitutions ;
- bras **imposé** : positions 2,3,13,14, 4 substitutions ;
- mêmes 6 parents, même chaîne, un parent en re-contrôle.

Plus les **3 meilleures du §31 en réplicats** (0,783 / 0,757 / 0,735), parce
qu'un maximum non répliqué ne vaut rien et que le §28 l'a déjà montré sur le
0,770 devenu 0,756.

Prédictions, écrites avant le lot :

1. si les positions imposées suppriment réellement l'effondrement : bras imposé
   ≈ 0 %, bras aléatoire ≈ 20 %, et un écart significatif ;
2. si c'était un effet de session : les deux bras au même niveau, quel que soit
   ce niveau — et alors le 0/22 était du bruit et 22 % du budget GPU restent
   perdus ;
3. le 0,783 **baisse** d'environ 0,015 (la baisse médiane des quatre re-mesures
   du §28) et atterrit vers 0,768. S'il se maintient au-dessus de 0,78, c'est le
   premier candidat soumissible sérieux de toute l'enquête.

## 33. Les trois prédictions tiennent, et un candidat répliqué à 0,786

43 séquences, 68,6 min de GPU. Témoin parent : 0,788 → 0,786.

### Prédictions 1 et 2 : c'était la session, pas le mode

| bras | n | médiane | max | effondrées |
|---|---|---|---|---|
| positions **tirées** | 21 | 0,619 | 0,712 | 1/21 = 5 % |
| positions **imposées** 2,3,13,14 | 18 | 0,621 | 0,743 | 1/18 = 6 % |
| | | | | Fisher **p = 1,000** |

Mann-Whitney sur les scores : **p = 0,944**. Les deux modes sont
indiscernables. **La prédiction 2 est celle qui se réalise** : le 0/22 du §31
était un effet de session, et les 22 % de budget perdus en effondrements le
restent.

Et le taux d'effondrement par session est maintenant 11 % / 31 % / 0 % / 5 %,
**χ² p = 0,0006**. Il ne dépend ni de la distance (§21), ni du nombre de
substitutions (§21), ni des positions (§31, §33), ni du mode de tirage — mais il
varie d'un facteur 6 selon la session. Une cause plausible, non testée : l'étape
de design de BoltzGen est stochastique et tout un lot partage sa session. Cela
n'infirme pas le §26 — les séquences effondrées restent effondrées d'une session
à l'autre — mais cela veut dire que les **séquences limites** basculent selon le
tirage, et que tout taux mesuré sur un seul lot est trompeur.

### Prédiction 3 : vérifiée au chiffre près, et elle livre un candidat

| mesure §31 | ici | écart |
|---|---|---|
| 0,78312 | **0,78889** | **+0,006** |
| 0,75679 | 0,74191 | −0,015 |
| 0,73504 | 0,70373 | −0,031 |
| | **médiane −0,015** | prédit −0,015 |

La médiane prédite était −0,015 ; elle est −0,015. Mais la meilleure des trois
**n'a pas baissé** : deux mesures indépendantes, 0,78312 et 0,78889, moyenne
**0,786**, écart 0,006 entre deux sessions dont les témoins parents étaient à
0,782 et 0,786.

### Le premier candidat sérieux de l'enquête

| | |
|---|---|
| estimation | **0,786 ± 0,010** (σ par mesure 0,014 sur 10 paires répétées) |
| gagnants récents à CDR3 = 18 | 0,759 / 0,762 / 0,762 / 0,772 / 0,772 / 0,785 |
| battus par la moyenne | **6/6** |

Le validateur ne note qu'une fois, donc en probabilité par époque : 97 %, 96 %,
95 %, 84 %, 83 %, 52 % — soit **5,1 époques sur 6** en espérance.

C'est le premier candidat de toute l'enquête qui ne soit pas un maximum non
répliqué, et c'est exactement ce que le §28 exigeait avant d'en parler.

### Ce qui manque, et c'est bloquant

**Sa PSH et son Aggrescan3D n'ont jamais été mesurées.** Le §20 avait inversé
l'ordre volontairement — liaison d'abord, développabilité ensuite, pour ne pas
dépenser des heures de CPU sur des séquences qui ne lient pas. La contrepartie
est que ce candidat **n'est pas soumissible en l'état** : au §17, seules 20
variantes sur 61 franchissaient PSH < 96 et les deux plafonds A3D, soit une sur
trois.

Le criblage des six meilleures de ce lot est lancé (3 replis chacune). Tant
qu'il n'est pas revenu, parler de soumission serait prématuré : une PSH à 100
rendrait ce 0,786 sans objet.

### Et la décision qui ne m'appartient pas

Si le criblage passe, nous aurons pour la première fois une séquence
soumissible dont l'espérance bat le champ. Soumettre engage du TAO, et
l'instruction debout depuis le début est qu'on n'en engage pas pendant la phase
de recherche. Je ne soumets rien et je ne pousse pas : je donnerai les chiffres,
la décision est à l'utilisateur.

### Le criblage du candidat (§33)

| `design_iiptm` | PSH | A3D moyenne | A3D pic | soumissible |
|---|---|---|---|---|
| **0,78889** | **88,9** | **−0,840** | **1,235** | **oui** |
| 0,74191 | 90,4 | −0,788 | 1,262 | oui |
| 0,70373 | 95,8 | −0,827 | 1,161 | oui |
| 0,74261 | 101,0 | −0,778 | 1,247 | non, PSH |
| 0,73982 | 110,4 | −0,737 | 1,238 | non, PSH |
| 0,73583 | 106,5 | −0,988 | 1,205 | non, PSH |

**3 sur 6**, conforme au une-sur-trois du §17. Le candidat passe avec de la
marge sur les trois portes : PSH 88,9 contre une rupture mesurée à 96, A3D
−0,840 et 1,235 contre des plafonds de −0,680 et 1,580.

Une mise en garde sur un chiffre que le journal de preflight affiche : il
annonce « 49 % de risque de refus » à PSH 88,9. Ce n'est pas un signal PSH, mais
le **taux de base des CDR3 longs**, dont 52 % seulement des soumissions sont
acceptées toutes causes confondues (§3ter). Le discriminant mesuré est la
rupture à 96.

État du candidat, porte par porte :

| | |
|---|---|
| `design_iiptm` | 0,786 ± 0,010, deux mesures, deux sessions |
| PSH | 88,9 ✓ |
| Aggrescan3D | −0,840 / 1,235 ✓ |
| portes déterministes (gène V, longueur, cystéines, marqueurs, identité germinale) | ✓ |
| **unicité contre l'archive** | **vérification en cours** |

L'unicité est la seule porte qui ait pu changer depuis la génération, et c'est
exactement la réserve 3 du §29 : l'archive grossit, un candidat généré il y a
quelques jours peut être devenu doublon de ce que quelqu'un a soumis depuis. Un
doublon est un rejet certain, donc cette vérification précède toute décision.

### L'unicité : le candidat passe, mais les marges sont minces par construction

Archive re-téléchargée, **15 803 séquences** (contre 14 424 dans notre jeu de
données — 1 379 sont arrivées depuis). Verdict complet des portes, archive
incluse : **aucun problème, soumissible**.

Le plus proche voisin de l'archive est **le parent lui-même** :

| | candidat contre parent | seuil de `is_duplicate` | marge |
|---|---|---|---|
| identité entière | 0,944 | ≥ 0,95 → doublon | **0,006** |
| CDR3 | 0,778 | ≥ 0,85 → doublon | 0,072 |
| CDR3 avec CDR1 et CDR2 ≥ 0,90 | 0,778 (CDR1 et CDR2 à 1,000) | ≥ 0,80 → doublon | **0,022** |

**Ces marges sont structurelles, pas accidentelles** : 7 substitutions sur 125
résidus, dont 4 dans le CDR3, c'est exactement la perturbation minimale que la
règle autorise. Nous sommes donc tout juste à l'intérieur de deux des trois
seuils, par conception.

Le risque que cela crée, et il faut le nommer : **tout écart entre notre
reproduction des fonctions de similarité et celle du validateur fait basculer le
verdict.** Six millièmes de marge sur l'identité entière, cela ne pardonne rien.
Notre `whole_identity` utilise `difflib.SequenceMatcher`, notre
`cdr_similarity` reproduit le moteur de recherche de NOVA-nanobody-filter ; les
deux ont été écrites d'après le code du validateur, mais jamais confrontées à
un verdict d'unicité réel sur une séquence à 0,944.

### L'atténuation possible, et son coût

Une ou deux substitutions de **charpente** en plus élargiraient la marge sans
toucher le CDR3 :

| substitutions totales | identité entière | marge sous 0,95 |
|---|---|---|
| 7 (le candidat) | 0,944 | 0,006 |
| 8 | 0,936 | 0,014 |
| 9 | 0,928 | 0,022 |

La charpente ne contacte pas la cible (§15) et nos déplacements vont vers le
consensus de l'archive, ce qui remonte l'identité germinale. Mais **changer la
séquence change le tirage** : le §26 a établi que le score est une propriété de
la séquence, donc une variante à 8 substitutions est un nouveau tirage dont la
liaison doit être re-mesurée. Coût : une séquence, ~2 min de GPU, et le risque
de retomber plus bas que 0,786.

### L'état final, et la décision

| porte | valeur | verdict |
|---|---|---|
| `design_iiptm` | 0,786 ± 0,010, deux sessions | **5,1 époques sur 6** en espérance |
| PSH | 88,9 | ✓ (rupture à 96) |
| Aggrescan3D | −0,840 / 1,235 | ✓ (plafonds −0,680 / 1,580) |
| portes déterministes | — | ✓ |
| unicité, archive de 15 803 | 0,944 / 0,778 | ✓ **à 0,006 du seuil** |

Pour la première fois en onze soumissions rejetées et sept hypothèses, toutes
les portes sont franchies et mesurées, et le score est répliqué. Soumettre
engage du TAO, ce que l'instruction debout exclut en phase de recherche : **je
ne soumets rien et je ne pousse pas.** Les deux options honnêtes sont de
soumettre ce candidat tel quel, en acceptant les 6 millièmes de marge sur une
fonction de similarité que nous n'avons jamais confrontée au validateur, ou de
dépenser 2 minutes de GPU sur une variante à 8 substitutions pour tripler cette
marge, au risque que son score retombe.

## 34. « Serons-nous capables de reproduire ce résultat régulièrement ? » — Non, et voici le taux

Sur **155 perturbations minimales distinctes** mesurées à CDR3 = 18, en six
sessions :

| seuil | atteint | taux |
|---|---|---|
| ≥ 0,772 (médiane des gagnants) | **1/155** | 1 sur 155 |
| ≥ 0,780 (notre candidat) | 1/155 | 1 sur 155 |
| ≥ 0,785 (meilleur gagnant récent) | **0/155** | — |

Médiane de nos tirages **0,602**, 90ᵉ percentile 0,705, maximum 0,783. **Le
candidat n'est pas un résultat typique de la méthode : c'est son maximum sur
155.**

### La nuance qui est tout le sujet

Nous savons comment nous l'avons **tiré**, pas comment le **produire**. La
procédure est connue — perturber au minimum légal un parent du haut de
l'archive — mais rien ne prédit lequel des 155 tirages sera le bon :

| prédicteur | mesure |
|---|---|
| PSH, A3D, identité V | ρ ≈ 0 (§19) |
| positions substituées | p = 0,894 et 0,944 (§31, §33) |
| distance de perturbation | le minimum est le meilleur, et il plafonne (§21) |
| classifieur sur la séquence seule | AUC 0,61 contre 0,52 au hasard (§29) |

Nous savons échantillonner une distribution dont le 99,4ᵉ percentile est ce
qu'il faut. C'est très différent de savoir fabriquer ce percentile.

### Le prix de « régulièrement »

| tirages par couple | P(au moins un ≥ 0,772) | GPU |
|---|---|---|
| 31 (ce que l'époque permettait) | 18 % | 0,8 h |
| 100 | 48 % | 2,7 h |
| 150 | 62 % | 4,0 h |
| **250** | **80 %** | **6,7 h** |

Le §29 tablait sur 100 candidats par couple avec un taux estimé sur un seul
lot ; le taux réel sur 155 tirages demande **~250 par couple**. Pour les 8
couples les plus fréquents — 32 % des époques — cela fait **~54 h de GPU**, soit
deux semaines de quota Kaggle. Faisable sans louer un seul GPU, mais ce n'est
pas « régulièrement » : c'est une campagne, et elle couvre un tiers des époques.

### Ce que cela ne dit pas

Que le candidat soit le maximum de 155 ne le rend pas suspect : il a été
**répliqué** (0,783 puis 0,789, dans deux sessions dont les témoins étaient à
0,782 et 0,786), ce qui est précisément le contrôle que le §28 exigeait et que
le 0,770 du §26 avait échoué. Sa valeur est solide ; c'est sa **fréquence
d'apparition** qui est faible.

## 35. Le lot v8 : élargir la marge d'unicité sans toucher au CDR3

La marge de 0,006 du §33 tient à la construction : 7 substitutions sur 125
résidus est le minimum légal. L'élargir demande une substitution de plus, et la
forme la plus conservatrice est de partir **du candidat lui-même** et d'ajouter
**un seul** mouvement de charpente vers le consensus de l'archive — le CDR3 ne
bouge pas d'un résidu, et la charpente ne contacte pas la cible (§15).

| | identité au parent | marge sous 0,95 |
|---|---|---|
| candidat, 7 substitutions | 0,944 | 0,006 |
| **variantes v8, 8 substitutions** | **0,936** | **0,014** |

Onze mouvements franchissent les portes déterministes. La vérification
d'unicité en a écarté un, et la façon dont elle l'a fait vaut d'être notée :
**`FR330P->A` remonte l'identité à 0,952 au lieu de la baisser**, parce que ce
mouvement ramenait vers le parent un résidu que le candidat avait déjà changé —
donc une variante à 6 substitutions, pas 8, et un doublon. Ajouter une
substitution n'éloigne pas toujours : il faut le mesurer séquence par séquence.
**10 sur 11** passent toutes les portes, archive de 15 803 incluse.

### Le lot

39 séquences, ~1 h de GPU (`kaggle_v8.ipynb`, run neuf `v8`) :

| bras | n | rôle |
|---|---|---|
| v8 | 10 | les variantes à 8 substitutions |
| candidat_7subs | 1 | **le candidat lui-même, témoin interne** |
| recontrol | 1 | le parent |
| frais | 27 | nouvelles perturbations minimales |

Le témoin interne n'est pas optionnel : sans lui, un 0,77 sur une variante v8 ne
se distingue pas d'une session médiocre. Il a deux mesures antérieures (0,783 et
0,789) et doit retomber près de 0,786.

Les 27 tirages frais ne sont pas du remplissage. Ils portent l'échantillon du
§34 de 155 à 182, et c'est sur ce taux — 1 sur 155 aujourd'hui — que repose la
décision d'engager ou non une campagne de deux semaines de quota Kaggle. Un
taux estimé sur un seul lot est exactement l'erreur du §26.

### Ce que le lot tranchera

1. **Une variante v8 garde-t-elle 0,786 ?** Si oui, nous aurons un candidat à
   0,014 de marge au lieu de 0,006, pour le même score. Si toutes retombent vers
   0,70, c'est que même un résidu de charpente compte, ce qui contredirait le
   §15 et serait en soi un résultat.
2. **Le taux de 1/155 tient-il ?** 27 tirages de plus ne le trancheront pas seuls,
   mais ils le resserrent, et un second tirage ≥ 0,772 dans ce lot le ferait
   passer de 1/155 à 2/182.

## 36. Le lot v8 : la marge double sans rien coûter, et la charpente n'est pas neutre partout

39 séquences, 58,8 min de GPU.

### Les deux témoins

| | avant | ici | écart |
|---|---|---|---|
| parent | 0,781 | 0,78333 | +0,002 |
| **candidat à 7 substitutions** | 0,786 (2 mesures) | **0,77907** | −0,007 |

La session est bonne, et le candidat tient sa valeur en **troisième mesure** :
0,78312 / 0,78889 / 0,77907, moyenne **0,78369**, écart-type 0,0049.

### Les dix variantes à 8 substitutions

| mouvement | `design_iiptm` | contre le candidat de cette session |
|---|---|---|
| **FR210Q→K** | **0,78405** | +0,005 |
| **FR320S→T** | 0,78354 | +0,004 |
| FR213D→E | 0,78125 | +0,002 |
| FR11Q→E | 0,77967 | +0,001 |
| FR12Q→V | 0,77805 | −0,001 |
| FR211Q→E | 0,76996 | −0,009 |
| FR112I→V | 0,75742 | −0,022 |
| FR123S→A | 0,74666 | −0,032 |
| FR31E→D | 0,73301 | −0,046 |
| FR217R→D | 0,62448 | **−0,155** |

**4 sur 10 égalent ou dépassent le candidat**, pour une marge d'unicité
**doublée** (0,014 contre 0,006). L'objectif du lot est atteint.

### Ce que cela nuance du §15

Le §15 posait que la charpente ne contacte pas la cible. **C'est vrai pour cinq
mouvements sur dix et faux pour deux** : `FR217R→D` coûte 0,155 et `FR31E→D`
0,046. Une position de charpente peut donc compter autant qu'une position de
CDR3 — et rien ne le laissait prévoir avant la mesure. C'est la raison pour
laquelle dix mouvements ont été mesurés plutôt qu'un seul choisi par
raisonnement.

### Le taux de reproductibilité se confirme, et il ne s'améliore pas

Les 27 tirages frais : médiane 0,547, maximum **0,736**, **0 au-dessus de
0,772**, 9 effondrées sur 27 (33 %).

**Le taux passe de 1/155 à 1 sur 182.** Le §34 tenait sur un seul tirage
exceptionnel ; 27 tirages de plus ne l'ont pas répété. La réponse à « saurons-
nous reproduire ce résultat régulièrement » reste non, et se durcit.

### L'état des deux candidats

| | mesures | valeur | marge d'unicité |
|---|---|---|---|
| 7 substitutions | **3** | 0,78369 ± 0,005 | 0,006 |
| 8 substitutions (FR210Q→K) | **1** | 0,78405 | **0,014** |

Les deux valent ~0,784 ; la différence est l'information qu'on en a. Le §28 a
montré qu'une mesure unique d'un maximum sélectionné redescend — c'est
exactement le cas du v8, choisi comme le meilleur de dix. **Avant de préférer le
v8, il faut le répliquer**, et mesurer sa PSH et son A3D, qui ne le sont pas
(criblage lancé sur les quatre meilleures).

En espérance, le candidat à 7 substitutions gagnerait **4,9 époques sur 6** à
CDR3 = 18.

## 37. Le lot v9 : répliquer, et tester si deux mouvements bon marché sont additifs

22 séquences, ~33 min (`kaggle_v9.ipynb`, run neuf `v9`).

| bras | n | rôle |
|---|---|---|
| v8 | 10 | **seconde mesure** des dix variantes à 8 substitutions |
| v9 | 10 | **première mesure** des combinaisons à 9 substitutions |
| candidat_7subs | 1 | quatrième mesure, témoin interne |
| recontrol | 1 | le parent |

### Pourquoi répliquer les dix et pas seulement les quatre meilleures

Le §36 classe dix mouvements de charpente sur **une mesure chacun**, et c'est
sur ce classement que repose la nuance au §15 (la charpente compte parfois :
`FR217R→D` coûte 0,155). Avec σ = 0,014 par mesure, un classement d'écarts de
0,005 entre les quatre premiers **ne veut rien dire**. Les dix repassent donc,
ce qui donne deux mesures à chacun et rend le classement lisible — y compris
l'écart de 0,155, qui lui est trop grand pour être du bruit mais mérite d'être
confirmé.

### Les variantes à 9 substitutions

Les cinq mouvements ayant coûté moins de 0,01 au §36 sont combinés deux à deux :
dix paires, **toutes** franchissant les portes et l'unicité.

| | identité au parent | marge sous 0,95 |
|---|---|---|
| 7 substitutions | 0,944 | 0,006 |
| 8 substitutions | 0,936 | 0,014 |
| **9 substitutions** | **0,928** | **0,022** |

La question est l'**additivité** : deux mouvements qui coûtent chacun moins de
0,01 coûtent-ils moins de 0,02 ensemble ? Rien ne le garantit — le §36 vient
justement de montrer qu'on ne peut pas raisonner sur la charpente sans mesurer.
Si c'est additif, nous aurons un candidat à **0,022 de marge**, presque quatre
fois celle du candidat actuel, pour le même score.

### Ce qui ne change pas

Le taux de reproductibilité reste 1 sur 182 (§36) : ce lot ne cherche pas un
meilleur candidat, il qualifie celui que nous avons. Et aucune de ces variantes
n'a de PSH ni d'A3D mesurées — le criblage des quatre meilleures v8 tourne, et
celui des v9 devra suivre avant toute décision.

## 38. Le classement du §36 était du bruit, et l'additivité tient

22 séquences, 34,7 min. Témoins : parent 0,77996, candidat à 7 substitutions
0,78077 en **quatrième** mesure (0,78312 / 0,78889 / 0,77907 / 0,78077, moyenne
**0,78296**).

### Les dix v8, deuxième mesure : retour à la moyenne

| mouvement | 1ʳᵉ | 2ᵉ | écart | moyenne |
|---|---|---|---|---|
| FR210Q→K | 0,78405 | 0,77538 | −0,009 | 0,77971 |
| FR320S→T | 0,78354 | 0,77268 | −0,011 | 0,77811 |
| FR213D→E | 0,78125 | **0,78659** | +0,005 | **0,78392** |
| FR11Q→E | 0,77967 | 0,77870 | −0,001 | 0,77919 |
| FR12Q→V | 0,77805 | 0,76949 | −0,009 | 0,77377 |
| FR211Q→E | 0,76996 | 0,78157 | +0,012 | 0,77577 |
| FR112I→V | 0,75742 | 0,77772 | +0,020 | 0,76757 |
| FR123S→A | 0,74666 | 0,76662 | +0,020 | 0,75664 |
| FR31E→D | 0,73301 | 0,76037 | +0,027 | 0,74669 |
| **FR217R→D** | 0,62448 | **0,61423** | −0,010 | **0,61936** |

Les quatre « meilleures » du §36 perdent 0,004 en moyenne, les cinq suivantes
gagnent **+0,014**. C'est un retour à la moyenne de manuel : **le classement du
§36 était du bruit**, et c'est exactement pourquoi les dix ont été répliquées
plutôt que les quatre.

**La nuance au §15 doit donc être corrigée à son tour.** Le §36 concluait qu'une
position de charpente « peut compter autant qu'une position de CDR3 » en
s'appuyant sur deux cas, `FR217R→D` (−0,155) et `FR31E→D` (−0,046). Seul le
premier survit : `FR31E→D` remonte à 0,760 et rejoint le peloton. L'énoncé juste
est plus étroit : **neuf mouvements de charpente sur dix sont gratuits, un seul
coûte, et il coûte beaucoup** (0,619 contre 0,783, confirmé sur deux mesures).

### Les dix combinaisons à 9 substitutions : l'additivité tient

| | |
|---|---|
| ≥ 0,75 | **9/10** |
| maximum | **0,78604** (FR210Q→K + FR320S→T) |
| exception | FR320S→T + FR213D→E à 0,69146 |

Deux mouvements gratuits restent gratuits ensemble, à une paire près sur dix.
**Un candidat à 9 substitutions garde donc le score avec une marge d'unicité de
0,022 — presque quatre fois celle du candidat actuel.**

### L'état des trois candidats

| | substitutions | mesures | moyenne | marge d'unicité |
|---|---|---|---|---|
| candidat initial | 7 | **4** | 0,78296 | 0,006 |
| FR213D→E | 8 | **2** | 0,78392 | 0,014 |
| FR210Q→K + FR320S→T | 9 | 1 | 0,78604 | **0,022** |

Les trois valent la même chose à 0,003 près, soit un cinquième de σ. **Le choix
ne se joue donc pas sur le score mais sur la marge et sur le nombre de
mesures.** Le v9 est le meilleur de dix et n'a qu'une mesure : le §38 vient de
montrer, sur les v8, ce que vaut un classement à une mesure.

### Ce qui manque toujours

Aucune PSH ni A3D sur les v8 et les v9. Le criblage des quatre meilleures v8 a
été lancé et doit être relevé ; celui des v9 reste à faire. Sans lui, aucun de
ces candidats n'est soumissible, quel que soit son `design_iiptm`.

### Le criblage des v8 : la charpente est gratuite pour la liaison, pas pour la PSH

| variante | `design_iiptm` (2 mesures) | PSH | A3D moy / pic | soumissible |
|---|---|---|---|---|
| **FR213D→E** | **0,78392** | **92,8** | −0,793 / 1,298 | **oui** |
| **FR210Q→K** | 0,77971 | **92,3** | −0,784 / 1,231 | **oui** |
| FR320S→T | 0,77811 | 97,0 | −0,827 / 1,302 | non, PSH |
| FR11Q→E | 0,77919 | 98,5 | −0,855 / 1,228 | non, PSH |

**2 sur 4 survivent**, et le résultat important n'est pas le compte : le
candidat à 7 substitutions a une PSH de **88,9**, et **un seul** mouvement de
charpente la fait monter à **92–98**.

Le §38 avait conclu que neuf mouvements de charpente sur dix sont gratuits. Ils
le sont **pour la liaison**. Ils ne le sont pas pour la développabilité : ajouter
un résidu de charpente coûte 3 à 10 points de PSH, et deux des quatre testés
franchissent la rupture de 96. L'arbitrage « élargir la marge d'unicité ne coûte
rien » du §35 est donc faux tel qu'il était écrit — **il ne coûte rien en
liaison, et il coûte en PSH.**

Conséquence directe pour les v9, à 9 substitutions : si un mouvement coûte
3 à 10 points de PSH, deux pourraient en coûter le double, et partir de 88,9
laisse peu de place avant 96. **Le criblage des v9 décidera donc s'il existe
seulement un candidat à 0,022 de marge**, et il est le premier travail de la
nuit.

## 39. Le tableau de décision : cinq candidats franchissent toutes les portes

### Le criblage des v9

| combinaison | `design_iiptm` | PSH | A3D moy / pic | soumissible |
|---|---|---|---|---|
| FR210Q→K + FR320S→T | **0,78604** | 98,1 | −0,813 / 1,236 | **non, PSH** |
| FR210Q→K + FR11Q→E | 0,78200 | 93,1 | −0,784 / 1,229 | **oui** |
| FR210Q→K + FR213D→E | 0,78071 | 93,4 | −0,825 / 1,272 | **oui** |
| FR213D→E + FR12Q→V | 0,77869 | 103,5 | −0,775 / 1,311 | non, PSH |
| FR210Q→K + FR12Q→V | 0,77814 | 96,1 | −0,851 / 1,242 | non, PSH |

**2 sur 5**, et la meilleure sur la liaison tombe sur la PSH. La prédiction
écrite avant le criblage — « deux mouvements devraient mener à 95–108, j'en
attends la majorité perdue » — tient dans sa direction : trois des cinq sont
perdues, et les trois recalées sont bien entre 96,1 et 103,5.

**Et la PSH n'est pas additive**, ce qui n'était pas prévisible : `FR11Q→E`
**seul** donne 98,5 et serait refusé, mais `FR210Q→K + FR11Q→E` donne **93,1**
et passe. La PSH est une mesure structurelle du repliement, pas une somme de
contributions par résidu — donc on ne peut pas prédire la PSH d'une combinaison
à partir de celle de ses parties, dans un sens comme dans l'autre.

### L'unicité, re-vérifiée sur l'archive à 15 803 séquences

Les cinq candidats passent, aux marges attendues par construction : 0,006 à
7 substitutions, 0,014 à 8, 0,022 à 9.

### Le tableau de décision

| candidat | subs | mesures | moyenne | PSH | A3D moy / pic | marge | époques (esp.) |
|---|---|---|---|---|---|---|---|
| 7 subs, candidat initial | 7 | **4** | 0,78296 | **88,9** | −0,840 / 1,235 | 0,006 | 4,8/6 |
| **8 subs FR213D→E** | 8 | 2 | **0,78392** | 92,8 | −0,793 / 1,298 | 0,014 | **4,9/6** |
| 8 subs FR210Q→K | 8 | 2 | 0,77971 | 92,3 | −0,784 / 1,231 | 0,014 | 4,5/6 |
| **9 subs FR210Q→K+FR11Q→E** | 9 | 1 | 0,78200 | 93,1 | −0,784 / 1,229 | **0,022** | 4,7/6 |
| 9 subs FR210Q→K+FR213D→E | 9 | 1 | 0,78071 | 93,4 | −0,825 / 1,272 | 0,022 | 4,6/6 |

### Ce qui départage, en trois phrases

Les cinq valent le même score : 0,780 à 0,784, soit un écart total de 0,004 pour
un σ de 0,014 par mesure — **l'espérance d'époques gagnées va de 4,5 à 4,9 sur
6, ce qui n'est pas une différence sur laquelle décider.** Le choix se joue donc
sur deux axes seulement : la **marge d'unicité**, où 0,006 ne pardonne aucun
écart entre notre reproduction des fonctions de similarité et celle du
validateur — une reproduction écrite d'après le code mais **jamais confrontée à
un verdict réel** —, et le **nombre de mesures**, où le candidat initial en a
quatre contre une seule pour les v9.

Le compromis se lit directement : plus de substitutions donne plus de marge
d'unicité **et** une PSH plus élevée (88,9 → 92,8 → 93,1), donc moins de marge
sur la développabilité. Les deux risques sont opposés, et aucun n'est nul.

### Ce que ce tableau ne dit pas

Le taux de production reste de l'ordre de **1 sur 180** (§36). Les cinq
candidats sont des descendants d'un **seul** tirage réussi — ce ne sont pas cinq
succès indépendants, mais une famille issue du même parent, par cinq chemins de
charpente différents. **Nous avons un candidat, pas une méthode.**

## 40. Le taux, consolidé — et son intervalle interdit de dimensionner la campagne

### Ce qui compte comme tirage

287 lignes de mesure sur huit sessions. Il faut en retirer trois choses :
les 25 **réplicats** (une re-mesure n'est pas un tirage), les 24 **parents**, et
les 20 variantes des bras **v8 et v9** — celles-là sont des dérivées du *même*
tirage réussi par cinq chemins de charpente, et les compter comme indépendantes
gonflerait le taux de 11 %.

**182 tirages indépendants** à CDR3 = 18, 4 substitutions.

| seuil | atteint | taux | IC95 |
|---|---|---|---|
| ≥ 0,750 | 3/182 | 1,65 % | 0,34 – 4,74 % |
| ≥ 0,772 (médiane des gagnants) | **1/182** | **0,55 %** | **0,01 – 3,02 %** |
| ≥ 0,780 | 1/182 | 0,55 % | 0,01 – 3,02 % |

Médiane 0,601, 90ᵉ percentile 0,701, maximum 0,783. Le chiffre du §36 est donc
confirmé après avoir proprement exclu les réplicats et les dérivées.

### L'intervalle est ce qui compte, et il est large

Pour 80 % de chance d'obtenir un candidat ≥ 0,772 :

| hypothèse de taux | candidats par couple | GPU |
|---|---|---|
| borne haute de l'IC (3,02 %) | **52** | 1,4 h |
| estimation ponctuelle (0,55 %) | **292** | 7,8 h |
| borne basse (0,01 %) | 16 000 | 429 h |

Un facteur **300** entre les bornes. **La campagne du §29 — 54 h pour huit
couples, soit ~250 candidats chacun — suppose implicitement la borne haute**, et
rien dans nos données ne la justifie plus que l'estimation ponctuelle, qui
demanderait 7,8 h par couple, soit 62 h pour huit.

Le seuil de 0,75 est mieux contraint (1,65 %, IC 0,34–4,74 %) mais il ne gagne
aucune des six époques récentes à CDR3 = 18 (§28) : il n'aide pas à décider.

**Conclusion sur la campagne : elle n'est pas dimensionnable avec ce que nous
savons.** Le seul moyen de resserrer l'intervalle est d'accumuler des tirages,
et il faudrait environ 400 tirages de plus pour distinguer 0,5 % de 2 %.

## 41. État de l'enquête : ce qui est établi, ce qui reste supposé

### Les sept hypothèses

| # | hypothèse | verdict | répliqué ? |
|---|---|---|---|
| 1 | le GPU comme puissance de calcul | réfutée (§12bis) | oui |
| 2 | les donneurs | réfutée (§10) | oui |
| 3 | la distance à l'archive | réfutée (ρ = −0,07) | oui |
| 4 | la surface de conception (3 CDR) | réfutée (§12) | une seule session |
| 5 | la recombinaison | réfutée (§19) | oui |
| 6 | la distance de perturbation | réfutée (§21) | oui |
| 7 | le volume de sélection | soutenue (§26), ampleur démentie (§28) | oui |

### Mesuré et reproduit

- Les portes du validateur, dont la rupture PSH à **96** et non 89 (§3ter).
- Le J-gene chimérique, 58 % de la production récupérée (§7ter).
- L'effondrement est une propriété de la séquence, pas du tirage (§26, réplicats
  3/3), mais son **taux** varie de 0 à 31 % selon la session (§33, χ² p = 0,0006).
- Le repliement se dégrade continûment, l'interface s'effondre d'un coup (§23).
- Neuf mouvements de charpente sur dix sont gratuits **pour la liaison** (§38)
  et aucun ne l'est pour la PSH (§39).
- Le taux de production : 1 sur 182, IC 0,01–3,02 % (§40).

### Ce qui ne repose encore que sur une mesure ou un lot

À traiter comme des hypothèses, le §36 et le §38 ayant montré **deux fois** qu'un
classement à une mesure est du bruit :

1. **Les deux candidats à 9 substitutions n'ont qu'une mesure de liaison.** Leur
   place dans le tableau du §39 n'est pas acquise.
2. **La non-additivité de la PSH** (§39) repose sur un seul couple :
   `FR11Q→E` seul à 98,5 contre 93,1 en combinaison. Un seul contre-exemple
   suffit à établir la non-additivité, mais pas à en dire l'ampleur.
3. **`FR217R→D` coûte 0,155** sur deux mesures — c'est le seul mouvement de
   charpente coûteux, et donc le seul point d'appui de l'idée qu'une position de
   charpente peut compter. Deux mesures d'une seule séquence.
4. **La conception *de novo* a été écartée sur une seule session** (§11, 20
   designs). Le verdict était net (0,21 contre 0,66 pour le champ) mais il n'a
   jamais été répété.
5. **Le gain d'époque est calculé contre six époques récentes** à CDR3 = 18.
   Six, c'est peu, et le champ progresse.

### Ce qui n'a jamais été testé

Les 46 autres couples (gène V, longueur) hors IGHV3-74/18, la campagne de
pré-calcul elle-même, et toute soumission réelle — nous n'avons donc **aucune
confrontation de notre reproduction des portes avec un verdict du validateur**,
ce qui est précisément le risque que porte la marge de 0,006.

## 42. Bilan de la nuit du 12 au 13 novembre

Quatre commits, aucune dépense, aucun GPU loué, aucune soumission, aucun
lancement Kaggle.

### Le tableau de décision

Cinq candidats franchissent **toutes** les portes — portes déterministes, PSH
sous la rupture mesurée de 96, Aggrescan3D sous les deux plafonds, et unicité
re-vérifiée sur l'archive à 15 803 séquences :

| candidat | subs | mesures | moyenne | PSH | A3D moy / pic | marge | époques (esp.) |
|---|---|---|---|---|---|---|---|
| candidat initial | 7 | **4** | 0,78296 | **88,9** | −0,840 / 1,235 | 0,006 | 4,8/6 |
| **FR213D→E** | 8 | 2 | **0,78392** | 92,8 | −0,793 / 1,298 | 0,014 | **4,9/6** |
| FR210Q→K | 8 | 2 | 0,77971 | 92,3 | −0,784 / 1,231 | 0,014 | 4,5/6 |
| FR210Q→K + FR11Q→E | 9 | **1** | 0,78200 | 93,1 | −0,784 / 1,229 | **0,022** | 4,7/6 |
| FR210Q→K + FR213D→E | 9 | **1** | 0,78071 | 93,4 | −0,825 / 1,272 | **0,022** | 4,6/6 |

Le candidat à 9 substitutions le mieux noté du §38 (0,78604) **est tombé** au
criblage, sur une PSH de 98,1.

### Ce qui départage

Les cinq valent le même score : 0,780 à 0,784, soit un écart total de 0,004 pour
un σ de 0,014 par mesure. **L'espérance d'époques gagnées va de 4,5 à 4,9 sur
6 : ce n'est pas une différence sur laquelle décider.**

Le choix se joue donc sur deux axes, et ils sont opposés. **La marge
d'unicité** : 0,006 pour le candidat initial ne pardonne aucun écart entre notre
reproduction des fonctions de similarité et celle du validateur — une
reproduction écrite d'après le code, **jamais confrontée à un verdict réel**, et
un doublon est un rejet certain. **La marge de développabilité** : plus on
ajoute de substitutions pour élargir la première, plus la PSH monte (88,9 → 92,8
→ 93,1) et plus on s'approche de la rupture de 96.

Mon avis, qui n'engage que mon jugement : **`FR213D→E`**. Il a deux mesures, la
meilleure moyenne, 0,014 de marge d'unicité et 3,2 points de marge sous la
rupture PSH — c'est le seul qui ne soit extrême sur aucun des deux axes. Les
deux candidats à 9 substitutions n'ont qu'une mesure, et le §36 comme le §38 ont
montré, deux fois, ce que vaut un classement à une mesure.

### Ce que ce tableau ne dit pas, et qu'il ne faut pas oublier

**Le taux de production est de 1 sur 182** (§40), IC95 de 0,01 % à 3,02 %. Les
cinq candidats ne sont pas cinq succès : ce sont cinq descendants d'un **seul**
tirage réussi, par cinq chemins de charpente différents. **Nous avons un
candidat, pas une méthode.**

Et l'intervalle sur ce taux est large d'un facteur 300, ce qui rend la campagne
de pré-calcul du §29 **non dimensionnable** en l'état : elle supposait
implicitement la borne haute. Il faudrait environ 400 tirages de plus pour
distinguer 0,5 % de 2 %.

### La décision

Soumettre engage du TAO, et l'instruction debout depuis le premier jour est
qu'on n'en engage pas en phase de recherche. **Rien n'a été soumis, et je ne
pousse pas à le faire.** Le tableau est là pour que l'arbitrage se fasse sur des
chiffres mesurés plutôt que sur une impression — c'est tout ce qu'il prétend
être.

## 43. La PSH n'était pas bruitée : notre DSSP ne tournait pas

### Ce que ce dépôt a soutenu pendant un mois, et qui est faux

Le post-mortem du 8 novembre concluait que la PSH ne se calibre pas : que notre
mesure et celle du validateur sont deux tirages de la même loi, que l'écart-type
d'une mesure vaut 4,01 pour un seuil à 89, et qu'il faut donc une courbe
empirique avec une rupture à 96 plutôt que le seuil lui-même. L'argument
reposait sur un fait apparemment solide : **44 à 62 % des soumissions que les
validateurs avaient *acceptées* dépassaient 89 sur notre échelle.** Un seuil que
la moitié des séquences acceptables franchissent n'est pas un seuil utilisable.

Ce fait était un artefact d'installation.

### Le mécanisme

TNP demande à Biopython une annotation DSSP du modèle NanoBodyBuilder2. Quand
elle échoue, TNP **n'arrête pas** : il imprime une ligne et continue avec une
PSH calculée autrement.

    ERROR: Failed to process output PDB. DSSP failed to produce an output

Le validateur construit DSSP depuis `cmbi/dssp`, qui est l'arbre autotools
**3.1.4**, et applique à TNP le correctif `process_pdb.py` de
NOVA-nanobody-filter. `mkdssp` 4.x — celui d'apt — refuse les modèles de
NanoBodyBuilder2, qui n'ont pas de ligne HEADER : il suppose du mmCIF et
échoue. La 3.1.4 les accepte dès que le correctif complète les lignes à
80 caractères, ce que son propre docstring dit en toutes lettres.

Sur un même candidat, la même structure lit **89,50 sans DSSP et 91,21 avec**,
et son drapeau de charge positive passe d'AMBER à GREEN. L'écart n'est pas un
décalage constant — c'est pourquoi aucune calibration ne reliait les deux
échelles, et c'est précisément ce que l'addendum du 8 novembre avait constaté
sans en trouver la cause.

### L'installation du validateur, reproduite depuis son propre code

Les épingles ne viennent pas de nos souvenirs mais de
`NOVA-nanobody-filter/install/requirements-minimal.txt` et de son `install.sh` :
TNP au commit `29dcac7`, ImmuneBuilder `0df4e2a`, pdbfixer v1.12,
biopython 1.86, le correctif `process_pdb`, DSSP depuis `cmbi/dssp`,
aggrescan3d-rs et sa FreeSASA 2.0.1. Plus le correctif du typo set/dict
d'ImmuneBuilder (`refine.py:346`), que `docs/DEPLOY.md` mentionnait déjà et qui
faisait revenir **3 séquences sur 37 sans aucun modèle**, donc sans PSH.

Détail qui a coûté une heure : TNP range ses modules dans un paquet de premier
niveau nommé `scripts`, comme nous. Un `PYTHONPATH` pointant sur ce dépôt masque
le sien, et TNP meurt sur `No module named 'scripts.region_definitions'` après
avoir profilé 0 séquence sur 100. `toolchain.tnp_env` retire l'entrée.

### Ce que ça mesure (results/psh_calibration_epoch26174.csv)

Les 36 soumissions acceptées à l'époque 26174, trois plis chacune :

| | |
|---|---|
| minimum | 70,3 |
| **médiane** | **83,5** |
| maximum | 88,9 |
| **au-dessus de 89** | **aucune, 0 sur 36** |

36 sur 36 sous le seuil que le validateur applique. Notre mesure est la leur.
Celle qu'il a refusée cette époque-là lit **88,0** ici contre **89,57** dans son
journal : 1,5 point plus bas, soit une raison de garder de la marge, pas de
douter de l'échelle.

**Le bruit, re-mesuré** sur 37 séquences × 3 plis : l'écart-type d'**une** mesure
est **1,74**, pas 4,01. Le 4,01 venait du chemin cassé et était deux fois trop
grand. La moyenne de k plis vaut 1,74/√k et le validateur tire une fois de son
côté, d'où

    P(le validateur lit au-dessus de 89) = Phi((m - 89) / (1,74 * racine(1 + 1/k)))

soit 1 % à m = 84,4, 18 % à 87,2 et 42 % à 88,6 pour k = 3.

### Validation externe : 296 soumissions dont on connaît le verdict

Les époques 26299-26306, un pli chacune
(`results/psh_field_26299_26306.csv`) :

| | n | médiane | > 89 |
|---|---|---|---|
| acceptées | 230 | 84,0 | 8 % |
| refusées | 66 | 87,1 | **36 %** |

AUC 0,671 (p = 1,1×10⁻⁵). Et le rapport des deux taux est maximal **au seuil
que le validateur utilise** : 1,7× à 85, 2,3× à 87, **4,6× à 89**. Les 8 %
d'acceptées au-dessus de 89 à un seul pli sont ce que le modèle de bruit prédit
(13 %, donc il est légèrement conservateur) : c'est exactement pourquoi on
moyenne des plis avant de décider.

Note : 36 % des refus seulement sont expliqués par la PSH. Le reste tombe sur
les autres portes, et nous ne savons pas lesquelles.

### Ce que ça change pour le mineur

Le plafond à 96 disparaît. La PSH n'est plus une porte mais un **coût** : on
mesure, on price, et on multiplie par la qualité prédite. Un plafond jette un
bon candidat qui est à un tirage de la ligne ; ignorer la PSH est ce qui a mis
onze soumissions à la poubelle.

Taille du vivier re-calibrée sur 300 candidats frais mesurés
(`results/psh_fresh_candidates_by_length.csv`) :

| CDR3 | PSH médiane | ≤ 87 | ≤ 89 | ≤ 96 |
|---|---|---|---|---|
| 8 | 85,0 | 62 % | 75 % | 95 % |
| 11 | 89,0 | 40 % | 50 % | 81 % |
| **18** | **100,6** | 10 % | **15 %** | 31 % |

L'ancienne table donnait 51 / 28 / 8,9 % pour la colonne ≤ 89. La forme survit —
le CDR3 long est là où le générateur souffre — mais tous les taux étaient
sous-estimés.

### Ce que ça change pour nos cinq finalistes du §42

Re-mesurés à 3 plis (`results/psh_candidates_v8v9.csv`), avec P(refus) calculée
pour un tirage du validateur :

| variante | PSH | P(refus) | design_iiptm | n | iiptm × P(passage) |
|---|---|---|---|---|---|
| **FR211Q→E** | **84,4** | **1 %** | 0,776 | 2 | **0,767** |
| FR12Q→V | 86,3 | 9 % | 0,774 | 2 | 0,707 |
| FR320S→T+FR213D→E | 85,0 | 2 % | 0,691 | 1 | 0,675 |
| FR210Q→K | 87,2 | 18 % | 0,780 | 2 | 0,638 |
| FR213D→E *(recommandé au §42)* | 88,6 | **42 %** | 0,784 | 2 | 0,453 |
| FR210Q→K+FR213D→E *(finaliste §42)* | 90,9 | **82 %** | 0,781 | 1 | 0,138 |
| FR217R→D | 93,7 | 99 % | 0,619 | 2 | 0,006 |

**La recommandation du §42 avait 42 % de chances d'être refusée avant d'être
notée, et l'un des cinq finalistes 82 %.** `FR211Q→E`, que le classement du §39
écartait — classement mesuré sur l'échelle cassée —, lit 84,4 pour 1 % de risque
et une liaison équivalente à 0,004 près, soit un quart d'un écart-type de mesure.

Le §39 concluait qu'« aucun mouvement de charpente n'est gratuit pour la PSH ».
C'est à refaire entièrement : les valeurs sur lesquelles il classait ne sont pas
celles du validateur.

### Une hypothèse testée et réfutée : biaiser la génération vers le polaire

Le docstring de `nova_lab/cdr3_generators.py` accuse depuis le début l'alphabet
uniforme des mutations — « 40 % des acides aminés sont hydrophobes alors que les
CDR3 naturels sont enrichis en résidus polaires ». `HydrophilicGenerator`
applique le biais des deux côtés : les positions à muter sont tirées vers les
résidus les plus hydrophobes du donneur, et le remplacement est tiré loin d'eux.

Sur la séquence il fait exactement ça : l'hydropathie médiane du CDR3 baisse de
plus de 2 unités. Sur la quantité qui décide l'époque, il ne fait presque rien.
100 candidats par bras à CDR3 18, mêmes donneurs, même graine
(`results/psh_generator_ab_cdr3_18.csv`) :

| bras | PSH médiane | ≤ 89 |
|---|---|---|
| chimeric | 100,2 | 14/100 |
| hydrophilic | 97,4 | 13/100 |

2,8 points de médiane (Mann-Whitney p = 0,064) et **aucun changement du taux qui
compte** (Fisher p = 1,00). La PSH est un score de *patch* sur une structure
prédite : ce qui compte est quels résidus se retrouvent exposés et voisins, pas
la somme d'hydropathie du CDR3.

La même conclusion sort d'une tentative de prédire la PSH depuis la séquence :
sur 596 séquences mesurées, un LightGBM sur les 147 descripteurs par région
atteint Spearman 0,62 toutes longueurs confondues — mais **0,25 à 0,31 à
longueur fixée**, et l'hydropathie brute des CDR 0,10 à 0,22. Le 0,62 ne dit que
« un CDR3 long est plus hydrophobe ». En croisant les distributions c'est pire :
entraîné sur le champ et appliqué à nos candidats, Spearman 0,25 — le même mode
d'échec que le modèle de substitution sur design_iiptm.

Le générateur est gardé, documenté comme réfuté, et n'est pas le défaut : sans
ça l'idée se retentera.

**Donc le levier à CDR3 long n'est pas un alphabet plus fin, c'est de mesurer
plus de candidats et de prendre la queue** : le meilleur de ces 100 lit 74,7. Le
vivier du mineur est dimensionné pour ça (67 candidats à CDR3 18 pour en trouver
une dizaine sous 89, soit 4 min de CPU sur un budget de 66). Si l'on y revient,
ce qu'il faut biaiser est l'exposition, qui demande la structure, c'est-à-dire
les 3,5 s qu'on cherchait à éviter.

### Six pannes silencieuses trouvées en chemin

Toutes de la même forme : un outil qui répond au lieu de refuser.

1. **Pas de `hmmscan`** → ANARCI rend `None` pour chaque séquence *sans lever*,
   et le jeu d'entraînement s'est construit sur **0 ligne sur 9 466** en
   affichant « saved ».
2. **DSSP absent ou 4.x** → une ligne d'erreur, puis une PSH qui n'est pas celle
   du validateur. Onze soumissions refusées sur une porte qu'on croyait cribler.
3. **Le typo set/dict d'ImmuneBuilder** → 3 séquences sur 37 sans modèle, donc
   sans PSH, ce qui ressemble à « rien à refuser ».
4. **`data/scores/` manquant** → `scripts.analyze_nanobodies` meurt, et c'est
   l'étape 2 de la mise en service.
5. **Le prébâti sans son `_X.npy`** → `bootstrap --build` jugeait le jeu de
   données complet, sautait sa reconstruction, et l'entraînement mourait. Le
   modèle de substitution était infabricable par le chemin documenté.
6. **`parse` sur un run vide** → KeyError trois cadres plus loin.

D'où `nova_lab/toolchain.py` (où sont les outils, et si ce sont les bons) et
`scripts/doctor.py`, qui les fait tourner sur un nanobody clinique et sort en
code 1 si un nombre manque. **Un outil présent et faux est pire qu'un outil
absent**, et c'est la leçon qui vaut pour tout ce dépôt : les six pannes
ci-dessus ont duré entre une semaine et un mois chacune, et aucune n'a jamais
produit de message d'erreur.

## 44. Qui refuse le champ, et ce que le PR 102 change vraiment

### Les portes, confrontées à 296 verdicts connus

`scripts.why_refused` passe les 296 soumissions des époques 26299-26306 — 230
acceptées, 66 refusées — dans toutes les portes que nous reproduisons, en
réutilisant les structures que TNP avait déjà construites
(`results/why_refused_26299_26306.csv`, un pli).

| porte | sur les refusées | sur les acceptées | lecture |
|---|---|---|---|
| déterministes (vérificateur hors ligne) | 1/66 | 4/230 | aucune séparation |
| **PSH > 89** | **24/66** | 18/230 | vraie cause de rejet |
| drapeau rouge TNP | 2/66 | 5/230 | aucune séparation |
| **Aggrescan3D au-delà d'un plafond** | **8/66** | 4/230 | vraie cause de rejet |

- **Nos portes expliquent 53 % des refus** (35 sur 66). La PSH en explique 36 %,
  l'A3D 12 %.
- **12 % des acceptées seraient refusées à tort** (28 sur 230) — sur un seul
  pli. C'est le chiffre qui justifie le moyennage : 18 des 28 sont des PSH, donc
  du bruit de tirage, et c'est exactement ce que `psh_screen` fait disparaître
  en repliant trois fois.
- Les portes déterministes et les drapeaux TNP **ne séparent rien** : ils tirent
  aussi souvent sur les acceptées. Ce ne sont pas les portes où le champ tombe.

### Les 47 % que rien n'explique, et une hypothèse

31 refus résistent à toutes nos portes. Leur PSH médiane est **83,4** — plus
basse que celle des acceptées (84,0). Ce ne sont donc pas des séquences
hydrophobes ni des séquences mauvaises sur ce que nous mesurons.

Et ils se groupent par époque : 9 à l'époque 26303, 8 à 26301, 6 à 26300, puis
4, 3, 1. Un refus qui dépend de l'époque et pas de la séquence ressemble à une
défaillance de la mesure, pas à un verdict. Le code amont le permet :

    rejected_sequences = [... if result.get("red_flags") or result.get("error")]

**Une séquence dont TNP *échoue* chez le validateur est refusée comme si elle
était mauvaise.** Et TNP échoue : le typo set/dict d'ImmuneBuilder nous a fait
perdre 3 séquences sur 37 ici (§43), et le filtre tourne avec
`max_concurrent_tnp = 16`, donc sous charge. Le groupement par époque est ce
qu'on attendrait d'un échec dépendant de la charge.

**Hypothèse, pas conclusion** : une part des refus du champ n'est pas un verdict
sur la séquence mais une panne chez le validateur. Si c'est vrai, notre taux
d'acceptation ne sera jamais de 100 % même avec une séquence parfaite, et il
n'y a rien à corriger de notre côté. Ce qui la testerait : les journaux W&B des
validateurs, où le refus PSH de l'époque 26209 avait déjà été lu en clair.

### Le PR 102 : la réaction change, les paramètres nanobody non

Annonce Discord du 13 novembre : fusion à l'époque 26331, et « pour les mineurs
qui ne font que du nanobody : ce changement affectera aussi la façon dont les
germlines et les longueurs de CDR3 valides sont choisies, du fait de la
consommation d'état du RNG ».

Lu dans `utils/challenge.py` sur la branche, c'est faux. Les deux sélections de
réaction sont de l'arithmétique sur `seed` et ne touchent jamais l'objet `rng` :

    ancienne : allowed_option = seed % total_reactions
    nouvelle : for group in rxn_groups: seed % len(group)

`rng` n'est consommé qu'ensuite, par le tirage du gène V puis celui de la
longueur. Vérifié sur 200 000 hashes : **zéro désaccord** sur le couple
(gène V, longueur CDR3).

**Nous n'avons pas à trancher qui a raison.** Une soumission nanobody n'est
acceptée que si son gène V et sa longueur correspondent à ceux de l'époque,
donc chaque soumission acceptée est un témoin des vrais paramètres.
`scripts.verify_challenge` les interroge : **8 époques sur 8, 230 soumissions
acceptées, accord parfait.** C'est un contrôle à retardement — l'époque doit
avoir été notée — donc il détecte une dérive en une époque, sans faire
confiance ni au Discord ni au dépôt amont. À lancer juste après 26331.

Les 11 filtres de molécules se resserrent en même temps (max_logp 7 → 5,5,
max_iodine 1 → 0, max_bromine 2 → 0, max_fluorine 9 → 3...). Nous ne minons pas
les molécules ; ils sont rangés dans `FILTERS_FROM_26331` à côté de `FILTERS`,
car `tests/test_parity.py` vérifie que les molécules archivées passent les
filtres d'aujourd'hui et l'archive a été construite sous les anciens.

## 45. Les deux faux refus de la calibration, et ce qu'ils ont révélé

Au §43 la calibration laissait **2 refus à tort sur 36** soumissions que les
validateurs avaient acceptées. Deux cas, deux causes différentes, et la seconde
a obligé à revoir la porte entière.

### Cas 1 : notre bug. Un seul pli décidait

La séquence refusée sur PPC lit, pli par pli : **1,2858 / 1,1078 / 1,1125**.
Le premier pli franchit la borne de 1,18, la moyenne des trois vaut **1,1687**
et passe.

Nous moyennions la PSH et l'Aggrescan3D sur les plis, et prenions **tous les
autres drapeaux du pli 1**. Sur une métrique dont l'écart-type est 0,037 — et
0,23 dans le pire cas — une lecture isolée à côté d'une borne ne dit rien.
Corrigé : toutes les métriques TNP sont moyennées, et une violation ne compte
qu'au-delà de deux écarts-types de la moyenne, la règle que l'A3D appliquait
déjà.

### Et au passage : ce ne sont pas les seuils de TNP qui décident

TNP imprime son propre verdict par métrique. **Le validateur ne l'utilise
pas** : son filtre prend les *valeurs* et les compare à sa propre
configuration, qui sont les seuils de TNP arrondis à deux décimales — dans les
deux sens.

| métrique | rouge de TNP | config du filtre | plus strict |
|---|---|---|---|
| C | 0,5634 – 1,6132 | 0,56 – 1,61 | le filtre, en haut |
| PSH | 73,4039 – 155,4699 | 73,4 – 155,47 | TNP, en haut |
| PPC | > 1,1853 | > 1,18 | le filtre |
| PNC | > 1,8838 | > 1,88 | le filtre |

Des centièmes, mais la porte est une falaise. `nova_lab/developability.py`
porte les seuils du filtre.

### Cas 2 : pas du bruit. Notre échelle lit bas

La séquence refusée sur la borne basse de PSH lit **69,50 / 71,59 / 69,72**,
soit 70,27 de moyenne, sous les 73,4 — aux trois plis. Le moyennage n'y change
rien. Le validateur l'a acceptée, donc il a lu au moins 73,4.

Trois estimations indépendantes disent la même chose :

| source | écart |
|---|---|
| la séquence refusée dont nous avons leur journal : nous 88,0, eux 89,57 | −1,6 |
| cette acceptée : nous 70,3, eux ≥ 73,4 | ≤ −3,1 |
| le seuil qui sépare au mieux leurs 296 verdicts sur notre échelle (Youden, époques 26299-26306) : 86,5 contre 89 | −2,5 |

**Notre PSH lit 2 à 3 points sous la leur.** Trois estimations faibles qui
s'accordent en signe et en ordre de grandeur, sur un axe où se tromper ne coûte
que dans un sens : sous-lire fait croire qu'un candidat est sûr alors que le
validateur le refuse.

Un biais de +2,0 est donc ajouté à notre lecture **partout où elle rencontre un
de leurs nombres** — la borne haute comme la rupture à 89. Et il faut l'ajouter
aux deux bornes : une échelle qui lit bas est fausse en haut (trop permissive)
*et* en bas (trop stricte). Ça explique aussi le cas 2 : 70,3 + 2,0 = 72,3, soit
*sur* la borne et non dessous.

**Résultat : 36 sur 36 des portes déterministes passent.** Zéro faux refus.

### Le bruit était sous-estimé, et de beaucoup

En relançant, `FR11Q→E` est passé de 88,1 à 83,2 — 4,8 points sur une moyenne
de trois plis, là où σ/√3 = 1,2 en prédisait trois fois moins. Cinq estimations
de l'écart-type d'UNE mesure :

| population | σ |
|---|---|
| les 36 acceptées, run 1 | 1,74 |
| les mêmes, run 2 | 2,49 |
| nos 20 candidats v8/v9, run 1 | 3,29 |
| les mêmes, run 2 | 2,70 |
| **entre ces deux runs de candidats** (apparié) | **2,85** |

La dernière ligne est celle qui tranche : une comparaison appariée entre runs,
donc indépendante des plis intra-run, **retrouve la même valeur**. Les plis sont
bien indépendants, et le 1,74 initial était le bas d'une estimation bruitée. Nos
candidats sont plus bruités que le champ accepté, ce qui se comprend — ils sont
plus haut sur la PSH, sur des structures moins stables.

σ est fixé à **2,7**. Et la conséquence mérite d'être dite franchement : **à
trois plis la moyenne a un écart-type de 1,6, donc cette porte sépare « sûr » de
« condamné » et ne sépare pas 84 de 86.** Deux runs indépendants des mêmes
20 candidats s'accordent à Spearman **0,64**. Il ne faut lire aucun classement
dans un écart d'un ou deux points.

### Le tableau de décision, sur 6 plis et deux runs

| variante | run 1 | run 2 | moyenne | P(refus) | design_iiptm | espérance |
|---|---|---|---|---|---|---|
| **FR211Q→E** | 84,4 | 83,4 | **83,9** | **14 %** | 0,776 | **0,664** |
| FR320S→T+FR213D→E | 85,0 | 84,7 | 84,9 | 23 % | 0,691 | 0,530 |
| FR11Q→E | 88,1 | 83,2 | 85,7 | 32 % | 0,779 | 0,528 |
| FR210Q→K | 87,2 | 86,8 | 87,0 | 50 % | 0,780 | 0,390 |
| FR213D→E *(recommandé au §42)* | 88,6 | 88,7 | 88,7 | 66 % | 0,784 | 0,267 |
| FR217R→D | 93,7 | 93,5 | 93,6 | 99 % | 0,619 | 0,008 |

`FR211Q→E` reste premier, et c'est le seul qui lise bas **aux deux runs**.
`FR11Q→E` a la meilleure liaison mais ses deux runs divergent de 4,8 points :
sa place n'est pas acquise.

Les chiffres de risque du §43 (1 % pour `FR211Q→E`, 42 % pour `FR213D→E`)
étaient optimistes de deux façons qui allaient dans le même sens : pas de
biais, et σ pris au plus bas. C'est 14 % et 66 %.

## 46. Le verdict du run réparé : le goulot n'est aucun de ceux qu'on surveillait

> **Révision du 15 novembre (§49) :** le « 0 victoire sur 8 » de cette section
> n'est pas mesurable — le bruit de BoltzGen (σ = 0,028) dépasse l'écart entre
> rivaux voisins (0,014). Les conclusions sur le *niveau* des candidats, qui
> reposent sur des médianes, tiennent ; le comptage de victoires est retiré.

64 séquences — 4 candidats et les 4 meilleures soumissions réelles de chacune
des époques 26299-26306 — notées dans un seul run BoltzGen, 92,2 min sur 2×T4,
64 sur 64 abouties (`results/scores_nanobody_repaired_26299_26306.csv`).

**Le résultat ne demande aucune interprétation :**

| | nos candidats | le champ |
|---|---|---|
| design_iiptm | min 0,141 · **médiane 0,416** · max 0,621 | min 0,615 · **médiane 0,718** · max 0,786 |
| interaction_pae | 10,9 – **15,6** – 23,9 | 8,8 – **9,2** – 10,9 |

**Huit époques sur huit : dernière place, zéro rival battu.** Pas une seule
époque où notre candidat passe devant ne serait-ce qu'une des quatre
soumissions portées dans le lot. Notre meilleur candidat des 64, à 0,621, ne
dépasse que le pire rival des 32, à 0,615.

Ce n'est pas une subtilité de classement : `design_iiptm` est une métrique
absolue, pas un rang, et l'écart entre les deux distributions ne se chevauche
presque pas. Aucune erreur d'ordonnancement n'explique un écart de 0,30.

### Ce que ce n'est pas

Trois suspects, écartés par les données du run lui-même.

**Ce n'est pas la porte PSH.** Sur nos 32 candidats, PSH contre design_iiptm
donne Spearman **+0,19 (p = 0,30)** — aucune relation, et si tendance il y a,
elle est inverse : la moitié la plus basse en PSH lie *mieux* (iiptm médiane
0,451 contre 0,307). Sélectionner sur l'acceptabilité ne détruit pas la
liaison. C'était l'hypothèse la plus inquiétante et elle est fausse.

**Ce n'est pas le sélecteur.** Un oracle parfait choisissant parmi nos
32 candidats perdrait quand même les 8 époques, puisque le meilleur de tout le
lot est sous la médiane du champ. Le classement ne peut pas réparer un vivier
dont le plafond est trop bas. (Le modèle de substitution fait d'ailleurs
Spearman +0,31 contre le rang réel, mieux que le −0,20 du §43 mais sans
importance ici.)

**Ce n'est pas le budget GPU.** À 1,4 min/séquence, une époque de 66 min
utiles permet d'en noter ~47. Nous en avons noté 4 par époque par choix de
protocole, pas par contrainte.

### Ce que c'est : le générateur ne lie pas

En approximation normale sur nos 32 tirages (moyenne 0,374, écart-type 0,171) :

| cible | p par candidat | candidats pour 80 % de chances |
|---|---|---|
| battre le pire des 4 rivaux (0,615) | 8,0 % | ~20 |
| atteindre la médiane du champ (0,718) | 2,3 % | ~71 |
| battre le meilleur rival (0,786) | 0,8 % | ~197 |

Des ordres de grandeur, pas des prévisions : 32 tirages ne contraignent pas une
queue. Mais ils disent la forme du problème, et **elle était déjà écrite au
§40** : le taux de production mesuré à IGHV3-74/18 était de **1 sur 182** pour
atteindre le niveau des gagnants. 32 tirages à ce taux donnent 0,18 gagnant
attendu. **Zéro est exactement le résultat prévu.**

### Donc les 14 à 29 % étaient faux, et voici ce qui diffère

Ce dépôt annonçait 14 à 29 % de victoires contre le champ complet, 48 % au
mieux de quatre. Cette mesure-ci les contredit directement, et c'est elle qu'il
faut croire : pipeline complet, champ actuel, même run pour les deux camps,
huit époques, sélection identique à celle du mineur.

Ce qui a changé depuis les balayages de 25996-26025 : **le champ**. Les
gagnants de cette époque-là lisaient 0,56 à 0,71 de `design_iiptm` ; ici le
*quatrième* rival de chaque époque lit 0,67 à 0,79, et QTTECH2023 occupe le
haut du tableau partout. Les anciens balayages comparaient en plus des
candidats issus de la lignée IGHV3-74/18 qu'on optimisait depuis des mois,
contre un champ plus faible. Rien de cela ne se généralise à huit couples
quelconques.

### Une réserve sur la mesure, qui ne sauve rien

L'accord de notre run avec le classement officiel des rivaux tombe à **33 paires
sur 47, soit 70 %**, contre 85-86 % sur les trois runs précédents. L'intervalle
à n = 47 est large (environ 55-82 %) et les quatre rivaux d'une époque sont très
serrés, ce qui rend l'ordonnancement plus dur — mais c'est plus faible, et il
faut le dire. La ligne « calibration: Spearman 0.33 » du journal, elle, est
l'artefact déjà documenté : elle compare des sommes de rangs d'époques
différentes, donc incomparables.

Cette réserve porte sur les *rangs*. Elle ne touche pas l'écart de 0,30 sur une
métrique absolue, qui est le fait central.

### Où cela laisse le projet

Le travail du 13 novembre valait : la porte de développabilité est désormais
fidèle, et les 8 candidats retenus seraient **acceptés** là où les onze
soumissions de production étaient toutes refusées. C'est un vrai acquis, et
c'était la condition nécessaire.

Elle n'est pas suffisante. Être accepté et finir dernier rapporte exactement
autant que d'être refusé : zéro. **Le goulot est la qualité de liaison du
générateur, et il est bien plus profond que ce dépôt ne le croyait.** Greffer un
CDR3 d'archive puis muter deux à quatre positions produit des séquences dont la
médiane est à 0,42 quand le champ est à 0,72. Aucun réglage de porte, de
sélecteur ou de vivier ne comble ça.

## 47. IgLM, et trois pistes mesurées à la place

La question posée après le verdict du §46 : IgLM ne serait-il pas la solution ?
La réponse courte est non, et la raison vaut d'être énoncée parce qu'elle
oriente tout le reste.

### IgLM répond à un problème que nous n'avons pas

`IgLMGenerator` est dans le dépôt depuis le début, avec ses tests de contrat
d'infilling — et **il n'a jamais produit une mesure** : aucune mention dans ce
fichier, aucun fichier de résultats. Ce qu'il ferait bien :

- il remplit un CDR3 **dans le contexte de la charpente**, là où notre
  générateur greffe un CDR3 issu d'une *autre* séquence sur une charpente
  étrangère. Le README accuse déjà cette chimère : « les CDR3 greffés ne créent
  pas assez de surface de contact — le transfert de charpente perd le contexte
  structural », delta_sasa au 14ᵉ percentile ;
- il produit des CDR3 **nouveaux**, donc la règle de nouveauté ne mord pas.

Et ce qui le disqualifie : **il est aveugle à la cible.** IgLM est entraîné sur
des répertoires d'anticorps naturels ; il ne sait rien de Q9NZQ7 ni de 8AOK. Il
produit un CDR3 plausible, pas un lieur. Or notre déficit est de 0,30 de
`design_iiptm` sur une cible précise. Aucun mécanisme ne fait qu'un modèle de
langage agnostique de la cible atteigne 0,72 là où notre greffage atteint 0,42 :
il changerait la loi du tirage, pas son espérance.

Il reste utile pour la développabilité — la campagne de production a été plombée
une fois par les *liabilities* (11ᵉ percentile) — mais ce n'est pas le goulot.

### Piste 1 : le plafond du vivier, et c'est lui qui décidait

Le fait central, mesuré sur le vivier que le rejeu a utilisé :

| | |
|---|---|
| séquences d'archive conformes **et** appelant le gène V de l'époque | médiane `design_iiptm` **0,22**, plafond 0,57 à 0,77 |
| le champ auquel on s'oppose | **0,72** |

On n'atteint pas 0,72 en partant de là, et la cause est connue : les bons
lieurs de l'archive précèdent la règle germinale. Sur 3 736 séquences qui lient
à 0,72 ou mieux, **20 sont conformes**. Les 3 716 autres échouent toutes sur le
plancher J, à une identité médiane de **64,3 pour un plancher à 86**.

**Et la longueur du CDR3 est le vrai couperet.** L'époque l'impose, et un
donneur de la mauvaise longueur ne sert à rien. Or les bons lieurs de l'archive
sont presque tous à CDR3 = 13 (6 209 séquences, 3 721 au-dessus de 0,72, jusqu'à
0,827), et deux longueurs seulement en avaient un au niveau du champ quand le
rejeu a tourné.

### Et une correction à mon verdict du §46

Le vivier de l'archive **bouge vite**, parce que ce sont les séquences gagnantes
des leaders qui l'alimentent :

| longueurs avec un donneur à 0,72 | |
|---|---|
| 13 novembre (le vivier du rejeu) | **2 sur 19** — 11 % des époques |
| 14 novembre, après rafraîchissement | **8 sur 19** — 42 % des époques |

460 séquences de plus en une nuit. Et pour trois des huit époques rejouées, le
plafond a bondi : CDR3 20 de 0,581 à 0,766 ; CDR3 11 de 0,692 à 0,785 ;
CDR3 18 de 0,649 à 0,766.

**Le 0 sur 8 reste le résultat du pipeline tel qu'il a tourné**, mais il a tourné
sur un vivier anormalement périmé pour les longueurs concernées. Avant de
conclure que le générateur est sans espoir, il faut le rejouer sur le vivier
frais. C'est une session Kaggle.

L'archive accusait par ailleurs un retard d'environ un jour sur l'API de partage
des scores : 196 des 230 soumissions notées en manquaient, à 0,694 de médiane et
0,779 de maximum, **exactement aux longueurs où nous n'avions rien**. Un jour,
c'est vingt époques des meilleurs donneurs qui existent.
`build_nanobody_dataset.scored_but_unarchived` fusionne désormais les deux
sources, ce qui supprime ce retard même quand l'archive rattrape ensuite.

### Piste 2 : humaniser les bons lieurs au lieu de greffer des CDR3

Si les bons lieurs sont exclus par la règle germinale, la question n'est pas
« quel CDR3 greffer » mais « comment rendre ce lieur conforme ». La moitié hors
ligne est **démontrée** (`scripts/humanise.py`) :

1. greffer les trois CDR sur la germline humaine que l'époque exige
   (`abnumber.graft_cdrs_onto_human_germline`, l'humanisation standard) →
   identité V de 91 à 95 pour un plancher à 82 ;
2. imposer au CDR3 la queue germinale `FDY` — 799 des CDR3 conformes de
   l'archive finissent ainsi → identité J de 93 à 100 pour un plancher à 86.
   Cette étape est nécessaire : la région J couvre la **queue du CDR3** autant
   que FR4, donc une charpente humanisée seule plafonne à 85,7, à un résidu du
   plancher ;
3. une ou deux substitutions de plus dans le CDR3, car le garder intact laisse
   l'identité CDR3 à 0,846 pour un seuil de rejet à 0,85.

Les apex des CDR, qui sont ce qui lie, ne sont jamais touchés ; la charpente
l'est entièrement. Vérifié : 8 donneurs sur 8 deviennent conformes, le gène V
visé est bien celui appelé (IGHV3-74, IGHV3-23, IGHV3-11 confirmés ; IGHV3-66
est revenu appelé IGHV3-23, et c'est vérifié par candidat), et la nouveauté
passe.

**Ce qui n'est pas démontré : que la liaison survive.** L'humanisation déplace
toute la charpente, et la perte d'affinité est son mode d'échec classique. Le
§38 a mesuré qu'une ou deux substitutions de charpente ne coûtent rien de
mesurable — c'est une raison d'essayer, pas un résultat. Le `design_iiptm` de
l'archive appartient à la séquence d'origine ; l'humanisée doit être notée.

### Piste 3 : RÉFUTÉE, et je l'avais présentée comme ouverte — correction

**Ce qui suit est faux et reste écrit pour la trace.** La conception à trois CDR
a été testée : §14 de ce même fichier, 60 designs, `design_iiptm` **médiane
0,169 et max 0,367** contre 0,214 pour la version à un CDR — donc *pire*, et
**0 %** au-dessus de notre propre greffage sélectionné (0,616). Le §14 réfute du
même coup l'explication causale du §11 : avec trois CDR la pose s'améliore
vraiment (5/60 passent `filter_rmsd` contre 0/20) et **les 5 designs à pose
correcte scorent 0,200 à 0,356** — corriger la pose ne remonte pas la liaison
d'un point.

J'ai lu la conclusion du §11 — « deux pistes, la 1 à tester d'abord » — sans
vérifier qu'une section ultérieure les avait déjà tranchées. Le tableau de
synthèse en tête de ce fichier le disait en une ligne : « trop peu de surface de
conception — **non**, trois CDR conçus font pire qu'un ». `inputs/design3/` et
`results/boltzgen_design_3cdr.csv` étaient là.

La leçon, pour ce dépôt comme pour moi : un §qui se termine par « piste à
tester » n'est pas une piste ouverte tant qu'on n'a pas cherché le § qui la
ferme.

### Ce que le §11 disait, et pourquoi je l'ai cru ouvert

Le §11 a écarté la conception de novo sur 20 designs à 0,21 de médiane. Mais sa
propre analyse d'échec désigne la contrainte : `filter_rmsd` 0/20 contre
`filter_rmsd_design` 7/20 — la boucle se replie comme prévu, mais **le nanobody
ne se repose pas où la conception l'avait placé**. Avec 11 résidus conçus et
CDR1, CDR2 et les deux bouts du CDR3 gelés, le modèle ne peut pas construire
d'interface.

L'exemple amont `nanobody_scaffolds/7eow.yaml` — caplacizumab — conçoit **les
trois CDR** (26..34, 52..59, 98..118) avec des insertions de longueur variable.
Le §11 listait déjà cette piste comme « à tester d'abord, elle ne coûte pas plus
cher ». Elle n'a jamais été testée. Son obstacle connu est le même que celui que
la piste 2 vient de résoudre : `inverse_folding` réassigne toute la chaîne et
faisait tomber l'identité J sous le plancher.

### Ce que font les leaders, et ce que ça dit

Les 69 soumissions acceptées de QTTECH2023 sur ces huit époques : **aucun de
leurs CDR3 n'est dans l'archive** (20 sur 20 nouveaux à moins de 85 %
d'identité), et ils lient à 0,727 de médiane. Ils ne greffent pas.

Ils ont en revanche leurs propres familles de motifs, adaptées à la longueur
imposée :

| époque | L | CDR3 | iiptm |
|---|---|---|---|
| 26302 | 11 | `ARFQWKHHFDY` | 0,761 |
| 26306 | 17 | `VRFQPRYPRTDPKHFDA` | 0,772 |
| 26303 | 18 | `ARFQPYHWDWRQENNFDV` | 0,768 |
| 26301 | 20 | `ARFQPRYPWKASRQTDNYDV` | 0,727 |

Le `FQPRYP` / `FQ[WY]` revient de la longueur 11 à la longueur 20. L'archive, de
son côté, porte un motif différent et très conservé à la longueur 13 :
`S[KR][FWY]P[YWF]R…`, 1 259 variantes sur le seul préfixe `SKFPYR`. **Deux
solutions de liaison distinctes pour la même cible, chacune réadaptée à la
longueur.** Nos 32 candidats ne portent ni l'une ni l'autre.

C'est le vrai enseignement : gagner ne demande pas de découvrir un lieur, mais
de **transposer un motif connu à la longueur imposée**. Les leaders le font.
Nous tirons au hasard autour de donneurs médiocres.

### Par ordre de coût

1. **Rejouer les 8 époques sur le vivier frais.** Et une partie de la réponse
   est gratuite : le substitut sur-promet de +0,137 (§13), et nos candidats du
   §46 étaient prédits à 0,57 de moyenne, ce qui donne 0,43 — soit exactement
   les 0,416 mesurés. Donc **il suffit de regarder si le vivier frais remonte
   la prédiction** pour savoir si une session GPU vaut la peine.
2. **Noter 32 humanisés** (`scripts.humanise`, une session). Tranche la piste 2,
   et le vivier passe de 20 donneurs conformes à des milliers si elle tient.
3. ~~La conception à trois CDR~~ — réfutée par le §14, voir ci-dessus.

### Et une vérification qui disculpe la règle en espérance

La porte PSH pourrait nous éloigner des meilleurs candidats prédits, puisqu'elle
escompte la qualité par le risque de refus. Mesuré sur les 8 époques du §46 :
les quatre candidats retenus occupent les rangs de qualité **1 à 9** dans la
plupart des époques, et le meilleur retenu est à **0,00 à 0,03** du meilleur
disponible dans le vivier mesuré. La règle ne coûte donc rien de mesurable en
qualité. Le plafond est bien celui du générateur, comme le §13 l'avait conclu.

## 48. Piste 1 répondue sans GPU : le vivier frais remonte la prédiction de 0,115

Les 8 époques rejouées sur le vivier rafraîchi (9 960 séquences contre 9 926
tronquées, et surtout les gagnantes des leaders archivées depuis). La mesure
utile est gratuite : le substitut sur-promet de **+0,137** (§13), donc sa
meilleure prédiction par époque suffit à savoir si une session GPU vaut la
peine.

| époque | L | vivier périmé (best / médiane) | vivier frais | gain |
|---|---|---|---|---|
| 26299 | 14 | 0,611 / 0,377 | **0,723** / 0,636 | +0,112 |
| 26300 | 8 | 0,415 / 0,367 | **0,600** / 0,566 | +0,185 |
| 26301 | 20 | 0,684 / 0,406 | 0,604 / 0,407 | −0,080 |
| 26302 | 11 | 0,589 / 0,515 | **0,677** / 0,557 | +0,088 |
| 26303 | 18 | 0,624 / 0,541 | **0,743** / 0,597 | +0,119 |
| 26304 | 14 | 0,590 / 0,504 | **0,732** / 0,640 | +0,142 |
| 26305 | 10 | 0,606 / 0,479 | **0,747** / 0,617 | +0,141 |
| 26306 | 17 | 0,623 / 0,463 | 0,644 / 0,446 | +0,021 |

**7 époques sur 8 en hausse, gain médian +0,115.** La meilleure prédiction du
lot passe de 0,684 à **0,747**, soit **0,610 attendu** après le biais — contre
0,718 pour la médiane du champ et 0,615 pour son quatrième rival.

### Ce que ça dit, et ce que ça ne dit pas

Le §46 a donc bien condamné le générateur sur un vivier périmé : à vivier égal,
la prédiction monte de 0,115 et l'écart au champ se réduit d'autant. Le
rafraîchissement du vivier est **gratuit** et c'est le levier le mieux rentable
qu'on ait trouvé — il suffit de lancer `fetch_submissions`,
`analyze_nanobodies` et `build_nanobody_dataset` avant chaque session.

Mais 0,610 attendu reste **sous** le 0,615 du quatrième rival, et très loin du
0,718 médian. **Le vivier frais ne suffit pas à gagner une époque.** Il fait
passer de « dernier de huit sur huit » à « au contact du dernier rival », ce qui
n'est pas la même chose que gagner.

Et une hausse de prédiction n'est pas une hausse de liaison : le biais de
+0,137 a été mesuré une fois, sur 24 candidats d'une seule tranche. C'est
l'hypothèse qu'il faut maintenant vérifier, et c'est à ça que sert le lot.

### Le lot à trois bras

`scripts.build_arms_batch` assemble en **un seul run** — puisque deux tirages de
BoltzGen diffèrent de 14,5 places et ne sont pas comparables :

| bras | n | ce qu'il teste |
|---|---|---|
| `grafted` | 32 | le générateur du mineur, vivier frais : la piste 1 |
| `humanised` | 31 | les bons lieurs d'archive humanisés : la piste 2 |
| `real` | 32 | les 4 meilleures soumissions de chaque époque, rang officiel connu |

95 séquences, ~2,2 h sur 2×T4, `RUN="arms_26299_26306"`. Les lignes sont rangées
en profondeur puis époque puis bras, donc un run interrompu couvre chaque bras
de chaque époque.

Les donneurs des humanisés vont de 0,537 à 0,785 de `design_iiptm` d'origine —
un étalement voulu : il permet de voir si l'humanisation **conserve** la liaison
du donneur, ce qui est la seule question de la piste 2. Si elle la conserve, le
vivier utilisable passe de 20 séquences conformes à des milliers, avec un
plafond à 0,827 au lieu de 0,65.

## 49. Le run à trois bras : la piste 1 confirmée, la piste 2 morte, et la mesure recalée sur ce qu'elle peut porter

Lot `arms_26299_26306`, 95 séquences notées dans le même run BoltzGen sur
2×T4 : 32 `grafted` (le générateur du mineur sur le vivier frais de §48),
31 `humanised` (la piste 2 de §47), 32 `real` (les 4 meilleures soumissions
réelles de chacune des 8 époques, rang officiel connu).
Preuves : `results/scores_nanobody_arms_26299_26306.csv`,
`results/boltzgen_arms_26299_26306.log`.

### Ce que le run dit

| bras | n | iiptm min / médiane / max | pae médian |
|---|---|---|---|
| `grafted` | 32 | 0,134 / **0,499** / 0,696 | 13,76 |
| `humanised` | 31 | 0,129 / **0,185** / 0,694 | 22,32 |
| `real` | 32 | 0,635 / **0,710** / 0,778 | 9,22 |

**Piste 1 (vivier frais) : confirmée.** Le bras `grafted` passe de 0,416 à
**0,499** de médiane et de 0,621 à **0,696** de maximum par rapport au même
générateur sur le vivier périmé du 14 novembre (§46). La prédiction de §48
était +0,115 sur le meilleur candidat ; le mesuré est +0,075. Le meilleur
candidat bat désormais 10 des 32 soumissions réelles sur l'iiptm, contre 1
sur 32. Rafraîchir le vivier est le seul levier gratuit qui ait jamais produit
un gain mesurable.

**Piste 2 (humanisation) : réfutée.** La corrélation entre l'iiptm du donneur
et celui du candidat humanisé est de **+0,027** (Spearman) : l'humanisation
ne conserve rien. La rétention médiane est de **28 %**. Une exception reste
inexpliquée — l'époque 26301 retient 88 % sur 4 candidats — et aucun
prédicteur ne la distingue : l'identité au donneur donne ρ = +0,007 (p = 0,97),
et à l'intérieur de 26301 même, quatre candidats traités identiquement
retiennent 93 / 89 / 86 / **28** %. Le vivier utilisable reste donc les
20 séquences conformes, pas les milliers qu'espérait §47.

### Ce que le run ne dit pas, et pourquoi le contrôle de `KAGGLE.md` était mal posé

L'accord avec le classement officiel des validateurs tombe à **27 paires
sur 47, soit 57 %**, après 70 % (§46) et 85-91 % sur les runs plus anciens.
Trois dégradations de suite : avant d'en tirer quoi que ce soit, il fallait
savoir si c'est notre chaîne qui dérive ou la question qui est insoluble.

Les 32 soumissions réelles ont été notées **deux fois**, dans le run du
14 novembre et dans celui-ci. Nos deux mesures des mêmes séquences ne
s'ordonnent entre elles que dans **60 % des paires** — à peine mieux que les
57 % d'accord avec les validateurs. Le bruit d'une mesure d'iiptm est
σ = **0,028** (écart absolu médian entre les deux runs : 0,024), alors que
l'écart médian entre deux rivaux voisins d'une même époque est de **0,014**.
Le rapport signal/bruit intra-époque vaut **0,51**, là où il faudrait > 2.

**Donc le classement intra-époque contre les 4 meilleures soumissions est
sous la résolution de l'instrument, et l'était déjà.** Les 85-91 % des runs
anciens ne mesuraient pas une chaîne plus saine : ils portaient sur des pools
contenant des rivaux faibles, donc largement étalés. Le contrôle tel que
`KAGGLE.md` le formulait — « si notre run les ordonne comme eux, la mesure
vaut quelque chose » — est invalide : il fait dépendre tout le protocole d'une
quantité que le bruit de BoltzGen interdit de mesurer sur 4 rivaux serrés.

Il faut donc retirer deux chiffres de ce dépôt, le mien compris :
**« 0 victoire sur 8 » et « place 3 sur 5 à l'époque 26304, premier rival
jamais battu » ne sont ni l'un ni l'autre mesurables.** Les deux reposent sur
l'ordre de candidats séparés par moins d'un demi-sigma.

### Ce qui reste mesurable, et ce que ça dit

Les médianes de bras, elles, moyennent le bruit sur 31-32 séquences
(erreur-type 0,005). L'écart `grafted` → `real` vaut **0,211 ± 0,007**, soit
**30 σ** : il est hors de tout doute. C'est la seule lecture que ce run
supporte, et elle est sévère — le champ nous devance de 0,21 d'iiptm, et le
vivier frais a comblé 0,08 de ce retard.

Conséquences pour le protocole :

- **Juger par médianes de bras, jamais par victoires par époque.** Un
  comptage de victoires sur 8 époques contre 4 rivaux serrés n'a aucune
  puissance ; 32 séquences par bras en ont.
- **Ne pas réclamer un meilleur accord avec les validateurs**, qui plafonne à
  ~60 % par construction. Le contrôle utile est la reproductibilité de notre
  propre mesure d'un run à l'autre, qui est désormais chiffrée (σ = 0,028).
- **Un classement intra-époque demanderait des répliques** : à σ = 0,028 et
  un écart réel de 0,014, séparer deux rivaux voisins demande ~16 répliques
  par séquence. Ce n'est pas le bon usage du quota.
- La porte PSH et les portes déterministes restent, elles, hors de ce débat :
  elles sont déterministes et toujours le préalable à toute soumission.

## 50. Pourquoi la liaison ne voyage pas avec le motif — et ce que ça règle pour IgLM

La piste 2 (§47, §49) a réfuté l'humanisation sans dire *où* la liaison se
perd. Les 32 humanisés notés ont été réappariés à leur donneur (`scripts.humanise`
étant déterministe à `--cutoff` fixé, 30 des 32 se réalignent sur un CDR3 de
même longueur), ce qui permet de compter les mutations et de corréler.

| variable | Spearman avec la rétention | p |
|---|---|---|
| mutations dans le CDR3 | **+0,58** | 0,001 |
| mutations dans la queue du CDR3 (3 résidus) | +0,58 | 0,001 |
| mutations dans la charpente (médiane 12) | **−0,13** | 0,49 |

| CDR3 conservé | n | rétention médiane |
|---|---|---|
| à l'identique (0 mutation) | 13 | **20 %** |
| 1-2 mutations | 15 | 31 % |
| 3 mutations | 2 | 42 % |

**Conserver le CDR3 exactement donne la pire rétention des trois groupes.**
Mon hypothèse était que la queue `FDY` forcée par le plancher d'identité J
coûtait la liaison ; elle est fausse, et dans le mauvais sens. Le signe positif
est à lire prudemment — n = 30 réparti sur 8 époques, et il suit la variation
entre époques plutôt qu'entre candidats — mais le **null sur la charpente**
(ρ = −0,13, p = 0,49, pour une médiane de 12 mutations) et les 20 % du groupe
« CDR3 intact » sont solides.

La lecture : **la liaison de ces séquences d'archive n'est pas portée par le
CDR3**. Elle dépend du contexte de charpente qu'on remplace, et 12 mutations de
charpente l'effacent quel que soit le soin apporté aux boucles.

### Ce que ça règle pour IgLM

IgLM est un modèle de langage à infilling entraîné sur ~558 M de séquences
d'anticorps, avec des jetons `[HEAVY]`/`[CAMEL]` et une longueur d'infill
variable. La seule variante que je m'étais engagé à défendre était
**l'infilling sous contrainte de motif** : garder le CDR3 d'un bon lieur de
l'archive et laisser IgLM régénérer la charpente conforme autour.

C'est **exactement l'opération que le bras `humanised` vient de mesurer**, et
dans sa version la plus favorable : la charpente n'était pas échantillonnée
par un modèle de langage mais prise sur un vrai germline humain, et les CDR
étaient greffés par `abnumber`. Résultat : 28 % de rétention médiane, 20 %
quand le CDR3 est conservé au résidu près, et aucune dépendance au nombre de
mutations de charpente. IgLM ferait la même chose avec deux handicaps en plus :

1. il est **agnostique à la cible** — rien dans son entraînement ne connaît
   Q9NZQ7, donc il ne peut pas retrouver le contexte de charpente perdu, il
   peut seulement en proposer un vraisemblable ;
2. la vraisemblance d'un modèle de langage d'anticorps est **décorrélée de la
   liaison** ; elle sélectionne des séquences d'allure naturelle, ce que les
   germlines humains de la piste 2 étaient déjà par construction.

**Donc non, IgLM n'est pas une solution**, et ce n'est plus un argument mais
une mesure : la piste a été testée par proxy sur 31 séquences GPU avec un
substitut strictement meilleur que lui. Ce qu'IgLM pourrait encore servir —
et seulement ça — est la **diversification de charpente à liaison déjà
acquise** : si un jour un candidat lie à 0,72 mais butte sur une porte
déterministe, IgLM propose des variantes de charpente conformes. C'est un
outil de conformité, pas un générateur de liaison.

### Ce qui reste comme pistes

Le retard est de 0,211 d'iiptm (§49) et aucune manipulation de séquence à
partir de l'archive ne l'a comblé. Les deux seules pistes que ce dépôt n'a pas
réfutées par mesure :

- **le vivier frais**, qui a rendu +0,075 et qui est gratuit, à refaire à
  chaque époque plutôt qu'une fois ;
- **le design conditionné par la structure** sur le vrai complexe 8AOK, qui
  est ce que fait le champ et que nous n'avons jamais eu les moyens de faire
  (§14 a réfuté la version naïve à 3 CDR de BoltzGen, pas la démarche).

## 51. Le vivier conforme n'est plus de 20 séquences mais de 256, et le mineur le rafraîchit seul

Le §48 et le §49 répètent un chiffre : « 3736 séquences de l'archive lient à
0,72 ou mieux mais seules 20 sont conformes », avec un plafond utilisable de
0,57-0,77 contre un champ à 0,72. **Les deux sont périmés**, et ce n'était pas
une erreur de raisonnement mais une mesure prise sur un vivier vieux de six
semaines.

Après un rafraîchissement complet (déchiffrement de 26345-26348, jointure des
métriques publiées, re-téléchargement de l'archive et re-caractérisation), le
vivier porte 10536 lignes au lieu de 9960, et :

| | avant (§48) | maintenant |
|---|---|---|
| séquences liant ≥ 0,72 | 3736 | 3972 |
| **dont conformes** (V ≥ 82, J ≥ 86) | **20** | **256** |
| longueurs de CDR3 (6-24) avec un donneur ≥ 0,72 | 2 puis 8 / 19 | **15 / 19** |
| meilleur donneur conforme | 0,65 | **0,788** |

Attribution, parce que la cause change ce qu'on en fait : 197 des 256 viennent
de **l'archive** re-téléchargée et 59 des soumissions du champ que nous avons
déchiffrées (époques 26306-26347). Les soumissions du champ ajoutent trois
longueurs de couverture (12/19 → 15/19) et le meilleur donneur de tous (0,788
contre 0,785). Ce n'est donc pas un seul des deux leviers : il faut les deux.

Ce que le générateur verrait, époque par époque, sur les huit que nous rejouons :

| époque | CDR3 | donneurs conformes ≥ 0,72 | meilleur |
|---|---|---|---|
| 26299 | 14 | 15 | 0,788 |
| 26300 | 8 | 1 | 0,721 |
| 26301 | 20 | 23 | 0,767 |
| 26302 | 11 | 41 | 0,785 |
| 26303 | 18 | 52 | 0,781 |
| 26304 | 14 | 15 | 0,788 |
| 26305 | 10 | 0 | — |
| 26306 | 17 | 20 | 0,772 |

Sept époques sur huit ont désormais un donneur conforme au niveau du champ, là
où le §48 n'en trouvait pour presque aucune. Le plafond du vivier (0,788) est
au-dessus du maximum du champ dans le run à trois bras (0,778).

**Attention à ce que cela ne dit pas.** Greffer depuis un donneur à 0,788 ne
produit pas un candidat à 0,788 : le bras `grafted` du §49 partait déjà de ce
vivier-là, en moins frais, et sa médiane était de 0,499. La perte à la greffe
reste le problème central, et le §50 a montré que la liaison ne voyage pas avec
le motif. Ce que le vivier frais change est le point de départ, pas le
rendement de l'opération.

### Le mineur ne refuse plus l'époque, il rafraîchit

Jusqu'ici le mineur détectait le vivier périmé et refusait l'époque en
imprimant trois commandes pour un humain (`status: donor_pool_stale`). Il
dépensait donc l'époque pour signaler un gain gratuit au lieu de le prendre.
Coûts mesurés sur cette machine :

| étape | coût |
|---|---|
| `scripts.fetch_submissions` | **22 s par époque** (86 s pour 4) |
| `scripts.analyze_nanobodies` | 7,5 s |
| `scripts.build_nanobody_dataset` | **4 min 14** (10536 lignes) |

Soit **7,6 min** pour le rattrapage par défaut de 8 époques, sur les 72 d'une
époque — contre 12 min pour 20 époques, et un quart d'heure si on laisse filer
une semaine. C'est l'argument décisif pour le faire **à chaque époque** plutôt
qu'une fois : tenir le rythme coûte une minute, rattraper coûte un quart
d'heure, et le vivier périmé coûte l'époque entière.

Deux points de conception valent d'être retenus :

- les trous sont comblés **du plus récent au plus ancien** dans une fenêtre de
  60 époques, pas seulement la queue. Ne déchiffrer que les dernières époques
  ferait retomber `donor_pool_lag` à zéro en laissant un trou permanent au
  milieu, et le mineur croirait son vivier frais alors qu'il lui manquerait la
  semaine ;
- les trois étapes tournent dans cet ordre (sans les deux premières les
  nouvelles séquences ne portent pas d'époque et le retard ne bouge pas),
  chacune dans son processus, et un rafraîchissement qui ne tient pas dans le
  budget restant est refusé **avant** d'être commencé plutôt qu'interrompu au
  milieu.

Le refus à `--max-donor-lag` reste, en dernier recours : si le
rafraîchissement échoue, mieux vaut ne pas soumettre que grafter depuis une
génération que la compétition a dépassée.

### Et PR 102 est en vigueur, sans conséquence pour nous

Nous avons passé l'époque 26331. `scripts.verify_challenge` compare désormais
notre dérivation au terrain sur **12 époques et 348 soumissions acceptées**,
dont quatre postérieures au changement (26345-26348, 118 soumissions) : **12
sur 12 d'accord**. L'avertissement Discord sur
`get_challenge_params_from_blockhash` ne touche pas les paramètres nanobody,
ce que le test sur 200 000 hashes prédisait et que le terrain confirme.

## 52. La conception BoltzGen est bornée à 0,39 — [la borne est retirée, voir §55] — et le goulot n'est pas le générateur, c'est le débit de notation

Le §14 avait réfuté la conception à trois CDR en la comparant au greffage. Il
laissait deux hypothèses ouvertes, dont « conception libre à grande échelle :
cohérent avec le `--num_designs` typiquement 10 000 du CLI — mais 10 000
designs font ~10 jours de Kaggle ». Cette hypothèse se teste sans GPU, par la
queue de la distribution déjà mesurée.

### La conception plafonne, et 10 000 designs n'y changeraient rien

Ajustement d'une loi de Pareto généralisée sur les dépassements des 151 designs
de `results/boltzgen_design_3cdr.csv` (protocole `nanobody-anything`, charpente
et CDR1/CDR2 fixes, CDR3 conçu) :

| seuil | ξ | max attendu sur 1 000 | sur 10 000 | sur 100 000 | borne |
|---|---|---|---|---|---|
| 0,235 (30 dépassements) | −0,53 | 0,378 | 0,384 | 0,386 | **0,387** |
| 0,260 (15 dépassements) | −1,89 | 0,367 | 0,367 | 0,367 | **0,367** |

**ξ < 0 aux deux seuils : la queue est bornée.** Concevoir le CDR3 avec
BoltzGen sur une charpente VHH conforme ne dépasse pas ~0,39, quel que soit le
nombre de designs. L'hypothèse de l'échelle est donc morte **dans cette forme**,
et le dixième de jour-Kaggle qu'elle réclamait est économisé.

Ce qui n'est pas réfuté : une conception libre, sans charpente imposée, puis
filtrage sur la conformité. Mais le champ, lui, a bien une charpente germinale
(identité V médiane 83,7), donc ce n'est pas ce qu'il fait.

### Le greffage, lui, n'est pas borné sous le champ

La même analyse sur le bras `grafted` du §49 donne le résultat inverse :

| | ξ | max attendu sur 32 | 128 | 256 | 512 | borne |
|---|---|---|---|---|---|---|
| `grafted` | −0,25 | 0,654 | 0,695 | **0,711** | 0,724 | **0,795** |
| `real` (le champ) | −1,27 | 0,775 | 0,778 | 0,778 | 0,778 | 0,778 |

La borne du greffage, **0,795**, est au-dessus du meilleur du champ dans ce run
(0,778) et au-dessus de son meilleur gagnant sur 12 époques (0,788). Ce qui
nous plafonne n'est pas le générateur : c'est que nous n'en notons que 32.

Bootstrap sur 400 rejeux, parce que n = 32 avec 10 dépassements est mince :
max attendu sur 256 candidats **0,702 [0,610 – 0,764]**. L'intervalle est
large et il faut le dire, mais son centre est au niveau du champ.

### Ce que cela coûte, et pourquoi c'est la vraie question

À 1,4 min par séquence sur 2×T4 :

| candidats notés | durée | tient dans une époque de 72 min ? |
|---|---|---|
| 47 | 1,1 h | tout juste |
| 128 | 3,0 h | non |
| 256 | 6,0 h | non |
| 512 | 11,9 h | non |

Le gain de 47 à 256 candidats vaut **+0,06 d'iiptm** sur l'espérance du
maximum. La question « faut-il louer un GPU » cesse donc d'être une question de
confort : elle a un prix et un rendement. Il faut environ **5× le débit de deux
T4** pour noter 256 candidats dans une époque.

### La réserve qui décide du prochain run

Le bras `grafted` du §49 a été produit sur le vivier **fuité** (§REPLAY_AUDIT).
Sa distribution est donc optimiste d'une quantité inconnue, et toute cette
section en hérite. Le run qui tranche n'est pas un run de conception : c'est un
lot **honnête et profond** — peu d'époques, beaucoup de candidats — qui mesure
la queue réelle du greffage. C'est le seul chiffre qui dise si augmenter le
débit de notation suffit, ou s'il faut changer de méthode malgré tout.

## 53. Le premier run honnête et profond : la queue tient, le témoin se retourne, et une époque s'effondre

Lot `deep_26304_26345_26347` : 170 séquences notées dans un seul run sur 2×T4,
158 candidats greffés sur trois époques et 12 soumissions réelles comme
ancrages. **Premier lot dont les deux défauts du rejeu sont corrigés** — la
coupure d'archive ne laisse plus passer les lignes non datées, et `observe`
donne au sélecteur le même vivier et les mêmes portes que le mineur.
Preuves : `results/scores_nanobody_deep_26304_26345_26347.csv`,
`results/boltzgen_deep_26304_26345_26347.log`.

### La mesure est fidèle

Les 12 ancrages portent leurs métriques publiées. Écart absolu médian sur
l'iiptm : **0,012**, maximum 0,042 — cohérent avec le σ = 0,028 d'une mesure
établi au §49. L'instrument lit juste ; ce qui suit peut être lu.

### La queue du greffage tient, et le §52 survit à la mesure honnête

| | ξ | borne | max attendu sur 256 |
|---|---|---|---|
| §52, vivier **fuité**, n = 32 | −0,25 | 0,795 | 0,711 |
| **ce run, vivier honnête, n = 158** | **−0,31** | **0,762** | **0,704** |

La fuite n'avait donc **pas** gonflé la queue de façon significative, et c'est
la première prédiction de cette série que la mesure confirme. L'espérance du
maximum passe de 0,651 sur 32 candidats à 0,704 sur 256 : noter plus de
candidats paie, et le prix reste celui du §52 (~5× le débit de deux T4 pour
tenir 256 candidats dans une époque).

Mais la borne honnête, **0,762**, est **sous** le meilleur gagnant du champ sur
douze époques (0,788) et sous deux des trois gagnants de ce run. Noter plus
nous amène au contact, pas devant.

### Le témoin se retourne : la disponibilité des donneurs ne commande pas

L'époque 26347 avait été mise dans le lot comme témoin négatif — **aucun**
donneur conforme à 0,72, plafond 0,675. La prédiction était que sa queue
s'effondrerait. C'est l'inverse :

| époque | donneurs ≥ 0,72 | plafond du vivier | médiane obtenue | max | ξ | borne |
|---|---|---|---|---|---|---|
| 26304 | 8 | 0,774 | **0,161** | 0,507 | +0,83 | — |
| 26345 | 16 | 0,779 | 0,490 | 0,726 | −0,24 | 0,802 |
| 26347 | **0** | 0,675 | **0,521** | 0,653 | −0,64 | 0,659 |

L'époque sans aucun bon donneur donne la **meilleure médiane des trois**, et
celle qui avait huit donneurs à 0,774 s'effondre. L'hypothèse que je voulais
tester est donc réfutée — mais rien de propre ne la remplace.

(Le ξ positif de 26304 est à jeter : sur 45 points dont 87 % sous 0,25,
l'ajustement extrapole un maximum de 1,56 pour une grandeur bornée par 1.)

### L'effondrement de 26304, et ce qu'il dit de notre variance

26304 et 26345 **imposent les mêmes paramètres** — IGHV3-74, CDR3 de 14 — et
leurs viviers honnêtes sont presque les mêmes : 201 contre 233 donneurs
conformes de longueur 14, 8 contre 16 à 0,72 ou mieux, plafonds 0,774 et 0,779.
Ce qui diffère est la graine du générateur et un crible PSH à une réplique.

| époque | p10 | p25 | p50 | p75 | p90 | max | sous 0,25 | au-dessus de 0,60 |
|---|---|---|---|---|---|---|---|---|
| 26304 | 0,131 | 0,141 | 0,161 | 0,190 | 0,254 | 0,507 | **87 %** | **0 %** |
| 26345 | 0,148 | 0,192 | 0,490 | 0,592 | 0,635 | 0,726 | 35 % | 20 % |
| 26347 | 0,181 | 0,368 | 0,521 | 0,572 | 0,601 | 0,653 | 17 % | 10 % |

**À entrées quasi identiques, la médiane varie d'un facteur trois et la part de
quasi non-lieurs de 17 % à 87 %.** C'est le résultat le plus dur de ce run : la
qualité d'un lot de candidats n'est pas une propriété du vivier ni de l'époque,
elle dépend d'un tirage que nous ne contrôlons pas et ne savons pas prédire. Un
mineur qui ne note que quelques candidats est donc exposé à tomber sur une
époque comme 26304, où **aucun** de ses 45 candidats n'atteint 0,60.

Cela renforce la conclusion du §52 par un autre chemin : la profondeur de
notation n'est pas seulement un gain d'espérance, c'est la seule protection
contre cette variance.

### Où nous en sommes contre le champ

| époque | notre meilleur | meilleur rival | écart |
|---|---|---|---|
| 26304 | 0,507 | 0,777 | 0,270 |
| 26345 | **0,726** | 0,788 | **0,063** |
| 26347 | 0,653 | 0,746 | 0,093 |

Sur l'époque 26345, notre meilleur candidat honnête est à **0,063** du meilleur
rival — l'écart le plus court jamais mesuré, contre 0,211 de médiane à médiane
au §49. Le classement du validateur porte sur dix métriques, donc cela ne dit
pas que nous aurions gagné ; cela dit que sur une époque favorable et avec
54 candidats notés, nous sommes entrés dans le bruit du champ.

## 54. Notre porte de nouveauté applique une règle qui n'existe pas, et elle jette la moitié du vivier

Le guide de stratégie amont (`docs/nanobody_stragegy_guide.md` dans
metanova-labs/nova) n'avait **jamais été lu par ce dépôt** — aucun document ne
le citait. Il nomme la vraie règle de diversité, et le code du filtre la
confirme.

### Ce que le validateur applique vraiment

`metanano/validators/diversity_validator.py` n'enforce que deux choses :

1. un regroupement MMseqs2 **intra-lot** (`global_cluster_identity`), que le
   guide déclare non activé pour l'instant ;
2. contre l'historique : **similarité de Jaccard sur k-mers de longueur 5**,
   seuil **0,9**, et la stratégie par défaut est `plan_b` — la comparaison
   porte sur les **50 meilleures** séquences seulement, pas sur l'archive.

Il n'y a **ni** règle d'identité de séquence entière à 95 %, **ni** règle de
CDR3 à 85 %. Notre `is_duplicate` les a tirées du `search_service` du même
dépôt — un moteur de recherche par similarité (Jaccard grossier puis
alignement fin) offert aux mineurs comme outil de consultation, que le
validateur de diversité n'appelle jamais.

### Ce que cela coûte, mesuré

Jaccard k=5 contre les 50 meilleures (approchées par notre `rank_sum` sur les
17 679 séquences notées) :

| | Jaccard max vs top 50 | refusés par la vraie règle |
|---|---|---|
| nos 158 candidats notés au §53 | médiane 0,172 · p90 0,239 · **max 0,286** | **0 sur 158** |
| 5064 lieurs d'archive à ≥ 0,72 (250 échantillonnés) | médiane 0,314 | 3 % |
| ceux d'entre eux conformes germline (71) | — | **0 %** |

Nos candidats sont à **0,29 au pire** pour un seuil à 0,90 : la contrainte
réelle n'est pas près de mordre. Pendant ce temps notre porte refusait, sur la
seule époque 26345, **1497 candidats sur 3000** au motif « CDR3 ≥ 85 % identique
à un CDR archivé », plus 6 des 60 retenus à ≥ 95 % d'identité entière.

### Les deux conséquences, par ordre d'importance

**1. `MAX_CDR3_IDENTITY = 0.85` est à supprimer.** Le générateur mute
délibérément le CDR3 greffé pour descendre sous 85 % d'identité avec tout CDR3
archivé. C'est-à-dire qu'il **dégrade exprès la seule chose qu'il emprunte à un
bon lieur**, pour satisfaire une règle inexistante. Greffer le CDR3 d'un donneur
**verbatim** est permis : l'unicité globale du Layer 1 porte sur la séquence
entière, et la charpente diffère.

**2. Le vivier double.** Ne plus jeter la moitié des candidats multiplie par
deux le nombre de candidats disponibles à coût de génération constant — et le
§52 comme le §53 ont montré que la profondeur de notation est le seul levier
mesuré (espérance du maximum de 0,651 sur 32 candidats à 0,704 sur 256).

### Ce qui reste à vérifier avant de changer la porte

La porte actuelle est notre filet avant soumission : la relâcher sur la foi
d'une lecture de code est exactement le genre de geste qui a coûté onze
soumissions à ce dépôt. Deux points à établir d'abord :

- le top 50 doit être calculé par `utils.minmax_weighted_rank`, l'algorithme que
  le guide nomme, et non approché par notre `rank_sum` ;
- le guide dit que la contrainte intra-lot (« au moins 2 mutations dans les CDR,
  dont 1 dans le CDR3 ») n'est pas active mais le sera. Un greffage verbatim
  doit donc rester compatible avec elle, ce qui est le cas si les séquences d'un
  même lot diffèrent entre elles.

## 55. Nos 210 designs n'avaient aucune cystéine : la réfutation du §14 et la borne du §52 sont des artefacts

Le guide amont avertit que « l'implémentation de BoltzGen incluse dans le dépôt
NOVA est adaptée pour servir de moteur de **notation** et non de pipeline
génératif. Vous utiliserez l'implémentation originale pour créer vos designs. »
Nos runs de conception (§11, §14) ont utilisé la version embarquée. Vérification
faite, l'avertissement est fondé — mais pas de la façon que je craignais.

### Ce qui est identique, et ce qui manque

| | embarquée (NOVA) | originale |
|---|---|---|
| version | 0.2.0 | 0.3.2 |
| poids de conception | `boltzgen1_diverse.ckpt`, `boltzgen1_adherence.ckpt`, `boltzgen1_ifold.ckpt` | **les mêmes** |
| `resources/config/design.yaml` | — | **identique à l'octet** |
| `model/modules/diffusion.py` | 2 lignes de différence, aucune touchant pas, bruit, sigma ou échantillonnage |  |
| `model/modules/inverse_fold.py` | **173 lignes absentes**, dont `build_constraint_logit_mask` et `_build_symmetric_groups` |  |

Le générateur n'est donc pas bridé : mêmes poids, même diffusion, même config.
Ce qui manque est la machinerie qui **contraint les résidus** que l'inverse
folding peut poser — précisément ce sur quoi notre spécification repose, puisque
elle fixe la charpente, CDR1, CDR2, et la tête et la queue du CDR3.

### La preuve, sur les designs que nous avons déjà

Le protocole `nanobody-anything` pose `filter_cysteine=true` et le CLI construit
une liste `inverse_fold_avoid`. Sans le masque de contrainte qui protège les
positions gelées, l'évitement de la cystéine s'est appliqué **partout**, y
compris aux positions que la spécification avait figées :

| | cystéines par séquence |
|---|---|
| conception 1 CDR (§11), n = 60 | min 0 · médiane 0 · max **0** |
| conception 3 CDR (§14), n = 150 | min 0 · médiane 0 · max **0** |
| les 348 soumissions acceptées du champ | min **2** · médiane 2 |

**210 designs sur 210 à zéro cystéine.** La charpente lit `...SLRLSEAASG` là où
le motif canonique est `...SLRLSCAASG` : la **Cys23**, que le guide décrit comme
formant avec Cys104 le pont disulfure « critique pour la stabilité du
repliement », a été remplacée par un glutamate dans tous les cas. FR1 dérive
aussi (quatre variantes des 26 premiers résidus là où la spécification en
imposait une).

Un nanobody sans son pont disulfure canonique ne se replie pas comme un
nanobody. Cela explique tout ce que le §11 et le §14 ont observé et n'ont pas
su expliquer : un `design_iiptm` plafonné à 0,21-0,37, des poses à 9,23 Å pour
un seuil de 2,5, et le fait que **corriger la pose ne remontait pas la liaison
d'un point**. Ce n'était pas une limite de la méthode, c'était une protéine
cassée.

Et ces designs n'étaient de toute façon **pas soumissibles** : la validité
Layer 1 exige `min_cysteines`, typiquement ≥ 1.

### Ce qui tombe

- **§14, « Trois CDR conçus au lieu d'un : réfutation nette »** : la réfutation
  ne porte pas. Elle compare le greffage à des séquences dépourvues du pont
  disulfure canonique.
- **§11** : même défaut, même population.
- **§52, la borne de 0,39** : l'ajustement de Pareto généralisé a été fait sur
  ces 151 designs. Il mesure la queue de distribution de protéines cassées, pas
  celle de la conception. **La borne est retirée.** Ce que le §52 disait de la
  conception à grande échelle — « 10 000 designs n'y changeraient rien » — n'est
  plus établi.
- Par conséquent, **la conception conditionnée par la structure n'est pas
  réfutée. Elle n'a jamais été testée.**

### Ce qui tient dans le §52

La moitié qui porte sur le **greffage** est intacte : elle a été ajustée sur le
bras `grafted` du §49, puis confirmée honnêtement au §53 (ξ = −0,31, borne
0,762, espérance du maximum de 0,651 sur 32 candidats à 0,704 sur 256). La
profondeur de notation reste le seul levier mesuré.

### Ce qu'il faut faire

Refaire la conception avec **BoltzGen 0.3.2 original**, pas la version
embarquée, et vérifier avant toute conclusion que les designs rendus portent
leurs deux cystéines canoniques et la charpente que la spécification impose.
C'est un contrôle de trois lignes qui aurait dû exister dès le §11 : compter les
cystéines et comparer les régions gelées à ce qu'on avait gelé.

## 56. La BoltzGen originale respecte les positions gelées : la conception est enfin mesurable

Run `design_cdr3_26345`, 120 designs, BoltzGen **0.3.2 originale** installée
depuis PyPI, 6 charpentes conformes IGHV3-74, CDR1 et CDR2 fixes, CDR3 conçu sur
9 résidus à la longueur 14. Preuves : `results/designs_design_cdr3_26345.csv`,
`results/design_design_cdr3_26345.log`.

### Le garde-fou passe, et il confirme le diagnostic du §55

| | cystéines | à zéro | sans la paire canonique | variantes de FR1 |
|---|---|---|---|---|
| §11 / §14, BoltzGen **embarquée** | médiane 0 | **210 / 210** | 210 / 210 | 7 (pour 4 charpentes) |
| ce run, BoltzGen **originale** | **médiane 2** | **0** | **0** | 6 (pour 6 charpentes) |

Le §55 avait attribué la destruction du pont disulfure à l'absence de
`build_constraint_logit_mask` dans la version embarquée. L'expérience le
confirme : à spécification identique, la version originale rend des nanobodies
avec leurs deux cystéines et une variante de FR1 par charpente, exactement ce
que la spécification impose.

**Et les 120 designs passent nos portes déterministes, 120 sur 120** : bon gène
V, CDR3 de longueur 14, planchers germinaux, marqueurs VHH, pas d'homopolymère,
pas de peptide signal. Ce sont des séquences soumissibles, là où les 210
précédentes auraient été refusées sur `min_cysteines`.

### Ce que les designs valent, sur l'échelle de BoltzGen 0.3.2

| | médiane | max | part > 0,50 | part > 0,60 |
|---|---|---|---|---|
| §11, 1 CDR (invalide) | 0,214 | 0,311 | — | — |
| §14, 3 CDR (invalide) | 0,171 | 0,367 | — | — |
| **ce run, CDR3 seul, valide** | **0,241** | **0,650** | **34 %** | **8 %** |

Le maximum passe de 0,367 à **0,650**. La médiane bouge peu, ce qui est
cohérent avec ce qu'on attend d'un générateur : la plupart des tirages ne lient
pas, et c'est la queue qui compte.

Queue ajustée sur les 271 lignes : ξ = **−0,59**, borne **0,658**, et le maximum
attendu ne monte qu'à 0,657 sur 10 000 designs. La borne retirée du §52 disait
0,387 ; la vraie borne de cette spécification est **0,658**. Concevoir plus ne
paiera donc pas beaucoup au-delà de quelques centaines — mais ce plafond vaut
pour *cette* spécification (6 charpentes, 9 résidus libres), pas pour la
méthode.

### La réserve qui interdit de conclure maintenant

Ces métriques sortent de **BoltzGen 0.3.2**, pas du moteur de notation des
validateurs. Les deux builds diffèrent sur 45 fichiers, et le §49 a mesuré
σ = 0,028 rien que sur le bruit d'une mesure. **Aucune comparaison au champ
n'est licite avant que les designs soient notés par la build des validateurs.**

C'est l'objet du lot `design3arms_26345` : 120 designs, les 54 candidats
greffés honnêtes du §53 et les 4 rivaux réels de l'époque, **notés dans le même
run** par le moteur des validateurs. 178 séquences, ~4,2 h sur 2×T4. Pour
mémoire, sur cette époque et sur l'échelle du validateur, notre meilleur greffé
honnête lit 0,726 et le meilleur rival 0,788.

## 57. Nous perdons la moitié « confiance » et gagnons la moitié « interaction » : place 3 sur 5, 2 rivaux battus

> **RETIRÉ le 18 novembre (§61).** Renoté dans un second run, le candidat de
> cette section lit 0,709 au lieu de 0,741 et ne bat plus aucun rival : les
> 170 candidats du lot `wide4arms` finissent tous derniers. Le résultat était
> un tirage favorable de la mesure.
>
> **Nuance apportée le 17 novembre (§58).** Le « place 3 sur 5 » de cette
> section est **un candidat sur 120** : confrontés un par un aux 4 rivaux,
> 98 % des designs et 98 % des greffés finissent derniers. La lecture
> métrique par métrique ci-dessous reste exacte pour ce candidat, mais ce
> n'est pas une propriété du bras.

Lot `design3arms_26345` : 178 séquences notées par le moteur des validateurs
dans un seul run — 120 designs du §56, les 54 greffés honnêtes du §53, et les
4 rivaux réels de l'époque. Preuves :
`results/scores_nanobody_design3arms_26345.csv`,
`results/boltzgen_design3arms_26345.log`.

Fidélité : écart absolu médian de **0,029** sur les 4 ancrages contre leurs
métriques publiées, cohérent avec σ = 0,028.

### Sur l'iiptm, le greffage reste meilleur que la conception

| bras | n | médiane | max | pae médian |
|---|---|---|---|---|
| `design` | 120 | 0,395 | **0,741** | 16,82 |
| `grafted` | 54 | **0,491** | 0,729 | 13,72 |
| `real` | 4 | 0,787 | 0,789 | 8,79 |

Les deux maxima sont séparés de 0,012, soit bien moins que σ : indiscernables.
La médiane, elle, favorise clairement le greffage. Queue sur l'échelle du
validateur, même run : `grafted` ξ = −0,29 et max attendu 0,747 sur
256 candidats, `design` ξ = −0,12 et 0,713. **Sur cette métrique, concevoir
n'apporte rien.**

### Mais le validateur ne classe pas sur l'iiptm, et là tout s'inverse

Notre meilleur **design** contre les 4 rivaux, dans un vivier de cinq séquences
— celui qui ressemble à une vraie époque, où nous soumettons quelques séquences
contre le champ :

| métrique | nous | rang | rivaux |
|---|---|---|---|
| `design_iiptm` | 0,637 | **5** | 0,755 – 0,789 |
| `design_ptm` | 0,808 | **5** | 0,828 – 0,837 |
| `design_to_target_iptm` | 0,596 | **5** | 0,668 – 0,689 |
| `min_design_to_target_pae` | 3,233 | **5** | 1,903 – 1,992 |
| `interaction_pae` | 10,643 | **5** | 8,750 – 8,902 |
| `plip_hbonds_refolded` | 10 | **1** | 5 – 9 |
| `plip_saltbridge_refolded` | 5 | **1** | 3 – 5 |
| `delta_sasa_refolded` | 873,6 | **1** | 751,5 – 799,7 |
| `liability_score` | 190 | **1** | 200 – 210 |
| `liability_num_violations` | 21 | **1** | 23 – 24 |
| **rang-somme** | | **30** | 21 – 31 |

**Dernier sur les cinq métriques de confiance, premier sur les cinq métriques
d'interaction et de développabilité.** Rang-somme 30 contre 21 à 31 :
**place 3 sur 5, 2 rivaux sur 4 battus**. C'est le meilleur résultat jamais
mesuré dans ce dépôt, et il vient du bras `design`, pas du greffage — dont le
meilleur candidat finit 5ᵉ sur 5 dans le même exercice.

### Ce que cela corrige dans ma façon de juger

J'ai passé toute cette série à juger sur `design_iiptm` seul. Le guide amont le
dit en clair : « vous êtes classé sur 10 dimensions simultanément », et « un
design 5ᵉ sur cinq métriques (total 25) bat un design 1ᵉʳ sur quatre et 100ᵉ sur
une (total 104) ». L'écart à combler n'est donc pas « 0,048 d'iiptm » mais
**cinq rangs sur le bloc confiance**, en conservant l'avance sur l'autre moitié.

Et cela réhabilite la conception : sur l'iiptm elle perd, sur le rang-somme elle
gagne. Les deux bras ne sont pas concurrents, ils sont complémentaires — le
greffage lie mieux, la conception est mieux équilibrée.

### Un piège de `scripts.review_observation`, à corriger

Le script annonce pour ce lot « place 1 sur 5 ; bat 4 rivaux sur 4 ». C'est
**vrai dans le vivier complet** (174 des nôtres plus 4 rivaux) et **faux dans
une vraie époque** : avec 174 de nos candidats dans le vivier, les rangs
par métrique des rivaux sont tirés vers le bas sur les métriques où nos
candidats les dominent en nombre, surtout les liabilities. Recalculé sur un
vivier de cinq séquences, le même candidat greffé passe de place 1 à **place 5
sur 5**. Le rang-somme n'est pas une propriété d'une séquence, c'est une
propriété d'une séquence **dans son vivier**, et le vivier réel comptait
31 soumissions acceptées, pas 4 ni 178.

## 58. Le mélange est pire que la conception seule, et le « place 3 sur 5 » est un candidat sur 120

### D'abord : il n'y a pas de mélange à soumettre

`config/config.yaml` du dépôt amont porte **`num_sequences: 1`**. Une soumission
ne contient qu'une séquence, donc l'agrégation multi-séquences que décrit le
guide (« moyennés puis sommés ») ne s'applique pas à la configuration actuelle.
Mélanger ne peut donc vouloir dire qu'une chose : réunir les deux **viviers de
candidats** et n'en soumettre que le meilleur.

### Et mélanger les viviers dégrade le résultat

Sélection à l'aveugle des rivaux — rang-somme calculé à l'intérieur de notre
propre vivier, ce qu'un mineur peut faire — puis confrontation aux 4 rivaux :

| vivier | candidat choisi | place |
|---|---|---|
| `design` seul (120) | iiptm 0,637 | **3 / 5**, bat 2 rivaux |
| `grafted` seul (54) | iiptm 0,604 | 5 / 5, bat 0 |
| **mélangé (174)** | iiptm 0,604 | **5 / 5, bat 0** |

Le vivier mélangé choisit le même candidat que le greffage seul, et finit
dernier. Trois critères de sélection ont été essayés — les 10 métriques, le
bloc confiance seul, le bloc interaction+développabilité seul — et **aucun ne
sauve le mélange** : tous choisissent un candidat greffé qui s'effondre face aux
rivaux. Le plafond d'un oracle parfait sur le vivier mélangé est la place 3,
exactement celui du bras `design` seul : mélanger n'ajoute aucun candidat
meilleur, il ajoute seulement des candidats qui trompent le sélecteur.

### Et la correction qui m'incombe sur le §57

J'ai présenté « place 3 sur 5, 2 rivaux battus » comme le résultat du bras
`design`. En confrontant **chacun** des 174 candidats aux 4 rivaux :

| bras | place 3 ou mieux | place 4 | place 5 |
|---|---|---|---|
| `design` (120) | **1 (1 %)** | 2 (2 %) | 117 (**98 %**) |
| `grafted` (54) | 0 (0 %) | 1 (2 %) | 53 (**98 %**) |

**98 % de nos candidats finissent derniers.** Le « place 3 sur 5 » est un
candidat sur 120, pas une propriété du bras. Ma formulation du §57 était
littéralement vraie et trompeuse : elle laissait entendre qu'un lot de designs
place un candidat à la 3ᵉ place, alors que c'est un événement à 1 %.

Ce qui reste notable, et qui demande confirmation : **le sélecteur a trouvé ce
candidat-là**, dans le bras `design`, avec deux critères différents, là où un
tirage au hasard avait 0,8 % de chances. C'est soit un bon signe sur le
sélecteur, soit de la chance, et un seul cas ne permet pas de trancher.

### Ce que la structure des métriques dit quand même

Le §57 tient sur un point : notre meilleur candidat est **dernier sur les cinq
métriques de confiance et premier sur les cinq autres**. Comme 98 % du vivier
finit dernier, l'avance sur le bloc interaction+développabilité ne suffit pas à
compenser cinq derniers rangs. Le levier reste le bloc confiance — `design_iiptm`,
`design_ptm`, `design_to_target_iptm` et les deux PAE — où les rivaux sont
devant sur la totalité du vivier, pas seulement sur le meilleur.

## 59. Le bloc confiance est une propriété de la charpente — et le cribler par la médiane coûterait notre meilleur candidat

Le §58 a désigné le bloc confiance (`design_iiptm`, `design_ptm`,
`design_to_target_iptm`, `min_design_to_target_pae`, `interaction_pae`) comme le
seul levier qui reste : nous y sommes derniers sur la totalité du vivier. Les
174 candidats du lot `design3arms_26345` permettent de chercher ce qui le
commande, sans nouveau GPU.

### C'est la charpente, et l'effet est énorme

Les 120 designs se répartissent sur les 6 charpentes de la spécification, qui
ne diffèrent que par leur FR1 :

| charpente | iiptm | ptm | iptm | min_pae | int_pae | n |
|---|---|---|---|---|---|---|
| 1 | **0,533** | 0,792 | 0,506 | 4,74 | 12,65 | 20 |
| 3 | 0,492 | 0,786 | 0,490 | 5,12 | 13,45 | 10 |
| 5 | 0,409 | 0,759 | 0,396 | 7,49 | 16,30 | 30 |
| 2 | 0,398 | 0,758 | 0,348 | 9,08 | 17,50 | 10 |
| 0 | 0,370 | 0,751 | 0,362 | 8,68 | 17,62 | 30 |
| 4 | **0,207** | 0,742 | 0,184 | 15,54 | 22,20 | 20 |

Un facteur **2,6** sur l'iiptm médian, et **le même ordre sur les cinq
métriques** : ce n'est pas du bruit, c'est la charpente qui porte le bloc
confiance. Dans le bras `grafted`, 20 charpentes distinctes donnent des
médianes de 0,172 à 0,604, soit un facteur **3,5**.

Le CDR3 conçu, lui, compte peu : la qualité de pose rendue par le pipeline de
conception ne prédit presque rien (`filter_rmsd` ρ = +0,01 sur les cinq
métriques, `filter_rmsd_design` ρ = +0,16 à +0,20). Cela confirme sur un run
valide ce que le §14 avait observé sur un run cassé.

### Mais la qualité d'une charpente ne se prédit pas à bas coût

| prédicteur | ρ avec l'iiptm médian de la charpente |
|---|---|
| identité V | +0,41 (n = 6, non concluant) |
| iiptm du donneur de CDR3 de la charpente | **+0,14** |

La pire charpente (0,207) portait le CDR3 du **deuxième meilleur** donneur
(0,764). Les CDR1 et CDR2 des six sont quasi identiques. Rien dans ce que nous
savons calculer sur CPU ne distingue la charpente 1 de la charpente 4.

### Le piège : cribler les charpentes par la médiane réduirait notre maximum

La conclusion qui vient naturellement — « cribler les 6 charpentes sur quelques
designs, puis concentrer le budget sur la meilleure » — est **fausse**, et c'est
mesurable :

| stratégie | max attendu sur 20 | 60 | 120 | 240 |
|---|---|---|---|---|
| budget concentré sur la meilleure médiane | 0,623 | 0,653 | 0,670 | 0,685 |
| budget étalé sur les six | **0,638** | **0,687** | **0,707** | **0,722** |

Et le fait brut qui l'explique : **notre meilleur design (0,741) vient de la
charpente de rang 4 sur 6** (médiane 0,398), pas de la première. Les cinq
meilleurs designs viennent des charpentes de rang 4, 1, 5, 5 et 1. Dans le bras
`grafted`, les trois meilleurs viennent des rangs 1, 5 et 2 sur 20.

Puisque nous ne soumettons qu'**une** séquence (`num_sequences: 1`), ce qui
compte est le maximum, pas la médiane — et la charpente à la meilleure médiane
a la queue la plus courte (borne 0,830 mais maximum observé 0,656, contre 0,741
pour le vivier étalé). **La diversité de charpentes est ce qui produit le
meilleur candidat**, et concentrer le budget l'aurait coûté.

### Ce que cela dit de la marche à suivre

- **Ne pas cribler les charpentes pour n'en garder qu'une.** Garder la
  diversité, et dépenser la profondeur de notation sur un vivier large — ce qui
  rejoint le §52 et le §53 par un troisième chemin.
- **Écarter les charpentes franchement mauvaises reste utile** : la charpente 4
  à 0,207 de médiane n'a donné aucun des cinq meilleurs designs. Il y a une
  différence entre éliminer le pire et ne garder que le meilleur.
- **Le bloc confiance reste ouvert.** Nous savons maintenant *où* il se décide
  (la charpente) mais pas comment choisir une bonne charpente autrement qu'en la
  mesurant. C'est la question à instruire ensuite, et elle est bon marché à
  poser : il suffit de noter quelques designs sur beaucoup de charpentes au lieu
  de beaucoup de designs sur six.

## 60. `--num_designs` échantillonne la liste de charpentes, il ne la parcourt pas : 9 sur 40

Lot `design_wide_26345` : 120 designs, spécification à **40 charpentes**,
BoltzGen 0.3.2. Preuves : `results/designs_design_wide_26345.csv`,
`results/design_design_wide_26345.log`.

### L'expérience que j'avais conçue n'a pas eu lieu

Le garde-fou passe (médiane 2 cystéines, aucun design à zéro, 120 sur 120
passent nos portes déterministes) mais il signale **9 variantes de FR1** pour
40 charpentes dans la spécification. Appariement fait contre `report.json` :

| charpentes dans la spec | effectivement utilisées | designs par charpente |
|---|---|---|
| 40 | **9** | 20 à 47 |

`--num_designs 120` ne distribue donc pas 3 designs sur chacune des 40 : il
**échantillonne** la liste. L'étape `design` affiche six lots dans les deux
runs et prend le même temps (2097 s contre 2175 s), quelle que soit la taille
de la liste. Mon plan « 40 × 3 » reposait sur une lecture fausse de l'outil ;
pour couvrir 40 charpentes il faudrait 40 travaux de conception, pas une
spécification à 40 entrées.

### Ce que le run dit quand même, et c'est en faveur de la diversité

Sur l'échelle de BoltzGen 0.3.2, donc comparable d'un run de conception à
l'autre mais **pas** aux métriques du validateur :

| | charpentes utilisées | médiane | max | designs > 0,60 | passent `filter_rmsd` |
|---|---|---|---|---|---|
| `design_cdr3_26345` | 6 | 0,241 | 0,650 | 23 | 3 / 120 |
| `design_wide_26345` | **9** | **0,349** | **0,708** | **43** | **11 / 120** |

Trois charpentes absentes de la première spécification valent mieux que tout
ce qu'elle contenait : la charpente 39 rend une médiane de **0,588** (max
0,708) et la 38 de 0,548, contre 0,533 pour la meilleure des six. La part de
designs au-dessus de 0,60 passe de 23 à 43, et la part qui passe le filtre de
pose du pipeline de 3 à 11 sur 120.

Donc **ouvrir le vivier de charpentes paie, même échantillonné** — neuf tirées
au hasard parmi quarante battent six choisies. C'est cohérent avec le §59 :
l'effet est porté par la charpente, sa qualité ne se prédit pas, et la seule
façon d'en trouver une bonne est d'en essayer.

### La prédiction que je m'étais engagé à lire, et pourquoi je ne la lis pas

J'avais écrit : « si le maximum dépasse 0,741, la diversité est le levier ».
Ce 0,741 est sur l'**échelle du validateur** et le 0,708 d'ici sur celle de
BoltzGen 0.3.2 — les comparer serait exactement l'erreur que le §56 interdit.
Le run précédent lisait 0,650 sur l'échelle BoltzGen et 0,741 sur celle du
validateur ; le rapport entre les deux n'est pas une constante qu'on peut
appliquer.

La prédiction reste donc ouverte, et le lot `wide4arms_26345` la tranche :
120 designs larges, les 30 meilleurs designs à 6 charpentes, les 20 meilleurs
greffés et les 4 rivaux réels, **notés dans le même run** par le moteur des
validateurs. 174 séquences, ~4,1 h sur 2×T4.

## 61. La prédiction se résout contre la diversité, et le « place 3 sur 5 » ne se reproduit pas

Lot `wide4arms_26345` : 174 séquences notées par le moteur des validateurs dans
un seul run — 120 designs à 9 charpentes, les 30 meilleurs designs du run à
6 charpentes, les 20 meilleurs greffés, et les 4 rivaux réels. Preuves :
`results/scores_nanobody_wide4arms_26345.csv`,
`results/boltzgen_wide4arms_26345.log`. Fidélité : écart absolu médian de
**0,012** sur les ancrages.

| bras | n | médiane | max | pae médian |
|---|---|---|---|---|
| `design_wide` (9 charpentes) | 120 | 0,393 | **0,699** | 17,64 |
| `design_6fw` (30 meilleurs du run précédent) | 30 | 0,574 | 0,709 | 12,07 |
| `grafted` (20 meilleurs) | 20 | 0,609 | **0,710** | 11,07 |
| `real` | 4 | 0,788 | 0,795 | 8,84 |

### La prédiction, lue comme je m'étais engagé à la lire

J'avais écrit, avant la mesure : « si le maximum dépasse **0,741**, la diversité
de charpentes est le levier et il faut en ouvrir davantage ; s'il plafonne vers
**0,70-0,74**, la queue est une propriété du générateur et non du choix de
charpente, et il faut revenir à la profondeur de notation du §53. »

`design_wide` lit **0,699**. Les trois bras plafonnent entre 0,699 et 0,710.
**La prédiction se résout contre la diversité de charpentes**, et la consigne
que j'avais inscrite s'applique : revenir à la profondeur de notation.

Les queues le confirment, sur l'échelle du validateur :

| bras | ξ | borne | max attendu sur 256 |
|---|---|---|---|
| `grafted` | −0,22 | **0,798** | **0,729** |
| `design_6fw` | −0,62 | 0,721 | 0,712 |
| `design_wide` | −0,54 | 0,707 | 0,691 |

Le greffage garde la meilleure queue, et élargir les charpentes la **raccourcit**
au lieu de l'allonger.

### Et le « place 3 sur 5 » du §57 ne se reproduit pas

Les 30 meilleurs designs du run précédent ont été **renotés ici**, et parmi eux
celui qui atteignait la place 3 sur 5 en battant 2 rivaux sur 4. Résultat :

| | place 3 ou mieux | place 5 |
|---|---|---|
| `design_wide` (120) | **0** | 120 |
| `design_6fw` (30, dont le candidat du §57) | **0** | 30 |
| `grafted` (20) | **0** | 20 |

**Les 170 candidats finissent derniers, sans exception.** Le candidat qui lisait
0,741 relit **0,709**, et il ne bat plus aucun rival. Sur les 50 séquences
notées dans les deux runs : écart absolu médian 0,015, σ d'une mesure 0,017,
Spearman +0,88 entre les deux lectures.

Donc le « meilleur résultat jamais mesuré dans ce dépôt » du §57 était **un
tirage favorable de la mesure**, et le §58 avait déjà montré qu'il s'agissait
d'un candidat sur 120. Il faut maintenant le retirer tout court : avec deux
mesures, il n'existe pas. C'est la troisième fois dans cette série qu'un
résultat annoncé sur une seule mesure ne survit pas à la seconde — après le
« 0 victoire sur 8 » et la borne de 0,39.

### Où cela laisse la compétition

Rien de ce que nous produisons ne bat un seul rival de cette époque, et le
meilleur rival lit 0,795 quand notre meilleure queue promet 0,729 sur
256 candidats notés. L'écart n'est pas comblé par la profondeur seule.

Ce qui reste acquis de tout cet arc : la conception est **légitime et
soumissible** (§55, §56 — 120 designs sur 120 passent les portes
déterministes), elle n'était pas réfutée mais cassée ; le bloc confiance est
porté par la charpente (§59) ; et la profondeur de notation reste le seul
levier qu'une mesure répétée confirme, sans suffire.

## 62. Reproduire la méthode du champ : les trois dernières différences testables sont éliminées, et le GPU ne suffit pas

Le §61 ayant refermé la diversité de charpentes, il restait à vérifier si la
méthode du champ diffère encore de la nôtre sur un point mesurable. Trois
candidats, tirés de sa propre signature dans `docs/FIELD_METHOD.md`. Tous les
trois sont éliminés **sans GPU**.

### 1. La bande d'identité V : un vrai défaut de notre côté, un petit effet

Notre sélection de charpentes importait `GERMLINE_MARGIN = 2.0` de
`propose_nanobody`. Cette marge existe pour le **greffage**, qui mute CDR1,
CDR2 et le corps du CDR3 à ~1 point d'identité V par mutation. Une conception
garde tout cela fixe, et la mesure sur les 120 designs du §60 le confirme :
elle coûte **0,00 point d'identité V en médiane et 1,02 au pire**.

Importer 2,0 nous interdisait donc la bande que le champ occupe. Nos
40 charpentes étaient à V 84,69-88,78 quand les soumissions acceptées sont à
**83,7 de médiane, 56 % à moins de deux points du plancher** — et 314 charpentes
conformes attendaient à V 82-84 sans jamais être utilisées. Corrigé :
`DESIGN_V_MARGIN = 1.1`, et les charpentes sont désormais tirées **près du
plancher**, ce qui place les 40 nouvelles à V 83,67, exactement la médiane du
champ.

Mais l'effet attendu est petit, et il faut le dire avant de dépenser quoi que
ce soit. Sur les 348 soumissions acceptées, l'iiptm médian par bande V vaut
0,698 (82-83), **0,719** (83-84), 0,680 (84-85), 0,696 (85-86), 0,682 (≥86) :
quasi plat. L'écart entre la bande 83-84 et la bande 85-86 est de **+0,023**,
contre **0,085** à combler. Le ρ intra-époque de −0,33 ordonne les candidats
d'une époque, il ne déplace pas le niveau.

### 2. La liberté sur le CDR3 : le champ est aussi contraint que nous

Nous ne concevons que 9 résidus sur 14 — la tête `AR` et la queue (`FDV`/`FDY`)
sont gelées. Et le gel **est respecté** par la 0.3.2 : tête conservée 60/60,
queue 59 à 60/60 sur les deux runs valides (l'observation contraire du
docstring de `scaffold_spec` venait du run cassé du §11).

Libérer la queue semblait donc être le levier. Il ne l'est pas, pour deux
raisons mesurées :

- la fenêtre J couvre 14 positions, donc **un pas vaut 7,14 points**. Depuis
  J = 100, **deux** écarts sur la queue suffisent à passer sous le plancher
  de 86. La marge n'existe pas ;
- surtout, les 348 soumissions acceptées du champ sont à J = **100** (43) ou
  **92,86** (305), c'est-à-dire **zéro ou un seul écart au germinal sur 14**.
  La queue de leur CDR3 est contrainte exactement comme la nôtre.

### 3. Ce qui reste, et pourquoi le GPU ne referme pas l'écart

Récapitulatif de ce qui a été éliminé comme explication des 0,085 :

| différence supposée | verdict |
|---|---|
| bande d'identité V | **+0,023** au mieux ; iiptm du champ plat selon la bande |
| diversité de charpentes | §61 : raccourcit la queue (borne 0,707 contre 0,798) |
| liberté sur la queue du CDR3 | le champ est à 0-1 écart germinal, comme nous |
| conception contre greffage | les deux plafonnent à 0,699-0,710 (§61) |
| profondeur de notation | seul levier confirmé, et insuffisant — voir ci-dessous |

Et le chiffrage de la profondeur, sur la queue du greffage mesurée au §61
(ξ = −0,22, borne 0,798) — nombre de candidats à noter pour en espérer **un**
au-dessus du seuil :

| seuil | candidats | durée sur 2×T4 |
|---|---|---|
| 0,710 (notre maximum actuel) | 82 | 1,9 h |
| 0,729 (max@256 prédit) | 248 | 5,8 h |
| 0,750 (au-dessus du 4ᵉ rival) | 1 313 | **30,6 h** |
| 0,769 (le rival le plus faible) | **13 384** | **312 h** |
| 0,788 (médiane des rivaux) | 2 millions | — |
| 0,795 (le meilleur rival) | 850 millions | — |

Le budget d'une époque est de **47 candidats** sur 2×T4, et le quota Kaggle
hebdomadaire d'environ 1 280. Battre le rival le **plus faible** de cette époque
demanderait 13 000 candidats, soit 312 h — et la médiane du champ est hors
d'atteinte par cette voie à **deux millions**.

**Conclusion.** La profondeur de notation est réelle mais ne change que la
deuxième décimale : elle nous mène de 0,710 à 0,729, pas à 0,769. Louer un GPU
mesurerait donc un débit, pas une victoire. Ce qui manque n'est pas du calcul,
c'est un générateur dont la **médiane** vaut 0,70 là où la nôtre vaut 0,61 — et
les cinq différences de méthode que ce dépôt savait tester sont maintenant
toutes éliminées. Ce qui reste n'est pas testable à notre échelle : le nombre de
designs que le champ génère réellement (le CLI suggère « typiquement 10 000 »),
les paramètres `--budget`, `--alpha`, `--metrics_override` du protocole, ou un
autre outil de conception — ODesign, RFdiffusion, ProteinMPNN — que le guide
amont cite et que nous n'avons jamais essayé.

## 63. Mille designs : le repliement de BoltzGen était du coût pur, et l'échantillonnage couvre la liste

Run `design_1k_26345` : **1000 designs** demandés sur la spécification à
40 charpentes, avec `--steps design inverse_folding`. Preuves :
`results/designs_design_1k_26345.csv`,
`results/design_design_1k_26345.log`.

### Ce que supprimer le repliement a rendu possible

| étape | 120 designs | 1000 designs |
|---|---|---|
| `design` | 2096 s | **19 347 s** (5,4 h) |
| `inverse_folding` | 67 s | 426 s |
| `folding` | 5466 s | **supprimée** |
| `analysis` + `filtering` | 103 s | supprimées |

Le temps de conception est bien linéaire (50 lots contre 6). Et les trois
étapes supprimées produisaient des métriques que nous jetions, puisque la
comparaison au champ exige la build des **validateurs** (§56). Mille designs
coûtent donc 5,5 h au lieu des ~21 h qu'aurait demandées le pipeline complet.

Le repli d'extraction a servi : sans filtrage BoltzGen n'écrit aucun CSV, et
les 992 séquences ont été sorties des structures.

### L'échantillonnage couvre la liste quand on demande assez de designs

Le §60 avait constaté que `--num_designs 120` n'utilisait que **9** des
40 charpentes. À 1000 designs :

| | charpentes utilisées | designs par charpente |
|---|---|---|
| 120 demandés | 9 / 40 | 20 à 47 |
| **1000 demandés** | **36 / 40** | **10 à 56** |

`--num_designs` échantillonne donc bien la liste, et le défaut du §60 n'était
pas structurel : il venait d'un budget trop petit pour couvrir 40 entrées. La
diversité de charpentes que j'avais essayé d'obtenir en élargissant la
spécification s'obtient en demandant plus de designs.

Et **992 designs sur 992 passent nos portes déterministes** — bon gène V, CDR3
de longueur 14, planchers germinaux, marqueurs VHH, médiane de 2 cystéines,
aucun à zéro. Mille séquences soumissibles pour 5,5 h de GPU gratuit.

### Ce que cela ne dit pas encore

Ces designs n'ont **aucune métrique de liaison** : nous avons supprimé l'étape
qui les produit. Les noter toutes coûterait 23 h sur 2×T4. Le lot
`k1_sample_26345` en prend un **échantillon stratifié** — 5 par charpente sur
36, soit 180 designs, plus les 4 rivaux réels comme ancrages, 184 séquences et
~4,3 h. Il répond à deux questions d'un coup :

1. **la queue s'allonge-t-elle avec mille tirages ?** Le §61 l'avait mesurée
   bornée à 0,707 sur 120 designs issus de 9 charpentes. Un échantillon de 180
   sur 36 charpentes permet d'ajuster la queue et d'extrapoler le maximum sur
   les 992 ;
2. **quelles charpentes portent ?** Cinq designs par charpente donnent un
   signal par charpente — de quoi diriger un second run concentré, ce que le
   §59 interdisait de faire sur la seule médiane mais qu'une mesure par
   charpente autorise.

## 64. Mille designs ne valent pas mieux que cent vingt : la conception est bornée vers 0,70

Lot `k1sample_26345` : échantillon stratifié de 180 designs du run à mille
(5 par charpente sur 36) et les 4 rivaux réels, notés par le moteur des
validateurs dans un seul run. Preuves :
`results/scores_nanobody_k1sample_26345.csv`,
`results/boltzgen_k1sample_26345.log`. Fidélité : écart absolu médian de
**0,014** sur les ancrages.

| | n | médiane | max |
|---|---|---|---|
| `design_1k` (échantillon des 992) | 180 | 0,353 | **0,692** |
| `real` | 4 | 0,783 | 0,785 |

### La queue ne s'allonge pas, et c'est une réponse définitive

| | ξ | borne | max attendu sur 992 | sur 5000 |
|---|---|---|---|---|
| §61 — 120 designs, **9** charpentes | −0,54 | 0,707 | — | — |
| ce run — 180 échantillonnés sur 992, **36** charpentes | **−0,61** | **0,698** | **0,692** | 0,696 |

Huit fois plus de designs, quatre fois plus de charpentes, et la borne est la
même à un centième près. Le maximum attendu sur **cinq mille** designs vaut
0,696. La question que le §52 avait laissée ouverte — et dont j'avais retiré la
borne parce qu'elle était ajustée sur des protéines cassées — est maintenant
tranchée sur un run valide : **concevoir le CDR3 avec BoltzGen sur une
charpente VHH conforme plafonne vers 0,70, et l'échelle n'y change rien.**

Le « typiquement 10 000 designs » du CLI amont ne nous amènerait donc pas au
niveau du champ. Ce n'est pas un problème de budget.

### Et aucun des 180 ne bat un rival

| place contre les 4 rivaux | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| designs | 0 | 0 | 0 | **1** | **179** |

Un seul atteint la 4ᵉ place, 179 finissent derniers. Le rival le plus faible de
cette époque lit 0,711 d'iiptm et notre meilleur design 0,692.

### Le signal par charpente existe, mais il ne se transfère pas au maximum

Les 36 charpentes mesurées à 5 designs chacune donnent des médianes de **0,142
à 0,619** — un facteur 4,4, ce qui confirme le §59 : le bloc confiance est porté
par la charpente. Mais la charpente qui porte notre **meilleur** design est au
**rang 25 sur 36** par médiane, et la pire par médiane (0,142) produit un design
à 0,587.

C'est exactement le piège du §59 sous une forme plus nette : la médiane d'une
charpente ne prédit pas son maximum, et comme nous ne soumettons qu'une
séquence, c'est le maximum qui compte. **Concentrer sur les bonnes médianes
reste une erreur**, même avec une mesure par charpente en main.

### Ce qui reste, après ce run

| voie | plafond mesuré | état |
|---|---|---|
| conception BoltzGen, charpente conforme | **0,698** | fermée par l'échelle (ce run) |
| greffage depuis l'archive | 0,798 | ouverte, mais 13 384 candidats pour 0,769 (§62) |
| le champ | 0,785-0,795 | — |

Le greffage garde la seule queue qui dépasse le niveau du champ, et le §62 en a
chiffré le prix. La conception, elle, est désormais bornée **sous** le rival le
plus faible : elle produit des séquences soumissibles en masse (992 sur 992) et
équilibrées sur les métriques d'interaction, mais elle ne franchira pas 0,70.

## 65. ODesign en *de novo* ne peut pas produire une soumission, et le mode qui pourrait est un autre test

### Pourquoi le *de novo* est fermé, par nature et non par rendement

Le validateur applique, **avant** toute notation, un filtre de nativité
(`metanano/validators/nativeness_validator.py`) qui exige que la séquence soit
numérotable en IMGT par `abnumber`, puis une nativité VHH ≥ 0,80 et une
humanité ≥ 0,75 par IgBLAST. Un binder conçu *de novo* n'est pas un anticorps :
il n'a ni gène V, ni gène J, ni architecture FR/CDR.

Démonstration sur une protéine réelle de la bonne taille — la chaîne A de 8AOK
elle-même, 116 résidus, dans la fenêtre 90-150 du validateur :

| | résultat |
|---|---|
| numérotation IMGT | **échec** |
| identité germinale V / J | **aucune** (pas de domaine à aligner) |
| nos portes déterministes | 3 refus, dont « ANARCI did not return exactly one heavy-chain domain » |

Ce n'est donc pas une question de taux de réussite qu'on pourrait compenser en
générant beaucoup : un repliement quelconque est rejeté **avant** d'être noté.
ODesign en *de novo* produirait des binders potentiellement excellents et
intégralement inéligibles.

### Le mode d'ODesign qui pourrait servir

Ses exemples montrent un `motif_scaffolding` dont la syntaxe mélange segments
repris d'une structure de référence et longueurs conçues :

```
"sequence": "1-80,A/1-11,25-35,B/1-11,15-25,C/1-11,20-30,D/1-10,1-80"
```

Appliqué à nous : garder les segments de charpente d'un VHH conforme et ne
laisser conçu que le CDR3. C'est la **même classe de contrainte** que notre
spécification BoltzGen — et c'est précisément ce qui en fait un test
intéressant, parce qu'il discrimine deux hypothèses que le §64 ne sépare pas :

- le plafond de **0,698** est celui de **BoltzGen** ; un autre modèle sous la
  même contrainte ferait mieux ;
- ou il est celui de la **contrainte** elle-même — une charpente VHH figée avec
  9 à 14 résidus libres ne permet pas de dépasser 0,70, quel que soit le modèle.

Aucune mesure de ce dépôt ne tranche entre les deux, et c'est la dernière
question ouverte qui ne soit pas une question de budget.

### Ce que monter ce test coûterait, mesuré

| élément | état |
|---|---|
| poids `odesign_base_prot_rigid` ou `_flex` | **4,32 Go** chacun, sur Hugging Face — accessible |
| `oinvfold_protein.ckpt` | sur Hugging Face — accessible |
| `components.v20240608.cif` + `.rdkit_mol.pkl` | **Google Drive**, dossier sans identifiants de fichiers, taille inconnue — et `src/utils/data/ccd.py` les charge **sans condition**, donc obligatoires même pour une conception protéine seule |
| pilotage | `inference_demo.sh` + hydra, et `conda install biotite` pour certains modes |

Les poids sont du même ordre que les ~6 Go de BoltzGen, qui passent sur Kaggle.
Le vrai risque est le dossier Google Drive : sans identifiant direct, un carnet
ne peut pas le récupérer de façon fiable, et le code refuse de démarrer sans.

C'est pourquoi ce test n'est pas lancé ici : la demande portait sur le *de
novo*, qui est fermé, et le mode qui reste suppose de résoudre d'abord un
blocage d'approvisionnement que rien dans ce dépôt ne contrôle.

## 66. 8AOK est un complexe nanobody–antigène, et nous n'avions jamais regardé la chaîne B

En creusant les points de méthode avant de monter ODesign, un fait resté invisible
depuis le début du projet : `data/structures/8AOK.cif`, que nous passons depuis
toujours comme « la cible », porte **deux** chaînes.

| chaîne | résidus | ce que c'est |
|---|---|---|
| A | 116 | la cible — PD-L1 (Q9NZQ7). Non numérotable IMGT, pas un anticorps |
| **B** | **127** | **un VHH cristallographique lié à elle** : numérotable IMGT, marqueurs VHH présents, V 78,6, J 85,7, CDR3 `AAEGTLYGSGGRTHQSAYDY` de longueur **20** |

Nous n'avons jamais utilisé que la chaîne A, en laissant BoltzGen inventer la
pose du nanobody — alors que le complexe de référence était dans le dépôt.

### Premier acquis : notre site de liaison est le bon

Contacts calculés sur les coordonnées, résidus de la chaîne A à 5 Å ou moins du
nanobody : **23 résidus**. Notre `BINDING_SITE`, repris de la configuration des
validateurs, en compte 21.

| | |
|---|---|
| de notre site, dans le vrai épitope | **21 / 21 = 100 %** |
| du vrai épitope, dans notre site | 21 / 23 = **91 %** (manquent 93 et 103) |

Les validateurs ont donc dérivé leur site du complexe lui-même, et notre
spécification désigne le bon épitope. Un point de méthode vérifié plutôt que
supposé — il restait possible que nous concevions depuis des mois contre la
mauvaise face de la protéine.

### Deuxième acquis : un CDR3 validé par cristallographie, et jamais soumis

Le CDR3 de la chaîne B lie la cible — c'est la structure qui le prouve. Et il
est **absent de notre vivier** : 0 occurrence sur les 474 CDR3 de longueur 20
que porte l'archive, le plus proche à **45 %** d'identité seulement.

Deux conséquences pour les portes du validateur :

- l'unicité globale du Layer 1 porte sur les séquences **déjà soumises**. Celle-ci
  vient d'une structure du PDB, pas de l'archive du concours : une greffe de ce
  CDR3 sur une charpente conforme est donc une séquence neuve ;
- la diversité (Jaccard k-mer < 0,9 contre le top 50, §54) est très loin de
  mordre à 45 % d'identité au plus proche.

L'époque **26301** impose justement une longueur de CDR3 de 20. Greffé sur les
charpentes conformes `IGHV3-11` disponibles à cette époque, ce CDR3 donne
**25 candidats qui passent toutes nos portes déterministes** (31 refusés sur
57 charpentes essayées), à V 84,7-92,9 et J 92,9.

### Ce que cela teste, et contre quoi

Le §50 avait conclu que « la liaison ne voyage pas avec le motif », sur des
donneurs d'archive déplacés d'une charpente à l'autre. Ici le motif vient du
**complexe cristallographique de la cible** : c'est la meilleure provenance
possible, et le test le plus favorable que cette hypothèse puisse recevoir. Si
la liaison ne voyage pas même dans ce cas, le §50 devient une loi de ce
problème plutôt qu'une observation.

Lot `xtal_26301` : **39 séquences, ~0,9 h** — 25 greffes du CDR3
cristallographique, 10 greffes ordinaires de la même époque comme témoin, et
les 4 meilleures soumissions réelles de 26301 comme ancrages. Le champ y lit
0,703 de médiane et 0,739 au mieux.

## 67. Nouvelle cible P05231, et deux règles que ce dépôt avait reconstruites de travers

La cible nanobody a changé le 20 novembre : **Q9NZQ7 (PD-L1, structure 8AOK)**
laisse la place à **P05231 (interleukine-6, structure 4O9H)**, avec
`clip_interval [[27, 212]]`, chaîne A, `res_index 21..186` et un site de liaison
de **8 résidus** au lieu de 21.

`nova_lab/config.py` porte désormais ces constantes, et
`scripts/boltzgen_score.py` comme `scripts/design_nanobody.py` les lisent au lieu
de littéraux : la bascule a montré que `8AOK` était en dur dans six scripts.
L'ancienne cible est gardée dans `PREVIOUS_NANOBODY`, puisque tout l'historique
de mesure de ce dépôt porte sur elle.

### Le site de liaison est juste, mais pas dans la numérotation qu'on croit

Contacts calculés sur 4O9H : **19 résidus** de la chaîne A à 5 Å ou moins du
Fab lié. Comparé naïvement, le site de la config n'en retrouve que **2 sur 8**.
Au décalage **−2**, il en retrouve **8 sur 8** :

| config | 24 | 77 | 80 | 82 | 131 | 184 | 185 | 186 |
|---|---|---|---|---|---|---|---|---|
| auteur | 22 | 75 | 78 | 80 | 129 | 182 | 183 | 184 |
| contact réel | oui | oui | oui | oui | oui | oui | oui | oui |

Les numéros de la config sont donc des `label_seq_id` de mmCIF — l'index
canonique séquentiel, qui est précisément ce que BoltzGen lit — et non la
numérotation auteur du cristal. La spécification peut les reprendre tels quels,
comme pour 8AOK. **Ce contrôle est à refaire à chaque changement de cible** : il
sépare « notre site est le bon » de « nous concevons contre la mauvaise face ».

Et 4O9H n'est pas un complexe nanobody comme 8AOK : chaîne A l'IL-6, chaîne B
une **chaîne lourde** d'anticorps conventionnel (200 résidus, V 88,8, CDR3
`ANRAGWGMGDY` de longueur 11), chaîne C une chaîne légère. Pas de charpente VHH
à réutiliser — mais le CDR3 de la chaîne B est un motif de liaison validé par
cristallographie, exploitable comme le §66 l'a fait pour 8AOK, sur les époques
imposant une longueur 11.

### L'archive de la nouvelle cible est incomparablement plus riche

| | Q9NZQ7 | **P05231** |
|---|---|---|
| séquences | 17 244 | **42 328** |
| `design_iiptm` médian | ~0,22 | **0,648** |
| à ≥ 0,72 | 5 064 | **17 489** |
| à ≥ 0,80 | — | 12 247 |
| maximum | 0,786 | **0,854** |

### Mais 0 sur 700 passe nos planchers germinaux — et nos planchers sont faux

Sur 700 lieurs à ≥ 0,72 tirés au hasard : **aucun** ne passe nos
`MIN_HUMAN_V_IDENTITY = 82` et `MIN_HUMAN_J_IDENTITY = 86`. Identité V médiane
**78,6**, identité J médiane **50,0**.

Avant d'en conclure que l'archive est inutilisable, j'ai lu le filtre réel.
`metanano/utils/igblast_nativeness.py` rejette durement sur :

| seuil réel | valeur |
|---|---|
| identité V | **≥ 0,70** |
| couverture V | ≥ 0,80 |
| longueur de CDR3 | 6 à 35 |
| identité J | **aucun seuil** |

Puis un score composite ≥ 0,80. Nos 82 et 86 sont donc une **reconstruction**,
plus stricte que la règle sur V (82 contre 70) et qui **invente** le plancher J.
C'est la deuxième fois après le §54 — où notre porte de nouveauté appliquait une
règle tirée du moteur de recherche et non du validateur — que ce dépôt s'interdit
du terrain que le concours autorise. Et la preuve empirique est dans l'archive
elle-même : 17 489 soumissions **notées** lient à ≥ 0,72 avec un J médian de 50.
Si un plancher J à 86 existait, elles n'auraient jamais été notées.

Réserve honnête : l'archive n'est pas datable, donc une partie peut précéder
l'introduction du filtre, et le score composite ≥ 0,80 reste exigeant même avec
V à 0,70. Nos planchers ne sont pas absurdes comme *approximation* de « passe la
nativité » — ils ne sont simplement pas la règle.

### Et la formule de nativité contient une contrainte que personne n'avait vue

```
nativité VHH = 0,25·A + 0,30·B + 0,20·C + 0,25·D     (seuil 0,80)
C = exp(−((longueur_CDR3 − 18) / 8)²)
```

**C culmine à une longueur de CDR3 de 18** et s'effondre loin d'elle :

| longueur | 6 | 8 | 10 | 11 | 14 | 18 | 20 | 23 |
|---|---|---|---|---|---|---|---|---|
| C | 0,105 | 0,210 | 0,368 | 0,465 | 0,779 | **1,000** | 0,939 | 0,677 |
| reste à trouver sur A, B, D | 0,779 | 0,758 | 0,726 | 0,707 | 0,644 | 0,600 | 0,612 | 0,665 |

Le maximum de A + B + D vaut 0,80. À une longueur de 18 il faut donc 75 % de ce
maximum ; à une longueur de **6**, il en faut **97 %**. Sur les 12 époques
connues, **4 imposent une longueur à C < 0,5** (6, 8, 10 et 11) : sur
celles-là, la nativité est quasi impossible à satisfaire, quelle que soit la
qualité de liaison.

C'est une explication structurelle que nous n'avions pas : les époques à CDR3
court ne sont pas seulement pauvres en donneurs, elles sont **pénalisées par le
filtre lui-même**. Et cela donne une consigne nette : sur ces époques, viser la
conformité germinale maximale (A, B, D) plutôt que la liaison.

## 68. Le CDR3 cristallographique greffé est le pire bras mesuré : le §50 devient une loi

Lot `xtal_26301` : 39 séquences notées dans un seul run — 25 greffes du CDR3
de la chaîne B de 8AOK (`AAEGTLYGSGGRTHQSAYDY`, lieur prouvé par la structure),
10 greffes ordinaires de la même époque en témoin, et les 4 meilleures
soumissions réelles de 26301 en ancrage. Preuves :
`results/scores_nanobody_xtal_26301.csv`,
`results/boltzgen_xtal_26301.log`.

| bras | n | min | médiane | max | pae médian |
|---|---|---|---|---|---|
| `real` | 4 | 0,691 | 0,700 | **0,723** | **9,19** |
| `grafted_ordinaire` | 10 | 0,168 | 0,250 | 0,595 | 20,10 |
| **`xtal_cdr3`** | 25 | 0,136 | **0,184** | **0,338** | 22,48 |

**Le motif cristallographique est le pire des trois.** Aucune de ses 25 greffes
n'atteint le meilleur greffage ordinaire (0,595), 3 sur 25 seulement dépassent
la *médiane* ordinaire, et aucune n'approche le rival le plus faible (0,691).
Au rang-somme : 93 au mieux, contre 56 pour les greffes ordinaires et 35 à 43
pour les rivaux.

### Pourquoi c'est un résultat fort et non un échec de plus

Le §66 avait posé le test comme **le plus favorable que l'hypothèse du §50
puisse recevoir** : le motif ne vient pas de l'archive mais du complexe
cristallographique de la cible, où la structure prouve qu'il lie. Si la liaison
voyageait avec le motif, c'était là.

Elle ne voyage pas. « La liaison ne voyage pas avec le motif » cesse d'être une
observation sur des donneurs d'archive et devient **une loi de ce problème**.

Et le mécanisme se lit dans le `pae` : 22,5 de médiane pour les greffes
cristallographiques contre 9,2 pour les rivaux. Le nanobody ne se pose pas. Le
CDR3 de 8AOK a été sélectionné par l'évolution dans le contexte d'une charpente
de lama à **V 78,6** ; posé sur une charpente germinale humaine conforme, il
perd la géométrie qui le faisait fonctionner. Les CDR3 de l'archive, eux, ont
co-évolué avec des charpentes du même vivier, ce qui explique qu'ils se
transplantent moins mal (0,595 contre 0,338).

### Ce que cela annule

**Le carnet ODesign ne doit pas être lancé.** Sa condition de lancement était
écrite au §65 et dans `docs/KAGGLE.md` : « ne lancer ce carnet que si
`xtal_26301` dit que la greffe conserve la liaison ». Elle dit non. ODesign
produirait des boucles conçues dans la vraie pose, donc inutilisables, puisque
sa charpente de référence (chaîne B de 8AOK, V 78,6 / J 85,7) est elle-même non
conforme et que transplanter la boucle détruit ce qu'elle a gagné. Une session
de quota économisée par une mesure d'une heure.

**Et la piste du motif validé sur la nouvelle cible tombe aussi.** Le §67 avait
noté que la chaîne B de 4O9H porte un CDR3 de longueur 11 validé par
cristallographie, exploitable « comme le §66 l'a fait ». Ce run montre que ce
n'était pas exploitable sur l'ancienne cible, et il n'y a aucune raison que
P05231 se comporte autrement — la charpente y est en plus une chaîne lourde
conventionnelle, appariée à une chaîne légère, donc encore plus éloignée d'un
VHH que ne l'était la chaîne B de 8AOK.

### Ce qui reste debout

| voie | plafond mesuré | état |
|---|---|---|
| greffage depuis l'archive | 0,798 (borne) | la seule qui dépasse le champ, au prix du §62 |
| conception BoltzGen | 0,698 | fermée, et l'échelle n'y fait rien (§64) |
| motif validé transplanté | **0,338** | fermée par ce run |
| ODesign | — | sans objet tant que la greffe ne conserve rien |

Reste donc, pour la nouvelle cible : une archive de 42 328 séquences dont
17 489 lient à ≥ 0,72, non conformes à nos planchers — mais nos planchers sont
faux (§67), et le vrai seuil est une identité V ≥ 0,70 sans plancher J. C'est
là qu'il faut regarder maintenant, et non du côté des motifs transplantés.

## 69. Le filtre de nativité, exécuté et non reproduit : le seuil de 0,80 n'est pas appliqué, celui d'humanité l'est

Plutôt que de reproduire la formule du §67 — ce qui aurait été une troisième
reconstruction après les deux fausses du §54 et du §67 — le filtre lui-même a
été rendu exécutable ici. `scripts/nativeness.py` appelle
`metanano.utils.igblast_nativeness.run`, sur les bases pré-construites du dépôt
du filtre.

### Ce qu'il a fallu, et le piège qui a failli produire un faux résultat

| pièce | provenance |
|---|---|
| binaires IgBLAST 1.22.0 | NCBI, `ncbi-igblast-1.22.0-x64-linux.tar.gz` |
| bases `camelid_V` et `human_V` | pré-construites dans NOVA-nanobody-filter, `custom_data/igblast_database/database` |
| `internal_data/camelid` | le même dépôt |
| **`internal_data/human`** | **l'archive NCBI** — absent du dépôt du filtre |

Sans cette dernière pièce, chaque appel mourait sur « Germline annotation
database .../internal_data/human/human_V could not be found ». Et
`NativenessFilter.analyze` **capture toute exception et rend `None`** : le
premier passage a donc noté les quatre soumissions **connues comme acceptées**
en « REFUSÉ ». C'était un échec d'installation lu comme un verdict — la même
classe d'erreur que le DSSP silencieux du §43, attrapée cette fois parce que le
contrôle portait sur des séquences dont le verdict est connu.

### Les trois seuils, mesurés contre des verdicts connus

| | rejets durs | nativité ≥ 0,80 | humanité ≥ 0,75 |
|---|---|---|---|
| 60 soumissions **acceptées** de Q9NZQ7 | **0 / 60** | **0 / 60** (médiane 0,507, max 0,596) | **40 / 40** |
| 60-80 lieurs de l'archive P05231 à ≥ 0,72 | **0 / 60** | 0 / 60 (médiane **0,646**) | **17 / 80 = 21 %** |

Trois lectures, et la deuxième contredit ce que j'avais écrit :

1. **Les rejets durs ne mordent pas** — V identité ≥ 0,70, couverture ≥ 0,80,
   CDR3 entre 6 et 35 : aucune séquence, d'aucune des deux cibles, n'y échoue.
2. **Le seuil de nativité de 0,80 n'est pas appliqué.** Aucune des
   60 soumissions acceptées ne l'atteint. S'il l'était, rien ne serait jamais
   accepté. Il est donc rapporté par `scripts/nativeness.py`, jamais enforcé.
3. **Le seuil d'humanité de 0,75, lui, discrimine** — et j'avais écrit l'inverse
   dans la sortie de ma première sonde, avant d'avoir le chiffre. Les 40
   acceptées le passent toutes ; seuls 21 % des lieurs de l'archive P05231 le
   passent. C'est le seul des trois qui sépare les deux populations.

Détail qui compte : les lieurs de P05231 ont une nativité **supérieure** à celle
des acceptées (0,646 contre 0,507) tout en échouant l'humanité. Les deux scores
ne mesurent pas la même chose, et c'est l'humanité qui coûte.

### Ce que cela change pour la nouvelle cible

| | Q9NZQ7 | P05231 |
|---|---|---|
| lieurs à ≥ 0,72 | 5 064 | 17 489 |
| utilisables selon **nos** planchers (§67) | 256 | **0** |
| utilisables selon l'**humanité** du vrai filtre | — | **~3 700** (21 %) |

Nos `MIN_HUMAN_V_IDENTITY = 82` et `MIN_HUMAN_J_IDENTITY = 86` rejetaient
**100 %** des lieurs de la nouvelle cible. Le vrai filtre en rejette 79 %, et
laisse donc de l'ordre de **3 700 donneurs** qui lient à 0,72 ou mieux — contre
256 sur l'ancienne cible, soit **quatorze fois plus**.

Réserve : « l'humanité discrimine » est une inférence sur 40 acceptées d'une
cible et 80 entrées d'archive de l'autre, pas une lecture du code de production.
Et le §68 vient de montrer que la liaison ne voyage pas avec le motif : disposer
de 3 700 donneurs ne dit pas encore qu'on sait s'en servir. Mais le vivier
cesse d'être le goulot, et c'est la première fois depuis le §48.

## 70. La porte germinale recalibrée, et le mur de la nouvelle cible : tout tient à la longueur 21

### Le crible calibré sur le vrai filtre

Le §69 a établi que la seule porte germinale qui mord est l'humanité ≥ 0,75 du
filtre, et qu'elle coûte une seconde d'IgBLAST par séquence — douze heures pour
les 42 328 de l'archive. Mesuré sur 190 séquences couvrant V 76,5 à 89,8,
l'identité V d'abnumber, que ce dépôt calcule en millisecondes, suit l'humanité
à **Spearman +0,93**. L'identité J, elle, à **−0,15** : notre plancher J n'était
pas seulement inventé, il était sans rapport avec quoi que ce soit.

| seuil V | précision | rappel | exactitude |
|---|---|---|---|
| 78,0 | 0,53 | 1,00 | 0,62 |
| **79,0** | **0,88** | **0,96** | **0,93** |
| 80,0 | 0,96 | 0,77 | 0,88 |
| **82,0 (le nôtre)** | 0,94 | **0,57** | 0,79 |

`nova_lab/nativeness_gate.py` retient **V ≥ 79**, optimisé sur le **rappel** :
un faux positif coûte un appel IgBLAST dans `verify`, un faux négatif perd un
donneur définitivement. Notre plancher à 82 écartait **43 %** des séquences que
la vraie porte accepte.

Sur un échantillon de 600 lieurs de P05231 à ≥ 0,72 :

| | retenus |
|---|---|
| crible calibré (V ≥ 79) | **154 / 600 = 26 %** → ~4 500 donneurs |
| notre ancien plancher (V ≥ 82) | **1 / 600** |
| avec le plancher J à 86 en plus | **0 / 600** |

### Mais les bons donneurs sont tous à une seule longueur de CDR3

Sur 900 séquences de l'archive tirées au hasard, toutes qualités confondues :

| longueur de CDR3 | n | iiptm médian | max |
|---|---|---|---|
| 13 | 529 | 0,618 | 0,683 |
| **21** | 348 | **0,814** | **0,854** |
| 10 | 6 | 0,602 | 0,639 |
| 17 | 4 | 0,222 | 0,296 |
| 20 | 3 | 0,205 | 0,236 |
| 19 | 3 | 0,164 | 0,168 |

**100 % des lieurs à ≥ 0,72 (340 sur 340) ont un CDR3 de longueur 21.**
L'archive est une monoculture bimodale : un sommet à 21 et un plateau à 13.
Toutes les autres longueurs sont rares *et* mauvaises.

C'est le même mur mécanique que sur l'ancienne cible — où 4 133 des 4 134
séquences à ≥ 0,70 avaient un CDR3 de longueur 13 — mais déplacé, et moins
étroit : **deux** longueurs exploitables au lieu d'une, et celle du sommet
porte des donneurs à **0,822 de médiane**, au-dessus du meilleur rival jamais
mesuré sur Q9NZQ7 (0,795).

### Et cela dit exactement quand le GPU sert

Sur 20 000 hashes simulés, la dérivation impose chaque longueur de 6 à 23 avec
une probabilité quasi uniforme de ~5,3 % :

| longueur imposée | fréquence | ce que nous avons |
|---|---|---|
| **21** | **5,5 %** | ~4 500 donneurs, médiane 0,822, max 0,850 |
| **13** | **5,5 %** | donneurs à 0,618 de médiane |
| les 16 autres | 89 % | rien d'exploitable |

**Une époque sur neuf** tombe sur une longueur où nous avons de quoi travailler,
et une sur dix-huit sur celle qui porte les très bons donneurs. Le reste du
temps, le vivier de cette cible ne vaut pas mieux que celui de l'ancienne.

Conséquence opérationnelle, et c'est la réponse à « quand relancer Kaggle » :
il ne sert à rien de noter des candidats pour une époque de longueur
quelconque. Le GPU doit être dépensé **sur une époque de longueur 21** — ou 13
en second choix — et c'est aussi là, et seulement là, que le §68 peut être
retesté utilement, puisque greffer depuis un donneur à 0,85 est le seul cas où
la perte à la greffe laisserait encore quelque chose.

## 71. Pourquoi nous ne gagnons pas : nous nous sommes mesurés sur les époques où le champ est le plus fort

Question posée directement : sur quoi bloque-t-on ? La réponse tient en deux
chiffres que ce dépôt n'avait jamais rapprochés.

### Le champ n'est pas également fort, et l'écart entre époques dépasse notre retard

| époque | L | gagnant réel | médiane du champ |
|---|---|---|---|
| 26348 | 6 | **0,589** | 0,458 |
| 26300 | 8 | **0,657** | 0,562 |
| 26305 | 10 | **0,663** | 0,699 |
| 26347 | 19 | 0,723 | 0,687 |
| 26301 | 20 | 0,727 | 0,703 |
| … | | | |
| 26346 | 13 | 0,787 | 0,734 |
| 26345 | 14 | **0,788** | 0,715 |

Le gagnant d'une époque vaut entre **0,589 et 0,788**. Cet écart de **0,199**
est plus grand que notre retard sur la médiane du champ. Et notre meilleur
candidat honnête jamais mesuré, 0,710, **aurait gagné 3 de ces 12 époques**.

### Or nos cinq runs GPU ont tous visé des époques où le champ est fort

Les lots `deep`, `design3arms`, `wide4arms`, `k1sample` et `xtal` ont porté sur
26304 (gagnant 0,763), 26345 (**0,788**), 26347 (0,723) et 26301 (0,727).
**Aucun** n'a porté sur 26348, 26300 ou 26305 — les trois époques où le gagnant
est sous 0,670. Nous avons donc passé une semaine à mesurer notre distance au
plus fort, et conclu que l'écart était structurel.

### Et sur ces époques, notre plafond de donneur dépasse déjà le gagnant

Avec la porte recalibrée du §70 (V ≥ 79, pas de plancher J), le meilleur donneur
honnête de la bonne longueur contre le gagnant réel :

| époque | L | donneurs | notre plafond | gagnant | écart |
|---|---|---|---|---|---|
| **26348** | 6 | 188 | **0,690** | 0,589 | **−0,101** |
| **26346** | 13 | 4767 | **0,827** | 0,787 | **−0,039** |
| **26300** | 8 | 204 | **0,675** | 0,657 | **−0,018** |
| **26304** | 14 | 406 | **0,774** | 0,763 | **−0,011** |
| 26345 | 14 | 438 | 0,779 | 0,788 | +0,010 |
| 26305 | 10 | 175 | 0,630 | 0,663 | +0,033 |
| 26347 | 19 | 122 | 0,675 | 0,723 | +0,047 |
| 26303 | 18 | 146 | 0,555 | 0,754 | +0,199 |

**Sur 4 époques des 12, notre vivier contient un donneur meilleur que le
gagnant.** L'écart médian n'est que de +0,040 — pas les 0,211 que le §49
annonçait, parce que ce 0,211 comparait des médianes de bras sur une époque
difficile, avec un vivier amputé par un plancher J inventé.

### Donc ce qui bloque, nommé précisément

Ce n'est pas la technique, et ce n'est pas le calcul :

1. **Ce n'est pas le générateur.** Le greffage a une borne de 0,798 et des
   donneurs à 0,827 sur certaines époques — au-dessus des gagnants réels.
2. **Ce n'est pas la conception.** Elle plafonne à 0,698 (§64) mais elle n'est
   pas nécessaire : le vivier suffit sur les bonnes époques.
3. **C'est la perte à la greffe.** Un donneur à 0,827 donne un candidat à 0,710
   au mieux. Le §68 a montré que la liaison ne voyage pas avec le motif, même
   cristallographique. **C'est le seul verrou.**
4. **Et c'est le choix des époques.** Nous avons systématiquement mesuré sur les
   plus dures. Sur 26348, le gagnant est à 0,589 et nous avons un donneur à
   0,690 : il suffirait de conserver 85 % de la liaison du donneur pour gagner,
   là où le §49 mesurait 28 % de rétention médiane sur l'humanisation.

La question n'est donc plus « comment atteindre 0,79 » mais « comment conserver
85 % d'un donneur à 0,69 », et elle se pose sur les époques à CDR3 court — que
le §67 a montrées structurellement pénalisées pour **tout le monde** par le
terme C du score de nativité.

### Ce que cela change pour le prochain run

Le §70 concluait qu'il fallait attendre une époque de longueur 21 sur la
nouvelle cible, parce que c'est là que sont les bons donneurs. Cette section le
corrige : il faut viser une époque où **l'écart entre notre plafond de donneur
et le gagnant est favorable**, ce qui n'est pas la même chose. Sur l'ancienne
cible, c'étaient les longueurs 6, 8, 13 et 14. Le critère n'est pas la qualité
absolue du vivier mais sa qualité **relative à ce que le champ parvient à faire
à cette longueur**.

## 72. L'hypothèse conformationnelle, telle que mesurable sur CPU : réfutée

Le §68 avait montré *où* la greffe casse : `pae` de 22,5 contre 9,2 pour les
rivaux, donc le nanobody ne se pose pas. L'hypothèse qui en découlait était que
ce n'est pas la séquence du CDR3 qui est perdue mais sa **conformation** — un
CDR3 est tenu par les résidus de charpente qui le bordent, et changer de
charpente le laisse flotter. Cela expliquerait aussi pourquoi aucun descripteur
de séquence ne prédit la rétention (ρ = −0,13 sur les mutations de charpente,
+0,027 sur l'iiptm du donneur, §49).

Test : 120 de nos candidats notés, échantillonnés par quantile d'iiptm pour
couvrir toute la plage mesurée (0,129 à 0,686), repliés par NanoBodyBuilder2, et
l'erreur prédite par résidu — que le modèle écrit dans la colonne B — moyennée
sur le CDR3. Preuve : `results/conformational_screen_120.csv`.

Seuils de lecture, écrits avant la mesure : ρ au-delà de ±0,35 sur le rapport
CDR3/charpente rendait le crible utile ; entre 0,20 et 0,35 il devenait un
descripteur de plus ; sous 0,20 l'hypothèse tombait.

| grandeur | ρ avec l'iiptm | p |
|---|---|---|
| erreur prédite sur le CDR3 | **−0,011** | 0,91 |
| erreur prédite sur la charpente (témoin) | +0,027 | 0,77 |
| maximum sur le CDR3 | −0,041 | 0,66 |
| **rapport CDR3 / charpente** | **−0,014** | 0,88 |

Et les médianes par quartile d'iiptm sont plates :

| quartile | bas | moyen− | moyen+ | haut |
|---|---|---|---|---|
| erreur CDR3 | 1,476 | 1,591 | 1,565 | 1,465 |
| rapport | 4,10 | 4,33 | 4,28 | 4,23 |

**Zéro signal**, très loin du seuil de 0,20. L'hypothèse est réfutée dans cette
forme.

### Ce qui a été testé, et ce qui ne l'a pas été

La grandeur mesurée est la **confiance du modèle en lui-même**, pas la
conformation. NanoBodyBuilder2 est incertain sur tout CDR3 — 1,5 Å contre
0,35 Å sur la charpente, un rapport de 4,2 quelle que soit la qualité du
candidat — et cette incertitude ne distingue pas un CDR3 contraint d'un CDR3
simplement long.

La version forte reste ouverte : superposer le candidat et **son donneur** sur
la charpente, et mesurer le RMSD du CDR3 entre les deux. Elle coûte deux
repliements par paire et l'appariement du donneur, que la mutation du CDR3
(`MAX_CDR3_IDENTITY = 0,85`) rend approximatif. Elle n'est pas lancée ici
parce que le témoin de ce run dit quelque chose de gênant pour elle : si
l'erreur prédite ne varie pas du tout avec la liaison, c'est que NanoBodyBuilder2
ne voit pas la différence entre une greffe qui marche et une qui ne marche pas —
et c'est le même modèle qui fournirait les deux conformations à comparer.

### Ce que cela laisse

Quatre voies ont été fermées par mesure : la conception BoltzGen (borne 0,698,
§64), le motif validé transplanté (0,338, §68), le ML de séquence comme
sélecteur suffisant (ρ = 0,31 groupé, trouve 0,594 là où le vivier contient
0,741), et le crible conformationnel accessible sur CPU (ici).

Reste ce que le §71 a mis au jour et que la stratégie de ce dépôt écarte depuis
le début : sur 4 époques des 12 connues, **notre vivier contient déjà un donneur
meilleur que le gagnant réel**, et les validateurs notent gratuitement ce que
nous leur soumettons. Vingt heures de quota Kaggle ont servi à approximer un
verdict que le subnet rend sans rien demander, et `scripts.mine` n'a journalisé
**qu'une époque, en dry-run, zéro soumission**.

## 73. La phase d'observation de la cible, enfin faite — et elle rend un levier

Le guide amont consacre son §5 à « Study the Target » : télécharger la
structure, **identifier le site de liaison**, chercher les liants connus. Le §66
l'avait fait pour l'ancienne cible. Pour `P05231` / `4O9H` je ne l'avais pas
fait : j'ai écrit la config et je suis parti tester des hypothèses de greffe.
C'était une erreur de méthode, et c'est le seul levier qui ne coûte ni GPU ni
TAO. Scripts : `results/{epitope,offset,label,converge,survey,vhh6}.py`.

### La config est juste, et nous la lisons dans la bonne convention

`4O9H` = **interleukine-6 en complexe avec le Fab camélidé 61H7** (chaînes H+L),
pas un VHH. Chaîne A = IL-6, 152 résidus observés.

La numérotation d'abord, parce qu'elle décide de tout : `label_seq_id = auth + 2`,
sans insertion. `res_index: 21..186` est donc **exactement l'empan observé** de la
chaîne A — rien ne tombe hors fenêtre — avec deux boucles désordonnées absentes
(label 51-62 et 107-108, hors épitope). Et le site de liaison, lu en index
mmcif canonique comme BoltzGen l'exige :

| config | résidu | auth | distance au Fab |
|---|---|---|---|
| 24 | S | 22 | 2,82 Å (H+L) |
| 77 | Q | 75 | 3,23 Å (H) |
| 80 | F | 78 | 3,40 Å (H) |
| 82 | E | 80 | 2,69 Å (H) |
| 131 | K | 129 | 4,00 Å (H) |
| 184 | R | 182 | 2,67 Å (H+L) |
| 185 | Q | 183 | 2,73 Å (H) |
| 186 | M | 184 | 3,00 Å (H) |

**8/8 au contact.** J'ai cherché un bug de décalage et il n'y en a pas : le −2
que j'avais d'abord vu est le passage label→auth, et `boltzgen_score.py` comme
`design_nanobody.py` passent bien la numérotation canonique.

Mieux : **la config est précisément l'empreinte de la chaîne lourde.** Six des
huit résidus ne sont touchés que par H, deux par H+L. Les sept résidus que seule
la chaîne légère touche (auth T20 S21 E23 D26 K27 R30 L178) sont **exclus** du
site. Le subnet a retiré la part qu'un domaine unique ne peut pas couvrir : la
géométrie demandée est atteignable par un VHH.

### Les liants connus ne convergent pas, sauf celui qui compte

17 entrées PDB contiennent IL-6, six avec un liant. Recouvrement avec les 8
résidus de la config :

| PDB | liant | épitope | recouvrement |
|---|---|---|---|
| 4O9H | Fab camélidé 61H7 | 19 res | **8/8** |
| 5FUC | **VHH6, domaine unique** | 11 res | **7/8** |
| 4ZS7 | Fab lama 68F2 | 23 res | 3/8 |
| 8YWQ | Fab I9H | 13 res | 2/8 |
| 8YWR | Fab 68F2 (2e cristal) | 12 res | 1/8 |
| 4CNI | Fab Olokizumab | — | 0/8 |

Le consensus à ≥4 liants sur 6 ne retient que 5 résidus (D28 R32 C75 F76 Q77),
dont un seul dans la config : **ces anticorps visent des sites différents**. Mais
le seul VHH authentique du PDB contre IL-6 couvre 7 des 8 résidus de la config.
Le site n'est donc pas « un épitope de Fab » : c'est le site validé pour un
domaine unique.

Réserve honnête sur VHH6 : il est **jonctionnel**. Il touche 15 résidus du
récepteur IL-6Rα contre 8 d'IL-6. Contre IL-6 seule — ce que le subnet demande —
sa liaison n'est pas démontrée.

### Le test décisif : les deux liants validés contre le vrai filtre

| | nativité | humanité | rejet dur | gène V | V id | CDR3 | L | Cys |
|---|---|---|---|---|---|---|---|---|
| VHH6 (5FUC) | **0,852** | 0,414 | non | IGHV3-66 | 69,1 | AARPSIRCASFSATEYKD | 18 | 4 |
| 61H7 VH (4O9H) | 0,577 | **0,992** | non | **IGHV3-23** | 88,8 | ANRAGWGMGDY | 11 | 2 |

Les deux sont aux antipodes, et le §69 dit lequel compte : le seuil 0,80 de
nativité **n'est pas appliqué**, l'humanité ≥0,75 **l'est**, et les soumissions
acceptées ont une nativité de 0,386 à 0,596. Donc :

- **VHH6 serait refusé** — trop camélidé, humanité 0,414. La solution de la
  nature à cet épitope ne franchit pas le portail du subnet.
- **61H7 VH serait accepté** — humanité 0,992, nativité 0,577, en plein dans la
  fourchette des acceptés, aucun rejet dur.

Et il franchit tout le reste : L=118 (90-150), 2 cystéines, IGHV3-23 avec
88,8 % d'identité V (bien au-dessus du plancher dur de 70 et de notre crible
calibré à 79), absent de l'archive, Jaccard k=5 maximal de 0,416 contre ses
35 358 séquences (règle : < 0,90).

**C'est la première séquence de ce projet qui est à la fois un liant prouvé par
cristallographie de la cible de la semaine, à l'épitope exact de la config, par
un seul domaine, et conforme à tous les portails que nous savons vérifier.**

### Pourquoi cela ne contredit pas le §68

Le §68 a mesuré **« CDR3 du donneur → notre charpente »** : 0,338, le pire bras.
Ce qui est en jeu ici est l'inverse — garder la charpente **du liant** et ne
redessiner que dans celle-ci. C'est la seule configuration que le §68 n'a jamais
testée, et c'est exactement le mode de rupture qu'il avait identifié (la perte au
transfert de charpente) que cette voie supprime.

Ce n'est pas une victoire mesurée. Les contraintes restent : l'époque impose le
gène V et la longueur du CDR3, et 61H7 est IGHV3-23 / L=11 ; l'unicité
s'évalue chez le validateur, pas dans notre archive ; et la liaison ne se mesure
pas sur CPU. C'est une hypothèse à tester en un lot GPU — mais c'est la
première depuis le §53 qui parte d'un liant réel conforme au filtre.

### Deux défauts de notre code, révélés en chemin

1. `scripts/nativeness.py` ne s'exécutait plus : le conteneur a été reconstruit
   et `biopython`, `abnumber`, `anarci` avaient disparu. Le filtre ne lève
   pas — il répond `hard_reject=True, reject_reason='abnumber_failed'` pour
   **tout**, y compris pour des soumissions acceptées. C'est la forme exacte du
   bug du §69 (`internal_data/human` manquant), une deuxième fois. `available()`
   vérifie désormais ces modules par leur nom, `requirements.txt` gagne
   `bittensor-drand`, et `tests/test_known_binders.py` fige les quatre nombres
   ci-dessus pour qu'une dérive devienne un test rouge.
2. La colonne `in_accepted_nativeness_range` ne testait qu'un plancher. Renommée
   `at_least_accepted_nativeness` : être **plus** camélidé que tout accepté n'est
   pas un indice de refus, et le §69 n'a établi que la borne haute à 0,80.

## 74. Le lot de calibration : 61H7 rendu soumettable, et la question jamais posée

Lot prêt : `inputs/kaggle_calibration.csv`, 85 séquences, `notebooks/kaggle_boltzgen.ipynb`
(517 kB, sous la limite d'import de 1 Mo), **2,0 h sur 2× T4**. Les lectures sont
fixées ci-dessous et dans le docstring de `scripts/build_calibration_batch.py`,
**avant** le run.

### Ce que le §73 permet, et ce qui le bloquait

61H7 lie IL-6 à l'épitope exact de la config (8/8) et son paratope n'emprunte
que **25 % au CDR3** — CDR1 31 %, CDR2 38 %, charpente 6 %. Les quatre autres
liants sont à 40-62 % de CDR3. C'est le seul pour lequel « garder la charpente,
redessiner le CDR3 » n'est pas absurde.

Mais il n'est **pas soumettable tel quel**, et j'avais tort de le croire. La
config en vigueur en amont pose trois règles qu'il manque :

- `enforce_vhh_hallmarks: true` (E/Q en IMGT 49, R en 50) — il a G49, L50 ;
- `min_human_v_identity: 82` — il passe à 88,8 ;
- `min_human_j_identity: 86` — il échoue à 85,7.

Ces deux planchers germinaux sont donc **réels**, lus dans `config/config.yaml` :
le §67 avait tort de les traiter comme notre reconstruction. Et la règle des
hallmarks est un hard gate du validateur, pas la composante D du score de
nativité du filtre — les deux existent et ne font pas la même chose.

### Trois mutations, aucune dans le paratope

| | hallmarks | V | J | nativité | humanité | couche 1 |
|---|---|---|---|---|---|---|
| 61H7 brut | non | 88,8 | 85,7 | 0,577 | 0,992 | 2 refus |
| + G49E, L50R | oui | 86,7 | 85,7 | 0,657 | 0,937 | 1 refus (J) |
| + FR4 Q→L seul | non | 88,8 | 92,9 | 0,577 | 0,992 | 1 refus (hallmarks) |
| **+ les trois** | **oui** | **86,7** | **92,9** | **0,657** | **0,937** | **OK** |

```
EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYRMYWVRQPPGKEREWVSAISAGGGSTYYGDSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCANRAGWGMGDYWGQGTLVTVSS
```

FR1, CDR1, CDR2, FR3 et CDR3 sortent **identiques au résidu** ; seuls FR2 et FR4
bougent, de trois positions, dont aucune n'est un contact. IGHV3-23 — le gène de
3 époques sur 12 — et 118 résidus.

### Le lot

| bras | n | ce qu'il demande |
|---|---|---|
| `ref_intact` | 6 | les six liants cristallographiques, recouvrement 8/8 à 0/8 |
| `camelised` | 3 | les trois mutations coûtent-elles la liaison ? |
| `scramble` | 5 | charpente camélisée, **intérieur** du CDR3 mélangé |
| `keep_fw` | 32 | charpente camélisée, CDR3 d'archive, L = 11/13/14/19 |
| `move_cdr3` | 9 | le bras du §68 : CDR3 de 61H7 sur charpente étrangère |
| `near_paratope` | 15 | archive la plus proche du paratope CDR1+CDR2 de 61H7 |
| `archive_far` | 15 | archive ordinaire, le mode à 4-5/11 |

`keep_fw` et `scramble` sont à **100 % conformes** : le filtrage sur le plancher
J était nécessaire, car échanger un CDR3 déplace l'alignement J et faisait
tomber les deux tiers du bras à 78,6 — un bras de séquences insoumettables
n'aurait pas mesuré ce qu'on veut. `ref_intact` et `move_cdr3` ne passent pas,
et c'est voulu : ce sont des ancres et des témoins, pas des soumissions.

**Pas de répliques, et pas par choix** : `boltzgen_score` nomme chaque entrée par
l'empreinte SHA-256 de la séquence, donc une séquence répétée est une seule
mesure payée deux fois. Mon premier jet en avait 3 par ancre et aurait gâché un
tiers de l'heure GPU ; `verify_boltzgen_notebook` l'a attrapé. La précision vient
donc du nombre de séquences **distinctes** par bras, et le plancher de l'échelle
est élargi par de vrais non-liants d'IL-6 : les 15 de `archive_far`, conçus
contre PD-L1, et le VHH de 8AOK qui lie PD-L1 par cristallographie.

### Lectures fixées avant le run

1. **Si 61H7 ne dépasse pas le plancher des 15 non-liants**, l'iiptm calculé à
   ce site de liaison ne lit pas la liaison, et chaque heure GPU depuis le §49
   mesurait autre chose. C'est le résultat le plus important que ce lot puisse
   rendre, et celui que j'attends le moins.
2. **Si `camelised` tombe sous 61H7 brut** de plus que l'écart entre séquences
   d'un même bras, les trois mutations qui le rendent soumettable coûtent la
   liaison, et la séquence soumettable n'est pas le liant. Cette voie se
   fermerait aussi nettement que le §68 a fermé la greffe. C'est le risque que
   ce bras est là pour chiffrer : les positions 49 et 50 sont celles où un VHH
   occupe la place qu'une chaîne légère laisse vide, et caméliser un VH apparié
   n'est pas gratuit dans la littérature non plus.
3. **Si `scramble` égale la charpente camélisée** dans le bruit, la mesure ignore
   le CDR3 et `keep_fw` contre `move_cdr3` ne s'y lit pas.
4. **Le §73 n'est soutenu que si** `keep_fw` dépasse `move_cdr3` de plus que
   l'écart intra-bras **et** atteint le 61H7 de `ref_intact`. Battre
   `archive_far` seul ne suffit pas : c'est la comparaison que le §59 a déjà
   montrée confondue par le choix d'échafaudage.

Ce lot ne mesure pas un classement — il n'y a aucun rival sur `P05231`, l'époque
n'a pas eu lieu. Il mesure si notre instrument lit la liaison, ce qu'aucun des
lots précédents n'a vérifié.

### Deux généralisations au passage

`build_boltzgen_notebook` et `verify_boltzgen_notebook` codaient `8AOK` en dur :
ils auraient embarqué l'ancienne structure et validé un carnet construit pour
l'ancienne cible. Les deux lisent désormais `nova_lab.config`, et le builder
lève si le `.cif` n'est pas en amont plutôt que de laisser
`scripts.boltzgen_score` le nommer dans un yaml sans qu'il existe.

## 75. L'analyse de cible en code, et deux corrections au §73

`scripts/analyze_target.py` + `nova_lab/structures.py` : le §5 du guide amont
(« Study the Target ») exécutable sur n'importe quelle cible, sans qu'un humain
lise des colonnes de mmCIF. Le §73 l'avait fait à la main pour `P05231`/`4O9H`
et ce n'était donc pas rejouable sur la cible suivante. Sortie :
`results/target_P05231_4O9H.md` et `..._binders.csv`.

Cinq questions, dans l'ordre : la numérotation de la config correspond-elle à
la chaîne que la structure contient (en `label_seq_id`, l'espace du `res_index`
de BoltzGen) ; le site de liaison est-il réellement au contact et atteignable
par un domaine unique ; quels liants connus porte cette entrée UniProt et
lesquels franchissent les portails du validateur ; quelle part du paratope
chaque liant emprunte à son CDR3 ; et quelles mutations hors CDR rendraient un
liant soumettable.

Le script redérive tous les chiffres du §73 — décalage label−auth de 2, fenêtre
21..186 couvrant exactement les 152 résidus observés, **8/8 au contact**, H
couvre 8/8 et L seulement 2/8, paratope de 61H7 à 25 % CDR3, nativité 0,577,
humanité 0,992, et les trois mutations G49E / L50R / FR4 Q→L qui le rendent
conforme sans toucher un CDR.

### Deux corrections

**Recouvrements de 4ZS7 et 8YWR.** Le §73 donnait 3/8 et 1/8 ; le script donne
**1/8 et 0/8**. Mon analyse à la main cherchait un décalage en numérotation
auteur ; le script aligne les séquences d'entité et travaille en `label_seq_id`,
qui est l'espace du `res_index`. C'est la bonne convention et elle corrige les
deux. Aucune conclusion ne change : les deux restent des liants d'un autre site.

**Deux faits neufs.** Olokizumab (4CNI, 0/8) a une humanité de **0,84** et
franchit donc le seul seuil appliqué — un liant d'IL-6 conforme au validateur
mais visant un autre épitope, ce qui en fait un témoin de spécificité meilleur
que prévu. Et I9H (8YWQ) est rendu soumettable en **2 mutations** (A49E, L50R) ;
à 1/8 de recouvrement il ne sert pas d'échafaudage, mais c'est une ancre de
calibration conforme.

### Ce que le script refuse de faire

Deux pièges que mon premier jet contenait et que la version finale ferme :

1. **Deviner les chaînes par préfixe de séquence.** Les premiers résidus d'une
   chaîne lourde de Fab sont souvent désordonnés, donc le préfixe observé ne
   correspond pas à la SEQRES. La correspondance entité→chaîne est dans le
   fichier (`_entity_poly.pdbx_strand_id`) et c'est elle qui est lue. Sans ça,
   4O9H sortait « aucun domaine variable lourd » et le récepteur de 5FUC était
   compté comme liant, gonflant son épitope de 15 résidus.
2. **Supposer que `label_seq_id` est comparable entre entrées.** IL-6 est déposée
   en constructions de 186, 183, 171, 168 et 166 résidus : comparer un épitope
   calculé sur 4ZS7 à un site exprimé dans la numérotation de 4O9H donnerait du
   bruit. `best_offset` recale chaque entité sur la référence (accord 0,988 à
   1,000 sur les six liants) avant toute comparaison.

Et un garde-fou : la recherche RCSB n'autorise que quelques requêtes par minute
et rend une limitation de débit sous la forme d'un **résultat vide**, pas d'une
erreur — ce qui se lirait « cette protéine n'a aucun liant connu », la chose la
plus trompeuse que ce script puisse dire. Le relevé est réessayé, et un vide
persistant est rapporté comme un échec réseau, jamais comme une absence.

188 tests passent, dont 10 nouveaux.

## 76. Le run a échoué sur une constante que le changement de cible avait épargnée

Premier lancement du lot de calibration : échec à 6,1 min sur
`missing /kaggle/temp/nova/data/msa_files/Q9NZQ7.a3m`. Le carnet avait pourtant
téléchargé la bonne MSA — `P05231.a3m`, vérifiée par empreinte, 14 s après le
démarrage.

La cause : `scripts/boltzgen_score.py` portait `TARGET = "Q9NZQ7"` **en dur**,
et la ligne juste en dessous était un commentaire affirmant que ces valeurs sont
lues dans `nova_lab.config`. Le commentaire était vrai pour la structure, la
chaîne, le `res_index` et le site de liaison — les quatre que le §67 avait
centralisés — et faux pour la cible, qui avait survécu à la bascule du
20 novembre parce qu'elle **ne sert qu'à construire un chemin** : nommer la MSA
et les fichiers yaml. Rien ne la contredisait en local, aucun test ne la
touchait, et elle n'a rien cassé jusqu'au moment où BoltzGen a rouvert la MSA
sous ce nom — après l'installation de torch et le téléchargement des poids.

C'est la forme de bug la plus chère de ce dépôt : invisible hors GPU, et
détectable en une seconde.

### Ce qui est corrigé

`TARGET = NANOBODY_TARGET[0]`. Et trois garde-fous, parce que corriger la
constante ne corrige pas la classe d'erreur :

1. **Dans le vérificateur local.** `verify_boltzgen_notebook` lit le `TARGET` de
   la copie de `boltzgen_score.py` **embarquée dans la charge utile** — pas la
   locale — et le compare au `MSA_NAME` du carnet. Reproduit sur le carnet qui
   avait échoué, il rend en une seconde : *« le carnet télécharge P05231.a3m
   mais le scoreur embarqué cherchera Q9NZQ7.a3m : régénérer le carnet »*.
2. **Dans le carnet, à l'étape 1.** Après le téléchargement de la MSA et avant
   l'installation de torch, le carnet importe le scoreur et vérifie la
   concordance. Un carnet mal construit s'arrête à 15 s au lieu de 6 min.
3. **Dans les tests.** `test_target_switch.py` vérifie que le scoreur lit sa
   cible dans la config, qu'aucun fichier de la charge utile ne fixe encore
   l'ancienne cible hors du dossier `PREVIOUS_NANOBODY` et des commentaires, et
   que les URL de MSA pointent vers la cible du scoreur. Les deux premiers ont
   été vérifiés en replantant la ligne : ils échouent et nomment
   `scripts/boltzgen_score.py:41`.

Le lot lui-même est inchangé — 85 séquences, les quatre lectures du §74 tiennent.
Carnet régénéré à 518 kB, vérificateur vert.

191 tests passent, dont 3 nouveaux.

## 77. L'instrument ne lit pas la liaison. Les quatre lectures du §74, telles qu'écrites.

85 séquences, 167 min sur 2× T4, zéro échec. Preuve :
`results/scores_nanobody_calibration.csv`. Les seuils étaient fixés au §74 avant
le run ; voici leur confrontation, dans l'ordre.

### Lecture 1 — l'échelle de calibration : ÉCHEC

| liant | recouvrement | iiptm | rang-somme |
|---|---|---|---|
| 61H7 (le liant de cet épitope) | **8/8** | 0,178 | 308 |
| VHH6 | 7/8 | 0,147 | 502 |
| 68F2 | 3/8 | 0,231 | **139** |
| I9H | 2/8 | 0,231 | 239 |
| **8AOK-VHH (lie PD-L1, pas IL-6)** | **0/8** | **0,299** | 164 |
| Olokizumab | 0/8 | 0,143 | 453 |

`Spearman(recouvrement, iiptm) = −0,232` (p = 0,66). Le seuil pré-enregistré
était **≥ +0,30**. Le signe est inversé.

Et la comparaison qui ne souffre pas du petit effectif : **9 des 15 non-liants
dépassent 61H7**, dont la médiane de rang-somme (308) est pire que la médiane
globale (269). Le liant prouvé par cristallographie de l'épitope exact de la
config est au-dessous du bruit de fond.

Deux ancres sortent « juste » — Olokizumab (0/8) et VHH6 (jonctionnel, a besoin
du récepteur) sont mal notés, ce qui est correct. Mais 8AOK-VHH, qui lie **une
autre protéine**, est l'ancre la mieux notée en iiptm, et 68F2 (3/8) est la
meilleure en rang-somme. Deux sur six dans le bon sens : le hasard.

### Lecture 2 — les trois mutations : sans objet

La camélisation **augmente** l'iiptm (0,178 → 0,214). Mais après la lecture 1,
cela ne dit rien sur la liaison.

### Lecture 3 — la mesure voit-elle le CDR3 ? NON, pas son ordre

Charpente identique, CDR3 natif contre cinq permutations de son intérieur
(composition conservée) :

- natif 0,2136 ; permutations 0,2011 à 0,2163, **le natif tombe dans l'étendue** ;
- écart au natif : 0,004, soit **0,65 écart-type intra-bras** et un cinquième
  du σ d'une mesure unique (0,017-0,028, §49).

Le CDR3 d'un liant prouvé est indiscernable d'un mélange de ses propres lettres.
Et le bras `scramble` est **le meilleur des sept** en rang-somme médian (152).

Mais le CDR3 n'est pas ignoré : à charpente fixe, 32 CDR3 différents étalent
l'iiptm sur **0,104 — 6,9 fois l'étendue des permutations**. Ce qui change le
score, c'est donc l'identité de la boucle, pas son ordre — et aucun descripteur
ne l'explique : sur ces 38 séquences à charpente identique, longueur,
hydrophobicité, charge, aromaticité et basicité du CDR3 donnent tous |ρ| ≤ 0,19.

### Lecture 4 — §73 contre §68 : non soutenu

| bras | n | rang-somme médian |
|---|---|---|
| scramble | 5 | **152** |
| camelised | 3 | 208 |
| move_cdr3 (§68) | 9 | 238 |
| keep_fw (§73) | 32 | 267 |
| ref_intact | 6 | 274 |
| archive_far | 15 | 349 |
| near_paratope | 15 | 362 |

`keep_fw` contre `move_cdr3` : p = 0,84. Contre `archive_far` : p = 0,39. Le §73
n'est pas soutenu — et le §68 n'est pas reproduit non plus, les deux bras étant
indiscernables. Les deux mesuraient du bruit.

### Le rang-somme n'est pas dix métriques

C'est **une** mesure en quatre copies :

| métrique | ρ avec le rang-somme |
|---|---|
| design_to_target_iptm | −0,941 |
| design_iiptm | −0,921 |
| interaction_pae | +0,919 |
| min_design_to_target_pae | +0,903 |
| delta_sasa_refolded | −0,294 |
| les cinq autres | ≤ 0,19 |

Les quatre premières sont mutuellement corrélées de 0,873 à 0,961, et
`liability_score` / `liability_num_violations` à 0,980. Sept paires sur 45
dépassent 0,60.

Cela règle une énigme du §49 : on ne peut pas moyenner le bruit sur dix
métriques quand quatre d'entre elles sont la même mesure. Le « rang-somme sur dix
métriques » est, en pratique, un classement sur la confiance d'appariement.

### Ce que cela explique rétroactivement

Quatre impasses deviennent une seule :

- §50 « la liaison ne voyage pas avec le motif » ;
- §68 « la greffe cristallographique est le pire bras » (0,338) ;
- §72 « la conformation ne prédit pas » (|ρ| ≤ 0,041) ;
- §73 « 61H7 est un liant prouvé conforme » — et il note sous les non-liants.

Elles portaient toutes sur la **liaison réelle**. La mesure ne la lit pas. Ce
n'est pas un défaut de notre montage : l'échec de l'ipTM/pAE sur les interfaces
anticorps-antigène est documenté dans la littérature. Vingt heures de GPU ont
servi à poser une question que l'instrument ne pouvait pas entendre.

### Ce qui reste, et c'est exploitable

La mesure est une **boîte noire** : déterministe, biologiquement muette, et
qu'aucun descripteur connu ne prédit. Mais à charpente fixe elle s'étale de
0,104 sur le choix du CDR3, soit 4 à 6 fois le bruit. Il y a donc du signal
réel à exploiter — seulement il ne s'obtient qu'en interrogeant l'oracle, à
2 min la requête.

Cela redéfinit le problème : **optimisation en boîte noire d'un objectif
coûteux**, pas biologie structurale. Les 1 065 séquences déjà notées sont l'actif
du dépôt, et le substitut (ρ = 0,31) est le seul levier qui restait — ce qui
explique aussi pourquoi il plafonnait : il était entraîné à prédire une mesure
qu'on croyait porteuse de sens.

### 77 bis. Correction du σ du §49, et ce que « reproduire le score » veut dire

Deux précisions, parce que le §77 était mal formulé.

**Nous reproduisons le score des validateurs.** Exactement : leur build de
BoltzGen (0.2.0 vendorisée), leur `design.yaml`, leur protocole
`nanobody-anything`, leurs dix métriques, leur structure, leur MSA, leur site de
liaison. 85 séquences sur 85, zéro échec. Les barrières amont sont également
reproduites — la couche 1 est lue dans `config/config.yaml`, et la couche 2
s'exécute avec le code du filtre lui-même (§69, §73). Il n'y a **aucun critère
du validateur que nous ne sachions calculer.** Ce que le §77 a montré est autre
chose : le score est clair, net et calculable, et il ne classe pas les liants
réels au-dessus des non-liants. L'examen est bien défini ; il ne porte pas sur
la biologie.

**Le σ du §49 était trois fois trop grand.** Le bras `scramble` est le témoin
qu'il nous manquait : cinq séquences à charpente identique et à composition de
CDR3 identique, seul l'ordre intérieur changeant. Sa dispersion est un
**majorant** du bruit de mesure.

| bras | n | ce qui varie | écart-type de l'iiptm |
|---|---|---|---|
| `scramble` | 5 | l'ordre du CDR3 seul | **0,0062** |
| `keep_fw` | 32 | l'identité du CDR3 | 0,0218 |
| tout le lot | 85 | tout | 0,0301 |

Le bruit d'une mesure est donc **≤ 0,006**, pas 0,017-0,028. Le §49 l'avait
estimé entre des runs séparés, où l'état du générateur aléatoire de BoltzGen
diffère ; à l'intérieur d'un même run il est bien plus faible. On ne peut pas le
mesurer directement, d'ailleurs : `boltzgen_score` nomme ses entrées par
l'empreinte de la séquence, donc une séquence ne peut pas être notée deux fois
dans un run (§74).

Conséquence, et elle est favorable : le choix du CDR3 déplace le score de
0,0218 d'écart-type, soit **3,5 fois le bruit**, et la variance de `keep_fw` est
12,5 fois celle de `scramble` (Bartlett p = 0,017). Le signal exploitable est
donc réel et bien au-dessus du bruit — contrairement à ce que le §49 laissait
croire, et contrairement à ce que le §77 a d'abord suggéré en parlant de « bruit ».

Ce qui reste hors de portée n'est pas le score : c'est sa **prédiction** sans
GPU. Aucun descripteur de séquence connu ne l'approche (|ρ| ≤ 0,19 sur 38
séquences à charpente identique), et c'est ce qui rend les 1 065 requêtes déjà
payées précieuses.

## 78. Le substitut réentraîné sur ce que le score est vraiment : ρ 0,58

`scripts/train_score_surrogate.py`, `tests/test_score_surrogate.py`. Le §77 a
montré que le rang-somme est une mesure en quatre copies et qu'il ne suit pas la
liaison. Le substitut utile ne prédit donc pas la biologie : il prédit
`design_iiptm`, d'où le rang-somme découle à ρ −0,92. Les substituts précédents
visaient une grandeur qu'on croyait porteuse de sens ; celui-ci vise l'oracle.

### Ce qui a été corrigé avant d'entraîner

**Le bruit, encore.** 248 séquences ont été notées dans plus d'un run : leur
écart-type intra-séquence médian est **0,0121**. C'est le bruit entre runs,
mesuré directement. Le §77 bis annonçait « ≤ 0,006 » d'après le bras `scramble`
et concluait que le §49 était trois fois trop haut : **c'était faux**, cela
confondait deux choses. À l'intérieur d'un run la mesure est bien plus propre
(≤ 0,006) ; entre runs elle vaut 0,012-0,016, donc le §49 (0,017-0,028) était à
peu près juste. Cela justifie rétroactivement la règle de
`scripts/build_arms_batch.py` : tout ce qu'on veut comparer doit tenir dans un
seul run.

**Une seule cible à la fois.** Un score est une fonction du couple
(séquence, cible). 1 490 de nos 1 575 requêtes sont contre `8AOK`, 85 contre
`4O9H`, et **exactement une séquence a été notée contre les deux**. Le transfert
entre cibles n'est donc pas mesurable avec ce que nous avons, et mutualiser les
deux serait l'inventer. Entraîné sur l'ancienne cible, où sont les données.

### La fuite que le test a attrapée

Mon premier regroupement utilisait Jaccard 5-mer ≥ 0,90, la règle de diversité
du validateur. Un test écrit pour vérifier que des jumelles tombent ensemble a
échoué : **deux permutations d'un même CDR3 sont à 0,854**, donc distinctes pour
le subnet et jumelles pour une validation croisée. Le seuil décide du résultat :

| seuil de regroupement | grappes | plus grosse | ρ hors-pli |
|---|---|---|---|
| 0,95 | 1157 | 1 | +0,634 |
| 0,90 (règle du validateur) | 1135 | 4 | +0,627 |
| **0,80 (défaut retenu)** | 659 | 45 | **+0,580** |
| 0,70 | 383 | 45 | +0,477 |

Une charpente partagée avec des CDR3 différents tombe vers 0,81, donc 0,80 est
le seuil qui fait d'une grappe retenue une **charpente jamais interrogée**.
C'est le défaut, et les deux lectures sont rapportées parce qu'elles répondent à
deux questions :

- **charpente jamais vue** : ρ **+0,580** ;
- **boucle choisie dans une charpente déjà interrogée** : ρ +0,627.

### Les modèles

| modèle | ρ (grappes 0,80) | sd plis |
|---|---|---|
| témoin (prédit la moyenne) | −0,109 | — |
| ridge sur descripteurs | +0,420 | 0,083 |
| gradient boosting sur descripteurs | +0,495 | 0,050 |
| ridge sur 3-mers | +0,488 | 0,064 |
| **gradient boosting sur descripteurs + 3-mers** | **+0,580** | 0,052 |

Contre ρ = 0,31 pour le substitut précédent. Le §77 disait qu'aucun descripteur
*isolé* ne dépasse |ρ| 0,19 : c'est toujours vrai, et un modèle qui les combine
atteint 0,58. L'oracle n'est donc pas imprédictible, seulement non résumable en
un descripteur.

### Ce qu'il vaut là où il servira

Un ρ global est flatté par l'amplitude : nos scores vont de 0,11 à 0,79 parce
qu'ils mélangent conceptions de novo et donneurs d'archive, alors que le vivier
d'une époque occupe un dixième de cette plage. Mesuré par bandes, sur la seule
question opérationnelle — tirer 30 candidats, en choisir un :

| bande | n | sd | ρ | hasard | substitut | oracle | part captée |
|---|---|---|---|---|---|---|---|
| tout le jeu | 1157 | 0,188 | 0,580 | 0,383 | 0,621 | 0,708 | **73 %** |
| réaliste 0,15-0,45 | 549 | 0,086 | 0,243 | 0,250 | 0,291 | 0,434 | 22 % |
| **haut du vivier 0,55-0,80** | 304 | 0,055 | 0,303 | 0,618 | **0,694** | 0,755 | **55 %** |

Le haut du vivier est la bande qui décide, puisque c'est là que se trouvent les
donneurs qui gagnent. Le substitut y récupère **55 % de l'écart hasard→oracle**,
soit **+0,076 d'iiptm**, ou **1,38 écart-type** du vivier.

### Comparaison avec du GPU, en unités comparables

L'extrapolation en valeurs extrêmes donnait, pour le choix du CDR3 à charpente
fixe : passer de 30 à 3 000 candidats interrogés vaut +0,030, soit **1,4
écart-type**. Le substitut en vaut **1,38**, gratuitement, et les deux se
composent — pré-classer 10 000 candidats pour n'en interroger que 300 cumule les
deux effets.

Autrement dit : **le substitut vaut à peu près un facteur 100 sur le GPU**, et
il ne coûte rien.

### Ce qui reste ouvert, et ce que ça coûte de le fermer

Le transfert entre cibles est non mesuré, et c'est le trou qui compte : le
substitut est entraîné sur `8AOK` et devra servir sur `4O9H`. Le combler est
bon marché — il suffit d'inclure une vingtaine de séquences déjà notées contre
l'ancienne cible dans le prochain lot de la nouvelle. Vingt séquences, 40 min de
GPU, et on saura si 1 490 requêtes payées sont réutilisables ou perdues.

196 tests passent, dont 5 nouveaux.

## 79. Le lot de transfert : 1,9 h pour savoir si 1 490 requêtes sont un actif

`scripts/build_transfer_batch.py`, `inputs/kaggle_transfer.csv`,
`tests/test_transfer_batch.py`. Carnet régénéré, vérificateur vert.

Le §78 a entraîné le substitut sur l'ancienne cible — 1 490 requêtes contre
`8AOK`, ρ 0,580 sur charpente jamais vue. Il devra servir sur `4O9H`, et
**exactement une séquence du dépôt a été notée contre les deux**. Le transfert
n'est donc pas mesuré, et mutualiser les deux jeux serait inventer ce qui donne
sa valeur au substitut. Ce lot mesure cela, et rien d'autre.

| bras | n | ce qu'il mesure |
|---|---|---|
| `transfer` | 48 | des séquences déjà notées contre `8AOK`, notées contre `4O9H` |
| `repeat` | 8 | des séquences du run de calibration, notées une seconde fois |

**56 séquences, 1,9 h sur 2× T4.** Erreur-type attendue du Spearman : ±0,15.

### Deux choix qui décident du résultat

**La stratification.** Les 48 couvrent `0,143 à 0,760` contre l'ancienne cible,
six par bande sur huit bandes. Ce n'est pas une élégance : le substitut a été
évalué dans une bande d'écart-type 0,055, et un échantillon tiré d'une seule
bande atténuerait par restriction d'étendue un transfert pourtant réel — et nous
dissuaderait de réutiliser les requêtes. Un test vérifie que l'étendue couverte
dépasse 85 % de celle du vivier et qu'aucune des huit bandes n'est vide.

**Les répétitions.** Nous n'avons **jamais** mesuré le bruit entre runs sur la
cible vivante. Sans lui, un ρ faible est ambigu : le transfert échoue, ou les
étiquettes sont bruitées. Huit séquences du run de calibration, étalées sur sa
plage, tranchent.

Les 48 sont tirées uniquement parmi les 248 séquences notées **plusieurs fois**,
donc l'abscisse de la comparaison porte moins de bruit que l'ordonnée — le bon
sens pour une corrélation sur laquelle on va agir.

### Lectures fixées avant le run

- **ρ ≥ 0,6** : les requêtes se transfèrent. Le substitut est utilisable sur la
  cible vivante dès maintenant, et le GPU dépensé sur l'ancienne n'est pas perdu.
- **0,3 ≤ ρ < 0,6** : transfert partiel. Utilisable comme a priori, à
  réentraîner au fur et à mesure.
- **ρ < 0,3** : le score est propre à la cible. Le substitut doit être
  reconstruit par cible, et les 1 490 requêtes sont de l'histoire, pas un actif.

Et un confondant à écarter sur les mêmes données, sans une séquence de plus :
`design_ptm` est la confiance de repliement du nanobody seul et n'implique pas
la cible. Si la corrélation survit au retrait de `design_ptm`, elle porte bien
sur l'interface ; sinon nous ne redécouvrons que « certaines séquences se
replient bien partout », ce qui n'aide pas à choisir un gagnant pour une époque.

La lecture est écrite dans le module (`analyse`, `verdict`) **avant** que les
données existent, et deux tests la font tourner sur un aller-retour synthétique
— une corrélation plantée doit être retrouvée, une absence doit être rapportée
comme telle. Le verdict ne pourra pas être choisi après coup.

206 tests passent, dont 10 nouveaux.

## 80. Zéro transfert entre cibles. Les 1 490 requêtes sont de l'histoire.

56 séquences, 121 min sur 2× T4, zéro échec. Preuve :
`results/scores_nanobody_transfer.csv`. Lecture exécutée par
`scripts.build_transfer_batch --analyse`, écrite au §79 avant que les données
existent.

```
48 séquences comparables
  Spearman(ancienne cible, cible vivante) = +0,003 ± 0,15   p = 0,98
  avec la confiance de repliement retirée  = +0,066
  bruit de la cible vivante (8 répétitions) : 0,0133, biais de run −0,0008
```

Le seuil pré-enregistré était ρ ≥ 0,6 pour « les requêtes se transfèrent » et
ρ < 0,3 pour « le score est propre à la cible ». **Le verdict est le troisième**,
et il n'est pas serré : avec ±0,15 sur 48 points, même un transfert partiel à
0,3 est exclu à deux erreurs-types.

Ce n'est pas un artefact de bruit, et le bras de répétitions le prouve : le
bruit de la cible vivante est **0,0133**, très petit devant l'étendue de 0,143
à 0,760 couverte sur l'abscisse. Le confondant nommé au §79 est écarté aussi —
retirer `design_ptm` ne fait pas apparaître de corrélation (+0,066).

Deux acquis annexes, tous deux utiles :

- **le bruit de la cible vivante est mesuré** pour la première fois : 0,0133,
  cohérent avec les 0,0121 de l'ancienne cible, donc la valeur ne dépend pas de
  la cible même si le score en dépend entièrement ;
- **les runs ne sont pas biaisés**, seulement bruités : le biais entre deux runs
  est de −0,0008 sur huit répétitions. Ce qui différait entre runs était du
  bruit, pas un décalage.

### Ce que cela coûte

Le substitut du §78 atteignait ρ 0,580 sur 1 157 séquences de l'ancienne cible.
Il ne vaut rien ici. Sur la cible vivante, avec les 132 séquences dont nous
disposons désormais (85 du §77 + 56 de ce lot), le même modèle rend :

| | ρ (grappes 0,80) | part captée sur 30 candidats |
|---|---|---|
| ancienne cible, 1 157 séquences | +0,580 | 73 % |
| **cible vivante, 132 séquences** | **+0,235 ± 0,14** | 20 % |

À peine distinguable de zéro. Et comme la cible change chaque semaine, **le
compteur repart de zéro à chaque fois**.

### La correction que ce résultat impose au §78

Le §78 appelait les 1 490 requêtes « l'actif du dépôt » et concluait que le
substitut valait « à peu près un facteur 100 sur le GPU ». Les deux affirmations
sont fausses telles qu'écrites. Ce qui reste vrai :

- la **méthode** est validée : un substitut atteint ρ 0,580 à partir d'environ
  un millier de requêtes **sur une même cible** ;
- la **valeur** mesurée (73 % de l'écart hasard→oracle, 1,38 écart-type) est
  réelle, mais elle n'est accessible qu'après avoir payé ces requêtes sur la
  cible en cours.

Autrement dit, le substitut n'est pas un substitut *au* GPU : c'est ce qu'on
obtient *avec* du GPU, et il faut le racheter à chaque cible. La question que
tu posais — « est-ce qu'un GPU permanent nous ferait avancer » — a donc une
réponse plus nette qu'au §78, et plus chère.

### Un défaut de code révélé au passage

`scripts/train_score_surrogate.py` reconnaissait les fichiers de la cible
vivante par un seul nom de fichier codé en dur, celui du run de calibration.
Le fichier de ce lot aurait donc été classé sous l'**ancienne** cible — et au vu
du ρ ci-dessus, ce n'est pas une dégradation du modèle, c'est l'entraîner sur
une fonction sans rapport. Remplacé par un registre explicite
(`LIVE_TARGET_FILES`) que chaque nouveau run doit alimenter.

### 80 bis. Courbe d'apprentissage : combien de requêtes par cible

Le §80 ferme le transfert, donc la seule question qui reste chiffrable est :
combien de requêtes faut-il sur une cible donnée. Mesuré en sous-échantillonnant
l'ancienne cible **par grappes** (pas par lignes, sinon on garde les jumelles),
trois tirages par point :

| requêtes | heures GPU à 2,2 min | ρ (grappes 0,80) |
|---|---|---|
| 100 | 3,7 h | +0,272 ± 0,224 |
| 200 | 7,3 h | **+0,422 ± 0,104** |
| 400 | 14,7 h | +0,394 ± 0,018 |
| 700 | 25,7 h | +0,474 ± 0,063 |
| 1 000 | 36,7 h | +0,540 ± 0,011 |
| 1 157 | 42,4 h | +0,577 |

C'est bien meilleur que ce que le §80 laissait craindre. **200 requêtes — 7,3 h —
donnent déjà ρ 0,42**, soit les trois quarts de ce que 1 157 requêtes
atteignent. Au-delà, le rendement est faible : tripler encore (700) ne gagne que
0,05, et il faut 1 000 requêtes pour franchir 0,54.

Nous avons déjà **132 requêtes** sur la cible vivante. Environ **70 de plus,
2,5 h de GPU**, nous placent au palier des 200.

Et une précision de structure qui change le plan : le substitut prédit le score
pour un couple (séquence, cible). La cible change chaque **semaine**, les
contraintes d'époque — gène V, longueur du CDR3 — changent chaque **époque**
(72 min). Le substitut est donc à racheter par semaine, pas par époque, et il
sert à toutes les époques de la semaine.

## 81. Le lot d'apprentissage : 70 requêtes pour atteindre le palier

`scripts/build_learning_batch.py`, `inputs/kaggle_learning.csv`,
`tests/test_learning_batch.py`. Carnet régénéré (519 kB), vérificateur vert,
run `learning`. **70 séquences, 2,6 h sur 2× T4.**

Le §80 bis a montré que 200 requêtes sur une même cible donnent ρ 0,42, les
trois quarts de ce que 1 157 atteignent. La cible vivante en a 141 et plafonne à
ρ 0,235 ± 0,14. Ces 70 l'amènent à 202.

### Pourquoi ces 70 et pas 70 au hasard

Le jeu de la cible vivante est déséquilibré : **39 de ses 132 séquences — 30 % —
partagent la seule charpente camélisée de 61H7**, parce que le lot de
calibration avait été construit pour répondre à une question sur cette
charpente (§74). Un substitut entraîné là-dessus extrapole partout ailleurs. Un
tirage au hasard, lui, reproduirait la forme de l'archive au lieu de corriger la
nôtre.

| bras | n | rôle | prédiction |
|---|---|---|---|
| `cover` | 45 | couverture par échantillonnage du point le plus éloigné | 0,149 à 0,315 (médiane 0,194) |
| `exploit` | 25 | les mieux notées par le substitut actuel | 0,294 à 0,320 (médiane 0,302) |

`cover` est choisi par max-min glouton dans un espace de 40 composantes
principales, amorcé par les 132 déjà interrogées : chaque pick est la séquence
du vivier la plus éloignée de tout ce qui est déjà tenu. Un test mesure que ces
picks sont **plus de 1,5 fois plus loin** du jeu connu qu'un tirage au hasard —
la couverture est vérifiée, pas supposée.

`exploit` sert deux fois : il peuple le haut de la plage, où se prennent toutes
les décisions, et c'est un **test prospectif** du substitut, le seul genre
qu'une validation croisée ne peut pas simuler.

Vivier : les 17 634 séquences d'archive jamais interrogées contre `4O9H`. C'est
le bon vivier de toute façon — c'est ce dans quoi notre générateur puise, donc
le substitut doit être entraîné sur la distribution qu'on lui demandera de
classer. 44 des 70 passent la couche 1 ; les 26 autres sont conservées, parce
que ce lot mesure l'oracle et qu'une séquence refusée par le validateur reste
une donnée valide. Un test vérifie d'ailleurs qu'elles n'ont **pas** été
filtrées, sans quoi le lot ne couvrirait plus l'espace voulu.

### Lectures fixées avant le run

1. **Réentraîné sur 202 requêtes, ρ doit approcher le palier** que l'ancienne
   cible montrait à 200 (+0,42 ± 0,10). S'il reste vers 0,235, la cible vivante
   est plus dure que l'ancienne et **la courbe d'apprentissage ne se transfère
   pas non plus** — ce qui vaudrait d'être su avant d'acheter du GPU.
2. **Les 25 `exploit` doivent noter au-dessus de la médiane des 45 `cover`.** À
   ρ 0,235 l'écart attendu est petit mais mesurable. Si le bras `exploit` ne
   bat pas `cover`, le ρ validé croisé ne se traduit pas en choix, et tout le
   plan du §80 bis est à revoir.
3. **Surveiller l'étendue.** Si les 70 reviennent dans 0,12-0,39 comme les 132
   existantes, la plage de score de la cible vivante est vraiment aussi étroite,
   et les 0,0133 de bruit en sont alors un cinquième utile et non un trentième.

### Un détail de coût

Le premier jet des tests mettait 6 min 49 : la passe max-min est quadratique.
La taille du vivier est désormais un paramètre, les tests travaillent sur un
dixième et passent en 18 s. Un fichier de tests lent est une taxe sur toutes
les exécutions futures.

213 tests passent, dont 7 nouveaux.

## 82. « Il suffit d'évaluer les gagnants précédents » : testé, et non

Question légitime après le §77 : si le problème est une optimisation en boîte
noire, les gagnants des époques passées devraient dire quoi viser. Elle a deux
formes et toutes deux sont maintenant mesurées.

### Forme faible : réutiliser les séquences qui ont bien noté

Le lot du §80 contenait précisément cela — des séquences d'archive, dont les
mieux notées contre l'ancienne cible, repassées contre la vivante :

| | médiane sur la cible vivante | rang-somme médian |
|---|---|---|
| les 8 **meilleures** de l'ancienne cible (0,677 à 0,760) | 0,1753 | 229 |
| tout le lot | 0,1865 | 194 |
| les 8 **pires** de l'ancienne cible | 0,1702 | — |

Les meilleures notent **en dessous** de la médiane du lot, et à peine au-dessus
des pires (0,1753 contre 0,1702). Savoir qu'une séquence a été excellente la
semaine dernière ne dit rien de cette semaine.

### Forme forte : un « style » transférable

C'est la version qui n'avait jamais été testée. Peut-être que les séquences qui
notent haut partagent des traits — longueur de CDR3, charge, hydrophobicité,
aromaticité — valables sur n'importe quelle cible. Corrélation de 21 traits avec
le score, sur chaque cible :

- **5 traits sur 21 pointent dans le même sens**, là où le hasard en donnerait
  ~10. Il n'y a donc pas d'accord, même pas partiel ;
- les deux profils de traits se ressemblent à **ρ = −0,021** (p = 0,93) ;
- le seul trait nettement prédictif sur l'ancienne cible — la charge du CDR3,
  ρ +0,287 — **change de signe** sur la vivante (−0,056).

Et un modèle entraîné sur les 1 157 séquences de l'ancienne cible, appliqué aux
132 de la vivante, rend **ρ = −0,075** (p = 0,39). Ce test est distinct du §80 :
le §80 comparait les scores des mêmes séquences, celui-ci demande si un modèle
emporte quoi que ce soit. Non.

### Pourquoi c'était à attendre

Le score est une fonction du **couple** (séquence, cible) : c'est une confiance
d'appariement d'interface. Rien dans une séquence seule ne détermine comment
elle se pose sur une protéine qu'elle n'a jamais vue. La surprise n'est pas que
cela ne transfère pas — c'est que nous l'ayons espéré.

### Ce qui transfère vraiment des gagnants

Non pas leur score, mais les **barrières** — et c'est déjà exploité. Les
barrières sont les mêmes d'une semaine à l'autre parce qu'elles ne regardent que
la séquence, pas la cible. Le §69 en a tiré que la nativité des soumissions
acceptées se situe entre 0,386 et 0,596 et que le seuil de 0,80 n'est pas
appliqué, et le §70 en a tiré le crible calibré à V ≥ 79. Ces deux résultats
viennent exactement de « regarder ce qui passe chez les autres », et ils
tiennent.

La ligne de démarcation est donc nette : **tout ce qui ne regarde que la
séquence s'apprend des autres mineurs une fois pour toutes ; tout ce qui regarde
l'interface se rachète à chaque cible.**
