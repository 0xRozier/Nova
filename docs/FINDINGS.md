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
