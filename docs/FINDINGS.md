# Ce que nous avons mesuré sur SN68 — synthèse

## Les données brutes derrière chaque chiffre

| Fichier | Contenu |
|---|---|
| `results/psh_by_generator.csv` | 540 candidats, 3 stratégies × 3 longueurs de CDR3, PSH mesurée (§3) |
| `results/place_vs_psh.csv` | 21 époques : place dans le champ complet, PSH, admissibilité (§3bis) |
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

## 8. Ce qui reste ouvert

**Bloquant pour un nouveau run :**

1. Confirmer la stratégie de génération sur 540 candidats et trois longueurs de
   CDR3 (en cours). Si `pssm` dégrade bien la PSH, il faut câbler
   `pssm_denovo`.
2. Un **preflight calibré** : dérouler chaque verrou dans l'ordre du validateur
   et le valider contre la vérité terrain — les 36 acceptées de l'époque 26174
   doivent passer, nos 11 rejetées doivent échouer. Sans ce calibrage on
   refait l'erreur du §6.
3. Mesurer le **temps réel de chaque étage** dans une époque : mon dernier
   essai à blanc est tombé sur `too_late`.

**Important :**

4. Persister les scores BoltzGen — sans eux, impossible de distinguer
   « rejeté par un verrou » de « noté et mal classé ».
5. Reprise sur HTTP 429 (3 époques et 5 heures perdues en production).
6. Vérifier le dépôt GitHub avant d'engager le pointeur (époque 26172 perdue
   sur un pointeur orphelin).

**Ce que le générateur ne sait pas faire :** nos candidats sont à 98,1 de PSH
médiane contre 90,0 pour les acceptées. Le criblage exploite la queue de
distribution ; il ne rend pas le générateur meilleur.
