# Observations de production — SN68 miner (2026-11-06)

## Contexte

- Machine : 2x RTX 4090 (Lium, 0.90$/h)
- UID 22 (hotkey sub68)
- Scorer : BoltzGen (GPU) depuis epoch 26174
- Générateur : greffe CDR3 + mutations légères CDR1/CDR2/CDR3

## Résultats des premières epochs

### Epoch 26174 (BoltzGen, V=IGHV3-23, CDR3=8)

| Métrique | Notre score | Percentile archive | Commentaire |
|---|---|---|---|
| design_iiptm | 0.406 | 30% | faible |
| design_ptm | 0.789 | 39% | sous médiane |
| design_to_target_iptm | 0.386 | 31% | faible |
| min_design_to_target_pae | 5.43 | 35% | |
| interaction_pae | 14.58 | 33% | loin du top (7.5) |
| plip_hbonds | 6 | 36% | |
| plip_saltbridge | 4 | 67% | correct |
| delta_sasa | 687.8 | 14% | **critique** |
| liability_score | 190 | 55% | correct |
| liability_violations | 22 | 41% | |

**Verdict** : non compétitif. Faible sur quasiment toutes les métriques.

### Epoch 26175 (BoltzGen, V=IGHV3-74, CDR3=11)

| Métrique | Notre score | Percentile archive | Commentaire |
|---|---|---|---|
| design_iiptm | 0.682 | 61% | correct |
| design_ptm | 0.824 | 67% | bon |
| design_to_target_iptm | 0.639 | 66% | bon |
| min_design_to_target_pae | 2.44 | 57% | correct |
| interaction_pae | 9.79 | 57% | correct |
| plip_hbonds | 10 | 84% | **excellent** |
| plip_saltbridge | 5 | 86% | **excellent** |
| delta_sasa | 718.9 | 20% | **critique** |
| liability_score | 220 | 16% | **critique** |
| liability_violations | 25 | 11% | **critique** |

**Verdict** : moyen. Bonnes interactions (hbonds, saltbridge) mais sabordé par
liability et delta_sasa. Le rank-sum souffre car un seul percentile < 20%
détruit le classement.

## Diagnostic

### Faiblesses identifiées

1. **delta_sasa systématiquement bas (14-20%)** : les CDR3 greffés ne créent
   pas assez de surface de contact avec la cible. L'archive contient des CDR3
   scorés sur d'autres frameworks — le transfert de framework perd le contexte
   structural.

2. **Liability catastrophique en 26175 (11-16%)** : trop de mutations de surface
   créent des sites de déamidation (NG, NS, ND) et d'oxydation (M exposé).
   Le générateur mute sans pénaliser les liabilities.

3. **Shortlist trop court** : 3000 générés → 1400-2100 valides → 8 shortlist →
   5-7 après développabilité. On ne score que 5-7 candidats sur GPU. Les gagnants
   en scorent probablement des dizaines à des centaines.

4. **CDR3 trop courts pour 26174** : CDR3=8 est très court, limitant les
   interactions possibles (d'où le delta_sasa bas).

### Axes d'amélioration

#### A. Volume de scoring (rapide)
- Augmenter `--shortlist` de 8 à 40-80
- 40 candidats × 0.5 min = ~20 min de GPU, dans le budget
- Plus de chance de trouver une perle rare

#### B. Filtre liability en amont (moyen)
- Compter les motifs NG/NS/ND/DG/DS + M exposé avant de scorer
- Éliminer les candidats à liability élevée avant le shortlist
- Le validator a un seuil à liability_score <= seuil : se caler dessus

#### C. Optimisation delta_sasa (difficile)
- Sélectionner les CDR3 dont l'archive montre un delta_sasa élevé
  **sur le même framework**, pas seulement en absolu
- Ou muter vers des résidus hydrophobes (W, F, Y, L, I) en position
  de contact connue

#### D. Augmenter la diversité de génération
- Plus de donors (40 → 100)
- Plus de CDR3 pool (60 → 200)
- Mutations structurées (pas random) : cibler les positions de contact

## Métriques de référence (archive, 12664 séquences)

| Stat | rank_sum | design_iiptm | interaction_pae | delta_sasa |
|---|---|---|---|---|
| best | 3000 | 0.827 | 7.49 | 1418 |
| p10 | 10602 | 0.782 | 7.91 | 862 |
| median | 33703 | 0.632 | 10.47 | 810 |
| p90 | 59268 | 0.346 | 15.08 | 743 |

Pour gagner une epoch, il faut un rank_sum dans le top-5 du pool de
l'epoch (~5-15 soumissions). Estimé : rank_sum < 5000 en pool global,
soit percentile > 85% sur la plupart des métriques.


---

# Analyse du 2026-11-06

## Ce que la production a réellement donné

**La mécanique fonctionne** : les engagements arrivent sur la chaîne, les
fichiers sont lus sur GitHub, le déchiffrement timelock réussit, une séquence
par époque. Les trois pièces jamais testées avant la mise en service marchent.

**Mais aucune soumission n'a été notée.** Les scores du tableau ci-dessus sont
ceux que notre mineur a calculés localement, pas ceux du validateur. Vérifié en
interrogeant directement l'API de partage de scores :

| Époque | Soumissions valides | UID 22 noté ? |
|---|---|---|
| 26173 | 19 sur 37 | **non** |
| 26174 | 36 sur 38 | **non** |
| 26175 | pas encore publié | — |

En 26173, 18 soumissions sur 37 ont été rejetées, dont celles de QTTECH et de
saon-miner : rejet massif, propre à cette époque. En 26174, seulement 2 rejets
sur 38 — nous et un UID de QTTECH.

## Pourquoi : un verrou jamais implémenté

Tout ce que notre vérificateur reproduit passe : règles de séquence, gène V,
identités germinales, VHH, cystéines canoniques, longueur de CDR3, unicité
(contre une archive rafraîchie à 12842 séquences), similarité (0,847 contre un
seuil à 0,95), TNP (tous verts) et Aggrescan3D.

Il restait un verrou dans `nanobody_validity.py` que je n'avais pas vu :

```
psh_strict_rejection_threshold: 89   # config.yaml, ajoute le 2026-10-20
```

Une soumission dont la PSH dépasse 89 est **refusée d'office**, soit 1,7 fois
plus strict que le seuil rouge de TNP (155,47). Nos trois séquences :

| Époque | PSH | Verdict de ce verrou |
|---|---|---|
| 26173 | 88,64 | passe |
| 26174 | **91,99** | **refusé** |
| 26175 | **103,80** | **refusé** |

Les deux époques où nous avons été rejetés sans rejet massif sont les deux dont
la PSH dépasse 89.

## La réserve qui compte

**Notre PSH n'est pas calibrée sur la leur.** 62 % des soumissions réelles que
les validateurs ont *acceptées* dépassent 89 sur notre calcul — impossible si le
verrou s'appliquait à cette valeur. Les cinq échelles d'hydrophobicité de TNP ont
été essayées : Kyte-Doolittle (défaut) est la plus basse à 62 %, les autres à
98-100 %. L'adaptateur du validateur lit pourtant le même champ
(`surface_hydrophobic_patches = seq_data["PSH"]`).

Donc le verrou existe et son seuil est connu, mais nous ne savons pas encore
reproduire la valeur qu'il teste. Il est donc implémenté comme un **classement**
et non comme un rejet : on garde les candidats les moins hydrophobes, sans
jamais soumettre le vide.

## Ce qui est vrai indépendamment de l'échelle

Nos candidats sont **systématiquement plus hydrophobes que ceux qui passent** :
PSH médiane 98,1 contre 90,0 pour les soumissions acceptées, même pipeline,
même échelle. C'est le signal noté « à surveiller » dans le README, et il a
mordu. La correction ne dépend pas de la calibration : générer et sélectionner
vers le bas de la distribution acceptée.

Si notre PSH était la bonne, appliquer le verrou laisserait 12 candidats sur 84
et 16 époques sur 21 sans aucun candidat — le taux mesuré tomberait de 29 % à
5 %. Ce chiffre n'est pas une mesure : c'est ce qui arriverait si l'échelle
correspondait.

## Sur les axes proposés plus haut

- **Liability (axe B)** : juste, et moins cher que prévu. `liability_score` et
  `liability_num_violations` sont **exactement déterministes** — écart 0,0000 sur
  84 séquences notées deux fois indépendamment, là où les 8 autres métriques
  varient. Ce sont 2 métriques sur 10, donc 20 % de la somme des rangs,
  optimisables **sans aucun GPU**.
- **Volume de scoring (axe A)** : à nuancer. `ranking.py` départage par somme de
  rangs → catégories → **bloc d'engagement** → **horodatage GitHub** → UID. Le
  plus tôt gagne, donc noter 40 candidats au lieu de 8 retarde la soumission et
  perd les égalités. Ce n'est pas gratuit.
- **delta_sasa (axe C)** : cohérent avec nos mesures hors ligne.
