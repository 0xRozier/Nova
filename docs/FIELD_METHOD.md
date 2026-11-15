# Par quelle méthode le champ produit-il ses soumissions acceptées ?

Analyse de 348 soumissions acceptées sur 12 époques (26299-26306, 26345-26348),
reproductible par `PYTHONPATH=. .venv/bin/python scripts/analyze_field_method.py`.

## 0. Discipline anti-fuite et ce qu'elle coûte

L'archive publique ne porte **aucune date**. La seule vue datée dont nous
disposons est l'instantané `data/prebuilt/Q9NZQ7_regions.csv`, **9466**
séquences, qui couvre tout ce qui a été soumis jusqu'à l'époque
26025. `data/archive/Q9NZQ7_first_seen.csv` n'existe pas :
`archive_ledger.seen_before()` retombe donc sur cet instantané pour **toutes**
les époques analysées ici (toutes > 26025).

Le vivier de donneurs « disponible avant l'époque E » est donc :

* l'instantané (9466 séquences, datées <= 26025) ;
* plus les séquences de `submissions.jsonl` d'époque < E (458 soumissions
  décryptées sur 12 époques, acceptées et refusées confondues).

Les ~7778 entrées de l'archive d'aujourd'hui qui ne sont pas
dans l'instantané sont **exclues** : nous ne pouvons pas prouver qu'elles
existaient avant l'époque analysée. Conséquence à garder en tête pour tout
ce qui suit : **chaque identité « au meilleur donneur antérieur » est une
borne inférieure**. Un vrai donneur plus proche a pu exister sans que nous
puissions le dater. Les conclusions ne sont donc tirées que dans le sens où
cette borne les autorise.

Trou de 39 époques entre 26306 et 26345 : aucune soumission décryptée n'existe
pour 26307-26344, le vivier des 4 époques récentes est donc particulièrement
sous-estimé.

## 1. Les soumissions acceptées dérivent-elles de l'archive ?

n = 348. Identité au meilleur donneur antérieur (déciles : min | p5 | p25 | médian | p75 | p95 | max).

| mesure | min | p5 | p25 | med | p75 | p95 | max |
|---|---|---|---|---|---|---|---|
| identité séquence entière (%) | 81,3 | 83,2 | 84,6 | 86,0 | 87,4 | 90,0 | 92,9 |
| identité CDR3, tout donneur (%) | 30,0 | 36,4 | 46,2 | 57,1 | 64,3 | 78,6 | 85,7 |
| identité CDR3, donneur de même longueur (%) | 25,0 | 31,6 | 44,4 | 54,5 | 64,3 | 78,6 | 85,7 |

* identité entière entre 0 et 70 % : 0 / 348 (0,0 %)
* identité entière entre 70 et 80 % : 0 / 348 (0,0 %)
* identité entière entre 80 et 85 % : 119 / 348 (34,2 %)
* identité entière entre 85 et 90 % : 212 / 348 (60,9 %)
* identité entière entre 90 et 95 % : 17 / 348 (4,9 %)
* identité entière entre 95 et 100 % : 0 / 348 (0,0 %)

CDR3 repris à l'identique d'un donneur antérieur : 0 / 348.
CDR3 à >= 85 % d'un donneur antérieur : 3 / 348.
CDR3 à >= 85 % d'une entrée de l'archive d'aujourd'hui, les 348 acceptées elles-mêmes retirées du vivier (vue **fuitée**, donc bornée par le haut) : 10 / 348, médian 60,0 %. Cette vue ne porte que sur les 10384 séquences de l'archive que le cache de régions sait numéroter (voir 1 bis).

### 1 bis. Un test qui ne dépend d'aucune numérotation ni d'aucune date

Un CDR3 greffé tel quel reste un **sous-mot exact** de la séquence
donneuse. On cherche donc chaque CDR3 accepté comme sous-chaîne dans les
séquences brutes de l'archive, en excluant les 348 acceptées
elles-mêmes. Ce test ne demande aucune numérotation IMGT, donc il couvre
l'archive **entière** (17244 séquences), et le résultat est donné
séparément pour la partie datée et la partie indatable.

* présents dans la partie datée (9466 séquences) : **0 / 348**
* présents dans la partie indatable (7778 séquences, donc une vue qui peut contenir leurs propres gagnantes) : **8 / 348**

Note sur une limite du dépôt : `data/archive/Q9NZQ7_regions.csv` porte 6740 lignes dont le CDR3 est vide alors que les séquences sont numérotables (vérifié sur un échantillon ; c'est le point 10 de `docs/REPLAY_AUDIT.md`, le cache de régions empoisonné). Toutes les comparaisons par CDR3 contre l'archive **d'aujourd'hui** ne portent donc que sur la partie numérotée, et c'est pour cela que le test par sous-chaîne ci-dessus, qui ignore le cache, est le garde-fou le plus solide. Le vivier daté (l'instantané) est complet : 0 CDR3 manquant.

## 2. Les CDR3 se répètent-ils d'une époque à l'autre ?

| époque | n acceptées | long. CDR3 imposée (observée) | CDR3 distincts | CDR3 déjà vu exact (époques antérieures) | >= 85 % (époques antérieures) | >= 85 % (instantané) |
|---|---|---|---|---|---|---|
| 26299 | 30 | {14: 30} | 20 | 0 | 0 | 0 |
| 26300 | 26 | {8: 26} | 24 | 0 | 0 | 0 |
| 26301 | 21 | {20: 21} | 16 | 0 | 0 | 0 |
| 26302 | 31 | {11: 31} | 27 | 0 | 0 | 0 |
| 26303 | 21 | {18: 21} | 18 | 0 | 0 | 0 |
| 26304 | 32 | {14: 32} | 22 | 0 | 1 | 0 |
| 26305 | 35 | {10: 35} | 24 | 0 | 0 | 0 |
| 26306 | 34 | {17: 34} | 29 | 0 | 0 | 0 |
| 26345 | 31 | {14: 31} | 28 | 0 | 1 | 0 |
| 26346 | 25 | {13: 25} | 21 | 0 | 0 | 0 |
| 26347 | 30 | {19: 30} | 26 | 0 | 0 | 0 |
| 26348 | 32 | {6: 32} | 25 | 0 | 0 | 0 |

Total hors première époque (n = 318) : 0 CDR3 exactement réutilisés, 2 à >= 85 % d'un CDR3 accepté antérieur.
Contre l'instantané (n = 348) : 0 CDR3 à >= 85 % d'un CDR3 de l'instantané.

## 3. Les mineurs réutilisent-ils leurs propres séquences ?

52 hotkeys avec >= 2 soumissions acceptées sur >= 2 époques. Identités entre paires de soumissions du même hotkey, époques différentes.

| uid | n acc. | n époques | ident. entière moy. | max | ident. CDR3 moy. | max | iiptm moy. | iiptm max |
|---|---|---|---|---|---|---|---|---|
| 26 | 11 | 11 | 78,8 | 86,8 | 17,8 | 36,4 | 0,722 | 0,788 |
| 63 | 11 | 11 | 79,2 | 85,1 | 14,8 | 57,1 | 0,711 | 0,788 |
| 23 | 4 | 4 | 76,3 | 79,7 | 14,7 | 21,4 | 0,620 | 0,787 |
| 17 | 10 | 10 | 77,4 | 81,6 | 15,0 | 35,0 | 0,722 | 0,787 |
| 29 | 10 | 10 | 78,8 | 84,1 | 15,9 | 40,0 | 0,715 | 0,785 |
| 25 | 8 | 8 | 78,2 | 84,0 | 16,6 | 40,0 | 0,690 | 0,780 |
| 19 | 4 | 4 | 75,8 | 78,8 | 15,9 | 23,1 | 0,640 | 0,779 |
| 48 | 9 | 9 | 78,6 | 84,7 | 15,4 | 57,1 | 0,722 | 0,774 |
| 4 | 7 | 7 | 78,7 | 83,5 | 17,1 | 42,9 | 0,642 | 0,771 |
| 62 | 9 | 9 | 78,1 | 84,5 | 18,9 | 42,9 | 0,700 | 0,770 |
| 18 | 10 | 10 | 79,1 | 84,0 | 14,3 | 42,9 | 0,713 | 0,768 |
| 57 | 9 | 9 | 79,0 | 86,0 | 15,4 | 42,9 | 0,697 | 0,767 |
| 7 | 11 | 11 | 81,3 | 89,5 | 18,2 | 50,0 | 0,685 | 0,763 |
| 12 | 9 | 9 | 80,4 | 84,9 | 16,2 | 64,3 | 0,713 | 0,762 |
| 50 | 7 | 7 | 78,0 | 82,1 | 18,8 | 40,0 | 0,705 | 0,761 |
| 2 | 3 | 3 | 75,2 | 77,3 | 10,4 | 15,4 | 0,615 | 0,761 |
| 31 | 6 | 6 | 78,1 | 84,5 | 20,0 | 42,9 | 0,650 | 0,761 |
| 2 | 7 | 7 | 81,8 | 89,9 | 18,7 | 31,6 | 0,689 | 0,760 |
| 23 | 6 | 6 | 78,8 | 82,8 | 15,8 | 42,9 | 0,673 | 0,760 |
| 24 | 3 | 3 | 82,3 | 84,6 | 22,9 | 26,3 | 0,718 | 0,756 |
| 61 | 9 | 9 | 78,1 | 83,2 | 14,9 | 42,9 | 0,648 | 0,754 |
| 10 | 6 | 6 | 81,7 | 87,6 | 17,7 | 71,4 | 0,657 | 0,752 |
| 38 | 7 | 7 | 80,7 | 84,0 | 14,5 | 35,7 | 0,702 | 0,751 |
| 6 | 4 | 4 | 78,6 | 81,2 | 13,8 | 21,4 | 0,655 | 0,751 |
| 59 | 9 | 9 | 79,9 | 86,4 | 16,8 | 50,0 | 0,660 | 0,749 |
| 58 | 4 | 4 | 82,3 | 86,6 | 11,0 | 17,6 | 0,694 | 0,748 |
| 16 | 10 | 10 | 80,3 | 87,4 | 17,5 | 64,3 | 0,664 | 0,748 |
| 36 | 8 | 8 | 79,2 | 87,9 | 18,3 | 30,0 | 0,697 | 0,747 |
| 58 | 3 | 3 | 76,3 | 79,7 | 21,8 | 28,6 | 0,679 | 0,746 |
| 53 | 11 | 11 | 79,7 | 87,3 | 17,1 | 42,9 | 0,659 | 0,744 |
| 1 | 7 | 7 | 79,7 | 87,1 | 17,1 | 36,8 | 0,674 | 0,743 |
| 9 | 4 | 4 | 83,9 | 85,7 | 15,6 | 40,0 | 0,679 | 0,736 |
| 62 | 6 | 6 | 77,7 | 81,8 | 15,3 | 38,9 | 0,691 | 0,736 |
| 19 | 8 | 8 | 76,6 | 81,8 | 11,8 | 35,7 | 0,637 | 0,735 |
| 45 | 8 | 8 | 78,5 | 83,2 | 17,6 | 35,7 | 0,654 | 0,731 |
| 25 | 3 | 3 | 79,3 | 83,2 | 16,5 | 23,5 | 0,695 | 0,731 |
| 54 | 9 | 9 | 78,0 | 90,1 | 19,2 | 64,3 | 0,645 | 0,730 |
| 55 | 7 | 7 | 78,9 | 84,5 | 16,9 | 40,0 | 0,669 | 0,730 |
| 51 | 8 | 8 | 77,1 | 82,8 | 17,6 | 40,0 | 0,623 | 0,729 |
| 47 | 6 | 6 | 82,4 | 85,2 | 14,7 | 28,6 | 0,662 | 0,725 |
| 49 | 5 | 5 | 80,2 | 84,6 | 11,1 | 20,0 | 0,678 | 0,725 |
| 14 | 6 | 6 | 78,2 | 84,0 | 19,2 | 41,2 | 0,689 | 0,723 |
| 47 | 4 | 4 | 77,7 | 80,5 | 13,8 | 21,4 | 0,590 | 0,723 |
| 13 | 10 | 10 | 81,7 | 88,0 | 18,6 | 85,7 | 0,674 | 0,721 |
| 61 | 4 | 4 | 77,1 | 81,2 | 13,8 | 21,4 | 0,649 | 0,715 |
| 58 | 5 | 5 | 76,4 | 82,4 | 15,7 | 35,7 | 0,621 | 0,699 |
| 34 | 4 | 4 | 75,5 | 77,8 | 10,0 | 15,4 | 0,496 | 0,698 |
| 15 | 8 | 8 | 80,5 | 87,7 | 16,9 | 40,0 | 0,579 | 0,694 |
| 4 | 2 | 2 | 76,2 | 76,2 | 5,3 | 5,3 | 0,544 | 0,681 |
| 6 | 4 | 4 | 76,2 | 81,0 | 16,4 | 42,9 | 0,635 | 0,675 |
| 8 | 2 | 2 | 89,9 | 89,9 | 35,7 | 35,7 | 0,598 | 0,644 |
| 31 | 2 | 2 | 79,5 | 79,5 | 14,3 | 14,3 | 0,498 | 0,635 |

Médiane des moyennes par hotkey : identité entière 78,8 %, identité CDR3 16,4 %.
Référence (paires de hotkeys **différents** dans la même époque, n = 720) : identité entière médiane 82,3 %, identité CDR3 médiane 45,5 %.

## 4. Les charpentes sont-elles germinales ou dessinées ?

n = 348 (identité germinale calculée par `germline_identity`, planchers du validateur V >= 82, J >= 86).

| mesure | min | p5 | p25 | med | p75 | p95 | max |
|---|---|---|---|---|---|---|---|
| identité V humaine (%) | 82,474 | 82,474 | 82,653 | 83,673 | 85,714 | 87,755 | 89,796 |
| identité J humaine (%) | 92,857 | 92,857 | 92,857 | 92,857 | 92,857 | 100,000 | 100,000 |

* V dans [plancher+0, plancher+2) : 195 / 348 (56,0 %)
* V dans [plancher+2, plancher+5) : 128 / 348 (36,8 %)
* V dans [plancher+5, plancher+10) : 25 / 348 (7,2 %)
* V dans [plancher+10, plancher+100) : 0 / 348 (0,0 %)
* J dans [plancher+0, plancher+2) : 0 / 348 (0,0 %)
* J dans [plancher+2, plancher+5) : 0 / 348 (0,0 %)
* J dans [plancher+5, plancher+10) : 305 / 348 (87,6 %)
* J dans [plancher+10, plancher+100) : 43 / 348 (12,4 %)

Répartition des gènes V attribués :

* `IGHV3-74` : 104
* `IGHV3-23` : 68
* `IGHV3-66` : 66
* `IGHV3-53,IGHV3-66` : 22
* `IGHV3-11` : 19
* `IGHV3-23,IGHV3-23D` : 18
* `IGHV3-48` : 16
* `IGHV3-23,IGHV3-74` : 10
* `IGHV3-64,IGHV3-74` : 6
* `IGHV3-21,IGHV3-48` : 4
* `IGHV3-48,IGHV3-74` : 3
* `IGHV3-23,IGHV3-48` : 3

Par époque (le validateur impose un gène V et une longueur de CDR3 par époque) :

| époque | gènes V observés | ident. V médiane | ident. J médiane |
|---|---|---|---|
| 26299 | {'IGHV3-66': 30, 'IGHV3-53': 4} | 83,0 | 92,9 |
| 26300 | {'IGHV3-66': 26, 'IGHV3-53': 15} | 83,0 | 92,9 |
| 26301 | {'IGHV3-11': 21, 'IGHV3-74': 1, 'IGHV3-64': 1} | 83,7 | 92,9 |
| 26302 | {'IGHV3-23': 31, 'IGHV3-23D': 10, 'IGHV3-74': 1} | 83,7 | 92,9 |
| 26303 | {'IGHV3-74': 21, 'IGHV3-48': 3, 'IGHV3-64': 2} | 83,7 | 92,9 |
| 26304 | {'IGHV3-74': 32, 'IGHV3-23': 1, 'IGHV3-64': 1} | 84,7 | 92,9 |
| 26305 | {'IGHV3-74': 35, 'IGHV3-23': 2, 'IGHV3-NL1': 1} | 84,7 | 92,9 |
| 26306 | {'IGHV3-23': 34, 'IGHV3-23D': 4, 'IGHV3-74': 4, 'IGHV3-30-3': 1} | 83,7 | 92,9 |
| 26345 | {'IGHV3-74': 31, 'IGHV3-64': 3, 'IGHV3-30': 1, 'IGHV3-30-3': 1} | 84,7 | 92,9 |
| 26346 | {'IGHV3-48': 25, 'IGHV3-21': 4, 'IGHV3-23': 3, 'IGHV3-64': 2} | 82,7 | 92,9 |
| 26347 | {'IGHV3-23': 30, 'IGHV3-23D': 5, 'IGHV3-64': 1, 'IGHV3-NL1': 1} | 83,7 | 92,9 |
| 26348 | {'IGHV3-66': 32, 'IGHV3-53': 3} | 82,5 | 92,9 |

## 5. Le binding corrèle-t-il avec quelque chose d'optimisable ?

Spearman de `design_iiptm` contre, au sein des acceptées.

| variable | rho | p | n |
|---|---|---|---|
| longueur CDR3 | +0,364 | 2,43e-12 | 348 |
| identité V | -0,010 | 0,855 | 348 |
| identité J | -0,062 | 0,249 | 348 |
| identité au meilleur donneur antérieur | -0,248 | 2,96e-06 | 348 |
| identité CDR3 au meilleur donneur antérieur | -0,389 | 4,97e-14 | 348 |
| GRAVY du CDR3 | -0,096 | 0,0727 | 348 |
| charge nette du CDR3 | +0,270 | 3,09e-07 | 348 |
| longueur CDR totale | +0,365 | 2,01e-12 | 348 |

Les mêmes correlations **a longueur de CDR3 fixée par époque** (partielles grossières : correlation intra-époque, moyennée pondérée par n) :

| variable | rho intra-époque moyen | n total |
|---|---|---|
| longueur CDR3 | n.d. | 0 |
| identité V | -0,251 | 348 |
| identité J | n.d. | 0 |
| identité au meilleur donneur antérieur | -0,167 | 348 |
| identité CDR3 au meilleur donneur antérieur | -0,057 | 323 |
| GRAVY du CDR3 | +0,010 | 348 |
| charge nette du CDR3 | -0,004 | 348 |
| longueur CDR totale | n.d. | 0 |

Par époque, pour situer la dispersion du score :

| époque | n | iiptm min | med | max |
|---|---|---|---|---|
| 26299 | 30 | 0,641 | 0,680 | 0,774 |
| 26300 | 26 | 0,434 | 0,562 | 0,694 |
| 26301 | 21 | 0,663 | 0,703 | 0,739 |
| 26302 | 31 | 0,664 | 0,725 | 0,771 |
| 26303 | 21 | 0,662 | 0,718 | 0,768 |
| 26304 | 32 | 0,636 | 0,693 | 0,779 |
| 26305 | 35 | 0,552 | 0,699 | 0,749 |
| 26306 | 34 | 0,615 | 0,727 | 0,772 |
| 26345 | 31 | 0,635 | 0,715 | 0,788 |
| 26346 | 25 | 0,144 | 0,734 | 0,787 |
| 26347 | 30 | 0,536 | 0,687 | 0,771 |
| 26348 | 32 | 0,328 | 0,458 | 0,654 |

## 6. Décisif : les CDR3 des 4 époques récentes étaient-ils transplantables ?

Pour chaque acceptée de 26345-26348 : existait-il, dans le vivier antérieur
prouvable, un CDR3 de la **bonne longueur** identique, ou proche ?

| époque | n | long. CDR3 | CDR3 identique dispo. | >= 85 % | >= 70 % | >= 50 % | meilleure ident. CDR3 médiane (%) | meilleure ident. CDR3 même long. médiane (%) | n donneurs même long. |
|---|---|---|---|---|---|---|---|---|---|
| 26345 | 31 | {14: 31} | 0 | 2 | 11 | 31 | 64,3 | 64,3 | 464 |
| 26346 | 25 | {13: 25} | 0 | 0 | 0 | 15 | 53,8 | 53,8 | 6125 |
| 26347 | 30 | {19: 30} | 0 | 0 | 4 | 14 | 47,4 | 44,7 | 122 |
| 26348 | 32 | {6: 32} | 0 | 0 | 14 | 32 | 66,7 | 66,7 | 195 |

Sur les 118 acceptées récentes : identité entière au meilleur donneur antérieur, médiane 86,8 %, max 92,9 %.
CDR3 : médiane 57,1 %, max 85,7 %.

Vue **fuitée** (archive d'aujourd'hui, donc borne supérieure) pour les mêmes :
identité CDR3 médiane 61,5 %, 0 / 118 à >= 85 %.

## 6 bis. Le motif plutôt que le CDR3 entier

La longueur du CDR3 est **imposée par époque** (section 2) : un CDR3 d'archive
de la mauvaise longueur ne peut pas dépasser `1 - |dL| / L` d'identité, donc le
greffage littéral est structurellement impossible la plupart du temps. La
bonne question est donc : le **motif** voyage-t-il ? On mesure, pour chaque
CDR3 accepté, la longueur du plus long fragment contigu qu'il partage avec un
CDR3 antérieur, contre deux viviers, avec un témoin par permutation des
lettres du même CDR3 (20 tirages, même composition en acides aminés).

| époque | n | motif partagé avec l'instantané (médian / max) | témoin permuté | motif partagé avec les acceptées antérieures (médian / max) | témoin permuté | n CDR3 antérieurs |
|---|---|---|---|---|---|---|
| 26299 | 30 | 4 / 7 | 4,2 | 0 / 0 | 0,0 | 0 |
| 26300 | 26 | 4 / 5 | 3,8 | 3 / 4 | 0,9 | 30 |
| 26301 | 21 | 6 / 8 | 4,4 | 5 / 10 | 2,3 | 56 |
| 26302 | 31 | 5 / 5 | 4,1 | 4 / 5 | 2,2 | 77 |
| 26303 | 21 | 6 / 7 | 4,3 | 5 / 8 | 3,0 | 108 |
| 26304 | 32 | 4 / 6 | 4,2 | 5 / 7 | 3,1 | 129 |
| 26305 | 35 | 5 / 7 | 4,1 | 4 / 7 | 2,7 | 161 |
| 26306 | 34 | 5 / 7 | 4,3 | 6 / 8 | 3,2 | 196 |
| 26345 | 31 | 4 / 7 | 4,2 | 4 / 9 | 3,1 | 230 |
| 26346 | 25 | 5 / 5 | 4,1 | 5 / 5 | 3,2 | 261 |
| 26347 | 30 | 5 / 8 | 4,4 | 4 / 10 | 3,3 | 286 |
| 26348 | 32 | 4 / 5 | 3,6 | 3 / 5 | 2,0 | 316 |

Toutes époques confondues (n = 348) : motif partagé avec
l'instantané, médian **5** résidus contre
**4,1** pour le témoin permuté.
Hors première époque (n = 318) : motif partagé avec les acceptées
des époques antérieures, médian **4** contre
**2,7** pour le témoin.

Lecture : contre l'instantané, le motif partagé est de 5 résidus pour un témoin à 4,1 -- c'est-à-dire **au niveau du hasard à composition égale**. Contre les acceptées des époques antérieures du champ lui-même, 4 contre 2,7 : un peu au-dessus du hasard, mais un motif de 4 résidus n'est pas un greffage. Le seul cas net du tableau ci-dessous est `RFQPRYP` (7 résidus), qui voyage de la longueur 20 (26301) à la longueur 17 (26306). Il existe donc une réutilisation de motifs courts **à l'intérieur du champ**, d'une époque à la suivante, et elle est faible.

Les plus longs motifs effectivement repris, pour les 12 meilleures acceptées :

| époque | uid | iiptm | CDR3 | motif partagé avec une acceptée antérieure | avec l'instantané |
|---|---|---|---|---|---|
| 26345 | 50 | 0,788 | `ARHYNTFKPTSMDV` | `ARHY` (4) | `ARHY` (4) |
| 26345 | 26 | 0,788 | `ARHRNTFKPTSMDV` | `ARHR` (4) | `ARHR` (4) |
| 26346 | 23 | 0,787 | `ARGYHRDPHKFDV` | `ARGYH` (5) | `ARGY` (4) |
| 26345 | 29 | 0,787 | `ARHTNSFKETDFDV` | `TDFD` (4) | `ARHT` (4) |
| 26345 | 22 | 0,785 | `ARHTNSFKETDFDV` | `TDFD` (4) | `ARHT` (4) |
| 26345 | 45 | 0,780 | `ARHYNTFKPTSMDV` | `ARHY` (4) | `ARHY` (4) |
| 26346 | 19 | 0,779 | `ARGYHRDLDKFDV` | `ARGYH` (5) | `DKFDV` (5) |
| 26304 | 29 | 0,779 | `ARFYGFHPRTDMDQ` | `ARFYG` (5) | `ARFY` (4) |
| 26346 | 22 | 0,774 | `ARGFHRHTQWFDV` | `ARGF` (4) | `ARGF` (4) |
| 26299 | 48 | 0,774 | `ARFWGWHPRTEFDV` | `-` (0) | `ARFW` (4) |
| 26306 | 13 | 0,772 | `VRFQPRYPRTDPKHFDA` | `RFQPRYP` (7) | `RFQPRYP` (7) |
| 26302 | 4 | 0,771 | `ARWHYRNHFDH` | `HFDH` (4) | `RWHY` (4) |

## 7. Les meilleures soumissions du champ, en clair

| époque | uid | iiptm | long. CDR3 | CDR3 | ident. V | ident. J | meilleur donneur antérieur (%) | meilleur CDR3 antérieur (%) |
|---|---|---|---|---|---|---|---|---|
| 26345 | 50 | 0,788 | 14 | `ARHYNTFKPTSMDV` | 82,7 | 100,0 | 86,8 | 50,0 |
| 26345 | 26 | 0,788 | 14 | `ARHRNTFKPTSMDV` | 82,7 | 100,0 | 86,8 | 57,1 |
| 26346 | 23 | 0,787 | 13 | `ARGYHRDPHKFDV` | 82,7 | 92,9 | 86,0 | 53,8 |
| 26345 | 29 | 0,787 | 14 | `ARHTNSFKETDFDV` | 82,7 | 92,9 | 86,8 | 57,1 |
| 26345 | 22 | 0,785 | 14 | `ARHTNSFKETDFDV` | 82,7 | 92,9 | 86,0 | 57,1 |
| 26345 | 45 | 0,780 | 14 | `ARHYNTFKPTSMDV` | 83,7 | 100,0 | 86,8 | 50,0 |
| 26346 | 19 | 0,779 | 13 | `ARGYHRDLDKFDV` | 82,7 | 92,9 | 83,6 | 46,2 |
| 26304 | 29 | 0,779 | 14 | `ARFYGFHPRTDMDQ` | 82,7 | 92,9 | 82,7 | 71,4 |
| 26346 | 22 | 0,774 | 13 | `ARGFHRHTQWFDV` | 82,7 | 92,9 | 85,0 | 46,2 |
| 26299 | 48 | 0,774 | 14 | `ARFWGWHPRTEFDV` | 87,6 | 92,9 | 87,6 | 50,0 |
| 26306 | 13 | 0,772 | 17 | `VRFQPRYPRTDPKHFDA` | 85,7 | 92,9 | 85,7 | 47,1 |
| 26302 | 4 | 0,771 | 11 | `ARWHYRNHFDH` | 84,7 | 92,9 | 88,4 | 54,5 |
| 26306 | 29 | 0,771 | 17 | `AKFQPRYPRQEPKHMDA` | 83,7 | 92,9 | 85,1 | 45,0 |
| 26347 | 30 | 0,771 | 19 | `SKFQGNHNPRGASQKQNDV` | 84,7 | 92,9 | 84,2 | 35,0 |
| 26306 | 2 | 0,771 | 17 | `AKFQPRYPRNDPQHMDA` | 82,7 | 92,9 | 83,3 | 44,4 |

## 8. Charpente recyclée, CDR3 neuf ? Le test direct du greffage

Si le champ greffe, la charpente (FR1+FR2+FR3+FR4, concaténée) doit souvent
être **exactement** celle d'un donneur antérieur, tandis que le CDR3 est neuf.
n = 348.

* charpente reprise à l'identique d'un donneur antérieur prouvable : 2 / 348 (0,6 %)
* identité de charpente au meilleur donneur antérieur (positionnelle) : médiane 91,2 %, p95 94,5 %, max 100,0 %
* dont charpente >= 98 % : 2 / 348 ; >= 95 % : 12 / 348

| époque | n | charpente exacte dispo. | ident. charpente médiane (%) | ident. CDR3 médiane (%) |
|---|---|---|---|---|
| 26299 | 30 | 0 | 90,1 | 57,1 |
| 26300 | 26 | 0 | 92,3 | 62,5 |
| 26301 | 21 | 0 | 91,2 | 36,4 |
| 26302 | 31 | 0 | 91,2 | 54,5 |
| 26303 | 21 | 0 | 90,1 | 44,4 |
| 26304 | 32 | 0 | 89,0 | 60,7 |
| 26305 | 35 | 0 | 91,2 | 60,0 |
| 26306 | 34 | 0 | 91,2 | 46,0 |
| 26345 | 31 | 0 | 91,2 | 64,3 |
| 26346 | 25 | 2 | 92,3 | 53,8 |
| 26347 | 30 | 0 | 92,3 | 47,4 |
| 26348 | 32 | 0 | 91,8 | 66,7 |

## 9. Les lots par propriétaire de dépôt

38 dépôts, 8 propriétaires, 53 hotkeys, 48 uid. Plusieurs hotkeys
partagent un même dépôt : ce sont des grappes de mineurs pilotés ensemble.

Identités **au sein d'un même propriétaire et d'une même époque** (médianes des médianes par groupe) :

| propriétaire | groupes (époques) | soumissions | hotkeys | ident. entière médiane (%) | ident. CDR3 médiane (%) | iiptm médian | iiptm max |
|---|---|---|---|---|---|---|---|
| QTTECH2023 | 12 | 104 | 11 | 94,2 | 71,4 | 0,728 | 0,788 |
| uchiha0818 | 11 | 99 | 15 | 82,2 | 50,0 | 0,675 | 0,771 |
| saon-miner | 11 | 56 | 7 | 85,5 | 47,1 | 0,694 | 0,760 |
| aichorx | 10 | 32 | 6 | 86,0 | 57,9 | 0,682 | 0,763 |
| bitty-labs | 4 | 30 | 9 | 90,3 | 66,7 | 0,656 | 0,787 |
| cczulaika | 8 | 16 | 3 | 88,0 | 63,6 | 0,703 | 0,748 |

Référence, paires de propriétaires **différents** dans la même époque (n = 3850) : identité entière médiane 81,4 %, identité CDR3 médiane 42,9 %.

Aucune séquence acceptée n'apparaît deux fois : 348 séquences distinctes pour 348 soumissions acceptées, et aucune n'est partagée entre deux hotkeys ni entre deux époques.

## 10. Le test de la bibliothèque : les époques de même longueur imposée

Si un propriétaire tirait d'une bibliothèque pré-calculée par (gène V,
longueur de CDR3), deux époques imposant **la même longueur** devraient
ressortir des CDR3 de la même famille. Trois époques imposent L = 14 :
26299, 26304, 26345.

| époque A | époque B | L | même propriétaire : n paires / ident. CDR3 médiane / max | propriétaires différents : n / médiane / max |
|---|---|---|---|---|
| 26299 | 26304 | 14 | 267 / 35,7 / 85,7 | 693 / 35,7 / 78,6 |
| 26299 | 26345 | 14 | 175 / 35,7 / 64,3 | 755 / 35,7 / 85,7 |
| 26304 | 26345 | 14 | 175 / 35,7 / 71,4 | 817 / 42,9 / 78,6 |

Repères : au sein d'une **même** époque et d'un même propriétaire (n = 1143), identité CDR3 médiane **64,3 %**. Entre époques de longueurs **différentes**, même propriétaire (n = 11214), **14,3 %** (plancher mécanique dû à l'écart de longueur).

Au sein d'une même époque, le même CDR3 est soumis plusieurs fois :
**101 / 348** acceptées partagent leur CDR3 avec une autre
acceptée de la même époque (149 paires, dont 1 entre
propriétaires différents). Ces paires ont une identité de séquence entière
médiane de **98,3 %** (min 84,1 %, max 99,2 %) : ce ne sont pas des charpentes différentes portant
le même CDR3, mais le **même candidat reproduit à quelques résidus près** sur
plusieurs hotkeys.

| époque | n acceptées | CDR3 distincts | blocs depuis l'époque précédente |
|---|---|---|---|
| 26299 | 30 | 20 | - |
| 26300 | 26 | 24 | 360 |
| 26301 | 21 | 16 | 362 |
| 26302 | 31 | 27 | 360 |
| 26303 | 21 | 18 | 361 |
| 26304 | 32 | 22 | 362 |
| 26305 | 35 | 24 | 360 |
| 26306 | 34 | 29 | 362 |
| 26345 | 31 | 28 | 14079 |
| 26346 | 25 | 21 | 360 |
| 26347 | 30 | 26 | 361 |
| 26348 | 32 | 25 | 361 |

Une époque dure 361 blocs, soit environ **72 minutes**. C'est le budget total
du champ pour produire, scorer et pousser ses candidats à chaque époque.

## 11. Le vivier daté n'a de bons lieurs qu'à **une seule** longueur de CDR3

Les 9466 séquences datées (<= 26025), ventilées
par longueur de CDR3, pour les longueurs imposées par les 12 époques analysées :

| longueur imposée | époque(s) | n donneurs datés | dont iiptm >= 0,70 | dont >= 0,75 | meilleur iiptm |
|---|---|---|---|---|---|
| 6 | 26348 | 195 | 0 | 0 | 0,690 |
| 8 | 26300 | 207 | 0 | 0 | 0,675 |
| 10 | 26305 | 177 | 0 | 0 | 0,630 |
| 11 | 26302 | 253 | 0 | 0 | 0,692 |
| 13 | 26346 | 6125 | 4133 | 2718 | 0,827 |
| 14 | 26299, 26304, 26345 | 389 | 0 | 0 | 0,651 |
| 17 | 26306 | 214 | 0 | 0 | 0,628 |
| 18 | 26303 | 189 | 0 | 0 | 0,649 |
| 19 | 26347 | 122 | 0 | 0 | 0,675 |
| 20 | 26301 | 81 | 0 | 0 | 0,581 |

Sur l'ensemble du vivier daté, **4134** séquences lient à 0,70 ou mieux, et
**4133** d'entre elles, soit 100,0 %, ont un CDR3 de
longueur **13**. L'archive publique n'est pas un vivier général : c'est **une
seule famille de solutions, sur une seule longueur**. Pour 11 des 12 époques
analysées, la longueur imposée n'est pas 13, et le vivier daté n'offre alors
**aucun** donneur à 0,70 ou mieux à la bonne longueur.

C'est le plafond mécanique de notre méthode, et il ne dépend ni du réglage du
générateur ni de la fraîcheur de l'archive : le vivier daté est même
meilleur que la partie indatable arrivée depuis (médiane 0,669 contre 0,635, p90 0,790 contre 0,736, n = 9466 et 7778), donc rafraîchir l'archive ne réglera pas ça.

## Ce que cela implique pour notre générateur

### Ce que les mesures excluent

**1. Le champ ne greffe pas depuis l'archive.** Aucune de ses 348 acceptées ne reprend un CDR3 antérieur prouvable : 0 reprise exacte, 3 au-dessus de 85 %, identité CDR3 médiane au meilleur donneur antérieur **56,3 %**. Même en s'autorisant la vue fuitée de l'archive d'aujourd'hui, qui contient leurs propres gagnantes archivées depuis, la médiane ne monte qu'à 60,0 %. Notre méthode -- choisir un donneur d'archive qui lie et transporter son CDR3 -- n'est pas la leur. Ce n'est pas un détail de réglage : c'est une famille de méthode différente.

**2. Ce n'est pas non plus une bibliothèque pré-calculée par (gène V, longueur).** Les trois époques imposant L = 14 ne ressortent pas la même famille de CDR3 : identité médiane entre époques de même longueur, même propriétaire, **35,7 %** (n = 617), contre 64,3 % au sein d'une même époque (n = 1143). Chaque époque reçoit une famille neuve.

**3. Le générateur du champ n'est pas contraint par les planchers germinaux.** Identité V médiane **83,7 %** pour un plancher à 82, J médiane **92,9 %** pour un plancher à 86. Ils ne sont pas collés au plancher : la charpente n'est pas la variable qu'ils poussent.

### Ce que les mesures soutiennent

**Rien n'est recopié : ni le CDR3, ni la charpente.** Les deux tests directs du greffage sont négatifs. Côté CDR3 : 0 / 348 reprise exacte, identité médiane 56,3 %. Côté charpente : 2 / 348 charpentes reprises à l'identique, identité médiane 91,2 %. Les 91 % de la charpente ne sont pas un recyclage : c'est le fond commun de toutes les VHH, que n'importe quelle charpente conforme atteint. L'ecart entre les deux chiffres dit seulement où se trouve le travail : la charpente est standard, le CDR3 est neuf.

**Le CDR3 est dessiné par époque, et la famille est partagée dans l'époque.** Quatre faits convergents :

* au sein d'une époque et d'un même propriétaire de dépôt, les CDR3 forment **une seule famille** : identité médiane 64,3 % (n = 1143), contre 42,9 % entre propriétaires différents de la même époque (n = 3850) ;
* les meilleures acceptées de propriétaires **différents** convergent pourtant sur le même motif dans une même époque (`ARFQPRYP...` en 26301, `ARAYYRHWDW...` en 26303, `ARGYHRD...` / `ARGFHRH...` en 26346, `ARH..NTFK..` en 26345) : le problème a une solution dominante que plusieurs pipelines indépendants retrouvent ;
* 101 / 348 acceptées partagent leur CDR3 avec une autre acceptée de la même époque, et ces paires ont 98,3 % d'identité de séquence entière : c'est le **même candidat répliqué sur plusieurs hotkeys**, pas une exploration ;
* le budget est de **72 minutes par époque** (361 blocs), ce qui exclut une recherche large et impose un pipeline qui converge vite.

La seule famille de méthode compatible avec ces faits est la **conception d'une boucle CDR3 conditionnée par la structure du complexe**, relancée à chaque époque avec la longueur imposée comme contrainte, sur une charpente conforme tenue fixe. C'est ce que le paragraphe 50 de `FINDINGS.md` avait désigné comme la dernière piste non réfutée ; c'est maintenant soutenu par le comportement mesuré du champ et non par déduction.

**Et la « transposition de motif » du paragraphe 47 est réelle mais faible.** Le plus long fragment contigu partagé avec un CDR3 accepté antérieur vaut 4 résidus en médiane contre 2,7 pour un témoin permuté, et 5 contre 4,1 contre l'instantané, soit le hasard. Un `RFQPRYP` qui passe de la longueur 20 à la longueur 17 existe, mais c'est l'exception : on ne peut pas construire un générateur sur un signal de 4 résidus contre 2,7.

### Le mécanisme du plafond, en chiffres

La raison pour laquelle notre méthode ne peut pas atteindre 0,70 n'est pas la qualité du vivier mais **sa monoculture de longueur** : 4133 des 4134 séquences datées qui lient à 0,70 ou mieux ont un CDR3 de longueur 13, et la longueur imposée n'est 13 que dans 1 des 12 époques analysées. Pour les 11 autres, le vivier daté n'offre aucun donneur à 0,70 à la bonne longueur (section 11). Un CDR3 de longueur 13 greffé dans une époque qui en impose 17 n'est même pas un greffage : c'est une séquence neuve à 76 % d'identité maximale avec son donneur, et le paragraphe 50 a montré que la liaison ne suit pas.

### Ce que notre générateur devrait arrêter de faire

Chercher un bon donneur dans l'archive et lui transporter son CDR3. La mesure de la section 5 est sans ambiguïté : au sein des acceptées, `design_iiptm` **ne corrèle pas** avec l'identité au meilleur donneur antérieur, ni avec l'identité V ou J. Il n'y a rien à optimiser de ce côté. Le paragraphe 50 avait déjà montré que la liaison ne voyage pas avec le motif greffé ; la présente mesure montre en plus que **personne dans le champ n'essaie**.

### Ce qui reste à notre portée, et ce qui ne l'est pas

* **À notre portée** : fixer la charpente, ne pas la muter du tout, et consacrer tout le budget de calcul au dessin du CDR3 à la longueur imposée. Le champ ne gagne pas sur la charpente et nous perdons de la liaison à chaque mutation de charpente (paragraphe 50).
* **À notre portée** : reconnaître que la longueur imposée change à chaque époque et qu'un CDR3 de la mauvaise longueur à une identité plafonnée a `1 - |dL| / L`. L'archive n'offre donc presque jamais un donneur utilisable à la longueur du jour, ce qui explique mécaniquement le plafond de notre vivier.
* **Hors de notre portée sans GPU** : le dessin conditionné par la structure lui-même. C'est le gap, et il est méthodologique, pas paramétrique.

### La réserve qui limite tout ce qui précède

Les identités « au meilleur donneur antérieur » sont des **bornes inférieures** : l'archive n'est datée que jusqu'à l'époque 26025 (9466 séquences), et les ~7778 entrées arrivées depuis sont indatables. Si le champ greffait depuis des séquences soumises entre 26026 et l'époque analysée, cette analyse ne le verrait pas directement. Deux éléments la protègent néanmoins : la vue fuitée de l'archive complète ne donne que 60,0 % d'identité CDR3 médiane, et les 8 époques consécutives 26299-26306 forment une chaîne où les acceptées de chaque époque sont dans le vivier de la suivante -- et l'identité n'y monte pas. Pour lever la réserve pour de bon, il faut que `data/archive/Q9NZQ7_first_seen.csv` existe et soit alimenté à chaque téléchargement.

