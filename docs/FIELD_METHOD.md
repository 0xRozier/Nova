# Par quelle methode le champ produit-il ses soumissions acceptees ?

Analyse de 348 soumissions acceptees sur 12 epoques (26299-26306, 26345-26348),
reproductible par `PYTHONPATH=. .venv/bin/python scripts/analyze_field_method.py`.

## 0. Discipline anti-fuite et ce qu'elle coute

L'archive publique ne porte **aucune date**. Le seul etat date dont nous
disposons est l'instantane `data/prebuilt/Q9NZQ7_regions.csv`, **9466**
sequences, qui couvre tout ce qui a ete soumis jusqu'a l'epoque
26025. `data/archive/Q9NZQ7_first_seen.csv` n'existe pas :
`archive_ledger.seen_before()` retombe donc sur cet instantane pour **toutes**
les epoques analysees ici (toutes > 26025).

Le vivier de donneurs « disponible avant l'epoque E » est donc :

* l'instantane (9466 sequences, datees <= 26025) ;
* plus les sequences de `submissions.jsonl` d'epoque < E (458 soumissions
  decryptees sur 12 epoques, acceptees et refusees confondues).

Les ~7778 entrees de l'archive d'aujourd'hui qui ne sont pas
dans l'instantane sont **exclues** : nous ne pouvons pas prouver qu'elles
existaient avant l'epoque analysee. Consequence a garder en tete pour tout
ce qui suit : **chaque identite « au meilleur donneur anterieur » est une
borne inferieure**. Un vrai donneur plus proche a pu exister sans que nous
puissions le dater. Les conclusions ne sont donc tirees que dans le sens ou
cette borne les autorise.

Trou de 39 epoques entre 26306 et 26345 : aucune soumission decryptee n'existe
pour 26307-26344, le vivier des 4 epoques recentes est donc particulierement
sous-estime.

## 1. Les soumissions acceptees derivent-elles de l'archive ?

n = 348. Identite au meilleur donneur anterieur (deciles : min | p5 | p25 | median | p75 | p95 | max).

| mesure | min | p5 | p25 | med | p75 | p95 | max |
|---|---|---|---|---|---|---|---|
| identite sequence entiere (%) | 81.3 | 83.2 | 84.6 | 86.0 | 87.4 | 90.0 | 92.9 |
| identite CDR3, tout donneur (%) | 30.0 | 36.4 | 46.2 | 57.1 | 64.3 | 78.6 | 85.7 |
| identite CDR3, donneur de meme longueur (%) | 25.0 | 31.6 | 44.4 | 54.5 | 64.3 | 78.6 | 85.7 |

* identite entiere dans [0.00, 0.70) : 0 / 348 (0.0 %)
* identite entiere dans [0.70, 0.80) : 0 / 348 (0.0 %)
* identite entiere dans [0.80, 0.85) : 119 / 348 (34.2 %)
* identite entiere dans [0.85, 0.90) : 212 / 348 (60.9 %)
* identite entiere dans [0.90, 0.95) : 17 / 348 (4.9 %)
* identite entiere dans [0.95, 1.01) : 0 / 348 (0.0 %)

CDR3 identique a l'identique a un donneur anterieur : 0 / 348.
CDR3 a >= 85 % d'un donneur anterieur : 3 / 348.
CDR3 a >= 85 % d'une entree de l'archive d'aujourd'hui (vue **fuitee**, bornee par le haut) : 60 / 348, median 71.4 %.

## 2. Les CDR3 se repetent-ils d'une epoque a l'autre ?

| epoque | n acceptees | long. CDR3 imposee (observee) | CDR3 distincts | CDR3 deja vu exact (epoques anterieures) | >= 85 % (epoques anterieures) | >= 85 % (instantane) |
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

Total hors premiere epoque (n = 318) : 0 CDR3 exactement reutilises, 2 a >= 85 % d'un CDR3 accepte anterieur.
Contre l'instantane (n = 348) : 0 CDR3 a >= 85 % d'un CDR3 de l'instantane.

## 3. Les mineurs reutilisent-ils leurs propres sequences ?

52 hotkeys avec >= 2 soumissions acceptees sur >= 2 epoques. Identites entre paires de soumissions du meme hotkey, epoques differentes.

| uid | n acc. | n epoques | ident. entiere moy. | max | ident. CDR3 moy. | max | iiptm moy. | iiptm max |
|---|---|---|---|---|---|---|---|---|
| 26 | 11 | 11 | 78.8 | 86.8 | 17.8 | 36.4 | 0.722 | 0.788 |
| 63 | 11 | 11 | 79.2 | 85.1 | 14.8 | 57.1 | 0.711 | 0.788 |
| 23 | 4 | 4 | 76.3 | 79.7 | 14.7 | 21.4 | 0.620 | 0.787 |
| 17 | 10 | 10 | 77.4 | 81.6 | 15.0 | 35.0 | 0.722 | 0.787 |
| 29 | 10 | 10 | 78.8 | 84.1 | 15.9 | 40.0 | 0.715 | 0.785 |
| 25 | 8 | 8 | 78.2 | 84.0 | 16.6 | 40.0 | 0.690 | 0.780 |
| 19 | 4 | 4 | 75.8 | 78.8 | 15.9 | 23.1 | 0.640 | 0.779 |
| 48 | 9 | 9 | 78.6 | 84.7 | 15.4 | 57.1 | 0.722 | 0.774 |
| 4 | 7 | 7 | 78.7 | 83.5 | 17.1 | 42.9 | 0.642 | 0.771 |
| 62 | 9 | 9 | 78.1 | 84.5 | 18.9 | 42.9 | 0.700 | 0.770 |
| 18 | 10 | 10 | 79.1 | 84.0 | 14.3 | 42.9 | 0.713 | 0.768 |
| 57 | 9 | 9 | 79.0 | 86.0 | 15.4 | 42.9 | 0.697 | 0.767 |
| 7 | 11 | 11 | 81.3 | 89.5 | 18.2 | 50.0 | 0.685 | 0.763 |
| 12 | 9 | 9 | 80.4 | 84.9 | 16.2 | 64.3 | 0.713 | 0.762 |
| 50 | 7 | 7 | 78.0 | 82.1 | 18.8 | 40.0 | 0.705 | 0.761 |
| 2 | 3 | 3 | 75.2 | 77.3 | 10.4 | 15.4 | 0.615 | 0.761 |
| 31 | 6 | 6 | 78.1 | 84.5 | 20.0 | 42.9 | 0.650 | 0.761 |
| 2 | 7 | 7 | 81.8 | 89.9 | 18.7 | 31.6 | 0.689 | 0.760 |
| 23 | 6 | 6 | 78.8 | 82.8 | 15.8 | 42.9 | 0.673 | 0.760 |
| 24 | 3 | 3 | 82.3 | 84.6 | 22.9 | 26.3 | 0.718 | 0.756 |
| 61 | 9 | 9 | 78.1 | 83.2 | 14.9 | 42.9 | 0.648 | 0.754 |
| 10 | 6 | 6 | 81.7 | 87.6 | 17.7 | 71.4 | 0.657 | 0.752 |
| 38 | 7 | 7 | 80.7 | 84.0 | 14.5 | 35.7 | 0.702 | 0.751 |
| 6 | 4 | 4 | 78.6 | 81.2 | 13.8 | 21.4 | 0.655 | 0.751 |
| 59 | 9 | 9 | 79.9 | 86.4 | 16.8 | 50.0 | 0.660 | 0.749 |
| 58 | 4 | 4 | 82.3 | 86.6 | 11.0 | 17.6 | 0.694 | 0.748 |
| 16 | 10 | 10 | 80.3 | 87.4 | 17.5 | 64.3 | 0.664 | 0.748 |
| 36 | 8 | 8 | 79.2 | 87.9 | 18.3 | 30.0 | 0.697 | 0.747 |
| 58 | 3 | 3 | 76.3 | 79.7 | 21.8 | 28.6 | 0.679 | 0.746 |
| 53 | 11 | 11 | 79.7 | 87.3 | 17.1 | 42.9 | 0.659 | 0.744 |
| 1 | 7 | 7 | 79.7 | 87.1 | 17.1 | 36.8 | 0.674 | 0.743 |
| 9 | 4 | 4 | 83.9 | 85.7 | 15.6 | 40.0 | 0.679 | 0.736 |
| 62 | 6 | 6 | 77.7 | 81.8 | 15.3 | 38.9 | 0.691 | 0.736 |
| 19 | 8 | 8 | 76.6 | 81.8 | 11.8 | 35.7 | 0.637 | 0.735 |
| 45 | 8 | 8 | 78.5 | 83.2 | 17.6 | 35.7 | 0.654 | 0.731 |
| 25 | 3 | 3 | 79.3 | 83.2 | 16.5 | 23.5 | 0.695 | 0.731 |
| 54 | 9 | 9 | 78.0 | 90.1 | 19.2 | 64.3 | 0.645 | 0.730 |
| 55 | 7 | 7 | 78.9 | 84.5 | 16.9 | 40.0 | 0.669 | 0.730 |
| 51 | 8 | 8 | 77.1 | 82.8 | 17.6 | 40.0 | 0.623 | 0.729 |
| 47 | 6 | 6 | 82.4 | 85.2 | 14.7 | 28.6 | 0.662 | 0.725 |
| 49 | 5 | 5 | 80.2 | 84.6 | 11.1 | 20.0 | 0.678 | 0.725 |
| 14 | 6 | 6 | 78.2 | 84.0 | 19.2 | 41.2 | 0.689 | 0.723 |
| 47 | 4 | 4 | 77.7 | 80.5 | 13.8 | 21.4 | 0.590 | 0.723 |
| 13 | 10 | 10 | 81.7 | 88.0 | 18.6 | 85.7 | 0.674 | 0.721 |
| 61 | 4 | 4 | 77.1 | 81.2 | 13.8 | 21.4 | 0.649 | 0.715 |
| 58 | 5 | 5 | 76.4 | 82.4 | 15.7 | 35.7 | 0.621 | 0.699 |
| 34 | 4 | 4 | 75.5 | 77.8 | 10.0 | 15.4 | 0.496 | 0.698 |
| 15 | 8 | 8 | 80.5 | 87.7 | 16.9 | 40.0 | 0.579 | 0.694 |
| 4 | 2 | 2 | 76.2 | 76.2 | 5.3 | 5.3 | 0.544 | 0.681 |
| 6 | 4 | 4 | 76.2 | 81.0 | 16.4 | 42.9 | 0.635 | 0.675 |
| 8 | 2 | 2 | 89.9 | 89.9 | 35.7 | 35.7 | 0.598 | 0.644 |
| 31 | 2 | 2 | 79.5 | 79.5 | 14.3 | 14.3 | 0.498 | 0.635 |

Mediane des moyennes par hotkey : identite entiere 78.8 %, identite CDR3 16.4 %.
Reference (paires de hotkeys **differents** dans la meme epoque, n = 720) : identite entiere mediane 82.3 %, identite CDR3 mediane 45.5 %.

## 4. Les charpentes sont-elles germinales ou dessinees ?

n = 348 (identite germinale calculee par `germline_identity`, planchers du validateur V >= 82, J >= 86).

| mesure | min | p5 | p25 | med | p75 | p95 | max |
|---|---|---|---|---|---|---|---|
| identite V humaine (%) | 82.474 | 82.474 | 82.653 | 83.673 | 85.714 | 87.755 | 89.796 |
| identite J humaine (%) | 92.857 | 92.857 | 92.857 | 92.857 | 92.857 | 100.000 | 100.000 |

* V dans [plancher+0, plancher+2) : 195 / 348 (56.0 %)
* V dans [plancher+2, plancher+5) : 128 / 348 (36.8 %)
* V dans [plancher+5, plancher+10) : 25 / 348 (7.2 %)
* V dans [plancher+10, plancher+100) : 0 / 348 (0.0 %)
* J dans [plancher+0, plancher+2) : 0 / 348 (0.0 %)
* J dans [plancher+2, plancher+5) : 0 / 348 (0.0 %)
* J dans [plancher+5, plancher+10) : 305 / 348 (87.6 %)
* J dans [plancher+10, plancher+100) : 43 / 348 (12.4 %)

Repartition des genes V attribues :

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

Par epoque (le validateur impose un gene V et une longueur de CDR3 par epoque) :

| epoque | genes V observes | ident. V mediane | ident. J mediane |
|---|---|---|---|
| 26299 | {'IGHV3-66': 30, 'IGHV3-53': 4} | 83.0 | 92.9 |
| 26300 | {'IGHV3-66': 26, 'IGHV3-53': 15} | 83.0 | 92.9 |
| 26301 | {'IGHV3-11': 21, 'IGHV3-74': 1, 'IGHV3-64': 1} | 83.7 | 92.9 |
| 26302 | {'IGHV3-23': 31, 'IGHV3-23D': 10, 'IGHV3-74': 1} | 83.7 | 92.9 |
| 26303 | {'IGHV3-74': 21, 'IGHV3-48': 3, 'IGHV3-64': 2} | 83.7 | 92.9 |
| 26304 | {'IGHV3-74': 32, 'IGHV3-23': 1, 'IGHV3-64': 1} | 84.7 | 92.9 |
| 26305 | {'IGHV3-74': 35, 'IGHV3-23': 2, 'IGHV3-NL1': 1} | 84.7 | 92.9 |
| 26306 | {'IGHV3-23': 34, 'IGHV3-23D': 4, 'IGHV3-74': 4, 'IGHV3-30-3': 1} | 83.7 | 92.9 |
| 26345 | {'IGHV3-74': 31, 'IGHV3-64': 3, 'IGHV3-30': 1, 'IGHV3-30-3': 1} | 84.7 | 92.9 |
| 26346 | {'IGHV3-48': 25, 'IGHV3-21': 4, 'IGHV3-23': 3, 'IGHV3-64': 2} | 82.7 | 92.9 |
| 26347 | {'IGHV3-23': 30, 'IGHV3-23D': 5, 'IGHV3-64': 1, 'IGHV3-NL1': 1} | 83.7 | 92.9 |
| 26348 | {'IGHV3-66': 32, 'IGHV3-53': 3} | 82.5 | 92.9 |

## 5. Le binding correle-t-il avec quelque chose d'optimisable ?

Spearman de `design_iiptm` contre, au sein des acceptees.

| variable | rho | p | n |
|---|---|---|---|
| longueur CDR3 | +0.364 | 2.43e-12 | 348 |
| identite V | -0.010 | 0.855 | 348 |
| identite J | -0.062 | 0.249 | 348 |
| identite au meilleur donneur anterieur | -0.248 | 2.96e-06 | 348 |
| identite CDR3 au meilleur donneur anterieur | -0.389 | 4.97e-14 | 348 |
| GRAVY du CDR3 | -0.096 | 0.0727 | 348 |
| charge nette du CDR3 | +0.270 | 3.09e-07 | 348 |
| longueur CDR totale | +0.365 | 2.01e-12 | 348 |

Les memes correlations **a longueur de CDR3 fixee par epoque** (partielles grossieres : correlation intra-epoque, moyennee ponderee par n) :

| variable | rho intra-epoque moyen | n total |
|---|---|---|
| longueur CDR3 | n.d. | 0 |
| identite V | -0.251 | 348 |
| identite J | n.d. | 0 |
| identite au meilleur donneur anterieur | -0.167 | 348 |
| identite CDR3 au meilleur donneur anterieur | -0.057 | 323 |
| GRAVY du CDR3 | +0.010 | 348 |
| charge nette du CDR3 | -0.004 | 348 |
| longueur CDR totale | n.d. | 0 |

Par epoque, pour situer la dispersion du score :

| epoque | n | iiptm min | med | max |
|---|---|---|---|---|
| 26299 | 30 | 0.641 | 0.680 | 0.774 |
| 26300 | 26 | 0.434 | 0.562 | 0.694 |
| 26301 | 21 | 0.663 | 0.703 | 0.739 |
| 26302 | 31 | 0.664 | 0.725 | 0.771 |
| 26303 | 21 | 0.662 | 0.718 | 0.768 |
| 26304 | 32 | 0.636 | 0.693 | 0.779 |
| 26305 | 35 | 0.552 | 0.699 | 0.749 |
| 26306 | 34 | 0.615 | 0.727 | 0.772 |
| 26345 | 31 | 0.635 | 0.715 | 0.788 |
| 26346 | 25 | 0.144 | 0.734 | 0.787 |
| 26347 | 30 | 0.536 | 0.687 | 0.771 |
| 26348 | 32 | 0.328 | 0.458 | 0.654 |

## 6. Decisif : les CDR3 des 4 epoques recentes etaient-ils transplantables ?

Pour chaque acceptee de 26345-26348 : existait-il, dans le vivier anterieur
prouvable, un CDR3 de la **bonne longueur** identique, ou proche ?

| epoque | n | long. CDR3 | CDR3 identique dispo. | >= 85 % | >= 70 % | >= 50 % | meilleure ident. CDR3 mediane (%) | meilleure ident. CDR3 meme long. mediane (%) | n donneurs meme long. |
|---|---|---|---|---|---|---|---|---|---|
| 26345 | 31 | {14: 31} | 0 | 2 | 11 | 31 | 64.3 | 64.3 | 464 |
| 26346 | 25 | {13: 25} | 0 | 0 | 0 | 15 | 53.8 | 53.8 | 6125 |
| 26347 | 30 | {19: 30} | 0 | 0 | 4 | 14 | 47.4 | 44.7 | 122 |
| 26348 | 32 | {6: 32} | 0 | 0 | 14 | 32 | 66.7 | 66.7 | 195 |

Sur les 118 acceptees recentes : identite entiere au meilleur donneur anterieur, mediane 86.8 %, max 92.9 %.
CDR3 : mediane 57.1 %, max 85.7 %.

Vue **fuitee** (archive d'aujourd'hui, donc borne superieure) pour les memes :
identite CDR3 mediane 78.9 %, 35 / 118 a >= 85 %.

## 7. Les meilleures soumissions du champ, en clair

| epoque | uid | iiptm | long. CDR3 | CDR3 | ident. V | ident. J | meilleur donneur anterieur (%) | meilleur CDR3 anterieur (%) |
|---|---|---|---|---|---|---|---|---|
| 26345 | 50 | 0.788 | 14 | `ARHYNTFKPTSMDV` | 82.7 | 100.0 | 86.8 | 50.0 |
| 26345 | 26 | 0.788 | 14 | `ARHRNTFKPTSMDV` | 82.7 | 100.0 | 86.8 | 57.1 |
| 26346 | 23 | 0.787 | 13 | `ARGYHRDPHKFDV` | 82.7 | 92.9 | 86.0 | 53.8 |
| 26345 | 29 | 0.787 | 14 | `ARHTNSFKETDFDV` | 82.7 | 92.9 | 86.8 | 57.1 |
| 26345 | 22 | 0.785 | 14 | `ARHTNSFKETDFDV` | 82.7 | 92.9 | 86.0 | 57.1 |
| 26345 | 45 | 0.780 | 14 | `ARHYNTFKPTSMDV` | 83.7 | 100.0 | 86.8 | 50.0 |
| 26346 | 19 | 0.779 | 13 | `ARGYHRDLDKFDV` | 82.7 | 92.9 | 83.6 | 46.2 |
| 26304 | 29 | 0.779 | 14 | `ARFYGFHPRTDMDQ` | 82.7 | 92.9 | 82.7 | 71.4 |
| 26346 | 22 | 0.774 | 13 | `ARGFHRHTQWFDV` | 82.7 | 92.9 | 85.0 | 46.2 |
| 26299 | 48 | 0.774 | 14 | `ARFWGWHPRTEFDV` | 87.6 | 92.9 | 87.6 | 50.0 |
| 26306 | 13 | 0.772 | 17 | `VRFQPRYPRTDPKHFDA` | 85.7 | 92.9 | 85.7 | 47.1 |
| 26302 | 4 | 0.771 | 11 | `ARWHYRNHFDH` | 84.7 | 92.9 | 88.4 | 54.5 |
| 26306 | 29 | 0.771 | 17 | `AKFQPRYPRQEPKHMDA` | 83.7 | 92.9 | 85.1 | 45.0 |
| 26347 | 30 | 0.771 | 19 | `SKFQGNHNPRGASQKQNDV` | 84.7 | 92.9 | 84.2 | 35.0 |
| 26306 | 2 | 0.771 | 17 | `AKFQPRYPRNDPQHMDA` | 82.7 | 92.9 | 83.3 | 44.4 |

## 8. Charpente recyclee, CDR3 neuf ? Le test direct du greffage

Si le champ greffe, la charpente (FR1+FR2+FR3+FR4, concatenee) doit souvent
etre **exactement** celle d'un donneur anterieur, tandis que le CDR3 est neuf.
n = 348.

* charpente identique a l'identique a un donneur anterieur prouvable : 2 / 348 (0.6 %)
* identite de charpente au meilleur donneur anterieur (positionnelle) : mediane 91.2 %, p95 94.5 %, max 100.0 %
* dont charpente >= 98 % : 2 / 348 ; >= 95 % : 12 / 348

| epoque | n | charpente exacte dispo. | ident. charpente mediane (%) | ident. CDR3 mediane (%) |
|---|---|---|---|---|
| 26299 | 30 | 0 | 90.1 | 57.1 |
| 26300 | 26 | 0 | 92.3 | 62.5 |
| 26301 | 21 | 0 | 91.2 | 36.4 |
| 26302 | 31 | 0 | 91.2 | 54.5 |
| 26303 | 21 | 0 | 90.1 | 44.4 |
| 26304 | 32 | 0 | 89.0 | 60.7 |
| 26305 | 35 | 0 | 91.2 | 60.0 |
| 26306 | 34 | 0 | 91.2 | 46.0 |
| 26345 | 31 | 0 | 91.2 | 64.3 |
| 26346 | 25 | 2 | 92.3 | 53.8 |
| 26347 | 30 | 0 | 92.3 | 47.4 |
| 26348 | 32 | 0 | 91.8 | 66.7 |

## 9. Les lots par proprietaire de depot

38 depots, 8 proprietaires, 53 hotkeys, 48 uid. Plusieurs hotkeys
partagent un meme depot : ce sont des grappes de mineurs pilotes ensemble.

Identites **au sein d'un meme proprietaire et d'une meme epoque** (medianes des medianes par groupe) :

| proprietaire | groupes (epoques) | soumissions | hotkeys | ident. entiere mediane (%) | ident. CDR3 mediane (%) | iiptm median | iiptm max |
|---|---|---|---|---|---|---|---|
| QTTECH2023 | 12 | 104 | 11 | 94.2 | 71.4 | 0.728 | 0.788 |
| uchiha0818 | 11 | 99 | 15 | 82.2 | 50.0 | 0.675 | 0.771 |
| saon-miner | 11 | 56 | 7 | 85.5 | 47.1 | 0.694 | 0.760 |
| aichorx | 10 | 32 | 6 | 86.0 | 57.9 | 0.682 | 0.763 |
| bitty-labs | 4 | 30 | 9 | 90.3 | 66.7 | 0.656 | 0.787 |
| cczulaika | 8 | 16 | 3 | 88.0 | 63.6 | 0.703 | 0.748 |

Reference, paires de proprietaires **differents** dans la meme epoque (n = 3850) : identite entiere mediane 81.4 %, identite CDR3 mediane 42.9 %.

Aucune sequence acceptee n'apparait deux fois : 348 sequences distinctes pour 348 soumissions acceptees, et aucune n'est partagee entre deux hotkeys ni entre deux epoques.

