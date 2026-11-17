# Le run Kaggle à lancer : refaire le taux de victoire sur le mineur réparé

## Pourquoi celui-ci avant tout autre

Les chiffres que ce dépôt annonce — 14 à 29 % contre le champ complet, 48 % au
mieux des quatre — ont été mesurés sur des candidats **choisis par la porte PSH
cassée** (§43 de `FINDINGS.md`). Ils décrivent un sélecteur qui n'existe plus.
Tant qu'ils ne sont pas refaits, on ne sait pas si le mineur d'aujourd'hui est
meilleur ou pire, et engager du TAO revient à parier sur une mesure prise à
travers un instrument faussé.

Ce run coûte **zéro euro** (GPU Kaggle gratuit, 2×T4, quota 30 h/semaine et 12 h
par session) et répond à la question.

## Ce que le lot contient, et pourquoi c'est mieux qu'avant

Le rejeu porte sur les **époques 26299-26306**, les huit dont nous avons
déchiffré les soumissions. Deux progrès par rapport aux balayages précédents :

1. **La sélection est celle du mineur.** `scripts.observe` appelle désormais
   `mine.select_by_expected_value` — vivier large, PSH sur trois plis,
   classement sur qualité × P(non refusé). Un rejeu qui choisit autrement
   mesure un autre mineur, et c'est précisément l'erreur passée.
2. **Les rivaux ont leurs rangs officiels.** Ces huit époques portent
   **230 soumissions acceptées** avec leurs métriques publiées, soit ~29 par
   époque : pour la première fois le champ est presque complet, là où les
   balayages précédents n'embarquaient que les 2 ou 4 meilleures.

BoltzGen étant stochastique (écart de ~14,5 places entre deux tirages des mêmes
séquences), nos candidats et leurs rivaux sont notés **dans le même run** : des
métriques venant de deux runs différents ne sont pas comparables.

## Le lot à lancer maintenant : profond et honnête (16 novembre)

Les lots précédents sont **périmés deux fois** : ils ont été produits avec la
coupure d'archive qui laissait passer les lignes non datées — donc des donneurs
venus du futur (§REPLAY_AUDIT, entrée 1) — et par un rejeu qui ne donnait pas
au sélecteur le même vivier que le mineur (entrée 5). Le lot courant est le
premier dont les deux défauts sont corrigés.

Le lot est `data/observation/batches/Q9NZQ7_deep.csv` : **170 séquences,
~4,0 h sur 2×T4**.

| époque | candidats | rivaux | vivier honnête |
|---|---|---|---|
| 26304 | 45 | 4 | 8 donneurs ≥ 0,72, plafond 0,774 (gagnant 0,763) |
| 26345 | 54 | 4 | 16 donneurs ≥ 0,72, plafond 0,779 (gagnant 0,788) |
| 26347 | 59 | 4 | aucun donneur ≥ 0,72, plafond 0,675 (gagnant 0,723) |

**La question de ce run, et elle n'est pas « gagnons-nous ».** Le §52 a montré
par analyse de valeurs extrêmes que la queue du greffage n'est pas bornée sous
le champ : borne 0,795, au-dessus du meilleur gagnant sur 12 époques (0,788),
et une espérance du maximum qui passe de 0,654 sur 32 candidats à 0,711 sur
256. Mais cette queue a été ajustée sur le bras `grafted` du §49, **produit
avec le vivier fuité**. Ce run mesure la queue **honnête**, à six fois la
profondeur du précédent.

Comment lire le résultat, par ordre :

1. **La queue.** Ajuster la même loi sur les 45-59 candidats de chaque époque.
   Si ξ reste négatif avec une borne au-dessus de 0,75, la voie « noter plus de
   candidats » est validée et le problème devient matériel, avec un prix connu
   (~5× le débit de deux T4 pour tenir 256 candidats dans une époque).
2. **L'époque 26347**, qui n'a aucun donneur au niveau du champ. Si sa queue
   s'effondre quand celles de 26304 et 26345 tiennent, c'est la disponibilité
   des donneurs qui commande tout, et aucun débit de notation ne sauvera les
   11 époques sur 12 où elle manque.
3. **Les 12 rivaux réels**, notés dans le même run, servent d'ancrage : leurs
   métriques publiées sont connues, donc ils disent si la mesure elle-même est
   fidèle. Le §49 a chiffré son bruit à σ = 0,028 sur l'iiptm, ce qui interdit
   tout classement intra-époque contre 4 rivaux serrés — on juge des
   distributions, pas des places.

**Une déviation à connaître** : le crible PSH a tourné à **une réplique au lieu
de trois**, faute de temps de calcul (12 min par réplique sur 160 candidats).
Cela ne touche pas la mesure de liaison, mais rend un peu plus bruité le choix
de *quels* candidats sont notés.

```sh
python -m scripts.build_arms_batch --depth 60 --epochs 26304,26345,26347 \
  --out data/observation/batches/Q9NZQ7_deep.csv
```

```sh
python -m scripts.build_boltzgen_notebook \
  --batch data/observation/batches/Q9NZQ7_deep.csv --run deep_26304_26345_26347
```

**Le carnet dans `notebooks/kaggle_boltzgen.ipynb` est déjà celui-là**
(`RUN="deep_26304_26345_26347"`, 170 séquences, vérifié). `--epochs` est
obligatoire : sans lui le script ramasse tous les lots journalisés sur disque,
y compris ceux des runs fuités.

## Le lot de conception (16 novembre) : la piste rouverte par le §55

Les deux runs de conception publiés (§11, §14) sont **retirés** : leurs 210
designs avaient **zéro cystéine**, le pont disulfure canonique Cys23-Cys104
détruit par l'inverse folding de la BoltzGen embarquée, qui manque la
machinerie de contrainte des résidus. Ils n'étaient même pas soumissibles
(`min_cysteines`). La conception conditionnée par la structure n'est donc pas
réfutée : elle n'a jamais été testée.

Ce lot la teste. **Deux étapes, délibérément, et non un seul carnet** — un DSSP
silencieusement cassé, une fuite d'archive et 210 designs sans cystéine ont
tous les trois survécu à des runs qui « marchaient ».

### Étape 1 — concevoir (`notebooks/kaggle_design.ipynb`)

`RUN = "design_cdr3_26345"`, **120 designs**, et il installe
**`boltzgen==0.3.2`, l'originale**, pas celle des validateurs.

La spécification est dans `data/design/26345_cdr3/` : **6 charpentes conformes
IGHV3-74** (V 84,7-85,7, J 92,9, **2 cystéines chacune**), CDR1 et CDR2 **fixes**,
le CDR3 conçu sur **9 résidus** à la longueur 14 que l'époque impose, tête `AR`
et queue (`FDV`/`FDY`) gelées. Les 6 structures sont repliées par
NanoBodyBuilder2 et les spans vérifiés résidu par résidu.

L'époque 26345 est choisie pour que ce bras soit **directement comparable** au
bras greffage honnête du §53 : même époque, mêmes rivaux, et notre meilleur
greffé honnête y lit 0,726 contre 0,788 pour le meilleur rival.

### Le lot courant (17 novembre) : 40 charpentes au lieu de 6

`notebooks/kaggle_design.ipynb`, `RUN = "design_wide_26345"`, **120 designs sur
40 charpentes** (environ 3 par charpente), spécification dans
`data/design/26345_wide/` — les 40 portent leurs 2 cystéines et 9 résidus
conçus, vérifié.

**Budget GPU identique** au run `design_cdr3_26345` (120 designs sur
6 charpentes), même époque, donc les deux sont directement comparables. C'est
le test que le §59 laisse ouvert : le bloc confiance est porté par la charpente
(facteur 2,6 entre la meilleure et la pire), la qualité d'une charpente ne se
prédit pas à bas coût, et concentrer le budget sur la meilleure médiane
**réduirait** le maximum — notre meilleur design venait de la charpente de
rang 4 sur 6.

**La prédiction, écrite avant la mesure.** Si le maximum dépasse **0,741** (le
meilleur à 6 charpentes), la diversité de charpentes est le levier et il faut en
ouvrir davantage. S'il plafonne vers 0,70-0,74, la queue est une propriété du
générateur et non du choix de charpente, et il faut revenir à la profondeur de
notation du §53.

Le repliement des 40 charpentes est groupé et reprenable depuis le
`524a3ff` : `fold_many` charge NanoBodyBuilder2 une fois et saute ce qui est
déjà sur disque.

### Étape 2 — vérifier, sur CPU, avant de dépenser du GPU de notation

```sh
python -m scripts.verify_designs --file results/designs_design_cdr3_26345.csv
```

**Ne pas lire les scores si ça échoue.** Le garde-fou compte les cystéines
contre le plancher du validateur *et* contre la paire canonique, et compare les
régions que la spécification a gelées. Passé sur les anciens runs il les refuse
tous les deux, et il y trouve au passage sept variantes distinctes de FR1 là où
la spécification en imposait une.

Si ça passe, construire le lot de notation avec le chemin éprouvé — celui des
170 séquences du §53 — pour obtenir des métriques comparables au champ :

```sh
python -m scripts.build_boltzgen_notebook \
  --batch data/observation/batches/Q9NZQ7_design.csv --run design_scored_26345
```

### Ce que le run tranche

Une seule question : **la BoltzGen originale respecte-t-elle les positions que
la spécification gèle ?** Si les 120 designs portent leurs deux cystéines et la
charpente imposée, la conception devient enfin mesurable et on lit ses scores.
S'ils ne les portent pas, le défaut n'est pas la version embarquée et la
spécification elle-même est en cause.

## La marche à suivre

### 1. Construire le lot et le notebook

```sh
python -m scripts.observe --combined-batch --depth 4
```

```sh
python -m scripts.build_boltzgen_notebook \
  --batch data/observation/batches/Q9NZQ7_combined.csv \
  --run repaired_26299_26306
```

Le `--run` doit être **neuf**. BoltzGen reprend les résultats d'un répertoire
existant au lieu de recalculer ; `boltzgen_rank` refuse maintenant un
`scores.csv` plus vieux que son lot, mais sur Kaggle c'est à vous de ne pas
réutiliser un nom.

À `--depth 4` les huit époques tiennent en **~64 séquences**, soit **~1 h 40**
à 1,5 min/séquence sur 2×T4. `--depth 8` donne ~128 séquences et ~3 h 20 ; le
quota le permet, et plus de rivaux par époque resserre le classement.

### 2. Sur Kaggle

1. Importer `notebooks/kaggle_boltzgen.ipynb`.
2. *Accelerator* → **GPU T4 x2**, *Internet* → **On**.
3. **Essai rapide d'abord** : mettre `TOP = 4` et **Run All** (~15 min). Ça
   vérifie que les poids se téléchargent (~6 Go) et que BoltzGen tourne. Ne
   sautez pas cette étape : trois obstacles d'environnement ont déjà été levés
   là (`MPLBACKEND` imposé par Kaggle, l'import de `bittensor` par la BoltzGen
   de NOVA, et le `TOP` qui ne prenait que les premières lignes).
4. Lot complet : remettre `TOP = None`, puis **Save Version → Save & Run All
   (Commit)**.

### 3. Lire le verdict

Récupérer `scores_nanobody_repaired_26299_26306.csv` dans les *Output* du run,
le placer dans `results/`, puis :

```sh
python -m scripts.review_observation results/scores_nanobody_repaired_26299_26306.csv
```

```sh
python -m scripts.review_observation --summary
```

`review_observation` applique le classement du validateur (rang dense par
métrique, somme des rangs) sur le pool réellement noté, époque par époque, et
dit à quelle place notre candidat serait arrivé.

## Comment juger le résultat, et le piège

**Juger par médianes de bras, pas par victoires par époque, et pas par
l'accord avec les validateurs.** §49 a mesuré le bruit de BoltzGen en notant
les mêmes 32 rivaux dans deux runs : σ = 0,028 sur l'iiptm, alors que deux
rivaux voisins d'une même époque ne sont séparés que de 0,014. Nos deux runs
ne s'ordonnent eux-mêmes que dans 60 % des paires. Un classement
intra-époque contre les 4 meilleures soumissions est donc **sous la
résolution de l'instrument** — l'ancien contrôle « 85-86 % d'accord avec les
validateurs » ne mesurait que des pools contenant des rivaux faibles, et
exiger de le retrouver est exiger l'impossible.

Ce qui porte : une médiane sur 31-32 séquences, dont l'erreur-type tombe à
0,005. C'est à cette échelle qu'il faut lire un bras, et comparer les bras
entre eux.

Ensuite :

- **8 époques, c'est peu.** Un taux de 3/8 a un intervalle de confiance à 95 %
  qui va de 9 % à 76 %. Ce run ne tranchera pas entre « compétitif » et
  « excellent » ; il tranchera entre « compétitif » et « non compétitif », ce
  qui est la question du moment.
- **Le rejeu reste un rejeu.** L'archive est tronquée à l'époque rejouée
  (`archive_cutoff`), donc le générateur ne grafte pas depuis le futur — c'est
  ce qui avait produit 4 des 5 victoires apparentes du premier balayage. Mais
  les rivaux, eux, ont joué avec l'information de leur temps, et nous rejouons
  avec un modèle de substitution entraîné sur tout.
- **Battre le champ noté n'est pas être accepté.** 47 % des refus du champ ne
  sont expliqués par aucune de nos portes
  (`results/why_refused_26299_26306.csv`), et 12 % des soumissions acceptées
  seraient refusées à tort par nos portes sur un seul pli. Une victoire en rejeu
  suppose que la soumission passe, ce qui n'est jamais arrivé en production.

## Ce qu'un GPU loué apporterait en plus, et que celui-ci n'apporte pas

Rien sur la qualité des candidats. Deux choses seulement :

- le **débit réel** sur la carte qu'on louerait (les facteurs ×3 à ×8 du
  `results/economics.txt` sont des estimations, jamais mesurées) ;
- le chemin `boltzgen_rank` **dans le mineur**, qui n'a jamais tourné faute de
  carte. Sa plomberie est testée (`tests/test_boltzgen_rank.py`), BoltzGen
  lui-même est éprouvé sur Kaggle ; ce qui reste à découvrir est que les deux
  marchent ensemble, dans le budget d'une époque.

Donc : ce run d'abord, la carte ensuite, et le `--live` seulement si les deux
disent oui.
