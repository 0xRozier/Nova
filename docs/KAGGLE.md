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

## Le lot à lancer maintenant : trois bras (14 novembre)

Le lot `repaired_26299_26306` est **périmé** : il a été noté (§46), et son
vivier de donneurs était anormalement vieux.
Le lot courant est `Q9NZQ7_arms.csv`, **95 séquences, ~2,2 h**, et il répond à
deux questions dans le même run plutôt qu'une :

| bras | n | question |
|---|---|---|
| `grafted` | 32 | le vivier frais suffit-il ? (prédiction +0,115, §48) |
| `humanised` | 31 | l'humanisation conserve-t-elle la liaison du donneur ? (§47) |
| `real` | 32 | les 4 meilleures soumissions réelles, rang officiel connu |

```sh
python -m scripts.build_arms_batch \
  --humanised-dir <dossier des sorties de scripts.humanise> --depth 4 \
  --out data/observation/batches/Q9NZQ7_arms.csv
```

```sh
python -m scripts.build_boltzgen_notebook \
  --batch data/observation/batches/Q9NZQ7_arms.csv --run arms_26299_26306
```

**Le carnet dans `notebooks/` est déjà celui-là** (`RUN="arms_26299_26306"`,
95 séquences, vérifié). Les deux commandes ci-dessus ne servent qu'à le
reconstruire si le vivier change d'ici là — et il faut le reconstruire si plus
d'une journée passe, puisque c'est précisément la fraîcheur du vivier qui a
valu +0,115.

**Comment lire ce run :** le bras `humanised` se juge contre le `design_iiptm`
de son propre donneur, pas contre le champ. S'il le conserve à 0,05 près,
la piste 2 tient et le vivier utilisable passe de 20 séquences à des milliers.
S'il s'effondre, l'humanisation coûte l'affinité et la piste est morte.

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
