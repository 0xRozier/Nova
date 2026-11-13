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

**Le contrôle à regarder d'abord n'est pas notre place, c'est l'accord avec les
validateurs.** Chaque époque porte ses rivaux avec leur rang officiel : si notre
run les ordonne comme eux (85-86 % des paires sur les trois runs précédents),
la mesure vaut quelque chose. Sinon, rien de ce qui suit ne compte.

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
