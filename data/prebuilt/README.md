# Artefacts pré-construits

Ce qui vaut la peine d'être transporté, c'est ce qui coûte cher à recalculer —
pas ce qui est gros. Mesuré :

| Artefact | Taille | Coût à froid | Ici ? |
|---|---|---|---|
| `P05231_regions.csv` (régions IMGT de l'archive, **cible vivante**) | 11,7 Mo | ~25 min | **oui** |
| `P05231_nanobodies.parquet` (jeu d'entraînement, **cible vivante**) | 6,4 Mo | ~30 min | **oui** |
| `Q9NZQ7_regions.csv` (régions IMGT de l'archive) | 2,3 Mo | ~10 min | **oui** |
| `Q9NZQ7_nanobodies.parquet` (jeu d'entraînement) | 1,6 Mo | ~10 min | **oui** |
| `P05231_nanobodies_X.npy` (matrice de traits) | 27,1 Mo | ~4 min | non : se refait depuis le parquet, et ne sert qu'au substitut |
| le modèle de substitution | 23,1 Mo | **1 min 45** | non : se réentraîne plus vite qu'il ne se télécharge |
| l'archive des soumissions | 12,0 Mo | secondes | non : elle grossit à chaque époque, elle doit être fraîche |

**Pourquoi les deux fichiers P05231 voyagent, alors que le principe est de ne
transporter que ce qui coûte cher.** Ils coûtent cher *et* ils exigent un outil
système que toutes les machines n'ont pas : `hmmscan`, le binaire de HMMER, sans
lequel `regions()` et `germline_identity()` rendent `None` et
`scripts.build_nanobody_dataset` ne peut pas démarrer. Sur Debian et Ubuntu,
`apt-get update` **puis** `apt-get install -y hmmer` (sans le `update`, le
paquet n'est pas localisé). Là où l'installer n'est pas possible, ces deux
fichiers suffisent à faire tourner le générateur :

```
python -m scripts.bootstrap --install-prebuilt --target P05231
```

qui les copie vers `data/archive/` et `data/dataset/`. **Ne pas enchaîner avec
`--build` dans ce cas** : `bootstrap` considère le jeu de données incomplet tant
que `P05231_nanobodies_X.npy` manque, et la reconstruction qu'il lancerait a
besoin de `hmmscan`. La matrice de traits ne sert qu'à `train_nanobody_surrogate`
— et le substitut vaut 0,145 en prospectif (§83), donc elle n'est pas sur le
chemin du générateur.

`python -m scripts.bootstrap --install-prebuilt` les met en place, puis
`python -m scripts.bootstrap --build` fabrique ce qui reste.

> **Le parquet livré ici est un instantané et il périme.** En production, le
> vivier de donneurs est resté 149 époques (7,5 jours) en retard, et aucune des
> 55 soumissions notées des époques 26173-26174 ne s'y trouvait. Reconstruisez-le
> avec `scripts.build_nanobody_dataset` avant de miner, puis quotidiennement.
> Le mineur refuse désormais de soumettre au-delà de `--max-donor-lag` époques.

Le cache des régions est indexé par séquence, donc il reste valable quand
l'archive grossit : seules les nouvelles soumissions sont numérotées (0,6 s
pour 50). Le parquet est un instantané du jeu d'entraînement ; reconstruisez-le
quotidiennement pour que le vivier de donneurs suive les soumissions récentes.
