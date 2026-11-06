# Artefacts pré-construits

Ce qui vaut la peine d'être transporté, c'est ce qui coûte cher à recalculer —
pas ce qui est gros. Mesuré :

| Artefact | Taille | Coût à froid | Ici ? |
|---|---|---|---|
| `Q9NZQ7_regions.csv` (régions IMGT de l'archive) | 2,3 Mo | ~10 min | **oui** |
| `Q9NZQ7_nanobodies.parquet` (jeu d'entraînement) | 1,6 Mo | ~10 min | **oui** |
| le modèle de substitution | 23,1 Mo | **1 min 45** | non : se réentraîne plus vite qu'il ne se télécharge |
| l'archive des soumissions | 2,5 Mo | secondes | non : elle grossit à chaque époque, elle doit être fraîche |

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
