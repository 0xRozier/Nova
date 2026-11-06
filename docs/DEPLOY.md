# Mise en service — ce qui est testé et ce qui ne l'est pas

## Testé de bout en bout, ici

- génération, contrôle de validité hors ligne, parité avec le validateur (`pytest tests`)
- porte de développabilité complète : TNP **et** Aggrescan3D, contrôle positif propre
- classement BoltzGen, et son accord avec les validateurs (86 % sur 118 paires)
- chiffrement timelock, format de fichier, et le verrou lui-même (drand HTTP 425)
- la boucle du mineur en essai à blanc, portes comprises

## Jamais exécuté — à découvrir en production

| Élément | Risque |
|---|---|
| `set_commitment` sur la chaîne | la chaîne limite la fréquence des engagements ; l'erreur est désormais journalisée au lieu de tuer la boucle |
| le PUT GitHub | un échec après l'engagement laisse un pointeur orphelin et l'époque est perdue |
| la recherche d'UID depuis un portefeuille | testé uniquement avec `--uid` en dur |
| `boltzgen_rank` dans le mineur | le chemin GPU du mineur n'a jamais tourné faute de carte |
| le budget temps réel d'une époque | estimé à ~20 min pour 8 candidats, jamais mesuré d'un bout à l'autre |

## À installer sur la machine louée

Dans cet ordre. **Si l'un manque, le mineur tourne quand même mais dégradé** —
et les deux premiers coûtent cher en silence.

1. **BoltzGen + ses poids** (~6 Go) et le clone amont (`--nova-dir`).
   Sans lui, le mineur se rabat sur le modèle de substitution : 14 % de
   victoires au lieu de 29 %, soit pire qu'un tirage au hasard.
2. **TNP + Aggrescan3D.** TNP dans son propre venv avec torch CPU uniquement
   (la version CUDA remplit le disque), **plus le correctif d'ImmuneBuilder 1.2**
   à `refine.py:346` où `{'Threads', str(n)}` doit devenir
   `{'Threads': str(n)}`. Puis `sh install/install-aggrescan3d.sh` dans le dépôt
   du filtre. Sans eux, 39 % des soumissions partent au rejet silencieux.
3. L'archive, le jeu de données et la MSA (`scripts.build_nanobody_dataset`).
4. Les variables GitHub et le portefeuille bittensor.

## Ordre de mise en service

```sh
python -m scripts.check_config          # les regles du subnet ont-elles bouge ?
python -m pytest tests -q               # parite + aller-retour du chiffrement
python -m scripts.mine --once --wallet mon --hotkey default   # a blanc, avec GPU
```

Le troisième doit montrer `ranked_by: boltzgen` dans le journal et tenir dans le
budget d'époque. S'il affiche `surrogate`, BoltzGen n'est pas joignable et il ne
faut pas passer en `--live`.

Puis `--live`, et surveiller les premières époques : vérifier que le pointeur
engagé correspond au fichier publié, et que `status` vaut `submitted`.

## Attentes chiffrées

- taux de victoire mesuré : **14 à 29 %** contre le champ complet, les deux
  portes appliquées (IC 95 % de 3 % à 52 %, sur 21 époques seulement)
- seuil de rentabilité : 0,3 à 1 %
- place médiane 2ᵉ sur 18, top 3 dans 81 % des époques
- point faible connu : les CDR3 courts (≤8), environ une époque sur cinq,
  place médiane 9
