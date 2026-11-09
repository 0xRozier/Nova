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

## Le dépôt GitHub des soumissions

**Rien n'a été créé.** Le `0xRozier/nova-submissions` visible dans les journaux
d'essai est un nom inventé pour l'essai à blanc : sans `--live`, le mineur
s'arrête avant tout appel à GitHub. Aucun dépôt, aucun jeton n'existe.

À faire, et c'est à vous : un jeton d'accès est un secret et ne doit jamais
passer par cette conversation.

1. **Créer un dépôt public.** Il doit être public : le validateur lit le fichier
   via `raw.githubusercontent.com` avec **son propre** jeton
   (`decrypt_submissions`, `github_headers`), qui n'a aucun accès à un dépôt
   privé — il recevrait un 404 et l'époque serait perdue.
2. **Garder les noms courts.** Le pointeur engagé est
   `owner/repo/branch[/chemin]/<hash20>.txt` et la chaîne **plafonne le
   commitment à 100 caractères**. Le mineur refuse de démarrer au-delà.
3. **Créer un jeton fin (fine-grained PAT)** limité à ce seul dépôt, avec la
   permission **Contents: write** et rien d'autre. Le mineur ne fait qu'un PUT
   sur l'API contents.
4. Renseigner sur la machine de minage : `GITHUB_REPO_OWNER`,
   `GITHUB_REPO_NAME`, `GITHUB_REPO_BRANCH`, `GITHUB_REPO_PATH` (vide convient),
   `GITHUB_TOKEN`.

Le validateur interroge aussi l'API des commits pour le chemin du fichier, afin
d'en lire l'horodatage. Le fichier doit donc être **commité sur la branche
annoncée**, pas seulement présent : un PUT par l'API contents crée bien un
commit, donc c'est satisfait, mais une branche ou un chemin qui ne correspond
pas au pointeur donne un horodatage vide.

### Le départage favorise le plus rapide

`ranking.py` départage les nanobodies dans cet ordre : somme de rangs →
sommes par catégorie → **bloc d'engagement** → **horodatage du push GitHub** →
UID. Le plus tôt gagne.

Cela crée un arbitrage réel : noter plus de candidats améliore notre rang mais
retarde la soumission, et nous fait perdre les égalités. Des égalités de somme
de rangs apparaissent dans nos données, donc ce n'est pas théorique. À budget
GPU donné, soumettre tôt n'est pas gratuit.

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
3. L'archive, le jeu de données et la MSA. **Ne reconstruisez pas à froid** :
   `python -m scripts.bootstrap` dit ce qui manque et ce que ça coûte.

   | Élément | Taille | Coût à froid | Transportable ? |
   |---|---|---|---|
   | archive des soumissions | 2,5 Mo | secondes (HuggingFace) | **non**, elle grossit à chaque époque |
   | cache des régions IMGT | 2,3 Mo | ~10 min | **oui**, et c'est le gain principal |
   | jeu de données | 1,6 Mo | ~10 min | oui, à refaire quotidiennement |
   | modèle de substitution | 23 Mo | minutes | oui |

   **Les deux premiers sont dans le dépôt**, sous `data/prebuilt/` — 3,9 Mo pour
   ~20 min économisées :

   ```sh
   python -m scripts.bootstrap --install-prebuilt
   python -m scripts.bootstrap --build        # le reste, dont le modele
   ```

   Le modèle de substitution n'y est **pas**, et c'est volontaire : 23 Mo pour
   1 min 45 d'entraînement, il se refait plus vite qu'il ne se télécharge.
   L'archive non plus : elle doit être fraîche. La MSA et la structure viennent
   avec le clone amont, les poids BoltzGen (~6 Go) avec HuggingFace.

   Le cache des régions est désormais **indexé par séquence** et non par
   position : un rafraîchissement de 50 nouvelles soumissions coûte **0,6 s**
   au lieu de renuméroter les 9466. C'est ce qui rend viable un rafraîchissement
   à chaque époque.
4. Les variables GitHub et le portefeuille bittensor.

## Ordre de mise en service

```sh
# toujours l'interprete du venv : /usr/bin/python porte une autre version de
# bittensor (11.3.0), ou l'API de soumission n'existe pas. Voir FINDINGS 5ter.
PY=/home/user/venv/bin/python

# Rafraichir le vivier AVANT tout le reste : le mineur refuse de soumettre
# au-dela de 40 epoques de retard, et les trois etapes sont necessaires. La
# troisieme seule rebatit le meme vivier perime (mesure le 2026-11-09 : retard
# de 63 epoques, inchange apres un build seul ; 3 apres les trois).
$PY -m scripts.fetch_submissions --last 60   # dechiffrer ce qui a ete soumis
$PY -m scripts.analyze_nanobodies            # joindre les metriques publiees
$PY -m scripts.build_nanobody_dataset        # re-telecharger l'archive, refeaturiser

$PY -m scripts.check_config          # les regles du subnet ont-elles bouge ?
$PY -m pytest tests -q               # parite + aller-retour du chiffrement
$PY -m scripts.mine --once --wallet mon --hotkey default   # a blanc, avec GPU
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
