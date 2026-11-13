# Mise en service — ce qui est testé et ce qui ne l'est pas

## Testé de bout en bout, ici

- génération, contrôle de validité hors ligne, parité avec le validateur (`pytest tests`)
- porte de développabilité complète : TNP **et** Aggrescan3D, contrôle positif
  propre, sur l'installation que le validateur fait tourner — les 36
  soumissions qu'il a acceptées à l'époque 26174 lisent toutes sous son seuil
  de 89 (`scripts.doctor`, `scripts.preflight --calibrate`)
- classement BoltzGen, et son accord avec les validateurs (86 % sur 118 paires)
- chiffrement timelock, format de fichier, et le verrou lui-même (drand HTTP 425)
- la boucle du mineur en essai à blanc, portes comprises

## Jamais exécuté — à découvrir en production

**Cette liste disait trois choses fausses jusqu'au 13 novembre.** La campagne
des époques 26172-26185 a bel et bien engagé 11 pointeurs sur la chaîne
(`publish_metadata_extrinsic`, `success: True` dans `journal_final.jsonl`) et
publié 11 fichiers sur GitHub. Ces deux chemins sont éprouvés en production ;
ce qui n'a jamais tourné est plus court, et c'est précisément ce qui décide le
rang.

| Élément | État |
|---|---|
| `set_commitment` sur la chaîne | **éprouvé** : 11 engagements réussis en production. La chaîne limite la fréquence ; l'erreur est journalisée au lieu de tuer la boucle |
| le PUT GitHub | **éprouvé** : 11 fichiers publiés, 1 échec en 404 à l'époque 26172. Le fichier part désormais AVANT l'engagement, donc un échec laisse un orphelin inoffensif et non un pointeur mort |
| le budget temps réel d'une époque | **mesuré** le 13 novembre en essai à blanc : 12,6 min jusqu'à BoltzGen (génération 4 min, PSH 7 min, archive et développabilité 1,5 min), 29 min de marge restantes |
| **`boltzgen_rank` dans le mineur** | **jamais exécuté, faute de carte.** C'est le bloquant : sans lui le mineur abandonne chaque époque, par choix, car se rabattre sur le modèle de substitution vaut 14 % de victoires contre 48 %. La structure (8AOK.cif) et la MSA sont dans le clone amont ; il manque le paquet BoltzGen, ses poids (~6 Go) et un GPU |
| la recherche d'UID depuis un portefeuille | testé uniquement avec `--uid` en dur |
| **le taux de victoire, depuis la réparation de la porte PSH** | les 14 à 29 % ont été mesurés sur des candidats choisis par une porte fausse. Le chiffre est à refaire avant d'engager quoi que ce soit |

## Le dépôt GitHub des soumissions

**Il existe depuis la campagne de novembre** : `0xRozier/sn68-submissions`, visible
dans les pointeurs engagés de `journal_final.jsonl`. Le paragraphe qui disait
ici que rien n'avait été créé datait d'avant cette campagne. Ce qui suit reste
la marche à suivre pour en remonter un, et la liste de ce qu'il faut vérifier
sur celui qui existe — notamment qu'il est **public**, puisque l'époque 26172 a
été perdue sur un HTTP 404 au PUT.

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
2. **La porte de développabilité, exactement comme le validateur l'installe.**
   C'est la partie la plus facile à rater en silence : elle a coûté les douze
   époques de production. Les épingles viennent de
   `NOVA-nanobody-filter/install/requirements-minimal.txt` et de
   `install/install.sh` du filtre, pas de nos souvenirs.

   Dans un venv à part, Python 3.12, torch **CPU uniquement** (la version CUDA
   remplit le disque) :

   ```sh
   python3.12 -m venv /home/user/venv-tnp
   /home/user/venv-tnp/bin/pip install torch --index-url https://download.pytorch.org/whl/cpu
   /home/user/venv-tnp/bin/pip install numpy scipy scikit-learn biopython==1.86 freesasa openmm anarci
   /home/user/venv-tnp/bin/pip install "git+https://github.com/openmm/pdbfixer.git@v1.12"
   /home/user/venv-tnp/bin/pip install "git+https://github.com/oxpig/ImmuneBuilder.git@0df4e2ad82a1aa60f37ea9dae335d1198159ef78"
   /home/user/venv-tnp/bin/pip install "git+https://github.com/oxpig/TNP.git@29dcac72f1380e8538e8870f45a699d3c6156162"
   ```

   Puis les trois choses sans lesquelles TNP répond quand même, et répond faux :

   ```sh
   # le correctif process_pdb du filtre : sans lui DSSP echoue sur TOUS les modeles
   cp /home/user/nanofilter/install/patches/TNP/process_pdb.py \
      /home/user/venv-tnp/lib/python3.12/site-packages/scripts/process_pdb.py
   ```

   ```sh
   # le typo set/dict d ImmuneBuilder 1.2 : sinon 3 sequences sur 37 reviennent sans modele
   sed -i "346s/{'Threads', str(n_threads)}/{'Threads': str(n_threads)}/" \
      /home/user/venv-tnp/lib/python3.12/site-packages/ImmuneBuilder/refine.py
   ```

   ```sh
   # DSSP 3.1.4, et PAS celui d apt : mkdssp 4.x refuse les modeles de
   # NanoBodyBuilder2 (pas de ligne HEADER) et TNP rend alors une PSH qui n est
   # pas celle du validateur. Meme structure : 89,50 sans DSSP, 91,21 avec.
   git clone --depth 1 https://github.com/cmbi/dssp.git && cd dssp \
     && ./autogen.sh && CXXFLAGS="-std=c++14 -Wno-error" ./configure --prefix=/home/user/dssp3 \
     && make -j"$(nproc)" CXXFLAGS="-std=c++14 -Wno-error -O2" && make install \
     && ln -sf /home/user/dssp3/bin/mkdssp /home/user/venv-tnp/bin/mkdssp \
     && ln -sf /home/user/dssp3/bin/mkdssp /home/user/venv-tnp/bin/dssp
   ```

   Enfin Aggrescan3D et sa FreeSASA épinglée, depuis le dépôt du filtre :

   ```sh
   sh /home/user/nanofilter/install/install-aggrescan3d.sh
   ```

   Les chemins par défaut sont `/home/user/venv-tnp/bin/TNP` et
   `/home/user/nanofilter` ; `NOVA_TNP`, `NOVA_DSSP` et `NOVA_NANOFILTER` les
   déplacent. Il faut aussi `hmmscan` (`apt install hmmer`) : sans lui ANARCI
   rend `None` **sans lever d'erreur**, et tout le reste se construit sur du
   vide.

   **Puis vérifier, parce qu'aucun de ces outils ne refuse de répondre :**

   ```sh
   python -m scripts.doctor
   ```

   Il fait tourner la chaîne entière sur un nanobody clinique et sort en code 1
   si un nombre manque. Ne démarrez pas le mineur avant qu'il soit vert.
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

$PY -m scripts.doctor                # les outils mesurent-ils ce que mesure le validateur ?
$PY -m scripts.preflight --calibrate # et le disent-ils des 36 acceptees de 26174 ?
$PY -m scripts.check_config          # les regles du subnet ont-elles bouge ?
$PY -m pytest tests -q               # parite + aller-retour du chiffrement
$PY -m scripts.mine --once --wallet mon --hotkey default   # a blanc, avec GPU
```

`scripts.doctor` et `scripts.preflight --calibrate` sont les deux qui
manquaient en novembre. Le second doit afficher **0 acceptee au-dessus de 89** :
s il en affiche davantage, notre echelle PSH a derive de la leur et la porte
mesure autre chose — c est exactement l etat dans lequel les douze epoques de
production ont ete perdues.

Le dernier doit montrer `ranked_by: boltzgen` dans le journal et tenir dans le
budget d'époque. S'il affiche `surrogate`, BoltzGen n'est pas joignable et il ne
faut pas passer en `--live`.

Puis `--live`, et surveiller les premières époques : vérifier que le pointeur
engagé correspond au fichier publié, et que `status` vaut `submitted`.

## Attentes chiffrées

**À lire en sachant que la seule campagne réelle a fait 0 sur 12.** Les chiffres
ci-dessous sont des rejeux hors ligne ; ils supposent que la soumission est
acceptée, ce qui n'était jamais arrivé avant que la porte PSH soit réparée.

- taux de victoire en rejeu : **14 à 29 %** contre le champ complet, les deux
  portes appliquées (IC 95 % de 3 % à 52 %, sur 21 époques seulement)
- campagne réelle, époques 26172-26185 : **11 soumissions, 0 notée**, toutes
  refusées sur PSH (`docs/POST_MORTEM.md`). La cause était l'installation, pas
  les séquences, et elle est corrigée ; rien ne le confirme encore en production
- seuil de rentabilité : 0,3 à 1 %
- place médiane 2ᵉ sur 18, top 3 dans 81 % des époques
- points faibles connus :
  - les CDR3 courts (≤8), environ une époque sur cinq, place médiane 9 ;
  - les CDR3 longs (≥17) pour la PSH : le candidat médian lit **100,6** là où
    le champ accepté lit 83,5, donc **15 %** seulement passent 89. C'est un
    problème de génération et aucun vivier ne le règle — il n'achète que la
    queue, où le meilleur de 100 lit 74,7. Le biais d'hydrophobicité a été
    essayé et mesuré sans effet (`results/psh_generator_ab_cdr3_18.csv`).
