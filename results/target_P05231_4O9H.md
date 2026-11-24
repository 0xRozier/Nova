# Analyse de cible : P05231 / 4O9H

Chaines du fichier : A, H, L. Chaine cible : A.

## 1. Numerotation

- residus observes : 152, label_seq_id 21..186
- fenetre de la config : 21..186 (166 positions)
- label - auth : [2] (decalage unique)
- trous (boucles desordonnees) : [(50, 63), (106, 109)]
- **la fenetre couvre tout l observe**

## 2. Site de liaison

Partenaires co-cristallises : H, L (H touche 12 residus, L touche 9 residus).
Epitope reel : 19 residus. Site de la config : 8 residus.

**8/8 du site sont au contact.**

| res_index | distance | atteint par |
|---|---|---|
| 24 | 2.82 A | H, L |
| 77 | 3.23 A | H |
| 80 | 3.4 A | H |
| 82 | 2.69 A | H |
| 131 | 4.0 A | H |
| 184 | 2.67 A | H, L |
| 185 | 2.73 A | H |
| 186 | 3.0 A | H |

Couverture par chaine partenaire : H atteint 8/8, L atteint 2/8.
Un site qu'une seule chaine couvre entierement est atteignable par un domaine unique.

## 3. Liants connus

17 entrees PDB portent P05231, 7 avec un liant, 6 domaines variables lisibles.

| PDB | chaine | recouvrement | hallmarks | V | J | nativite | humanite | couche 1 | % CDR3 |
|---|---|---|---|---|---|---|---|---|---|
| 4O9H | H | **8/8** | non | 88.8 | 85.7 | 0.577 | 0.992 | human J identity 85.7 < 86,  | 25% |
| 5FUC | E | **7/8** | oui | 69.1 | 71.4 | 0.852 | 0.414 | human V identity 69.1 < 82,  | 62% |
| 4ZS7 | H | **1/8** | non | 81.8 | 85.7 | 0.662 | 0.782 | human V identity 81.8 < 82,  | 40% |
| 8YWQ | H | **1/8** | non | 92.0 | 92.9 | 0.0 | 1.0 | VHH hallmarks absent (H49=A, | 43% |
| 4CNI | A | **0/8** | non | 84.0 | 85.7 | 0.387 | 0.84 | human J identity 85.7 < 86,  | 38% |
| 8YWR | H | **0/8** | non | 81.8 | 85.7 | 0.662 | 0.782 | human V identity 81.8 < 82,  | 38% |

Planchers en vigueur : V >= 82, J >= 86, hallmarks VHH imposes.

### Paratopes

- **4O9H H** (8/8), 16 contacts, CDR3 = `ANRAGWGMGDY` (L=11) : CDR1 31%, CDR2 38%, CDR3 25%, framework 6%
- **5FUC E** (7/8), 8 contacts, CDR3 = `AARPSIRCASFSATEYKD` (L=18) : CDR1 38%, CDR2 0%, CDR3 62%, framework 0%
- **4ZS7 H** (1/8), 15 contacts, CDR3 = `ARDPDVVTGFHYDY` (L=14) : CDR1 13%, CDR2 40%, CDR3 40%, framework 7%
- **8YWQ H** (1/8), 7 contacts, CDR3 = `ARRANYGTSYDYGMDY` (L=16) : CDR1 0%, CDR2 43%, CDR3 43%, framework 14%
- **4CNI A** (0/8), 13 contacts, CDR3 = `ARESYYGFTSY` (L=11) : CDR1 15%, CDR2 38%, CDR3 38%, framework 8%
- **8YWR H** (0/8), 8 contacts, CDR3 = `ARDPDVVTGFHYDY` (L=14) : CDR1 25%, CDR2 38%, CDR3 38%, framework 0%

Un liant dont le CDR3 porte peu du paratope est le seul pour lequel « garder la charpente, redessiner le CDR3 » a un sens (FINDINGS 73).

### Rendus soumettables

- **4O9H H** : 3 mutations (G49E, L50R, FR4 WGQGTQVTVSS->WGQGTLVTVSS), tous les CDR identiques, couche 1 ok
  ```
  EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYRMYWVRQPPGKEREWVSAISAGGGSTYYGDSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCANRAGWGMGDYWGQGTLVTVSS
  ```
- **8YWQ H** : 2 mutations (A49E, L50R), tous les CDR identiques, couche 1 ok
  ```
  QITLKESGPALVKPTQTLTLTCTFSGFSLSTSGMGVGWIRQPPGKEREWLAHIWWDDDKYYNPSLKSRLTISKDTSKNQVVLTMTNMDPVDTATYYCARRANYGTSYDYGMDYWGQGTLVTVSS
  ```

---
Produit par `scripts.analyze_target`. Chaque nombre vient des fichiers de structure ou du filtre de nativite du validateur.