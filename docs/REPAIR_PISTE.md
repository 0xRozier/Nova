# Réparer un lieur de l'archive au lieu de le regreffer

*Mesuré le 15 novembre 2026, sur CPU, cible Q9NZQ7. Tout ce qui est écrit ici
porte sur la **conformité** et la **nouveauté** — les deux murs que le
validateur oppose. **Rien ici ne mesure la liaison**, qui demande un GPU. Les
candidats produits sont étiquetés non testés, et ils le sont.*

Code : `scripts/repair_binder.py`, tests `tests/test_repair_binder.py`.

## La question

L'archive publique porte 16 700 séquences, dont **3972 lient la cible à
design_iiptm ≥ 0,72**, le niveau du champ. Notre générateur plafonne à 0,50.
La liaison dont nous avons besoin est donc déjà sur le disque, et deux murs
l'en séparent :

- **conformité** : identité germinale humaine V ≥ 82 et J ≥ 86. **256** des
  3972 passent ; les autres échouent surtout sur J, médiane 64,3 ;
- **nouveauté** : `nova_lab.nanobodies.is_duplicate` refuse une soumission trop
  proche de l'archive — et un donneur d'archive y est présent à 100 %.

Le §50 a réfuté la greffe (`scripts/humanise.py`) sur GPU : 28 % de rétention
médiane, 20 % quand le CDR3 est conservé au résidu près, et aucune corrélation
avec le nombre de mutations de charpente (ρ = −0,13, p = 0,49, médiane 12). La
question posée ici est donc : **existe-t-il une transformation simultanément
conforme, nouvelle, et minimale et confinée loin du paratope ?**

## 1. L'arithmétique de la nouveauté, avant tout code

Lue sur `is_duplicate` et `whole_identity`, pas sur une formule supposée.

`whole_identity` est le `ratio()` de `difflib.SequenceMatcher`, soit
2M/(2L) = (L − k)/L pour k substitutions sur une séquence de longueur L
(les substitutions ne changent pas la longueur, donc le CDR3 de l'époque est
préservé). Mesure sur 40 donneurs réels × 20 tirages aléatoires par k :

| k substitutions | identité min | médiane | max | fraction < 0,95 |
|---|---|---|---|---|
| 5 | 0,9583 | 0,9583 | 0,9667 | **0 %** |
| 6 | 0,9500 | 0,9500 | 0,9590 | **0 %** |
| 7 | 0,9333 | 0,9417 | 0,9512 | 96,6 % |
| 8 | 0,9167 | 0,9333 | 0,9431 | 100 % |

**Le plancher est donc k = ⌊0,05 L⌋ + 1 : 7 substitutions pour un 120-mère, 6
pour un 114-mère, et 8 pour être sûr.** C'est la branche à 0,95, qui est un
refus certain.

Mais 7 à 8 substitutions laissent l'identité entière vers **0,94, donc au-dessus
de 0,85** : la branche CDR3 d'`is_duplicate` devient alors elle aussi un refus
certain (`certain = identity >= 0.85`) à moins que la similarité CDR3 ne tombe
sous 0,85. `cdr_similarity` valant 1 − mismatches/longueur à longueur égale, cela
**force des substitutions dans le CDR3 lui-même** :

| longueur CDR3 | pour CDR3 < 0,85 (refus) | pour CDR3 < 0,80 (avertissement) |
|---|---|---|
| 6 | 1 | 2 |
| 7–9 | 2 | 2 |
| 10–13 | 2 | 3 |
| 14–19 | 3 | 4 |
| 20–24 | 4 | 5 |

Descendre l'identité entière sous 0,85 pour échapper au « certain » coûterait
19 substitutions : hors de question. **La nouveauté n'est donc pas confinable à
la charpente ; elle passe obligatoirement par le CDR3.**

Contre l'**archive entière**, et non contre le seul donneur : après réparation,
le plus proche voisin non-soi mesuré sur les candidats retenus est à
**0,9333–0,9474** (n = 5), c'est-à-dire que **le donneur lui-même est la
contrainte active** — aucun voisin n'oblige à aller plus loin. À noter aussi
que 6740 des 17 244 entrées de l'archive ne sont pas numérotables : pour elles
`is_duplicate` bascule sur `elif identity >= 0.90`, un avertissement ; le plus
proche voisin non numérotable mesuré est à 0,83.

**Minimums réellement obtenus sur donneurs réels**, toutes mutations dirigées
vers un germinal humain et confinées au pied du CDR3 et à FR4 :

| donneur (iiptm) | longueur CDR3 | k total | dont CDR3 | dont apex | dont FR1/FR2/FR3 | voisin le plus proche |
|---|---|---|---|---|---|---|
| 0,7972 | 13 | **7** | 4 | 0 | 0 | 0,9417 |
| 0,7975 | 13 | 8 | 5 | 0 | 0 | 0,9417 |
| 0,8064 | 13 | 8 | 5 | 0 | 0 | 0,9417 |
| 0,7988 | 13 | 9 | 3 | 0 | 0 | 0,9333 |
| 0,7254 | 7 | 6 | 5 | 0 | 0 | 0,9474 |

## 2. La réparation du J

`germline_identity` note le J sur les **14 positions IMGT 115 à 128** — les 14
allèles IGHJ humains couvrent exactement cette fenêtre, sans trou. L'identité J
est donc quantifiée par pas de 7,14 points, et le plancher de 86 signifie
**13/14 au minimum, soit au plus un seul mésappariement**. (12/14 = 85,7 échoue.)

Distribution du nombre minimal de substitutions dans la fenêtre J, sur les
**3972** lieurs à ≥ 0,72 (minimum sur les 14 allèles de mismatches − 1) :

| substitutions nécessaires | n | part |
|---|---|---|
| 0 (déjà conforme) | 256 | 6,4 % |
| 1 | **0** | — |
| 2 | 102 | 2,6 % |
| 3 | 751 | 18,9 % |
| 4 | 1888 | 47,5 % |
| 5 | 798 | 20,1 % |
| 6 | 139 | 3,5 % |
| 7 | 38 | 1,0 % |
| inatteignable | **0** | — |

Aucune n'est inatteignable : réparer toute la fenêtre vers un allèle donne
J = 100 %. Restreint aux donneurs qui passent déjà V ≥ 82 (**n = 400**) :
0 → 256, 2 → 22, 3 → 64, 4 → 43, 5 → 14, 6 → 1.

Vérification de la méthode : ma lecture de la fenêtre (3 derniers résidus du
CDR3 + FR4) concorde avec `germline_identity` sur **3970 des 3972** ; les deux
écarts sont des séquences où le découpage d'ANARCI diffère de celui d'abnumber,
et mon compte y est pessimiste.

**La réparation n'est presque jamais confinable à FR4.** Sur les 144 donneurs
V ≥ 82 qui ont besoin d'une réparation J :

| mésappariements dans la queue du CDR3 (115-117) | n |
|---|---|
| 1 | 1 |
| 2 | 66 |
| 3 | 77 |

**1 seul des 144** peut être réparé en ne touchant que FR4 (118-128) ; les 143
autres doivent déplacer au moins un des trois derniers résidus du CDR3. FR4
seul — l'espoir le mieux placé structurellement — ne suffit pas.

## 3. Les deux exigences se heurtent-elles ?

Non pour les donneurs qui échouent sur J, oui pour ceux qui passent déjà. C'est
le résultat central, et il est l'inverse de ce que j'attendais.

Les mutations admissibles ici sont celles qui vont **vers un résidu germinal
humain** et qui tombent dans le **pied du CDR3 (IMGT 105-106, couverts par la
référence V, et 115-117, couverts par la référence J)** ou dans **FR4
(118-128)**. L'apex du CDR3 (107-114), CDR1, CDR2, FR1, FR2 et FR3 ne sont
jamais touchés.

- **La réparation du J finance la nouveauté.** Un donneur dont le J est à 64 %
  a 4 à 9 mésappariements dans sa fenêtre : les corriger, c'est à la fois le
  rendre conforme et le rendre nouveau. Sur les 3972 lieurs, **3730 disposent
  de ≥ 7 substitutions dirigées vers un germinal à l'intérieur de
  {105,106} ∪ {115…128}**, et **3666 de ≥ 8**. La queue du CDR3 fournit les
  substitutions que la branche CDR3 exige, FR4 fournit le reste.
- **Un donneur déjà conforme n'a plus rien à déplacer.** Sa fenêtre J est déjà
  germinale (médiane : 1 mésappariement). Sur les 256 conformes, le quota CDR3
  est atteignable pour 167 au seuil 0,85, mais le budget total de 8 est
  **inatteignable dans la fenêtre** : il faut aller chercher FR1/FR3 (médiane
  14 positions dirigeables vers le germinal). C'est exactement le territoire
  que le §50 condamne.

Par longueur de CDR3, le quota (tableau du §1) est couvert par le pied et la
queue jusqu'à 19 ; à 20-24 il faut 4 ou 5 substitutions dans le CDR3 et le pied
n'en offre au mieux que 5.

