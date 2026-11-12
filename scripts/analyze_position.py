"""Lire le test de position -- ecrit AVANT que les resultats existent.

Le §21 avait mesure des ecarts par position allant jusqu a 0,43 et avait refuse
de les reporter, parce que les positions etaient tirees au hasard et que
l appartenance a une position etait confondue avec "quelle variante s est
effondree". Le dispositif corrige impose l ensemble des positions a nombre de
substitutions constant, et il n a de valeur que si sa lecture est fixee d avance
-- sinon on choisira, apres coup, la statistique qui dit quelque chose.

Donc cette lecture, et rien d autre :

1. le temoin parent. S il ne retrouve pas sa valeur a 0,03 pres, la session est
   suspecte et le reste ne se lit pas.
2. les deux medianes, les deux maximums, les deux taux d effondrement.
3. Mann-Whitney sur les scores des deux bras (bilateral), et la taille d effet
   en delta de Cliff, qui ne suppose rien de la forme des distributions.
4. la meme comparaison sur les seules survivantes, parce que le §21 a montre
   que l effondrement et la degradation graduelle sont deux mecanismes distincts
   et qu un bras peut differer sur l un sans differer sur l autre.

Ce que l on NE fera PAS : chercher un effet position par position (le
dispositif n a que deux ensembles), ni tester d autres seuils que 0,30 pour
l effondrement, ni rapporter un p apres avoir essaye plusieurs decoupages.

Verdict, decide d avance :
  p < 0,05 et |delta| > 0,4  -> la position est un levier, le generateur doit
                                cesser de tirer au hasard et tout le §29 gagne
                                le meme facteur.
  sinon                      -> pas d effet LARGE. Ce n est pas "pas d effet" :
                                avec 13 contre 9 un effet modere passe
                                inapercu, et il faudra le dire ainsi.

Usage : python scripts/analyze_position.py <scores.csv>
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parent.parent
COLLAPSE = 0.30
CONTROL_TOLERANCE = 0.03


def cliffs_delta(a, b) -> float:
    """P(a > b) - P(a < b). Sans hypothese sur la forme des distributions."""
    gt = sum(x > y for x in a for y in b)
    lt = sum(x < y for x in a for y in b)
    return (gt - lt) / (len(a) * len(b))


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit(__doc__.strip().splitlines()[-1])
    d = pd.read_csv(sys.argv[1])

    print("=== 1. le temoin ===")
    ctl = d[d.role == "parent"]
    ok = True
    for r in ctl.itertuples():
        prev = r.previous_design_iiptm
        if pd.isna(prev) or prev == "":
            print(f"  parent jamais note : {r.design_iiptm:.3f}")
            continue
        delta = r.design_iiptm - float(prev)
        flag = "" if abs(delta) <= CONTROL_TOLERANCE else "  <-- HORS TOLERANCE"
        ok &= abs(delta) <= CONTROL_TOLERANCE
        print(f"  {float(prev):.3f} -> {r.design_iiptm:.3f}  ({delta:+.3f}){flag}")
    if not ok:
        print(f"\n  Le temoin derive de plus de {CONTROL_TOLERANCE}. La session a mesure")
        print("  autre chose que les sequences : ne pas lire la suite comme un resultat.")

    v = d[d.role == "variant"]
    arms = {a: g.design_iiptm.values for a, g in v.groupby("arm")}
    if set(arms) != {"periph", "milieu"}:
        raise SystemExit(f"bras inattendus : {sorted(arms)}")

    print("\n=== 2. les deux bras ===")
    print("  bras        n   median    max     effondrees   survivantes: median")
    for name in ("periph", "milieu"):
        s = pd.Series(arms[name])
        surv = s[s >= COLLAPSE]
        print(f"  {name:<10} {len(s):>2}   {s.median():.3f}   {s.max():.3f}   "
              f"{(s < COLLAPSE).sum():>2}/{len(s)} = {(s < COLLAPSE).mean():>4.0%}   "
              f"{surv.median() if len(surv) else float('nan'):.3f} (n={len(surv)})")

    print("\n=== 3. tous les tirages ===")
    u = mannwhitneyu(arms["periph"], arms["milieu"], alternative="two-sided")
    delta = cliffs_delta(arms["periph"], arms["milieu"])
    print(f"  Mann-Whitney bilateral : p = {u.pvalue:.3f}")
    print(f"  delta de Cliff (periph vs milieu) : {delta:+.2f}")

    print("\n=== 4. les survivantes seules ===")
    sp = [x for x in arms["periph"] if x >= COLLAPSE]
    sm = [x for x in arms["milieu"] if x >= COLLAPSE]
    if len(sp) >= 3 and len(sm) >= 3:
        us = mannwhitneyu(sp, sm, alternative="two-sided")
        print(f"  n = {len(sp)} contre {len(sm)}   p = {us.pvalue:.3f}   "
              f"delta de Cliff {cliffs_delta(sp, sm):+.2f}")
    else:
        print(f"  n = {len(sp)} contre {len(sm)} : trop peu pour comparer, non teste")

    print("\n=== VERDICT, selon la regle ecrite avant les donnees ===")
    if u.pvalue < 0.05 and abs(delta) > 0.4:
        better = "la peripherie" if delta > 0 else "le milieu"
        print(f"  Effet LARGE : {better} preserve mieux (p={u.pvalue:.3f}, "
              f"delta={delta:+.2f}).")
        print("  La position est un levier. Le generateur doit cesser de tirer ses")
        print("  positions au hasard, et le plan du §29 gagne le meme facteur.")
    else:
        print(f"  Pas d effet LARGE (p={u.pvalue:.3f}, delta={delta:+.2f}).")
        print("  Ce n est PAS 'pas d effet' : avec ces effectifs un effet modere")
        print("  passe inapercu. Deux suites possibles, au choix de l utilisateur :")
        print("    - doubler les effectifs (les variantes sont deja generees) ;")
        print("    - considerer la position comme close et passer au §29, test 2.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
