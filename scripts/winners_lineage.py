import sys, random
sys.path.insert(0, "/home/user/Nova")
import pandas as pd, statistics
from nova_lab.nanobodies import cdr_similarity

RECENT = 26170
df = pd.read_parquet("data/dataset/Q9NZQ7_nanobodies.parquet").dropna(subset=["epoch","design_iiptm"])
df["L"] = df["cdr3"].str.len()
w = pd.read_csv("data/measurements/winners_distance.csv")

print("=== CONTROLE : deux soumissions quelconques de MEME longueur ===")
rng = random.Random(0)
base = {}
for L in sorted(set(w.L)):
    sub = df[(df.L == L)].reset_index(drop=True)
    if len(sub) < 10: continue
    sims = []
    for _ in range(400):
        i, j = rng.sample(range(len(sub)), 2)
        a, b = sub.iloc[i], sub.iloc[j]
        sims.append(cdr_similarity({"CDR1":a.cdr1,"CDR2":a.cdr2,"CDR3":a.cdr3},
                                   {"CDR1":b.cdr1,"CDR2":b.cdr2,"CDR3":b.cdr3})["CDR3"])
    base[L] = statistics.median(sims)
m = w.merge(pd.Series(base, name="base_sim").rename_axis("L").reset_index(), on="L")
print(f"  similarite CDR3 entre deux soumissions au hasard : median {statistics.median(base.values()):.3f}")
print(f"  similarite CDR3 gagnant -> son plus proche anterieur : median {m.cdr3_sim.median():.3f}")
m["excess"] = m.cdr3_sim - m.base_sim
print(f"  ecart (gagnant - hasard), par longueur : median {m.excess.median():+.3f}, "
      f"{(m.excess > 0).sum()}/{len(m)} au-dessus du hasard")

print("\n=== LIGNEE : le plus proche anterieur d un gagnant est-il lui-meme un gagnant ? ===")
winner_epochs = set(w.epoch)
w["nn_is_winner"] = w.nn_epoch.isin(winner_epochs)
print(f"  {w.nn_is_winner.sum()}/{len(w)} des gagnants ont pour plus proche voisin un gagnant anterieur")
w["gain"] = w.iiptm - w.nn_iiptm
print(f"  gain sur le plus proche voisin : median {w.gain.median():+.3f}, "
      f"{(w.gain > 0).sum()}/{len(w)} en hausse")
chains = w[w.nn_is_winner]
print(f"  parmi les chaines gagnant->gagnant (n={len(chains)}) : gain median {chains.gain.median():+.3f}")

print("\n=== LA LIGNEE A CDR3 = 18, epoque par epoque ===")
for L in (18, 16, 15):
    sub = w[w.L == L].sort_values("epoch")
    if sub.empty: continue
    print(f"  L={L} : " + "  ".join(f"{int(r.epoch)}:{r.iiptm:.3f}" for r in sub.itertuples()))

print("\n=== le champ monte-t-il ? max par epoque, par tranche ===")
g = df.groupby("epoch")["design_iiptm"].max()
for lo, hi in [(26170,26190),(26190,26210),(26210,26235)]:
    s = g[(g.index>=lo)&(g.index<hi)]
    print(f"  epoques {lo}-{hi} : max median {s.median():.3f} (n={len(s)})")
