"""Comment les gagnants satisfont-ils la regle d unicite ?

Pour chaque gagnant d une epoque recente, sa distance a l archive TELLE QU ELLE
ETAIT : uniquement les soumissions d epoques strictement anterieures. Si les
gagnants se tiennent juste au-dessus du seuil, ils perturbent comme nous et le
§19 doit etre nuance. S ils sont loin de tout, ils concoivent vraiment.
"""
import sys
sys.path.insert(0, "/home/user/Nova")
import pandas as pd
from nova_lab.nanobodies import cdr_similarity, whole_identity

RECENT = 26170
df = pd.read_parquet("data/dataset/Q9NZQ7_nanobodies.parquet").dropna(subset=["epoch", "design_iiptm"])
df["L"] = df["cdr3"].str.len()
df = df.sort_values("epoch")

winners = df.loc[df.groupby("epoch")["design_iiptm"].idxmax()]
winners = winners[winners.epoch >= RECENT]
print(f"{len(winners)} gagnants d epoques >= {RECENT}\n")

rows = []
for w in winners.itertuples():
    prior = df[df.epoch < w.epoch]
    a = {"CDR1": w.cdr1, "CDR2": w.cdr2, "CDR3": w.cdr3}
    best = None
    for p in prior.itertuples():
        s = cdr_similarity(a, {"CDR1": p.cdr1, "CDR2": p.cdr2, "CDR3": p.cdr3})
        if best is None or s["CDR3"] > best[0]["CDR3"]:
            best = (s, p)
    s, p = best
    ident = whole_identity(w.sequence, p.sequence)
    tier3 = s["CDR3"] >= 0.80 and s["CDR1"] >= 0.90 and s["CDR2"] >= 0.90
    # Le voisin LUI-MEME, pas le gagnant de son epoque : une premiere version
    # avait compare chaque gagnant au gagnant de l epoque de son voisin, ce qui
    # ne mesure pas son pas. La sequence est donc stockee.
    rows.append(dict(epoch=int(w.epoch), L=w.L, iiptm=w.design_iiptm,
                     nn_sequence=p.sequence, nn_cdr3=p.cdr3, cdr3=w.cdr3,
                     cdr3_sim=s["CDR3"], cdr1=s["CDR1"], cdr2=s["CDR2"],
                     whole=ident, nn_iiptm=p.design_iiptm, nn_epoch=int(p.epoch),
                     tier3=tier3))
out = pd.DataFrame(rows)
out.to_csv("data/measurements/winners_distance.csv", index=False)
pd.set_option("display.width", 200)
print(out.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
print("\n=== resume ===")
print(f"similarite CDR3 au plus proche voisin anterieur : median {out.cdr3_sim.median():.3f}, "
      f"max {out.cdr3_sim.max():.3f}, min {out.cdr3_sim.min():.3f}")
print(f"identite entiere : median {out.whole.median():.3f}, max {out.whole.max():.3f}")
print(f"gagnants a >= 0,80 de similarite CDR3 (la bande ou nous perturbons) : "
      f"{(out.cdr3_sim >= 0.80).sum()}/{len(out)}")
print(f"gagnants a >= 0,70 : {(out.cdr3_sim >= 0.70).sum()}/{len(out)}")
print(f"gagnants a < 0,60 (loin de tout) : {(out.cdr3_sim < 0.60).sum()}/{len(out)}")
print(f"score du plus proche voisin : median {out.nn_iiptm.median():.3f}")
