import sys, numpy as np
from collections import defaultdict
AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
     'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
     'TYR':'Y','VAL':'V'}
REF = ("MAPVPPGEDSKDVAAPHRQPLTSSERIDKQIRYILDGISALRKETCNKSNMCESSKEALAENNLNLPKMAEKDGCFQSGF"
       "NEETCLVKIITGLLEFEVYLEYLQNRFESSEEQARAVQMSTKVLIQFLQKKAKNLDAITTPDPTTNASLLTKLQAQNQWL"
       "QDMTTHLILRSFKEFLQSSLRALRQM")           # P05231 mature-construct, index i = res_index i

def load(path):
    cols, rows, inloop = [], [], False
    for line in open(path):
        if line.startswith('_atom_site.'):
            cols.append(line.strip().split('.')[1]); inloop = True; continue
        if inloop:
            if line.startswith('#') or line.startswith('loop_'): break
            if line.startswith(('ATOM','HETATM')): rows.append(line.split())
    i = {c: n for n, c in enumerate(cols)}
    out = defaultdict(list)
    for r in rows:
        if r[i['group_PDB']] != 'ATOM': continue
        try: num = int(r[i['auth_seq_id']])
        except ValueError: continue
        out[(r[i['auth_asym_id']], num)].append(
            ([float(r[i[k]]) for k in ('Cartn_x','Cartn_y','Cartn_z')], r[i['label_comp_id']]))
    return out

def offset(atoms, chain):
    """auth -> res_index shift, chosen by sequence agreement with REF."""
    obs = {n: AA3.get(v[0][1], 'X') for (c, n), v in atoms.items() if c == chain}
    best = (0, -1)
    for off in range(-30, 31):
        ok = sum(1 for n, aa in obs.items()
                 if 1 <= n + off <= len(REF) and REF[n + off - 1] == aa)
        if ok > best[1]: best = (off, ok)
    return best[0], best[1], len(obs)

CASES = [("4O9H", "A", "HL",  "Fab camelide 61H7"),
         ("4CNI", "C", "AB",  "Fab Olokizumab (humanise)"),
         ("4ZS7", "A", "HL",  "Fab lama 68F2"),
         ("5FUC", "A", "EV",  "VHH6 (domaine unique)"),
         ("8YWQ", "A", "HL",  "Fab I9H"),
         ("8YWR", "A", "HL",  "Fab 68F2 (2e cristal)")]
CFG = {24, 77, 80, 82, 131, 184, 185, 186}
table = {}
for pdb, tgt, binders, label in CASES:
    atoms = load(f"target4o9h/{pdb}.cif")
    off, ok, n = offset(atoms, tgt)
    pts = np.array([xyz for (c, _), v in atoms.items() if c in binders for xyz, _ in v])
    if not len(pts): print(f"{pdb}: aucune chaine liante {binders}"); continue
    epi = set()
    for (c, num), v in atoms.items():
        if c != tgt: continue
        p = np.array([xyz for xyz, _ in v])
        if np.linalg.norm(p[:, None, :] - pts[None, :, :], axis=2).min() <= 5.0:
            epi.add(num + off)
    table[pdb] = (epi, label, off, ok, n)
    inter = sorted(CFG & epi)
    print(f"{pdb}  {label:<28} auth{off:+d} ({ok}/{n} seq ok)  "
          f"epitope {len(epi):2d} res  |  recouvre {len(inter)}/8 de la config {inter}")

print("\n--- residus de la config, par structure (res_index) ---")
print(f"{'res':>5} " + " ".join(f"{p:>5}" for p in table))
for r in sorted(CFG):
    print(f"{REF[r-1]}{r:<4} " + " ".join(f"{'  X  ' if r in table[p][0] else '  .  '}" for p in table))
print("\n--- consensus : residus touches par >=4 des 6 liants ---")
cnt = defaultdict(int)
for epi, *_ in table.values():
    for r in epi: cnt[r] += 1
cons = sorted(r for r, c in cnt.items() if c >= 4)
print(f"{len(cons)} residus : " + " ".join(f"{REF[r-1]}{r}" for r in cons))
print(f"  dont dans la config : {sorted(set(cons) & CFG)}")
print(f"  absents de la config : {sorted(set(cons) - CFG)}")
