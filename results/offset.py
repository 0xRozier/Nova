import sys, numpy as np
from collections import defaultdict
path = sys.argv[1]
cols, rows, inloop = [], [], False
for line in open(path):
    if line.startswith('_atom_site.'):
        cols.append(line.strip().split('.')[1]); inloop = True; continue
    if inloop:
        if line.startswith('#') or line.startswith('loop_'): break
        if line.startswith(('ATOM','HETATM')): rows.append(line.split())
i = {c: n for n, c in enumerate(cols)}
atoms = defaultdict(list)
for r in rows:
    if r[i['group_PDB']] != 'ATOM': continue
    try: num = int(r[i['auth_seq_id']])
    except ValueError: continue
    atoms[(r[i['auth_asym_id']], num)].append(
        ([float(r[i[k]]) for k in ('Cartn_x','Cartn_y','Cartn_z')], r[i['label_comp_id']]))
AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V'}
A = {k[1]: v for k, v in atoms.items() if k[0] == 'A'}
by = {c: np.array([xyz for k, v in atoms.items() if k[0] == c for xyz, _ in v]) for c in 'HL'}
both = np.vstack([by['H'], by['L']])

def mind(num, pts):
    p = np.array([xyz for xyz, _ in A[num]])
    return float(np.linalg.norm(p[:, None, :] - pts[None, :, :], axis=2).min())

epi = sorted(n for n in A if mind(n, both) <= 5.0)
cfg = [24, 77, 80, 82, 131, 184, 185, 186]
print("decalage  recouvrement config -> epitope reel")
for off in range(-5, 6):
    hit = [n for n in cfg if (n + off) in epi]
    bar = '#' * len(hit)
    print(f"  {off:+d}   {len(hit)}/8  {bar}")
print()
off = -2
print(f"avec decalage {off:+d} (config = numerotation UniProt P05231, PDB = UniProt-2) :")
for n in cfg:
    m = n + off
    dh, dl = mind(m, by['H']), mind(m, by['L'])
    who = 'H seule' if dh <= 5 < dl else ('L seule' if dl <= 5 < dh else ('H + L' if max(dh,dl) <= 5 else '-'))
    print(f"   config {n:<4} -> PDB {AA3.get(A[m][0][1],'X')}{m:<4} "
          f"dH={dh:5.2f} dL={dl:5.2f}  {who}")
print(f"\nepitope complet, repartition des {len(epi)} residus :")
tally = {'H seule': [], 'L seule': [], 'H + L': []}
for n in epi:
    dh, dl = mind(n, by['H']), mind(n, by['L'])
    k = 'H seule' if dh <= 5 < dl else ('L seule' if dl <= 5 < dh else 'H + L')
    tally[k].append(f"{AA3.get(A[n][0][1],'X')}{n}")
for k, v in tally.items():
    print(f"   {k:<8} {len(v):2d}  {' '.join(v)}")
