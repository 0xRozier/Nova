import numpy as np
from collections import defaultdict
AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
     'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
     'TYR':'Y','VAL':'V'}
cols, rows, inloop = [], [], False
for line in open("target4o9h/5FUC.cif"):
    if line.startswith('_atom_site.'):
        cols.append(line.strip().split('.')[1]); inloop = True; continue
    if inloop:
        if line.startswith('#') or line.startswith('loop_'): break
        if line.startswith(('ATOM','HETATM')): rows.append(line.split())
i = {c: n for n, c in enumerate(cols)}
at = defaultdict(list)
for r in rows:
    if r[i['group_PDB']] != 'ATOM': continue
    try: num = int(r[i['auth_seq_id']])
    except ValueError: continue
    at[(r[i['auth_asym_id']], num)].append(
        ([float(r[i[k]]) for k in ('Cartn_x','Cartn_y','Cartn_z')], r[i['label_comp_id']]))
chains = sorted({c for c, _ in at})
print("chaines 5FUC :", " ".join(f"{c}({len([1 for x,_ in at if x==c])})" for c in chains))
seq = "".join(AA3.get(at[('E', n)][0][1], 'X') for n in sorted(n for c, n in at if c == 'E'))
print(f"\nVHH6 (chaine E), {len(seq)} residus observes :\n{seq}")
pts = lambda cs: np.array([xyz for (c, _), v in at.items() if c in cs for xyz, _ in v])
E = pts('E')
for cs, name in (('AB', 'IL-6'), ('CD', 'IL-6R alpha')):
    p = pts(cs)
    if not len(p): continue
    d = np.linalg.norm(E[:, None, :] - p[None, :, :], axis=2)
    n = sum(1 for (c, num), v in at.items() if c == 'E' and
            np.linalg.norm(np.array([x for x, _ in v])[:, None, :] - p[None, :, :], axis=2).min() <= 5)
    print(f"  VHH6 <-> {name:<12} contact minimal {d.min():.2f} A, "
          f"{n} residus de VHH6 a l'interface")
