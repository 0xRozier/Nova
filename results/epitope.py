import sys, numpy as np
from collections import defaultdict
path = sys.argv[1]
cols, rows, inloop = [], [], False
for line in open(path):
    if line.startswith('_atom_site.'):
        cols.append(line.strip().split('.')[1]); inloop = True; continue
    if inloop:
        if line.startswith('#') or line.startswith('loop_'): break
        if line.startswith('ATOM') or line.startswith('HETATM'):
            rows.append(line.split())
i = {c: n for n, c in enumerate(cols)}
atoms = defaultdict(list)
for r in rows:
    if r[i['group_PDB']] != 'ATOM': continue
    ch = r[i['auth_asym_id']]
    try: num = int(r[i['auth_seq_id']])
    except ValueError: continue
    xyz = [float(r[i[k]]) for k in ('Cartn_x','Cartn_y','Cartn_z')]
    atoms[(ch, num)].append((xyz, r[i['label_comp_id']]))

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
A  = {k: v for k, v in atoms.items() if k[0] == 'A'}
HL = {k: v for k, v in atoms.items() if k[0] in ('H', 'L')}
print(f"chaine A : {len(A)} residus, numerotation {min(k[1] for k in A)}..{max(k[1] for k in A)}")
for c in ('H','L'):
    n = [k[1] for k in atoms if k[0]==c]
    print(f"chaine {c} : {len(n)} residus, {min(n)}..{max(n)}")

bind = np.array([xyz for v in HL.values() for xyz, _ in v])
contacts = {}
for (ch, num), v in sorted(A.items(), key=lambda x: x[0][1]):
    pts = np.array([xyz for xyz, _ in v])
    d = np.linalg.norm(pts[:, None, :] - bind[None, :, :], axis=2).min()
    if d <= 5.0:
        contacts[num] = (AA3.get(v[0][1], 'X'), round(float(d), 2))

cfg = [24, 77, 80, 82, 131, 184, 185, 186]
print(f"\nepitope reel (<=5 A de H ou L) : {len(contacts)} residus")
print('  ' + ' '.join(f"{aa}{n}" for n, (aa, _) in contacts.items()))
print(f"\nsite de liaison de la config : {cfg}")
hit = [n for n in cfg if n in contacts]
print(f"  {len(hit)}/{len(cfg)} dans l'epitope reel : {hit}")
for n in cfg:
    if n in contacts:
        aa, d = contacts[n]; print(f"   {aa}{n:<4} contact  {d} A")
    elif ('A', n) in A:
        pts = np.array([xyz for xyz, _ in A[('A', n)]])
        d = np.linalg.norm(pts[:, None, :] - bind[None, :, :], axis=2).min()
        print(f"   {AA3.get(A[('A',n)][0][1],'X')}{n:<4} HORS     {d:.2f} A du Fab")
    else:
        print(f"   ?{n:<4} ABSENT de la chaine A")
