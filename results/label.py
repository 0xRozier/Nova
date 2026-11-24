import sys
path = sys.argv[1]
cols, rows, inloop = [], [], False
for line in open(path):
    if line.startswith('_atom_site.'):
        cols.append(line.strip().split('.')[1]); inloop = True; continue
    if inloop:
        if line.startswith('#') or line.startswith('loop_'): break
        if line.startswith(('ATOM','HETATM')): rows.append(line.split())
i = {c: n for n, c in enumerate(cols)}
seen = {}
for r in rows:
    if r[i['group_PDB']] != 'ATOM' or r[i['auth_asym_id']] != 'A': continue
    try: lab, auth = int(r[i['label_seq_id']]), int(r[i['auth_seq_id']])
    except ValueError: continue
    seen[lab] = (auth, r[i['label_comp_id']])
labs = sorted(seen)
print(f"chaine A : label_seq_id {labs[0]}..{labs[-1]} ({len(labs)} residus observes)")
offs = {lab - auth for lab, (auth, _) in seen.items()}
print(f"label - auth = {offs}   (un seul decalage = pas d'insertion)")
gaps = [(a, b) for a, b in zip(labs, labs[1:]) if b != a + 1]
print(f"trous dans label_seq_id : {gaps if gaps else 'aucun'}")
AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
     'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
     'TYR':'Y','VAL':'V'}
print("\nsite de liaison de la config, lu en label_seq_id :")
for n in (24, 77, 80, 82, 131, 184, 185, 186):
    if n in seen:
        auth, comp = seen[n]
        print(f"   res_index {n:<4} = {AA3.get(comp,'X')} (auth {auth})")
    else:
        print(f"   res_index {n:<4} NON OBSERVE")
lo, hi = 21, 186
inside = [n for n in labs if lo <= n <= hi]
print(f"\nres_index 21..186 : {len(inside)} residus observes sur {hi-lo+1} positions ; "
      f"hors de la fenetre : {[n for n in labs if not lo <= n <= hi]}")
