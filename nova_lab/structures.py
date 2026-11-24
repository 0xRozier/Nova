"""Reading mmCIF files: numbering, chains, contacts.

FINDINGS 73 was assembled from throwaway scripts, which meant the target
analysis it produced could not be replayed on the next target. This is that
work made reusable. It deliberately does no network and no antibody reasoning:
those belong to scripts.analyze_target.

The one thing to understand here is the numbering, because it is what decides
whether a binding site points at the right residues. An mmCIF carries two
numberings per polymer chain:

  label_seq_id   sequential from 1 over the entity's full SEQRES, gaps and all.
                 This is BoltzGen's res_index, and its README is explicit that
                 the author numbering is "incorrect" for the purpose.
  auth_seq_id    the depositors' numbering, usually UniProt-based, which can
                 start anywhere and skip.

On 4O9H they differ by exactly 2, which is why a binding site that looks 2/8
correct in author numbering is 8/8 correct in canonical numbering. Reading the
wrong one does not fail loudly -- it aims the design at a different patch of
the protein.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}
CONTACT_CUTOFF = 5.0      # angstroms, heavy-atom to heavy-atom
FR4_END = re.compile(r"W[GA][QKR]G[TA][A-Z]{0,2}VTVSS")


@dataclass(frozen=True)
class Residue:
    chain: str
    label: int          # label_seq_id: BoltzGen's res_index
    auth: int           # auth_seq_id: the depositors' numbering
    aa: str
    xyz: np.ndarray     # (n_atoms, 3)


class Structure:
    """The polymer residues of one mmCIF, both numberings kept."""

    def __init__(self, residues: list[Residue], path: Path | None = None):
        self.residues = residues
        self.path = path

    @classmethod
    def read(cls, path: str | Path) -> "Structure":
        """Parsed with Biopython, which handles mmCIF quoting we would not."""
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
        path = Path(path)
        d = MMCIF2Dict(str(path))
        need = ("_atom_site.group_PDB", "_atom_site.auth_asym_id",
                "_atom_site.auth_seq_id", "_atom_site.label_seq_id",
                "_atom_site.label_comp_id", "_atom_site.Cartn_x",
                "_atom_site.Cartn_y", "_atom_site.Cartn_z")
        missing = [k for k in need if k not in d]
        if missing:
            raise ValueError(f"{path.name}: champs atom_site absents {missing}")
        grouped: dict[tuple[str, int], list] = {}
        meta: dict[tuple[str, int], tuple[int, str]] = {}
        rows = zip(d["_atom_site.group_PDB"], d["_atom_site.auth_asym_id"],
                   d["_atom_site.auth_seq_id"], d["_atom_site.label_seq_id"],
                   d["_atom_site.label_comp_id"], d["_atom_site.Cartn_x"],
                   d["_atom_site.Cartn_y"], d["_atom_site.Cartn_z"])
        for group, ch, auth, label, comp, x, y, z in rows:
            # HETATM is ligands, waters and modified residues: not the polymer
            # whose numbering a res_index names.
            if group != "ATOM" or label in (".", "?"):
                continue
            try:
                key = (ch, int(label))
                auth_i = int(auth)
            except ValueError:
                continue
            grouped.setdefault(key, []).append((float(x), float(y), float(z)))
            meta.setdefault(key, (auth_i, AA3.get(comp, "X")))
        residues = [Residue(ch, label, *meta[(ch, label)], xyz=np.array(pts))
                    for (ch, label), pts in grouped.items()]
        residues.sort(key=lambda r: (r.chain, r.label))
        return cls(residues, path)

    # -- chains -------------------------------------------------------------

    @property
    def chains(self) -> list[str]:
        return sorted({r.chain for r in self.residues})

    def chain(self, name: str) -> list[Residue]:
        return [r for r in self.residues if r.chain == name]

    def sequence(self, name: str) -> str:
        """The OBSERVED sequence: disordered residues are simply absent."""
        return "".join(r.aa for r in self.chain(name))

    def labels(self, name: str) -> list[int]:
        return [r.label for r in self.chain(name)]

    # -- numbering ----------------------------------------------------------

    def label_auth_offsets(self, name: str) -> set[int]:
        """{label - auth}. More than one value means an insertion or a renumber."""
        return {r.label - r.auth for r in self.chain(name)}

    def gaps(self, name: str) -> list[tuple[int, int]]:
        """Breaks in label_seq_id: disordered loops the file does not contain."""
        labels = self.labels(name)
        return [(a, b) for a, b in zip(labels, labels[1:]) if b != a + 1]

    # -- geometry -----------------------------------------------------------

    def contacts(self, target: str, partners: list[str] | str,
                 cutoff: float = CONTACT_CUTOFF) -> dict[int, float]:
        """{label_seq_id of `target`: nearest heavy-atom distance to `partners`}.

        Keyed by label_seq_id, i.e. in res_index space, so the result compares
        directly against a config binding site.
        """
        partners = list(partners)
        cloud = np.vstack([r.xyz for r in self.residues if r.chain in partners]) \
            if any(r.chain in partners for r in self.residues) else np.empty((0, 3))
        if not len(cloud):
            return {}
        out = {}
        for r in self.chain(target):
            d = float(np.linalg.norm(r.xyz[:, None, :] - cloud[None, :, :], axis=2).min())
            if d <= cutoff:
                out[r.label] = round(d, 2)
        return out


def parse_res_index(spec: str) -> list[int]:
    """"21..186" or "24,77,80,184..186" -> the residue numbers it names.

    This is the config's own syntax, shared by structure_res_index and
    structure_binding_site, and it is read as label_seq_id.
    """
    out: list[int] = []
    for piece in str(spec).split(","):
        piece = piece.strip()
        if not piece:
            continue
        if ".." in piece:
            lo, hi = (int(x) for x in piece.split("..", 1))
            out.extend(range(lo, hi + 1))
        else:
            out.append(int(piece))
    return out


def entity_sequences(path: str | Path) -> list[str]:
    """Every polymer SEQRES in the file, construct tags and all.

    SEQRES rather than the observed sequence, because a binder's framework must
    be complete to be numbered or mutated -- a disordered loop would silently
    shorten a CDR.
    """
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    d = MMCIF2Dict(str(path))
    key = "_entity_poly.pdbx_seq_one_letter_code_can"
    if key not in d:
        return []
    values = d[key] if isinstance(d[key], list) else [d[key]]
    return [re.sub(r"[^A-Z]", "", v.upper()) for v in values]


def polymers(path: str | Path) -> list[dict]:
    """[{"sequence": SEQRES, "chains": [auth chain ids]}] for each polymer entity.

    The entity-to-chain map is in the file; guessing it from observed sequence
    prefixes fails whenever a chain's first residues are disordered, which is
    most Fab heavy chains. label_seq_id counts from 1 over the SEQRES returned
    here, which is what makes a cross-entry numbering offset computable.
    """
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    d = MMCIF2Dict(str(path))
    seq_key = "_entity_poly.pdbx_seq_one_letter_code_can"
    chain_key = "_entity_poly.pdbx_strand_id"
    if seq_key not in d or chain_key not in d:
        return []
    seqs = d[seq_key] if isinstance(d[seq_key], list) else [d[seq_key]]
    strands = d[chain_key] if isinstance(d[chain_key], list) else [d[chain_key]]
    out = []
    for seq, strand in zip(seqs, strands):
        out.append({"sequence": re.sub(r"[^A-Z]", "", seq.upper()),
                    "chains": [c.strip() for c in str(strand).split(",") if c.strip()]})
    return out


def best_offset(sequence: str, reference: str, span: int = 60) -> tuple[int, float]:
    """(k, agreement) such that `sequence[i]` lines up with `reference[i + k]`.

    The same protein is deposited as different constructs -- IL-6 appears at
    186, 183, 171, 168 and 166 residues -- so label_seq_id is not comparable
    across entries and an epitope must be shifted before it can be compared
    against a config site. Plain offset search, no gaps: two constructs of one
    protein differ at their ends, not in the middle.
    """
    best = (0, -1.0)
    if not sequence or not reference:
        return best
    for k in range(-span, span + 1):
        pairs = [(a, reference[i + k]) for i, a in enumerate(sequence)
                 if 0 <= i + k < len(reference)]
        if len(pairs) < min(len(sequence), 40):
            continue
        frac = sum(1 for a, b in pairs if a == b) / len(pairs)
        if frac > best[1]:
            best = (k, frac)
    return best


def vh_domain(sequence: str) -> str | None:
    """A heavy variable domain cut at the end of FR4, or None.

    A Fab heavy chain in the PDB runs on into CH1 for another hundred residues;
    everything downstream here wants the variable domain alone.
    """
    end = FR4_END.search(sequence)
    if not end or not sequence[:1] in set("EQDGA"):
        return None
    vh = sequence[:end.end()]
    return vh if 90 <= len(vh) <= 150 else None
