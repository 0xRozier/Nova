"""The target analysis the strategy guide puts first, as code.

The upstream guide's section 5, "Study the Target", asks for three things in
the ~24 hours before an epoch: download and inspect the structure, identify
the binding site, and look for known binders. FINDINGS 73 did that for
P05231/4O9H by hand and it produced the only usable lever in weeks -- after
four GPU-measured dead ends. It was done by hand, so it was not replayable on
the next target. This is it as a script.

What it answers, in order:

  1. NUMBERING. Does the config's structure_res_index match the chain the
     structure actually contains, read as label_seq_id -- BoltzGen's res_index?
     A config read in author numbering does not fail loudly, it aims the design
     at a different patch. On 4O9H the two differ by 2.
  2. BINDING SITE. Is the config's structure_binding_site really in contact
     with whatever the structure co-crystallises? And is it reachable by a
     single domain, or does it include a light chain's share?
  3. KNOWN BINDERS. What else in the PDB binds this UniProt entry, how much of
     the config site does each cover, and -- the part that matters -- does any
     of them clear the validator's gates? FINDINGS 73's answer was that the one
     true VHH against IL-6 is refused on humanness while a camelid Fab's heavy
     chain passes, which is why the strategy points at humanised frameworks.
  4. PARATOPE. For each binder, how much of its contact surface is CDR3? That
     decides whether "keep the framework, redesign the CDR3" is even coherent:
     61H7 is 25% CDR3, the other four 40-62%.
  5. SUBMITTABILITY. Three rules in the live config refuse raw 61H7 --
     enforce_vhh_hallmarks, min_human_v_identity, min_human_j_identity. Where
     a binder fails only the hallmark rule or the J floor, this proposes the
     minimal fix and reports whether it leaves every CDR untouched.

It writes a markdown report and a CSV, and makes no claim it has not measured:
every number below comes from the structure files or from the validator's own
nativeness filter, never from a model of them.

Usage:
    python -m scripts.analyze_target                      # the config's target
    python -m scripts.analyze_target --no-survey          # skip the PDB search
    python -m scripts.analyze_target --target P05231 --structure 4O9H
"""
from __future__ import annotations

import argparse
import json
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

from nova_lab import REPO_ROOT
from nova_lab.config import (MIN_HUMAN_J_IDENTITY, MIN_HUMAN_V_IDENTITY,
                             NANOBODY_BINDING_SITE, NANOBODY_STRUCTURE,
                             NANOBODY_STRUCTURE_CHAIN, NANOBODY_STRUCTURE_RES_INDEX,
                             NANOBODY_TARGET)
from nova_lab.nanobodies import check_sequence, germline_identity, is_vhh, regions
from nova_lab.structures import (Structure, best_offset, entity_sequences,
                                 parse_res_index, polymers, vh_domain)

UPSTREAM_STRUCTURES = Path("/home/user/nova-upstream/data/structures")
RCSB_FILE = "https://files.rcsb.org/download/{pdb}.cif"
RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query?json="
RCSB_GRAPHQL = "https://data.rcsb.org/graphql"
# The FR4 of human IGHJ4. Substituted wholesale when a binder fails only the J
# floor: it is the same length as the FR4 it replaces in every case seen, and
# it touches no CDR.
HUMAN_FR4 = "WGQGTLVTVSS"
BINDER_WORDS = ("nanobody", "vhh", "single-domain", "single domain", "fab",
                "antibody", "immunoglobulin", "heavy chain", "light chain",
                "variable domain", "variable lymphocyte")


# --------------------------------------------------------------------------- 1

def structure_path(pdb: str, cache: Path) -> Path:
    """Upstream's copy if it has one, else fetched from the PDB into `cache`.

    Upstream is preferred because it is the file the validators score against;
    a fetched copy is a different revision of the same entry in principle.
    """
    local = UPSTREAM_STRUCTURES / f"{pdb}.cif"
    if local.exists():
        return local
    cache.mkdir(parents=True, exist_ok=True)
    out = cache / f"{pdb}.cif"
    if not out.exists():
        with urllib.request.urlopen(RCSB_FILE.format(pdb=pdb), timeout=120) as r:
            out.write_bytes(r.read())
    return out


def numbering_audit(st: Structure, chain: str, res_index: str) -> dict:
    """Does the config's window match the chain the file contains?"""
    labels = st.labels(chain)
    if not labels:
        return {"ok": False, "why": f"la chaine {chain} n'est pas dans le fichier"}
    window = parse_res_index(res_index)
    outside = [n for n in labels if n not in set(window)]
    offsets = st.label_auth_offsets(chain)
    return {"ok": not outside,
            "why": "" if not outside else
                   f"{len(outside)} residus observes hors de la fenetre : {outside[:8]}",
            "observed": len(labels), "span": (labels[0], labels[-1]),
            "window": (min(window), max(window)), "window_size": len(window),
            "label_minus_auth": sorted(offsets),
            "single_offset": len(offsets) == 1,
            "gaps": st.gaps(chain)}


# --------------------------------------------------------------------------- 2

def epitope_check(st: Structure, chain: str, binding_site: str) -> dict:
    """The config site against every contact the structure actually makes."""
    others = [c for c in st.chains if c != chain]
    site = parse_res_index(binding_site)
    epitope = st.contacts(chain, others)
    per_chain = {c: set(st.contacts(chain, [c])) for c in others}
    inside = [n for n in site if n in epitope]
    missing = [n for n in site if n not in epitope]
    absent = [n for n in site if n not in set(st.labels(chain))]
    # Which partner chains reach each site residue: a residue only the light
    # chain of a Fab touches is not reachable by a single domain.
    reach = {n: sorted(c for c, s in per_chain.items() if n in s) for n in inside}
    return {"partners": others, "epitope": epitope, "site": site,
            "inside": inside, "missing": missing, "absent": absent,
            "per_chain": {c: len(s) for c, s in per_chain.items()},
            "reach": reach,
            "single_chain_covers": {c: sum(1 for n in inside if c in reach[n])
                                    for c in others}}


# --------------------------------------------------------------------------- 3

def survey(accession: str, limit: int = 100, attempts: int = 3) -> list[dict]:
    """Every PDB entry carrying this UniProt accession, with its polymers.

    Retried, because the RCSB search allows only a few queries a minute and a
    throttled call comes back as an empty result set rather than an error --
    which would be reported as "this protein has no known binders", the single
    most misleading thing this script could say.
    """
    import time
    for attempt in range(attempts):
        found = _survey_once(accession, limit)
        if found:
            return found
        if attempt + 1 < attempts:
            time.sleep(20)
    return []


def _survey_once(accession: str, limit: int) -> list[dict]:
    query = {"query": {"type": "terminal", "service": "text", "parameters": {
        "attribute": "rcsb_polymer_entity_container_identifiers."
                     "reference_sequence_identifiers.database_accession",
        "operator": "exact_match", "value": accession}},
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": limit}}}
    url = RCSB_SEARCH + urllib.parse.quote(json.dumps(query))
    with urllib.request.urlopen(url, timeout=120) as r:
        ids = [x["identifier"] for x in json.load(r).get("result_set", [])]
    if not ids:
        return []
    gql = ("{entries(entry_ids:[%s]){rcsb_id struct{title} polymer_entities{"
           "rcsb_polymer_entity{pdbx_description} "
           "rcsb_polymer_entity_container_identifiers{auth_asym_ids}}}}"
           % ",".join(f'"{i}"' for i in ids))
    req = urllib.request.Request(RCSB_GRAPHQL, data=json.dumps({"query": gql}).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        entries = json.load(r)["data"]["entries"]
    out = []
    for e in entries:
        polymers = [{"description": (p["rcsb_polymer_entity"]["pdbx_description"] or "?"),
                     "chains": p["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"]}
                    for p in e["polymer_entities"]]
        blob = (e["struct"]["title"] + " " + " ".join(p["description"] for p in polymers)).lower()
        out.append({"pdb": e["rcsb_id"], "title": e["struct"]["title"],
                    "polymers": polymers,
                    "has_binder": any(w in blob for w in BINDER_WORDS)})
    return sorted(out, key=lambda e: e["pdb"])


def classify_chains(path: Path, reference: str) -> dict:
    """Which chains are the target, which are heavy variable domains.

    Both come from the file's own entity-to-chain map, so nothing depends on a
    chain's first residues being present. `offset` shifts this entry's
    label_seq_id into the reference's numbering.
    """
    target_chains: list[str] = []
    binder: list[str] = []
    offset, agreement = 0, 0.0
    for entity in polymers(path):
        seq = entity["sequence"]
        k, frac = best_offset(seq, reference)
        if frac >= 0.8 and len(seq) >= 40:
            if frac > agreement:
                offset, agreement = k, frac
            target_chains += entity["chains"]
        elif vh_domain(seq):
            binder += entity["chains"]
    return {"target": sorted(set(target_chains)), "binder": sorted(set(binder)),
            "offset": offset, "agreement": round(agreement, 3)}


# --------------------------------------------------------------------------- 4

def paratope_split(st: Structure, binder: str, targets: list[str], vh: str) -> dict:
    """How a binder's contacts divide between CDR1, CDR2, CDR3 and framework.

    The regions come from the complete SEQRES-derived VH, then are located in
    the OBSERVED chain by substring search. A CDR that is partly disordered
    will not be found, and this says so rather than guessing: a wrong region
    assignment would mislead exactly the decision this is here to inform.
    """
    r = regions(vh)
    if not r:
        return {"ok": False, "why": "VH non numerotable"}
    observed = st.sequence(binder)
    labels = st.labels(binder)
    assign: dict[int, str] = {}
    cursor, unplaced = 0, []
    for name in ("CDR1", "CDR2", "CDR3"):
        at = observed.find(r[name], cursor)
        if at < 0:
            unplaced.append(name)
            continue
        for j in range(at, at + len(r[name])):
            assign[labels[j]] = name
        cursor = at + len(r[name])
    contacts = st.contacts(binder, targets)
    split: dict[str, list[int]] = {"CDR1": [], "CDR2": [], "CDR3": [], "framework": []}
    for n in contacts:
        split[assign.get(n, "framework")].append(n)
    total = len(contacts)
    return {"ok": not unplaced, "why": f"regions introuvables dans l'observe : {unplaced}"
            if unplaced else "", "total": total,
            "counts": {k: len(v) for k, v in split.items()},
            "share": {k: (len(v) / total if total else 0.0) for k, v in split.items()},
            "residues": split, "cdr3": r["CDR3"]}


# --------------------------------------------------------------------------- 5

def assess(vh: str) -> dict:
    """Every gate we can check offline, plus the validator's own filter."""
    g = germline_identity(vh)
    hall, why = is_vhh(vh)
    row = {"length": len(vh), "cysteines": vh.count("C"),
           "hallmarks": hall, "hallmark_why": why,
           "v_gene": g["V"]["genes"][0] if g else None,
           "v_identity": round(g["V"]["identity_pct"], 1) if g else None,
           "j_identity": round(g["J"]["identity_pct"], 1) if g else None,
           "layer1": ", ".join(check_sequence(vh)) or "ok"}
    from scripts import nativeness
    usable, reason = nativeness.available()
    if usable:
        f = nativeness.verdict(nativeness.score([vh])).iloc[0]
        row |= {"nativeness": round(float(f.nativeness), 3),
                "humanness": round(float(f.humanness), 3) if pd.notna(f.humanness) else None,
                "hard_reject": bool(f.hard_reject),
                "passes_humanness": bool(f.passes_humanness)}
    else:
        # Never scored as a refusal: FINDINGS 69 and 73 were both caused by a
        # missing piece being read as one.
        row |= {"nativeness": None, "humanness": None, "hard_reject": None,
                "passes_humanness": None, "filter_unavailable": reason}
    return row


def make_submittable(vh: str) -> dict | None:
    """The minimal edits that clear the hallmark rule and the J floor, or None.

    Only FR2 positions 49/50 and FR4 are touched, never a CDR, and the result
    is returned only if check_sequence then passes and every CDR comes out
    byte-identical. On 61H7 this is three mutations: G49E, L50R, FR4 Q->L.
    """
    from abnumber import Chain
    before = regions(vh)
    if not before:
        return None
    seq, muts = vh, []
    if not is_vhh(seq)[0]:
        c = Chain(seq, scheme="imgt")
        want = {"49": "E", "50": "R"}
        rebuilt, changed = [], []
        for p, aa in zip(c.positions, c.seq):
            target = want.get(str(p.number))
            if target and aa != target and not (str(p.number) == "49" and aa in "EQ"):
                rebuilt.append(target)
                changed.append(f"{aa}{p.number}{target}")
            else:
                rebuilt.append(aa)
        seq, muts = "".join(rebuilt), muts + changed
    g = germline_identity(seq)
    if g and g["J"]["identity_pct"] < MIN_HUMAN_J_IDENTITY:
        fr4 = regions(seq)["FR4"]
        if len(fr4) == len(HUMAN_FR4) and seq.count(fr4) == 1 and fr4 != HUMAN_FR4:
            candidate = seq.replace(fr4, HUMAN_FR4)
            after = germline_identity(candidate)
            if after and after["J"]["identity_pct"] >= MIN_HUMAN_J_IDENTITY:
                seq = candidate
                muts.append(f"FR4 {fr4}->{HUMAN_FR4}")
    if seq == vh or not muts:
        return None
    after = regions(seq)
    if not after or any(before[k] != after[k] for k in ("CDR1", "CDR2", "CDR3")):
        return None
    problems = check_sequence(seq)
    return {"sequence": seq, "mutations": muts, "n_mutations":
            sum(1 for a, b in zip(vh, seq) if a != b) if len(vh) == len(seq) else None,
            "layer1": ", ".join(problems) or "ok", "submittable": not problems,
            "cdrs_untouched": True}


# --------------------------------------------------------------------------- report

def analyse(target: str, pdb: str, chain: str, res_index: str, binding_site: str,
            cache: Path, do_survey: bool = True) -> dict:
    st = Structure.read(structure_path(pdb, cache))
    numbering = numbering_audit(st, chain, res_index)
    epitope = epitope_check(st, chain, binding_site)
    reference = st.sequence(chain)
    # The reference is the entity that carries the config's chain -- not the
    # longest in the file, which on 4O9H is the Fab heavy chain and would make
    # every cross-entry alignment fail.
    own = structure_path(pdb, cache)
    full_reference = next((e["sequence"] for e in polymers(own)
                           if chain in e["chains"]), reference)

    binders: list[dict] = []
    entries = survey(target) if do_survey else []
    survey_note = ("" if entries or not do_survey else
                   "le releve PDB est revenu vide apres plusieurs essais : "
                   "a lire comme un echec reseau ou une limitation de debit, "
                   "PAS comme une absence de liant connu")
    site = set(parse_res_index(binding_site))
    for entry in entries:
        if not entry["has_binder"]:
            continue
        try:
            path = structure_path(entry["pdb"], cache)
            other = Structure.read(path)
        except Exception as exc:
            binders.append({"pdb": entry["pdb"], "title": entry["title"],
                            "error": f"{type(exc).__name__}: {exc}"})
            continue
        kinds = classify_chains(path, full_reference)
        if not kinds["target"]:
            binders.append({"pdb": entry["pdb"], "title": entry["title"],
                            "error": "chaine cible non identifiee"})
            continue
        if not kinds["binder"]:
            binders.append({"pdb": entry["pdb"], "title": entry["title"],
                            "error": "aucun domaine variable lourd"})
            continue
        tgt = kinds["target"][0]
        # Contacts of the whole binder, shifted into the reference numbering.
        reached = {n + kinds["offset"] for n in other.contacts(tgt, kinds["binder"])}
        overlap = sorted(site & reached)
        for entity in polymers(path):
            vh = vh_domain(entity["sequence"])
            if not vh or not set(entity["chains"]) & set(kinds["binder"]):
                continue
            bch = sorted(set(entity["chains"]) & set(kinds["binder"]))[0]
            row = {"pdb": entry["pdb"], "title": entry["title"], "chain": bch,
                   "sequence": vh, "overlap": len(overlap),
                   "overlap_residues": overlap, "site_size": len(site),
                   "label_offset": kinds["offset"], "alignment": kinds["agreement"]}
            row |= assess(vh)
            row["paratope"] = paratope_split(other, bch, kinds["target"], vh)
            row["fix"] = make_submittable(vh)
            binders.append(row)
    return {"target": target, "pdb": pdb, "chain": chain,
            "numbering": numbering, "epitope": epitope,
            "entries": entries, "binders": binders, "survey_note": survey_note,
            "surveyed": do_survey, "chains": st.chains, "observed": reference}


def report(a: dict) -> str:
    n, e = a["numbering"], a["epitope"]
    L = [f"# Analyse de cible : {a['target']} / {a['pdb']}", "",
         f"Chaines du fichier : {', '.join(a['chains'])}. "
         f"Chaine cible : {a['chain']}.", "",
         "## 1. Numerotation", ""]
    L += [f"- residus observes : {n['observed']}, label_seq_id "
          f"{n['span'][0]}..{n['span'][1]}",
          f"- fenetre de la config : {n['window'][0]}..{n['window'][1]} "
          f"({n['window_size']} positions)",
          f"- label - auth : {n['label_minus_auth']}"
          f"{' (decalage unique)' if n['single_offset'] else ' (PLUSIEURS : insertion ou renumerotation)'}",
          f"- trous (boucles desordonnees) : {n['gaps'] or 'aucun'}",
          f"- **{'la fenetre couvre tout l observe' if n['ok'] else 'ECART : ' + n['why']}**", ""]
    L += ["## 2. Site de liaison", "",
          f"Partenaires co-cristallises : {', '.join(e['partners'])} "
          f"({', '.join(f'{c} touche {k} residus' for c, k in e['per_chain'].items())}).",
          f"Epitope reel : {len(e['epitope'])} residus. "
          f"Site de la config : {len(e['site'])} residus.", "",
          f"**{len(e['inside'])}/{len(e['site'])} du site sont au contact.**", ""]
    if e["missing"]:
        L += [f"Hors contact : {e['missing']}. "
              f"{'Absents de la chaine : ' + str(e['absent']) if e['absent'] else ''}", ""]
    L += ["| res_index | distance | atteint par |", "|---|---|---|"]
    for res in e["site"]:
        d = e["epitope"].get(res)
        L.append(f"| {res} | {str(d) + ' A' if d else 'hors contact'} | "
                 f"{', '.join(e['reach'].get(res, [])) or '-'} |")
    covers = ", ".join(f"{c} atteint {k}/{len(e['inside'])}"
                       for c, k in e["single_chain_covers"].items())
    L += ["", f"Couverture par chaine partenaire : {covers}.",
          "Un site qu'une seule chaine couvre entierement est atteignable par un domaine unique.", ""]
    usable = [b for b in a["binders"] if "error" not in b]
    L += ["## 3. Liants connus", ""]
    if a.get("survey_note"):
        L += [f"> **{a['survey_note']}**", ""]
    elif not a.get("surveyed"):
        L += ["> Releve PDB non effectue (`--no-survey`).", ""]
    L += [
          f"{len(a['entries'])} entrees PDB portent {a['target']}, "
          f"{sum(1 for x in a['entries'] if x['has_binder'])} avec un liant, "
          f"{len(usable)} domaines variables lisibles.", ""]
    if usable:
        L += ["| PDB | chaine | recouvrement | hallmarks | V | J | nativite | humanite | couche 1 | % CDR3 |",
              "|---|---|---|---|---|---|---|---|---|---|"]
        for b in sorted(usable, key=lambda x: -x["overlap"]):
            p = b["paratope"]
            L.append(f"| {b['pdb']} | {b['chain']} | **{b['overlap']}/{b['site_size']}** | "
                     f"{'oui' if b['hallmarks'] else 'non'} | {b['v_identity']} | "
                     f"{b['j_identity']} | {b['nativeness']} | {b['humanness']} | "
                     f"{b['layer1'][:28]} | "
                     f"{100 * p['share']['CDR3']:.0f}% |" if p["ok"] else
                     f"{b['layer1'][:28]} | ? |")
        L += ["", f"Planchers en vigueur : V >= {MIN_HUMAN_V_IDENTITY}, "
                  f"J >= {MIN_HUMAN_J_IDENTITY}, hallmarks VHH imposes.", ""]
        L += ["### Paratopes", ""]
        for b in sorted(usable, key=lambda x: -x["overlap"]):
            p = b["paratope"]
            if not p["ok"]:
                L.append(f"- {b['pdb']} {b['chain']} : {p['why']}")
                continue
            L.append(f"- **{b['pdb']} {b['chain']}** ({b['overlap']}/{b['site_size']}), "
                     f"{p['total']} contacts, CDR3 = `{p['cdr3']}` (L={len(p['cdr3'])}) : "
                     + ", ".join(f"{k} {100 * p['share'][k]:.0f}%" for k in
                                 ("CDR1", "CDR2", "CDR3", "framework")))
        L += ["", "Un liant dont le CDR3 porte peu du paratope est le seul pour lequel "
                  "« garder la charpente, redessiner le CDR3 » a un sens (FINDINGS 73).", ""]
        fixes = [b for b in usable if b["fix"] and b["fix"]["submittable"]]
        L += ["### Rendus soumettables", ""]
        if not fixes:
            L.append("Aucun liant n'est rendu conforme par des mutations hors CDR.")
        for b in fixes:
            f = b["fix"]
            L += [f"- **{b['pdb']} {b['chain']}** : {f['n_mutations']} mutations "
                  f"({', '.join(f['mutations'])}), tous les CDR identiques, couche 1 "
                  f"{f['layer1']}", f"  ```", f"  {f['sequence']}", f"  ```"]
    L += ["", "---", "Produit par `scripts.analyze_target`. Chaque nombre vient des "
          "fichiers de structure ou du filtre de nativite du validateur."]
    return "\n".join(L)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default=NANOBODY_TARGET[0])
    parser.add_argument("--structure", default=NANOBODY_STRUCTURE)
    parser.add_argument("--chain", default=NANOBODY_STRUCTURE_CHAIN)
    parser.add_argument("--res-index", default=NANOBODY_STRUCTURE_RES_INDEX)
    parser.add_argument("--binding-site", default=NANOBODY_BINDING_SITE)
    parser.add_argument("--no-survey", action="store_true",
                        help="skip the PDB search and the known-binder work")
    parser.add_argument("--cache", type=Path,
                        default=REPO_ROOT / "data" / "structures")
    parser.add_argument("--out", type=Path, default=REPO_ROOT / "results")
    args = parser.parse_args()

    a = analyse(args.target, args.structure, args.chain, args.res_index,
                args.binding_site, args.cache, not args.no_survey)
    args.out.mkdir(parents=True, exist_ok=True)
    md = args.out / f"target_{args.target}_{args.structure}.md"
    md.write_text(report(a))
    rows = [{k: v for k, v in b.items()
             if k not in ("paratope", "fix", "overlap_residues")}
            | {"cdr3_share": b["paratope"]["share"]["CDR3"] if b.get("paratope", {}).get("ok") else None,
               "fix_mutations": ", ".join(b["fix"]["mutations"]) if b.get("fix") else None,
               "fix_sequence": b["fix"]["sequence"] if b.get("fix") else None}
            for b in a["binders"] if "error" not in b]
    csv = args.out / f"target_{args.target}_{args.structure}_binders.csv"
    if rows:
        pd.DataFrame(rows).to_csv(csv, index=False)

    n, e = a["numbering"], a["epitope"]
    print(f"{a['target']} / {a['pdb']} chaine {a['chain']}")
    print(f"  numerotation : {'ok' if n['ok'] else 'ECART -- ' + n['why']}"
          f"  (label - auth = {n['label_minus_auth']})")
    print(f"  site de liaison : {len(e['inside'])}/{len(e['site'])} au contact")
    print(f"  liants lisibles : {sum(1 for b in a['binders'] if 'error' not in b)}")
    if a.get("survey_note"):
        print(f"  ATTENTION : {a['survey_note']}")
    print(f"  -> {md}")
    if rows:
        print(f"  -> {csv}")
    return 0 if n["ok"] and len(e["inside"]) == len(e["site"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
