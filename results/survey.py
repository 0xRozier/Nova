import json, urllib.request
Q = {"query": {"type": "terminal", "service": "text", "parameters": {
        "attribute": "rcsb_polymer_entity_container_identifiers."
                     "reference_sequence_identifiers.database_accession",
        "operator": "exact_match", "value": "P05231"}},
     "return_type": "entry", "request_options": {"paginate": {"start": 0, "rows": 100}}}
url = "https://search.rcsb.org/rcsbsearch/v2/query?json=" + urllib.parse.quote(json.dumps(Q))
res = json.load(urllib.request.urlopen(url, timeout=60))
ids = [r["identifier"] for r in res["result_set"]]
print(f"{res['total_count']} entrees PDB contiennent P05231 (IL-6) : {' '.join(ids)}\n")
G = """{entries(entry_ids:[%s]){rcsb_id struct{title}
 polymer_entities{rcsb_polymer_entity{pdbx_description}
 entity_poly{rcsb_sample_sequence_length}
 rcsb_polymer_entity_container_identifiers{auth_asym_ids}}}}""" % ",".join(f'"{i}"' for i in ids)
r = urllib.request.Request("https://data.rcsb.org/graphql",
        data=json.dumps({"query": G}).encode(), headers={"Content-Type": "application/json"})
d = json.load(urllib.request.urlopen(r, timeout=60))["data"]["entries"]
KEY = ("nanobody", "vhh", "single-domain", "single domain", "fab", "antibody",
       "immunoglobulin", "heavy chain", "light chain", "variable domain")
for e in sorted(d, key=lambda x: x["rcsb_id"]):
    parts = []
    for p in e["polymer_entities"]:
        desc = (p["rcsb_polymer_entity"]["pdbx_description"] or "?")
        ln = p["entity_poly"]["rcsb_sample_sequence_length"]
        ch = ",".join(p["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"])
        parts.append(f"{ch}:{desc[:46]}({ln})")
    hit = any(k in (e["struct"]["title"] + " ".join(parts)).lower() for k in KEY)
    print(f"{'>>' if hit else '  '} {e['rcsb_id']}  {e['struct']['title'][:70]}")
    for p in parts: print(f"      {p}")
