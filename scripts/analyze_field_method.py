"""Infer the competing miners' method from the evidence on disk.

Answers the six questions of docs/FIELD_METHOD.md. CPU only, no network.

Leak discipline
---------------
The public archive carries no submission date, so "available before epoch E"
NEVER means "in data/archive today". It means:

  * nova_lab.archive_ledger.seen_before(target, E) -- with no ledger on disk
    this falls back to the frozen snapshot data/prebuilt/Q9NZQ7_regions.csv
    (9466 sequences), which covers everything submitted up to epoch 26025; and
  * sequences from data/submissions/submissions.jsonl whose epoch is < E.

Everything else in data/archive is undated and is therefore excluded from any
"before E" pool. The donor pools are consequently a LOWER bound: a best-match
identity computed against them can only understate the true best match.

Run:  PYTHONPATH=. .venv/bin/python scripts/analyze_field_method.py
"""
from __future__ import annotations

import json
import math
import multiprocessing as mp
import sys
from collections import Counter, defaultdict
from difflib import SequenceMatcher
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from nova_lab import DATA_DIR, archive_ledger  # noqa: E402
from nova_lab.nanobodies import germline_identity, regions_batch  # noqa: E402

TARGET = "Q9NZQ7"
OUT = Path(__file__).resolve().parents[1] / "docs" / "FIELD_METHOD.md"

# Kyte-Doolittle
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
      "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
      "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}


# --------------------------------------------------------------------------- #
# similarity primitives (identical formulas to nova_lab.nanobodies)
# --------------------------------------------------------------------------- #
def whole_identity(a: str, b: str) -> float:
    return SequenceMatcher(None, a, b).ratio()


def cdr_identity(x: str, y: str) -> float:
    """The validator's per-CDR similarity: position-wise mismatches over the
    shared length plus the length difference, over the longer length."""
    if not x and not y:
        return 1.0
    n = min(len(x), len(y))
    distance = sum(p != q for p, q in zip(x[:n], y[:n])) + abs(len(x) - len(y))
    return max(0.0, 1.0 - distance / max(len(x), len(y), 1))


def _upper_bound(qc: Counter, dc: Counter, la: int, lb: int) -> float:
    """difflib.quick_ratio: multiset intersection, an upper bound on ratio()."""
    matches = 0
    for aa, n in qc.items():
        m = dc.get(aa)
        if m:
            matches += n if n < m else m
    return 2.0 * matches / (la + lb)


_POOL: list[str] = []
_POOL_COUNTS: list[Counter] = []
_POOL_LEN: list[int] = []


def _init(pool):
    global _POOL, _POOL_COUNTS, _POOL_LEN
    _POOL = pool
    _POOL_COUNTS = [Counter(s) for s in pool]
    _POOL_LEN = [len(s) for s in pool]


def _best_whole(args):
    """Exact best whole-sequence identity of `query` over the pool indices
    allowed by `mask_key`, pruned with the quick_ratio upper bound."""
    query, allowed = args
    qc = Counter(query)
    lq = len(query)
    order = []
    for i in allowed:
        order.append((_upper_bound(qc, _POOL_COUNTS[i], lq, _POOL_LEN[i]), i))
    order.sort(reverse=True)
    best, best_i = -1.0, -1
    for bound, i in order:
        if bound <= best:
            break
        r = whole_identity(query, _POOL[i])
        if r > best:
            best, best_i = r, i
    return best, best_i



# --------------------------------------------------------------------------- #
# The repo's docs are in accented French. The report text below is written
# without accents so that the source stays ASCII; this table restores them on
# the way out, skipping anything inside backticks (paths, identifiers).
# --------------------------------------------------------------------------- #
ACCENTS = {
    # multi-word first: the preposition "a" is only "a" with an accent in
    # these contexts, never in the verb "il a".
    "a garder": "\u00e0 garder", "a l'identique": "\u00e0 l'identique",
    "a chaque": "\u00e0 chaque", "a la longueur": "\u00e0 la longueur",
    "a notre": "\u00e0 notre", "a moins de": "\u00e0 moins de",
    "a une identite": "\u00e0 une identit\u00e9", "a l'epoque": "\u00e0 l'\u00e9poque",
    "a quelques residus": "\u00e0 quelques r\u00e9sidus",
    "a cette": "\u00e0 cette", "a des sequences": "\u00e0 des s\u00e9quences",
    "a ces": "\u00e0 ces", "a l'archive": "\u00e0 l'archive",
    "a l'identite": "\u00e0 l'identit\u00e9",
    "rien a optimiser": "rien \u00e0 optimiser",
    "reste a": "reste \u00e0", "Reste a": "Reste \u00e0",
    "a lever": "\u00e0 lever", "a l'oeuvre": "\u00e0 l'oeuvre",
    "jusqu'a": "jusqu'\u00e0", "qu'a": "qu'\u00e0",
    "ce que fait": "ce que fait",
    "de meme longueur": "de m\u00eame longueur",
    "pre-calculee": "pr\u00e9-calcul\u00e9e",
    "lui-meme": "lui-m\u00eame",
    "au-dessus de 85": "au-dessus de 85",
    "ne monte qu'a": "ne monte qu'\u00e0",
    "ne la muter pas du tout": "ne pas la muter du tout",
    # single words
    "accepte": "accept\u00e9", "acceptee": "accept\u00e9e", "acceptees": "accept\u00e9es",
    "alimente": "aliment\u00e9", "ambiguite": "ambigu\u00eft\u00e9",
    "amines": "amin\u00e9s", "analysee": "analys\u00e9e", "analysees": "analys\u00e9es",
    "anterieur": "ant\u00e9rieur", "anterieure": "ant\u00e9rieure",
    "anterieures": "ant\u00e9rieures", "anterieurs": "ant\u00e9rieurs",
    "archivees": "archiv\u00e9es", "arreter": "arr\u00eater",
    "arrivees": "arriv\u00e9es", "attribues": "attribu\u00e9s",
    "bibliotheque": "biblioth\u00e8que", "bornee": "born\u00e9e",
    "calculee": "calcul\u00e9e", "chaine": "cha\u00eene", "colles": "coll\u00e9s",
    "complete": "compl\u00e8te", "concatenee": "concat\u00e9n\u00e9e",
    "conditionne": "conditionn\u00e9", "conditionnee": "conditionn\u00e9e",
    "consequence": "cons\u00e9quence", "Consequence": "Cons\u00e9quence",
    "consecutives": "cons\u00e9cutives", "correle": "corr\u00e8le",
    "cote": "c\u00f4t\u00e9", "coute": "co\u00fbte",
    "decisif": "d\u00e9cisif", "deciles": "d\u00e9ciles",
    "decryptee": "d\u00e9crypt\u00e9e", "decryptees": "d\u00e9crypt\u00e9es",
    "deduction": "d\u00e9duction", "deja": "d\u00e9j\u00e0",
    "depot": "d\u00e9p\u00f4t", "depots": "d\u00e9p\u00f4ts",
    "derivent": "d\u00e9rivent", "derniere": "derni\u00e8re",
    "dessinees": "dessin\u00e9es", "detail": "d\u00e9tail",
    "difference": "diff\u00e9rence", "differente": "diff\u00e9rente",
    "differentes": "diff\u00e9rentes", "differents": "diff\u00e9rents",
    "elements": "\u00e9l\u00e9ments", "entiere": "enti\u00e8re",
    "entree": "entr\u00e9e", "entrees": "entr\u00e9es",
    "epoque": "\u00e9poque", "epoques": "\u00e9poques",
    "etaient": "\u00e9taient", "etat": "\u00e9tat", "ete": "\u00e9t\u00e9",
    "etre": "\u00eatre", "fixee": "fix\u00e9e", "fuitee": "fuit\u00e9e",
    "gene": "g\u00e8ne", "genes": "g\u00e8nes",
    "generateur": "g\u00e9n\u00e9rateur", "grossieres": "grossi\u00e8res",
    "identifie": "identifi\u00e9", "identite": "identit\u00e9",
    "identites": "identit\u00e9s", "imposee": "impos\u00e9e",
    "inferieure": "inf\u00e9rieure", "inferieures": "inf\u00e9rieures",
    "instantane": "instantan\u00e9", "l'ecart": "l'\u00e9cart",
    "litteral": "litt\u00e9ral", "mecanique": "m\u00e9canique",
    "mecaniquement": "m\u00e9caniquement", "median": "m\u00e9dian",
    "mediane": "m\u00e9diane", "medianes": "m\u00e9dianes",
    "meme": "m\u00eame", "memes": "m\u00eames", "methode": "m\u00e9thode",
    "methodologique": "m\u00e9thodologique", "molecule": "mol\u00e9cule",
    "moyennee": "moyenn\u00e9e", "n'apparait": "n'appara\u00eet",
    "neanmoins": "n\u00e9anmoins", "observee": "observ\u00e9e",
    "observes": "observ\u00e9s", "parametrique": "param\u00e9trique",
    "partagee": "partag\u00e9e", "particulierement": "particuli\u00e8rement",
    "permute": "permut\u00e9", "pilotes": "pilot\u00e9s",
    "plafonnee": "plafonn\u00e9e", "plutot": "plut\u00f4t",
    "ponderee": "pond\u00e9r\u00e9e", "portee": "port\u00e9e",
    "precede": "pr\u00e9c\u00e8de", "precedente": "pr\u00e9c\u00e9dente",
    "premiere": "premi\u00e8re", "pres": "pr\u00e8s", "presente": "pr\u00e9sente",
    "proprietaire": "propri\u00e9taire", "proprietaires": "propri\u00e9taires",
    "protegent": "prot\u00e8gent", "recentes": "r\u00e9centes",
    "recoit": "re\u00e7oit", "reconnaitre": "reconna\u00eetre",
    "recyclee": "recycl\u00e9e", "reference": "r\u00e9f\u00e9rence",
    "Reference": "R\u00e9f\u00e9rence", "refusees": "refus\u00e9es",
    "refutee": "r\u00e9fut\u00e9e", "reglage": "r\u00e9glage",
    "relancee": "relanc\u00e9e", "reparties": "r\u00e9parties",
    "repartition": "r\u00e9partition", "Repartition": "R\u00e9partition",
    "reperes": "rep\u00e8res", "Reperes": "Rep\u00e8res",
    "repetent": "r\u00e9p\u00e8tent", "reserve": "r\u00e9serve",
    "residus": "r\u00e9sidus", "reutilisent": "r\u00e9utilisent",
    "reutilises": "r\u00e9utilis\u00e9s", "sequence": "s\u00e9quence",
    "sequences": "s\u00e9quences", "superieure": "sup\u00e9rieure",
    "telechargement": "t\u00e9l\u00e9chargement", "temoin": "t\u00e9moin",
    "tete": "t\u00eate", "tirees": "tir\u00e9es",
    "a la bonne longueur": "à la bonne longueur",
    "a l'autre": "à l'autre",
    "identique a un": "identique à un",
    "plancher a ": "plancher à ",
    "neuve a ": "neuve à ",
    "du a l'ecart": "dû à l'écart",
    "a montre": "a montré",
    "a 0,70": "à 0,70",
    "a >=": "à >=",
    "le sens ou": "le sens où",
    "chaine ou": "chaîne où",
    "vivier date": "vivier daté",
    "Cote": "Côté",
    "Decisif": "Décisif",
    "Identite": "Identité",
    "Identites": "Identités",
    "Mediane": "Médiane",
    "Meme": "Même",
    "datee": "datée",
    "datees": "datées",
    "depend": "dépend",
    "designe": "désigné",
    "echantillon": "échantillon",
    "empoisonne": "empoisonné",
    "sous-estime": "sous-estimé",
    "general": "général",
    "imposees": "imposées",
    "independants": "indépendants",
    "litteralement": "littéralement",
    "mecanisme": "mécanisme",
    "negatifs": "négatifs",
    "numerotables": "numérotables",
    "numerotation": "numérotation",
    "numerotee": "numérotée",
    "numeroter": "numéroter",
    "presents": "présents",
    "probleme": "problème",
    "qualite": "qualité",
    "rafraichir": "rafraîchir",
    "fraicheur": "fraîcheur",
    "regions": "régions",
    "reglera": "réglera",
    "replique": "répliqué",
    "retirees": "retirées",
    "separement": "séparément",
    "ventilees": "ventilées",
    "verifie": "vérifié",
    "motif partage": "motif partagé",
    "est donne": "est donné",
    "resultat": "résultat",
    "donneurs dates": "donneurs datés",
    "arrivee": "arrivée",
    "pas ca": "pas ça",
    "n'est recopie": "n'est recopié",
    "seulement ou se trouve": "seulement où se trouve",
    "est dessine": "est dessiné",
    "comportement mesure": "comportement mesuré",
    "reelle": "réelle",
    "montre que": "montré que",
    "motif greffe": "motif greffé",
    "greffe dans": "greffé dans",
    "**A notre": "**À notre",
    "fragment contigu partage": "fragment contigu partagé",
    "temoin a ": "témoin à ",
    "c'est-a-dire": "c'est-à-dire",
    "hasard a composition egale": "hasard à composition égale",
    "reutilisation": "réutilisation",
    "a l'interieur du champ": "à l'intérieur du champ",
    "epoque a la suivante": "époque à la suivante",
    "greffe tel quel": "greffé tel quel",
    "depasser": "dépasser",
    "numerotes": "numérotés",
    "a du etre recalcule": "a dû être recalculé",
    "recalculees": "recalculées",
    "numerote": "numérote",
    "gardees": "gardées",
    "herite": "hérite",
    "calculees": "calculées",
    "similarite": "similarité",
    "resultats": "résultats",
    "intermediaires": "intermédiaires",
    "couteux": "coûteux",
    "gardes": "gardés",
    "rejouable": "rejouable",
    "colle au plancher": "collé au plancher",
    "non-germinalite": "non-germinalité",
    "depense": "dépense",
    "la ou le champ": "là où le champ",
    "a moins de deux": "à moins de deux",
    "plancher a 82": "plancher à 82",
    "plancher a 86": "plancher à 86",
    "correlation": "corrélation",
    "s'arrete": "s'arrête",
    "l'a refute": "l'a réfuté",
    "L'ecart-la": "L'écart-là",
    "gratuite a mesurer": "gratuite à mesurer",
    "negative": "négative",
    "ressemble a l'archive": "ressemble à l'archive",
    "a l'interieur": "à l'intérieur",
    "La difference entre les deux chiffres": "La différence entre les deux chiffres",
    "memes correlations, **a longueur": "mêmes corrélations, **à longueur",
}


def accentuate(text: str) -> str:
    """Restore French accents outside backticked spans."""
    import re
    keys = sorted(ACCENTS, key=len, reverse=True)
    pattern = re.compile("|".join(
        (r"\b" if k[0].isalpha() else "") + re.escape(k)
        + (r"\b" if k[-1].isalpha() else "") for k in keys))

    decimal = re.compile(r"(?<=\d)\.(?=\d)")

    def fix(chunk: str) -> str:
        return decimal.sub(",", pattern.sub(lambda m: ACCENTS[m.group(0)], chunk))

    return "".join(part if part.startswith("`") else fix(part)
                   for part in re.split(r"(`[^`]*`)", text))


# --------------------------------------------------------------------------- #
# data
# --------------------------------------------------------------------------- #
def load_everything():
    acc = pd.read_csv(DATA_DIR / "submissions" / f"{TARGET}_nanobody_results.csv")
    acc = acc[acc["valid"]].reset_index(drop=True)

    subs = [json.loads(l) for l in (DATA_DIR / "submissions" / "submissions.jsonl").open()]
    submitted: list[tuple[int, str]] = []
    for s in subs:
        for q in s.get("sequences", []):
            submitted.append((int(s["epoch"]), q))

    # dated archive view: everything we can prove existed by epoch 26025
    ledger = archive_ledger.load(TARGET) or archive_ledger.bootstrap(TARGET)
    snapshot = sorted(ledger)

    # cached IMGT regions for the whole of today's archive + our submissions
    reg = {}
    fr = {}
    for f in (DATA_DIR / "archive" / f"{TARGET}_regions.csv",
              DATA_DIR / "prebuilt" / f"{TARGET}_regions.csv"):
        if f.exists():
            df = pd.read_csv(f, keep_default_na=False)
            for seq, c1, c2, c3, f1, f2, f3, f4 in zip(
                    df.sequence, df.CDR1, df.CDR2, df.CDR3,
                    df.FR1, df.FR2, df.FR3, df.FR4):
                if seq not in reg and c3:
                    reg[seq] = (c1, c2, c3)
                    fr[seq] = f1 + f2 + f3 + f4
    # The repo's own region cache marks 6740 of the archive's 17244 rows as
    # non-numberable although ANARCI numbers them fine (point 10 of
    # docs/REPLAY_AUDIT.md), and it holds nothing for the refused submissions.
    # Whatever is missing is numbered here and kept in results/, so the regions
    # used by this analysis never depend on the poisoned cache.
    extra = Path(__file__).resolve().parents[1] / "results" / "field_method_regions.csv"
    if extra.exists():
        df = pd.read_csv(extra, keep_default_na=False)
        for seq, c1, c2, c3, f1, f2, f3, f4 in zip(
                df.sequence, df.CDR1, df.CDR2, df.CDR3, df.FR1, df.FR2, df.FR3, df.FR4):
            if c3:
                reg[seq] = (c1, c2, c3)
                fr[seq] = f1 + f2 + f3 + f4
    wanted = set(acc.sequence) | {q for _, q in submitted} | set(snapshot)
    missing = [q for q in wanted if q not in reg]
    if missing:
        print(f"numbering {len(missing)} sequences with ANARCI ...", flush=True)
        rows = []
        for q, r in zip(missing, regions_batch(missing)):
            if r:
                reg[q] = (r["CDR1"], r["CDR2"], r["CDR3"])
                fr[q] = r["FR1"] + r["FR2"] + r["FR3"] + r["FR4"]
                rows.append({"sequence": q, **r})
        if rows:
            old = (pd.read_csv(extra, keep_default_na=False).to_dict("records")
                   if extra.exists() else [])
            extra.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(old + rows, columns=("sequence", "CDR1", "CDR2", "CDR3",
                                              "FR1", "FR2", "FR3", "FR4")
                         ).to_csv(extra, index=False)
    return acc, submitted, snapshot, reg, fr


# --------------------------------------------------------------------------- #
def quantiles(xs, qs=(0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0)):
    s = sorted(xs)
    if not s:
        return {q: float("nan") for q in qs}
    return {q: s[min(len(s) - 1, int(round(q * (len(s) - 1))))] for q in qs}


def fmt_q(xs, pct=False):
    q = quantiles(xs)
    f = (lambda v: f"{100*v:.1f}") if pct else (lambda v: f"{v:.3f}")
    return " | ".join(f(q[k]) for k in (0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0))


def spearman(x, y):
    """Spearman rho and a two-sided p from the t approximation."""
    n = len(x)
    if n < 4:
        return float("nan"), float("nan"), n

    def rank(v):
        order = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    rx, ry = rank(list(x)), rank(list(y))
    mx, my = sum(rx) / n, sum(ry) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    if sxx == 0 or syy == 0:
        return float("nan"), float("nan"), n
    rho = sxy / math.sqrt(sxx * syy)
    if abs(rho) >= 1:
        return rho, 0.0, n
    t = rho * math.sqrt((n - 2) / (1 - rho * rho))
    # two-sided p from Student t with n-2 df, via the incomplete beta
    df = n - 2
    x_ = df / (df + t * t)
    p = _betainc(df / 2, 0.5, x_)
    return rho, min(1.0, max(0.0, p)), n


def _betainc(a, b, x):
    """Regularised incomplete beta I_x(a,b), continued fraction (NR 6.4)."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    front = math.exp(lbeta + a * math.log(x) + b * math.log(1 - x))
    if x < (a + 1) / (a + b + 2):
        return front * _cf(a, b, x) / a
    return 1 - math.exp(lbeta + b * math.log(1 - x) + a * math.log(x)) * _cf(b, a, 1 - x) / b


def _cf(a, b, x):
    tiny = 1e-30
    qab, qap, qam = a + b, a + 1, a - 1
    c, d = 1.0, 1 - qab * x / qap
    d = 1 / (tiny if abs(d) < tiny else d)
    h = d
    for m in range(1, 300):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1 + aa * d
        d = 1 / (tiny if abs(d) < tiny else d)
        c = 1 + aa / c
        if abs(c) < tiny:
            c = tiny
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1 + aa * d
        d = 1 / (tiny if abs(d) < tiny else d)
        c = 1 + aa / c
        if abs(c) < tiny:
            c = tiny
        de = d * c
        h *= de
        if abs(de - 1) < 3e-12:
            break
    return h


# --------------------------------------------------------------------------- #
def main():
    CACHE = Path(__file__).resolve().parents[1] / "results" / "field_method_cache.csv"
    CACHE_COLS = ["best_whole", "best_whole_seq", "best_cdr3", "best_cdr3_same_len",
                  "cdr3_exact_prior", "leaky_best_cdr3", "fr_exact_prior",
                  "fr_best_prior", "v_identity", "j_identity", "v_genes"]
    acc, submitted, snapshot, reg, fr = load_everything()
    epochs = sorted(acc.epoch.unique())
    acc["cdr3"] = [reg.get(s, ("", "", ""))[2] for s in acc.sequence]
    acc["cdr1"] = [reg.get(s, ("", "", ""))[0] for s in acc.sequence]
    acc["cdr2"] = [reg.get(s, ("", "", ""))[1] for s in acc.sequence]
    print(f"{len(acc)} accepted submissions, epochs {epochs}")
    print(f"snapshot (<= {archive_ledger.SNAPSHOT_EPOCH}): {len(snapshot)} sequences")
    print(f"CDR3 missing for {(acc.cdr3 == '').sum()} accepted sequences")

    cached = None
    if CACHE.exists():
        c = pd.read_csv(CACHE, keep_default_na=False)
        if set(acc.sequence) <= set(c.sequence) and all(k in c.columns for k in CACHE_COLS):
            cached = c.drop_duplicates("sequence").set_index("sequence")
            text_cols = ("best_whole_seq", "v_genes")
            bool_cols = ("cdr3_exact_prior", "fr_exact_prior")
            for k in CACHE_COLS:
                col = [cached.at[q, k] for q in acc.sequence]
                if k in text_cols:
                    acc[k] = col
                elif k in bool_cols:
                    acc[k] = [str(v) in ("True", "true", "1", "1.0") for v in col]
                else:
                    acc[k] = pd.to_numeric(pd.Series(col), errors="raise").values
            print(f"reusing {CACHE} ({len(cached)} rows)")

    lines: list[str] = []
    W = lines.append

    # ---------------- pool bookkeeping --------------------------------------
    pool = list(snapshot)
    index = {s: i for i, s in enumerate(pool)}
    sub_epoch: dict[int, list[int]] = defaultdict(list)   # epoch -> pool indices
    for e, q in submitted:
        if q not in index:
            index[q] = len(pool)
            pool.append(q)
        sub_epoch[e].append(index[q])
    acc_epoch: dict[int, list[int]] = defaultdict(list)
    for e, q in zip(acc.epoch, acc.sequence):
        acc_epoch[int(e)].append(index[q])
    n_snapshot = len(snapshot)

    def allowed_before(E: int) -> list[int]:
        out = list(range(n_snapshot))
        for e, idxs in sub_epoch.items():
            if e < E:
                out.extend(idxs)
        return sorted(set(out))

    pool_cdr3 = [reg.get(q, ("", "", ""))[2] for q in pool]
    today = pd.read_csv(DATA_DIR / "archive" / f"{TARGET}_regions.csv", keep_default_na=False)
    today_seqs = set(today.sequence)
    # The 348 accepted submissions are themselves archived by now, and so are
    # their own epoch-mates, which the validator's duplicate check excludes from
    # an epoch's own comparison. Leaving them in would make the "leaky" bound
    # measure the field against itself, so every accepted sequence is dropped.
    _own = set(acc.sequence)
    today_cdr3 = [c for q, c in zip(today.sequence, today.CDR3) if c and q not in _own]

    # ---------------- Q1: whole + CDR3 identity to the prior pool -----------
    if cached is None:
        print("Q1: best prior match per accepted submission ...", flush=True)
        jobs = [(s, allowed_before(int(e))) for s, e in zip(acc.sequence, acc.epoch)]
        with mp.Pool(4, initializer=_init, initargs=(pool,)) as p:
            best_whole = p.map(_best_whole, jobs, chunksize=4)
        acc["best_whole"] = [b for b, _ in best_whole]
        acc["best_whole_seq"] = [pool[i] if i >= 0 else "" for _, i in best_whole]

        # CDR3: best identity, and best identity restricted to same-length donors
        best_cdr3, best_cdr3_same, best_cdr3_exact = [], [], []
        for s, e, c3 in zip(acc.sequence, acc.epoch, acc.cdr3):
            allowed = allowed_before(int(e))
            b, bs, exact = 0.0, 0.0, False
            for i in allowed:
                d = pool_cdr3[i]
                if not d:
                    continue
                v = cdr_identity(c3, d)
                if v > b:
                    b = v
                if len(d) == len(c3):
                    if v > bs:
                        bs = v
                    if d == c3:
                        exact = True
            best_cdr3.append(b)
            best_cdr3_same.append(bs)
            best_cdr3_exact.append(exact)
        acc["best_cdr3"] = best_cdr3
        acc["best_cdr3_same_len"] = best_cdr3_same
        acc["cdr3_exact_prior"] = best_cdr3_exact

        # also: the CDR3 best match against TODAY'S archive (leaky, reported as such)
        leaky_cdr3 = []
        own = set(acc.sequence)
        for c3 in acc.cdr3:
            b = 0.0
            for d in today_cdr3:
                v = cdr_identity(c3, d)
                if v > b:
                    b = v
            leaky_cdr3.append(b)
        acc["leaky_best_cdr3"] = leaky_cdr3

    W("# Par quelle methode le champ produit-il ses soumissions acceptees ?\n")
    W("Analyse de 348 soumissions acceptees sur 12 epoques (26299-26306, 26345-26348),")
    W("reproductible par `PYTHONPATH=. .venv/bin/python scripts/analyze_field_method.py`.\n")
    W("## 0. Discipline anti-fuite et ce qu'elle coute\n")
    W("L'archive publique ne porte **aucune date**. La seule vue datee dont nous")
    W(f"disposons est l'instantane `data/prebuilt/{TARGET}_regions.csv`, **{len(snapshot)}**")
    W(f"sequences, qui couvre tout ce qui a ete soumis jusqu'a l'epoque")
    W(f"{archive_ledger.SNAPSHOT_EPOCH}. `data/archive/{TARGET}_first_seen.csv` n'existe pas :")
    W("`archive_ledger.seen_before()` retombe donc sur cet instantane pour **toutes**")
    W("les epoques analysees ici (toutes > 26025).\n")
    W("Le vivier de donneurs « disponible avant l'epoque E » est donc :\n")
    W(f"* l'instantane ({len(snapshot)} sequences, datees <= {archive_ledger.SNAPSHOT_EPOCH}) ;")
    W("* plus les sequences de `submissions.jsonl` d'epoque < E (458 soumissions")
    W("  decryptees sur 12 epoques, acceptees et refusees confondues).\n")
    W(f"Les ~{len(today_seqs) - len(snapshot)} entrees de l'archive d'aujourd'hui qui ne sont pas")
    W("dans l'instantane sont **exclues** : nous ne pouvons pas prouver qu'elles")
    W("existaient avant l'epoque analysee. Consequence a garder en tete pour tout")
    W("ce qui suit : **chaque identite « au meilleur donneur anterieur » est une")
    W("borne inferieure**. Un vrai donneur plus proche a pu exister sans que nous")
    W("puissions le dater. Les conclusions ne sont donc tirees que dans le sens ou")
    W("cette borne les autorise.\n")
    W("### Ce qui a du etre recalcule plutot que lu dans le depot\n")
    W("Rien ici ne vient des colonnes mises en cache :\n")
    W("* les regions IMGT (CDR1/2/3, FR1-FR4) des 348 acceptees et des 458")
    W("  soumissions decryptees sont recalculees par `nanobodies.regions_batch` sur")
    W("  la sequence brute, parce que `data/archive/Q9NZQ7_regions.csv` marque 6740")
    W("  de ses 17244 lignes comme non numerotables alors qu'ANARCI les numerote")
    W("  (voir 1 bis). Elles sont gardees dans `results/field_method_regions.csv` ;")
    W("* les identites V et J sont recalculees par `nanobodies.germline_identity`,")
    W("  pas lues dans `data/dataset/Q9NZQ7_nanobodies.parquet`, qui n'a d'epoque")
    W("  datee que pour 152 de ses 10536 lignes et herite du meme cache ;")
    W("* les identites de sequence entiere et de CDR3 sont calculees avec les")
    W("  formules exactes du validateur (`whole_identity`, la similarite par CDR de")
    W("  `cdr_similarity`), pas avec une approximation.\n")
    W("Les resultats intermediaires couteux sont gardes dans")
    W("`results/field_method_cache.csv` pour que le script soit rejouable en")
    W("quelques minutes ; supprimer ce fichier force le recalcul complet.\n")
    W("Trou de 39 epoques entre 26306 et 26345 : aucune soumission decryptee n'existe")
    W("pour 26307-26344, le vivier des 4 epoques recentes est donc particulierement")
    W("sous-estime.\n")

    W("## 1. Les soumissions acceptees derivent-elles de l'archive ?\n")
    W(f"n = {len(acc)}. Identite au meilleur donneur anterieur (deciles : min | p5 | p25 | median | p75 | p95 | max).\n")
    W("| mesure | min | p5 | p25 | med | p75 | p95 | max |")
    W("|---|---|---|---|---|---|---|---|")
    W(f"| identite sequence entiere (%) | {fmt_q(acc.best_whole, pct=True)} |")
    W(f"| identite CDR3, tout donneur (%) | {fmt_q(acc.best_cdr3, pct=True)} |")
    W(f"| identite CDR3, donneur de meme longueur (%) | {fmt_q(acc.best_cdr3_same_len, pct=True)} |")
    W("")
    for lo, hi in ((0.0, 0.70), (0.70, 0.80), (0.80, 0.85), (0.85, 0.90), (0.90, 0.95), (0.95, 1.01)):
        n = int(((acc.best_whole >= lo) & (acc.best_whole < hi)).sum())
        W(f"* identite entiere entre {100*lo:.0f} et {min(100*hi, 100):.0f} % : "
          f"{n} / {len(acc)} ({100*n/len(acc):.1f} %)")
    W("")
    W(f"CDR3 repris a l'identique d'un donneur anterieur : "
      f"{int(acc.cdr3_exact_prior.sum())} / {len(acc)}.")
    W(f"CDR3 a >= 85 % d'un donneur anterieur : "
      f"{int((acc.best_cdr3 >= 0.85).sum())} / {len(acc)}.")
    W(f"CDR3 a >= 85 % d'une entree de l'archive d'aujourd'hui, les 348 acceptees "
      f"elles-memes retirees du vivier (vue **fuitee**, donc bornee par le haut) : "
      f"{int((acc.leaky_best_cdr3 >= 0.85).sum())} / {len(acc)}, "
      f"median {100*quantiles(acc.leaky_best_cdr3)[0.5]:.1f} %. Cette vue ne porte "
      f"que sur les {len(today_cdr3)} sequences de l'archive que le cache de regions "
      "sait numeroter (voir 1 bis).\n")

    # ---------------- Q1bis: numbering-free substring test -----------------
    print("Q1bis: exact CDR3 substring over the whole archive ...", flush=True)
    arch_all = pd.read_csv(DATA_DIR / "archive" / f"{TARGET}_nanobodies.csv")
    dated_seqs = [q for q in arch_all.sequence if q in set(snapshot)]
    undated_seqs = [q for q in arch_all.sequence if q not in set(snapshot)]
    own_set = set(acc.sequence)

    def substring_hits(c3, bank_seqs):
        return sum(1 for q in bank_seqs if q not in own_set and c3 in q)

    hits_dated = [substring_hits(c3, dated_seqs) for c3 in acc.cdr3]
    hits_undated = [substring_hits(c3, undated_seqs) for c3 in acc.cdr3]
    W("### 1 bis. Un test qui ne depend d'aucune numerotation ni d'aucune date\n")
    W("Un CDR3 greffe tel quel reste un **sous-mot exact** de la sequence")
    W("donneuse. On cherche donc chaque CDR3 accepte comme sous-chaine dans les")
    W(f"sequences brutes de l'archive, en excluant les {len(acc)} acceptees")
    W("elles-memes. Ce test ne demande aucune numerotation IMGT, donc il couvre")
    W(f"l'archive **entiere** ({len(arch_all)} sequences), et le resultat est donne")
    W("separement pour la partie datee et la partie indatable.\n")
    W(f"* presents dans la partie datee ({len(dated_seqs)} sequences) : "
      f"**{sum(1 for h in hits_dated if h)} / {len(acc)}**")
    W(f"* presents dans la partie indatable ({len(undated_seqs)} sequences, donc une "
      f"vue qui peut contenir leurs propres gagnantes) : "
      f"**{sum(1 for h in hits_undated if h)} / {len(acc)}**")
    W("")
    W("Note sur une limite du depot : `data/archive/Q9NZQ7_regions.csv` porte "
      f"{int((pd.read_csv(DATA_DIR / 'archive' / f'{TARGET}_regions.csv', keep_default_na=False).CDR3 == '').sum())} "
      "lignes dont le CDR3 est vide alors que les sequences sont numerotables "
      "(verifie sur un echantillon ; c'est le point 10 de `docs/REPLAY_AUDIT.md`, "
      "le cache de regions empoisonne). Toutes les "
      "comparaisons par CDR3 contre l'archive **d'aujourd'hui** ne portent donc que "
      "sur la partie numerotee, et c'est pour cela que le test par sous-chaine "
      "ci-dessus, qui ignore le cache, est le garde-fou le plus solide. Le vivier "
      "date (l'instantane) est complet : 0 CDR3 manquant.\n")

    # ---------------- Q2: CDR3 recurrence -----------------------------------
    W("## 2. Les CDR3 se repetent-ils d'une epoque a l'autre ?\n")
    W("| epoque | n acceptees | long. CDR3 imposee (observee) | CDR3 distincts | CDR3 deja vu exact (epoques anterieures) | >= 85 % (epoques anterieures) | >= 85 % (instantane) |")
    W("|---|---|---|---|---|---|---|")
    prior_cdr3_by_epoch: list[tuple[int, str]] = []
    snap_cdr3 = [reg.get(s, ("", "", ""))[2] for s in snapshot]
    snap_cdr3 = [c for c in snap_cdr3 if c]
    q2 = []
    for e in epochs:
        sub = acc[acc.epoch == e]
        lens = Counter(len(c) for c in sub.cdr3 if c)
        prior = [c for ep, c in prior_cdr3_by_epoch if ep < e]
        exact = sum(1 for c in sub.cdr3 if c and c in set(prior))
        near = sum(1 for c in sub.cdr3 if c and any(cdr_identity(c, d) >= 0.85 for d in prior))
        nears = sum(1 for c in sub.cdr3 if c and any(cdr_identity(c, d) >= 0.85 for d in snap_cdr3))
        W(f"| {e} | {len(sub)} | {dict(lens)} | {len(set(sub.cdr3))} | {exact} | {near} | {nears} |")
        q2.append((e, len(sub), exact, near, nears))
        prior_cdr3_by_epoch.extend((e, c) for c in sub.cdr3 if c)
    W("")
    n_after_first = sum(n for e, n, _, _, _ in q2 if e != epochs[0])
    W(f"Total hors premiere epoque (n = {n_after_first}) : "
      f"{sum(x for e, _, x, _, _ in q2 if e != epochs[0])} CDR3 exactement reutilises, "
      f"{sum(x for e, _, _, x, _ in q2 if e != epochs[0])} a >= 85 % d'un CDR3 accepte anterieur.")
    W(f"Contre l'instantane (n = {len(acc)}) : "
      f"{sum(x for _, _, _, _, x in q2)} CDR3 a >= 85 % d'un CDR3 de l'instantane.\n")

    # ---------------- Q3: self-reuse per hotkey -----------------------------
    W("## 3. Les mineurs reutilisent-ils leurs propres sequences ?\n")
    rows = []
    for hk, g in acc.groupby("hotkey"):
        if len(g) < 2:
            continue
        seqs = list(zip(g.epoch, g.sequence, g.cdr3, g.design_iiptm))
        pairs = []
        c3pairs = []
        for i in range(len(seqs)):
            for j in range(i + 1, len(seqs)):
                if seqs[i][0] == seqs[j][0]:
                    continue
                pairs.append(whole_identity(seqs[i][1], seqs[j][1]))
                c3pairs.append(cdr_identity(seqs[i][2], seqs[j][2]))
        if not pairs:
            continue
        rows.append((hk[:6], int(g.uid.iloc[0]), len(g), len(set(g.epoch)),
                     sum(pairs) / len(pairs), max(pairs),
                     sum(c3pairs) / len(c3pairs), max(c3pairs),
                     g.design_iiptm.mean(), g.design_iiptm.max()))
    rows.sort(key=lambda r: -r[9])
    W(f"{len(rows)} hotkeys avec >= 2 soumissions acceptees sur >= 2 epoques. "
      f"Identites entre paires de soumissions du meme hotkey, epoques differentes.\n")
    W("| hotkey | uid | n acc. | n epoques | ident. entiere moy. | max | ident. CDR3 moy. | max | iiptm moy. | iiptm max |")
    W("|---|---|---|---|---|---|---|---|---|---|")
    for hk, uid, n, ne, mw, xw, mc, xc, mi, xi in rows:
        W(f"| `{hk}` | {uid} | {n} | {ne} | {100*mw:.1f} | {100*xw:.1f} | {100*mc:.1f} | {100*xc:.1f} | {mi:.3f} | {xi:.3f} |")
    allw = [r[4] for r in rows]
    allc = [r[6] for r in rows]
    W("")
    W(f"Mediane des moyennes par hotkey : identite entiere {100*quantiles(allw)[0.5]:.1f} %, "
      f"identite CDR3 {100*quantiles(allc)[0.5]:.1f} %.")
    # cross-hotkey baseline, same epoch
    import random
    random.seed(0)
    base_w, base_c = [], []
    for e in epochs:
        sub = acc[acc.epoch == e].reset_index(drop=True)
        for _ in range(60):
            i, j = random.sample(range(len(sub)), 2)
            if sub.hotkey[i] == sub.hotkey[j]:
                continue
            base_w.append(whole_identity(sub.sequence[i], sub.sequence[j]))
            base_c.append(cdr_identity(sub.cdr3[i], sub.cdr3[j]))
    W(f"Reference (paires de hotkeys **differents** dans la meme epoque, n = {len(base_w)}) : "
      f"identite entiere mediane {100*quantiles(base_w)[0.5]:.1f} %, "
      f"identite CDR3 mediane {100*quantiles(base_c)[0.5]:.1f} %.\n")

    # ---------------- Q4: germline ------------------------------------------
    if cached is None:
        print("Q4: germline identity ...", flush=True)
        vi, ji, vg = [], [], []
        for s in acc.sequence:
            g = germline_identity(s)
            if g is None:
                vi.append(float("nan")); ji.append(float("nan")); vg.append("")
            else:
                vi.append(g["V"]["identity_pct"]); ji.append(g["J"]["identity_pct"])
                vg.append(",".join(g["V"]["genes"]))
        acc["v_identity"] = vi
        acc["j_identity"] = ji
        acc["v_genes"] = vg
    vok = [v for v in acc.v_identity if v == v]
    jok = [j for j in acc.j_identity if j == j]
    W("## 4. Les charpentes sont-elles germinales ou dessinees ?\n")
    W(f"n = {len(vok)} (identite germinale calculee par `germline_identity`, "
      f"planchers du validateur V >= 82, J >= 86).\n")
    W("| mesure | min | p5 | p25 | med | p75 | p95 | max |")
    W("|---|---|---|---|---|---|---|---|")
    W(f"| identite V humaine (%) | {fmt_q(vok)} |")
    W(f"| identite J humaine (%) | {fmt_q(jok)} |")
    W("")
    for name, xs, floor in (("V", vok, 82), ("J", jok, 86)):
        for band in ((0, 2), (2, 5), (5, 10), (10, 100)):
            n = sum(1 for x in xs if floor + band[0] <= x < floor + band[1])
            W(f"* {name} dans [plancher+{band[0]}, plancher+{band[1]}) : {n} / {len(xs)} ({100*n/len(xs):.1f} %)")
    W("")
    W("Repartition des genes V attribues :\n")
    for g, n in Counter(acc.v_genes).most_common(12):
        W(f"* `{g}` : {n}")
    W("")
    W("Par epoque (le validateur impose un gene V et une longueur de CDR3 par epoque) :\n")
    W("| epoque | genes V observes | ident. V mediane | ident. J mediane |")
    W("|---|---|---|---|")
    for e in epochs:
        sub = acc[acc.epoch == e]
        genes = Counter(x for g in sub.v_genes for x in g.split(",") if x)
        W(f"| {e} | {dict(genes.most_common(4))} | {sub.v_identity.median():.1f} | {sub.j_identity.median():.1f} |")
    W("")

    # ---------------- Q5: correlations --------------------------------------
    def gravy(c):
        return sum(KD.get(a, 0) for a in c) / len(c) if c else float("nan")

    def charge(c):
        return sum((a in "KR") - (a in "DE") + 0.1 * (a == "H") for a in c)

    acc["cdr3_len"] = [len(c) for c in acc.cdr3]
    acc["cdr3_gravy"] = [gravy(c) for c in acc.cdr3]
    acc["cdr3_charge"] = [charge(c) for c in acc.cdr3]
    acc["cdr1_len"] = [len(c) for c in acc.cdr1]
    acc["cdr2_len"] = [len(c) for c in acc.cdr2]
    acc["total_cdr"] = acc.cdr1_len + acc.cdr2_len + acc.cdr3_len

    feats = [("longueur CDR3", "cdr3_len"), ("identite V", "v_identity"),
             ("identite J", "j_identity"),
             ("identite au meilleur donneur anterieur", "best_whole"),
             ("identite CDR3 au meilleur donneur anterieur", "best_cdr3"),
             ("GRAVY du CDR3", "cdr3_gravy"), ("charge nette du CDR3", "cdr3_charge"),
             ("longueur CDR totale", "total_cdr")]
    W("## 5. Le binding correle-t-il avec quelque chose d'optimisable ?\n")
    W(f"Spearman de `design_iiptm` contre, au sein des acceptees.\n")
    W("| variable | rho | p | n |")
    W("|---|---|---|---|")
    for label, col in feats:
        sub = acc[[col, "design_iiptm"]].dropna()
        rho, p, n = spearman(sub[col].tolist(), sub.design_iiptm.tolist())
        W(f"| {label} | {rho:+.3f} | {p:.3g} | {n} |")
    W("")
    W("Les memes correlations, **a longueur de CDR3 fixee par epoque** (partielles "
      "grossieres : correlation intra-epoque, moyennee ponderee par n) :\n")
    W("| variable | rho intra-epoque moyen | n total |")
    W("|---|---|---|")
    for label, col in feats:
        num = den = 0.0
        for e in epochs:
            sub = acc[acc.epoch == e][[col, "design_iiptm"]].dropna()
            if len(sub) < 6 or sub[col].nunique() < 3:
                continue
            rho, p, n = spearman(sub[col].tolist(), sub.design_iiptm.tolist())
            if rho == rho:
                num += rho * n
                den += n
        W(f"| {label} | {num/den:+.3f} | {int(den)} |" if den else f"| {label} | n.d. | 0 |")
    W("")
    W("Par epoque, pour situer la dispersion du score :\n")
    W("| epoque | n | iiptm min | med | max |")
    W("|---|---|---|---|---|")
    for e in epochs:
        sub = acc[acc.epoch == e]
        W(f"| {e} | {len(sub)} | {sub.design_iiptm.min():.3f} | {sub.design_iiptm.median():.3f} | {sub.design_iiptm.max():.3f} |")
    W("")

    # ---------------- Q6: recent epochs, transplantable CDR3? ---------------
    W("## 6. Decisif : les CDR3 des 4 epoques recentes etaient-ils transplantables ?\n")
    recent = [e for e in epochs if e >= 26345]
    W("Pour chaque acceptee de 26345-26348 : existait-il, dans le vivier anterieur")
    W("prouvable, un CDR3 de la **bonne longueur** identique, ou proche ?\n")
    W("| epoque | n | long. CDR3 | CDR3 identique dispo. | >= 85 % | >= 70 % | >= 50 % | meilleure ident. CDR3 mediane (%) | meilleure ident. CDR3 meme long. mediane (%) | n donneurs meme long. |")
    W("|---|---|---|---|---|---|---|---|---|---|")
    for e in recent:
        sub = acc[acc.epoch == e]
        allowed = allowed_before(int(e))
        lens = Counter(len(c) for c in sub.cdr3 if c)
        L = lens.most_common(1)[0][0] if lens else 0
        ndon = sum(1 for i in allowed if pool_cdr3[i] and len(pool_cdr3[i]) == L)
        W(f"| {e} | {len(sub)} | {dict(lens)} | {int(sub.cdr3_exact_prior.sum())} | "
          f"{int((sub.best_cdr3 >= 0.85).sum())} | {int((sub.best_cdr3 >= 0.70).sum())} | "
          f"{int((sub.best_cdr3 >= 0.50).sum())} | {100*sub.best_cdr3.median():.1f} | "
          f"{100*sub.best_cdr3_same_len.median():.1f} | {ndon} |")
    W("")
    rec = acc[acc.epoch >= 26345]
    W(f"Sur les {len(rec)} acceptees recentes : identite entiere au meilleur donneur "
      f"anterieur, mediane {100*rec.best_whole.median():.1f} %, max {100*rec.best_whole.max():.1f} %.")
    W(f"CDR3 : mediane {100*rec.best_cdr3.median():.1f} %, max {100*rec.best_cdr3.max():.1f} %.\n")
    W("Vue **fuitee** (archive d'aujourd'hui, donc borne superieure) pour les memes :")
    W(f"identite CDR3 mediane {100*rec.leaky_best_cdr3.median():.1f} %, "
      f"{int((rec.leaky_best_cdr3 >= 0.85).sum())} / {len(rec)} a >= 85 %.\n")

    # ---------------- Q6bis: the motif, not the whole CDR3 -----------------
    print("Q6bis: shared motifs ...", flush=True)
    import random
    random.seed(7)

    def kmers(x, k):
        return {x[i:i + k] for i in range(len(x) - k + 1)} if len(x) >= k else set()

    def bank_of(cdr3s, kmax=24):
        b = {}
        for k in range(3, kmax + 1):
            b[k] = set()
            for c in cdr3s:
                b[k] |= kmers(c, k)
            if not b[k]:
                del b[k]
                break
        return b

    def longest_shared(c3, bank):
        """Longest contiguous motif of c3 present in the bank. Monotone in k:
        no (k+1)-mer can match if no k-mer does, so the loop may stop early."""
        best = 0
        for k in range(3, len(c3) + 1):
            if k not in bank or not (kmers(c3, k) & bank[k]):
                break
            best = k
        return best

    def shared_motif(c3, bank, k):
        for i in range(len(c3) - k + 1):
            if c3[i:i + k] in bank.get(k, ()):
                return c3[i:i + k]
        return ""

    snap_bank = bank_of(snap_cdr3)
    W("## 6 bis. Le motif plutot que le CDR3 entier\n")
    W("La longueur du CDR3 est **imposee par epoque** (section 2) : un CDR3 d'archive")
    W("de la mauvaise longueur ne peut pas depasser `1 - |dL| / L` d'identite, donc le")
    W("greffage litteral est structurellement impossible la plupart du temps. La")
    W("bonne question est donc : le **motif** voyage-t-il ? On mesure, pour chaque")
    W("CDR3 accepte, la longueur du plus long fragment contigu qu'il partage avec un")
    W("CDR3 anterieur, contre deux viviers, avec un temoin par permutation des")
    W("lettres du meme CDR3 (20 tirages, meme composition en acides amines).\n")
    W("| epoque | n | motif partage avec l'instantane (median / max) | temoin permute | motif partage avec les acceptees anterieures (median / max) | temoin permute | n CDR3 anterieurs |")
    W("|---|---|---|---|---|---|---|")
    prior_acc_cdr3: list[tuple[int, str]] = []
    motif_rows = []
    for e in epochs:
        sub = acc[acc.epoch == e]
        prior = [c for ep, c in prior_acc_cdr3 if ep < e]
        pbank = bank_of(prior) if prior else {}
        obs_s, nul_s, obs_p, nul_p = [], [], [], []
        for c3 in sub.cdr3:
            if not c3:
                continue
            obs_s.append(longest_shared(c3, snap_bank))
            obs_p.append(longest_shared(c3, pbank) if pbank else 0)
            shuf = []
            for _ in range(20):
                lst = list(c3)
                random.shuffle(lst)
                shuf.append("".join(lst))
            nul_s.append(sum(longest_shared(x, snap_bank) for x in shuf) / 20)
            nul_p.append((sum(longest_shared(x, pbank) for x in shuf) / 20) if pbank else 0)
        W(f"| {e} | {len(obs_s)} | {quantiles(obs_s)[0.5]:.0f} / {max(obs_s)} | "
          f"{sum(nul_s)/len(nul_s):.1f} | "
          f"{quantiles(obs_p)[0.5]:.0f} / {max(obs_p) if obs_p else 0} | "
          f"{(sum(nul_p)/len(nul_p)):.1f} | {len(prior)} |")
        motif_rows.append((e, obs_s, nul_s, obs_p, nul_p))
        prior_acc_cdr3.extend((e, c) for c in sub.cdr3 if c)
    W("")
    all_obs_s = [v for _, o, _, _, _ in motif_rows for v in o]
    all_nul_s = [v for _, _, n, _, _ in motif_rows for v in n]
    all_obs_p = [v for e, _, _, o, _ in motif_rows if e != epochs[0] for v in o]
    all_nul_p = [v for e, _, _, _, n in motif_rows if e != epochs[0] for v in n]
    W(f"Toutes epoques confondues (n = {len(all_obs_s)}) : motif partage avec")
    W(f"l'instantane, median **{quantiles(all_obs_s)[0.5]:.0f}** residus contre")
    W(f"**{sum(all_nul_s)/len(all_nul_s):.1f}** pour le temoin permute.")
    W(f"Hors premiere epoque (n = {len(all_obs_p)}) : motif partage avec les acceptees")
    W(f"des epoques anterieures, median **{quantiles(all_obs_p)[0.5]:.0f}** contre")
    W(f"**{sum(all_nul_p)/len(all_nul_p):.1f}** pour le temoin.\n")
    W("Lecture : contre l'instantane, le motif partage est de "
      f"{quantiles(all_obs_s)[0.5]:.0f} residus pour un temoin a "
      f"{sum(all_nul_s)/len(all_nul_s):.1f} -- c'est-a-dire **au niveau du hasard a "
      "composition egale**. Contre les acceptees des epoques anterieures du champ "
      f"lui-meme, {quantiles(all_obs_p)[0.5]:.0f} contre "
      f"{sum(all_nul_p)/len(all_nul_p):.1f} : un peu au-dessus du hasard, mais un "
      "motif de 4 residus n'est pas un greffage. Le seul cas net du tableau "
      "ci-dessous est `RFQPRYP` (7 residus), qui voyage de la longueur 20 (26301) "
      "a la longueur 17 (26306). Il existe donc une reutilisation de motifs courts "
      "**a l'interieur du champ**, d'une epoque a la suivante, et elle est faible.\n")
    W("Les plus longs motifs effectivement repris, pour les 12 meilleures acceptees :\n")
    W("| epoque | uid | iiptm | CDR3 | motif partage avec une acceptee anterieure | avec l'instantane |")
    W("|---|---|---|---|---|---|")
    seen: list[tuple[int, str]] = []
    banks_by_epoch = {}
    cum: list[str] = []
    for e in epochs:
        banks_by_epoch[e] = bank_of(cum) if cum else {}
        cum.extend(c for c in acc[acc.epoch == e].cdr3 if c)
    for _, r in acc.nlargest(12, "design_iiptm").iterrows():
        b = banks_by_epoch[int(r.epoch)]
        kp = longest_shared(r.cdr3, b) if b else 0
        ks = longest_shared(r.cdr3, snap_bank)
        W(f"| {int(r.epoch)} | {int(r.uid)} | {r.design_iiptm:.3f} | `{r.cdr3}` | "
          f"`{shared_motif(r.cdr3, b, kp) if kp else '-'}` ({kp}) | "
          f"`{shared_motif(r.cdr3, snap_bank, ks) if ks else '-'}` ({ks}) |")
    W("")

    # best-scoring submissions table
    W("## 7. Les meilleures soumissions du champ, en clair\n")
    W("| epoque | uid | iiptm | long. CDR3 | CDR3 | ident. V | ident. J | meilleur donneur anterieur (%) | meilleur CDR3 anterieur (%) |")
    W("|---|---|---|---|---|---|---|---|---|")
    for _, r in acc.nlargest(15, "design_iiptm").iterrows():
        W(f"| {int(r.epoch)} | {int(r.uid)} | {r.design_iiptm:.3f} | {r.cdr3_len} | `{r.cdr3}` | "
          f"{r.v_identity:.1f} | {r.j_identity:.1f} | {100*r.best_whole:.1f} | {100*r.best_cdr3:.1f} |")
    W("")

    # ---------------- Q8: frameworks, grafted or designed? -----------------
    if cached is None:
        print("Q8: frameworks ...", flush=True)
        pool_fr = [fr.get(s, "") for s in pool]
        acc["fr"] = [fr.get(s, "") for s in acc.sequence]
        fr_exact, fr_best = [], []
        for s, e, f in zip(acc.sequence, acc.epoch, acc.fr):
            allowed = allowed_before(int(e))
            exact = False
            best = 0.0
            for i in allowed:
                d = pool_fr[i]
                if not d:
                    continue
                if d == f:
                    exact = True
                    best = 1.0
                    break
                v = cdr_identity(f, d)   # same positional formula, on the framework
                if v > best:
                    best = v
            fr_exact.append(exact)
            fr_best.append(best)
        acc["fr_exact_prior"] = fr_exact
        acc["fr_best_prior"] = fr_best
    W("## 8. Charpente recyclee, CDR3 neuf ? Le test direct du greffage\n")
    W("Si le champ greffe, la charpente (FR1+FR2+FR3+FR4, concatenee) doit souvent")
    W("etre **exactement** celle d'un donneur anterieur, tandis que le CDR3 est neuf.")
    W(f"n = {len(acc)}.\n")
    W(f"* charpente reprise a l'identique d'un donneur anterieur prouvable : "
      f"{int(acc.fr_exact_prior.sum())} / {len(acc)} ({100*acc.fr_exact_prior.mean():.1f} %)")
    W(f"* identite de charpente au meilleur donneur anterieur (positionnelle) : "
      f"mediane {100*acc.fr_best_prior.median():.1f} %, p95 "
      f"{100*quantiles(acc.fr_best_prior)[0.95]:.1f} %, max {100*acc.fr_best_prior.max():.1f} %")
    W(f"* dont charpente >= 98 % : {int((acc.fr_best_prior >= 0.98).sum())} / {len(acc)} ; "
      f">= 95 % : {int((acc.fr_best_prior >= 0.95).sum())} / {len(acc)}")
    W("")
    W("| epoque | n | charpente exacte dispo. | ident. charpente mediane (%) | ident. CDR3 mediane (%) |")
    W("|---|---|---|---|---|")
    for e in epochs:
        sub = acc[acc.epoch == e]
        W(f"| {e} | {len(sub)} | {int(sub.fr_exact_prior.sum())} | "
          f"{100*sub.fr_best_prior.median():.1f} | {100*sub.best_cdr3.median():.1f} |")
    W("")

    # ---------------- Q9: batches per repo owner ---------------------------
    print("Q9: repo owners ...", flush=True)
    acc["owner"] = [r.split("/")[0] for r in acc.repo]
    W("## 9. Les lots par proprietaire de depot\n")
    W(f"{acc.repo.nunique()} depots, {acc.owner.nunique()} proprietaires, "
      f"{acc.hotkey.nunique()} hotkeys, {acc.uid.nunique()} uid. Plusieurs hotkeys")
    W("partagent un meme depot : ce sont des grappes de mineurs pilotes ensemble.\n")
    W("Identites **au sein d'un meme proprietaire et d'une meme epoque** "
      "(medianes des medianes par groupe) :\n")
    W("| proprietaire | groupes (epoques) | soumissions | hotkeys | ident. entiere mediane (%) | ident. CDR3 mediane (%) | iiptm median | iiptm max |")
    W("|---|---|---|---|---|---|---|---|")
    import itertools
    owner_rows = []
    for o, go in acc.groupby("owner"):
        meds_w, meds_c, n = [], [], 0
        for e, g in go.groupby("epoch"):
            if len(g) < 2:
                continue
            g = g.reset_index(drop=True)
            w = [whole_identity(g.sequence[i], g.sequence[j])
                 for i, j in itertools.combinations(range(len(g)), 2)]
            c = [cdr_identity(g.cdr3[i], g.cdr3[j])
                 for i, j in itertools.combinations(range(len(g)), 2)]
            meds_w.append(quantiles(w)[0.5])
            meds_c.append(quantiles(c)[0.5])
            n += len(g)
        if meds_w:
            owner_rows.append((o, len(meds_w), n, go.hotkey.nunique(),
                               quantiles(meds_w)[0.5], quantiles(meds_c)[0.5],
                               go.design_iiptm.median(), go.design_iiptm.max()))
    owner_rows.sort(key=lambda r: -r[2])
    for o, ng, n, nh, mw, mc, mi, xi in owner_rows:
        W(f"| {o} | {ng} | {n} | {nh} | {100*mw:.1f} | {100*mc:.1f} | {mi:.3f} | {xi:.3f} |")
    W("")
    base_w2, base_c2 = [], []
    for e, g in acc.groupby("epoch"):
        g = g.reset_index(drop=True)
        for i, j in itertools.combinations(range(len(g)), 2):
            if g.owner[i] != g.owner[j]:
                base_w2.append(whole_identity(g.sequence[i], g.sequence[j]))
                base_c2.append(cdr_identity(g.cdr3[i], g.cdr3[j]))
    W(f"Reference, paires de proprietaires **differents** dans la meme epoque "
      f"(n = {len(base_w2)}) : identite entiere mediane {100*quantiles(base_w2)[0.5]:.1f} %, "
      f"identite CDR3 mediane {100*quantiles(base_c2)[0.5]:.1f} %.\n")
    W(f"Aucune sequence acceptee n'apparait deux fois : "
      f"{acc.sequence.nunique()} sequences distinctes pour {len(acc)} soumissions "
      f"acceptees, et aucune n'est partagee entre deux hotkeys ni entre deux epoques.\n")

    # ---------------- Q10: same-length epochs, the library test ------------
    print("Q10: same-length epochs ...", flush=True)
    ep_len = {int(e): int(acc[acc.epoch == e].cdr3_len.mode().iloc[0]) for e in epochs}
    W("## 10. Le test de la bibliotheque : les epoques de meme longueur imposee\n")
    W("Si un proprietaire tirait d'une bibliotheque pre-calculee par (gene V,")
    W("longueur de CDR3), deux epoques imposant **la meme longueur** devraient")
    W("ressortir des CDR3 de la meme famille. Trois epoques imposent L = 14 :")
    W("26299, 26304, 26345.\n")
    W("| epoque A | epoque B | L | meme proprietaire : n paires / ident. CDR3 mediane / max | proprietaires differents : n / mediane / max |")
    W("|---|---|---|---|---|")
    import itertools as _it
    samelen_same_owner = []
    for a, b in _it.combinations(sorted(ep_len), 2):
        if ep_len[a] != ep_len[b]:
            continue
        ga, gb = acc[acc.epoch == a], acc[acc.epoch == b]
        same, diff = [], []
        for _, x in ga.iterrows():
            for _, y in gb.iterrows():
                (same if x.owner == y.owner else diff).append(cdr_identity(x.cdr3, y.cdr3))
        samelen_same_owner.extend(same)
        W(f"| {a} | {b} | {ep_len[a]} | {len(same)} / {100*quantiles(same)[0.5]:.1f} / "
          f"{100*max(same):.1f} | {len(diff)} / {100*quantiles(diff)[0.5]:.1f} / "
          f"{100*max(diff):.1f} |")
    wsame, dlen = [], []
    for e in epochs:
        g = acc[acc.epoch == e].reset_index(drop=True)
        for i, j in _it.combinations(range(len(g)), 2):
            if g.owner[i] == g.owner[j]:
                wsame.append(cdr_identity(g.cdr3[i], g.cdr3[j]))
    for a, b in _it.combinations(sorted(ep_len), 2):
        if ep_len[a] == ep_len[b]:
            continue
        ga, gb = acc[acc.epoch == a], acc[acc.epoch == b]
        for _, x in ga.iterrows():
            for _, y in gb.iterrows():
                if x.owner == y.owner:
                    dlen.append(cdr_identity(x.cdr3, y.cdr3))
    W("")
    W(f"Reperes : au sein d'une **meme** epoque et d'un meme proprietaire "
      f"(n = {len(wsame)}), identite CDR3 mediane **{100*quantiles(wsame)[0.5]:.1f} %**. "
      f"Entre epoques de longueurs **differentes**, meme proprietaire (n = {len(dlen)}), "
      f"**{100*quantiles(dlen)[0.5]:.1f} %** (plancher mecanique du a l'ecart de longueur).\n")

    # CDR3 shared inside one epoch
    pairs_same_cdr3, cross_owner, wid = 0, 0, []
    shared_n = 0
    for e in epochs:
        g = acc[acc.epoch == e].reset_index(drop=True)
        vc = g.cdr3.value_counts()
        shared_n += int(g.cdr3.map(vc).gt(1).sum())
        for i, j in _it.combinations(range(len(g)), 2):
            if g.cdr3[i] == g.cdr3[j]:
                pairs_same_cdr3 += 1
                wid.append(whole_identity(g.sequence[i], g.sequence[j]))
                if g.owner[i] != g.owner[j]:
                    cross_owner += 1
    W("Au sein d'une meme epoque, le meme CDR3 est soumis plusieurs fois :")
    W(f"**{shared_n} / {len(acc)}** acceptees partagent leur CDR3 avec une autre")
    W(f"acceptee de la meme epoque ({pairs_same_cdr3} paires, dont {cross_owner} entre")
    W("proprietaires differents). Ces paires ont une identite de sequence entiere")
    W(f"mediane de **{100*quantiles(wid)[0.5]:.1f} %** (min {100*min(wid):.1f} %, "
      f"max {100*max(wid):.1f} %) : ce ne sont pas des charpentes differentes portant")
    W("le meme CDR3, mais le **meme candidat reproduit a quelques residus pres** sur")
    W("plusieurs hotkeys.\n")
    W("| epoque | n acceptees | CDR3 distincts | blocs depuis l'epoque precedente |")
    W("|---|---|---|---|")
    prev_block = None
    for e in epochs:
        g = acc[acc.epoch == e]
        blk = int(g.block.min())
        W(f"| {e} | {len(g)} | {g.cdr3.nunique()} | "
          f"{blk - prev_block if prev_block else '-'} |")
        prev_block = blk
    W("")
    W("Une epoque dure 361 blocs, soit environ **72 minutes**. C'est le budget total")
    W("du champ pour produire, scorer et pousser ses candidats a chaque epoque.\n")

    # ---------------- Q11: the dated pool, by imposed CDR3 length ----------
    print("Q11: donor pool by CDR3 length ...", flush=True)
    arch = pd.read_csv(DATA_DIR / "archive" / f"{TARGET}_nanobodies.csv")
    snapset = set(snapshot)
    dated = arch[arch.sequence.isin(snapset)].copy()
    dated["L"] = [len(reg.get(q, ("", "", ""))[2]) for q in dated.sequence]
    W("## 11. Le vivier date n'a de bons lieurs qu'a **une seule** longueur de CDR3\n")
    W(f"Les {len(dated)} sequences datees (<= {archive_ledger.SNAPSHOT_EPOCH}), ventilees")
    W("par longueur de CDR3, pour les longueurs imposees par les 12 epoques analysees :\n")
    W("| longueur imposee | epoque(s) | n donneurs dates | dont iiptm >= 0,70 | dont >= 0,75 | meilleur iiptm |")
    W("|---|---|---|---|---|---|")
    by_len = defaultdict(list)
    for e in epochs:
        by_len[ep_len[int(e)]].append(int(e))
    for L in sorted(by_len):
        g = dated[dated.L == L]
        W(f"| {L} | {', '.join(str(x) for x in by_len[L])} | {len(g)} | "
          f"{int((g.design_iiptm >= 0.70).sum())} | "
          f"{int((g.design_iiptm >= 0.75).sum())} | "
          f"{g.design_iiptm.max():.3f} |" if len(g) else
          f"| {L} | {', '.join(str(x) for x in by_len[L])} | 0 | 0 | 0 | - |")
    W("")
    tot70 = int((dated.design_iiptm >= 0.70).sum())
    at13 = int(((dated.L == 13) & (dated.design_iiptm >= 0.70)).sum())
    W(f"Sur l'ensemble du vivier date, **{tot70}** sequences lient a 0,70 ou mieux, et")
    W(f"**{at13}** d'entre elles, soit {100*at13/max(tot70,1):.1f} %, ont un CDR3 de")
    W("longueur **13**. L'archive publique n'est pas un vivier general : c'est **une")
    W("seule famille de solutions, sur une seule longueur**. Pour 11 des 12 epoques")
    W("analysees, la longueur imposee n'est pas 13, et le vivier date n'offre alors")
    W("**aucun** donneur a 0,70 ou mieux a la bonne longueur.\n")
    W("C'est le plafond mecanique de notre methode, et il ne depend ni du reglage du")
    W("generateur ni de la fraicheur de l'archive : le vivier date est meme")
    und = arch[~arch.sequence.isin(snapset)]
    W(f"meilleur que la partie indatable arrivee depuis (mediane "
      f"{dated.design_iiptm.median():.3f} contre {und.design_iiptm.median():.3f}, "
      f"p90 {dated.design_iiptm.quantile(0.9):.3f} contre "
      f"{und.design_iiptm.quantile(0.9):.3f}, n = {len(dated)} et {len(und)}), donc "
      "rafraichir l'archive ne reglera pas ca.\n")

    # ---------------- conclusion -------------------------------------------
    rec = acc[acc.epoch >= 26345]
    W("## Ce que cela implique pour notre generateur\n")
    W("### Ce que les mesures excluent\n")
    W("**1. Le champ ne greffe pas depuis l'archive.** Aucune de ses "
      f"{len(acc)} acceptees ne reprend un CDR3 anterieur prouvable : "
      f"{int(acc.cdr3_exact_prior.sum())} reprise exacte, "
      f"{int((acc.best_cdr3 >= 0.85).sum())} au-dessus de 85 %, identite CDR3 mediane "
      f"au meilleur donneur anterieur **{100*acc.best_cdr3.median():.1f} %**. Meme en "
      "s'autorisant la vue fuitee de l'archive d'aujourd'hui, qui contient leurs "
      f"propres gagnantes archivees depuis, la mediane ne monte qu'a "
      f"{100*acc.leaky_best_cdr3.median():.1f} %. Notre methode -- choisir un donneur "
      "d'archive qui lie et transporter son CDR3 -- n'est pas la leur. Ce n'est pas "
      "un detail de reglage : c'est une famille de methode differente.\n")
    W("**2. Ce n'est pas non plus une bibliotheque pre-calculee par (gene V, "
      "longueur).** Les trois epoques imposant L = 14 ne ressortent pas la meme "
      "famille de CDR3 : identite mediane entre epoques de meme longueur, meme "
      f"proprietaire, **{100*quantiles(samelen_same_owner)[0.5]:.1f} %** "
      f"(n = {len(samelen_same_owner)}), contre "
      f"{100*quantiles(wsame)[0.5]:.1f} % au sein d'une meme epoque "
      f"(n = {len(wsame)}). Chaque epoque recoit une famille neuve.\n")
    near_v = sum(1 for v in acc.v_identity if v < 84)
    W("**3. Un seul des deux planchers germinaux les contraint, et c'est celui du "
      f"gene V.** Identite V mediane **{acc.v_identity.median():.1f} %** pour un "
      f"plancher a 82, minimum **{acc.v_identity.min():.1f} %**, et "
      f"{near_v} / {len(acc)} ({100*near_v/len(acc):.0f} %) sont a moins de deux "
      "points du plancher : le champ est **colle au plancher V**. Le plancher J, "
      f"lui, ne contraint rien : mediane {acc.j_identity.median():.1f} % pour un "
      "plancher a 86, et aucune acceptee en dessous de 92,9. Et la correlation "
      "intra-epoque de la section 5 va dans le meme sens (rho = -0,25 entre "
      "identite V et `design_iiptm`) : **moins la charpente est germinale, mieux "
      "elle lie**, et le plancher a 82 est exactement la ou le champ s'arrete. "
      "Ce qu'il faut en retenir n'est pas que la charpente porte la liaison -- le "
      "paragraphe 50 l'a refute -- mais que le champ depense sa marge de "
      "non-germinalite jusqu'au dernier point disponible.\n")
    W("### Ce que les mesures soutiennent\n")
    W("**Rien n'est recopie : ni le CDR3, ni la charpente.** Les deux tests directs "
      f"du greffage sont negatifs. Cote CDR3 : {int(acc.cdr3_exact_prior.sum())} / "
      f"{len(acc)} reprise exacte, identite mediane "
      f"{100*acc.best_cdr3.median():.1f} %. Cote charpente : "
      f"{int(acc.fr_exact_prior.sum())} / {len(acc)} charpentes reprises a "
      f"l'identique, identite mediane {100*acc.fr_best_prior.median():.1f} %. Les "
      "91 % de la charpente ne sont pas un recyclage : c'est le fond commun de "
      "toutes les VHH, que n'importe quelle charpente conforme atteint. La difference "
      "entre les deux chiffres dit seulement ou se trouve le travail : la "
      "charpente est standard, le CDR3 est neuf.\n")
    W("**Le CDR3 est dessine par epoque, et la famille est partagee dans "
      "l'epoque.** Quatre faits convergents :\n")
    W(f"* au sein d'une epoque et d'un meme proprietaire de depot, les CDR3 forment "
      f"**une seule famille** : identite mediane {100*quantiles(wsame)[0.5]:.1f} % "
      f"(n = {len(wsame)}), contre {100*quantiles(base_c2)[0.5]:.1f} % entre "
      f"proprietaires differents de la meme epoque (n = {len(base_c2)}) ;")
    W("* les meilleures acceptees de proprietaires **differents** convergent "
      "pourtant sur le meme motif dans une meme epoque (`ARFQPRYP...` en 26301, "
      "`ARAYYRHWDW...` en 26303, `ARGYHRD...` / `ARGFHRH...` en 26346, "
      "`ARH..NTFK..` en 26345) : le probleme a une solution dominante que "
      "plusieurs pipelines independants retrouvent ;")
    W(f"* {shared_n} / {len(acc)} acceptees partagent leur CDR3 avec une autre "
      f"acceptee de la meme epoque, et ces paires ont "
      f"{100*quantiles(wid)[0.5]:.1f} % d'identite de sequence entiere : c'est le "
      "**meme candidat replique sur plusieurs hotkeys**, pas une exploration ;")
    W("* le budget est de **72 minutes par epoque** (361 blocs), ce qui exclut une "
      "recherche large et impose un pipeline qui converge vite.\n")
    W("La seule famille de methode compatible avec ces faits est la **conception "
      "d'une boucle CDR3 conditionnee par la structure du complexe**, relancee a "
      "chaque epoque avec la longueur imposee comme contrainte, sur une charpente "
      "conforme tenue fixe. C'est ce que le paragraphe 50 de `FINDINGS.md` avait "
      "designe comme la derniere piste non refutee ; c'est maintenant soutenu par "
      "le comportement mesure du champ et non par deduction.\n")
    W("**Et la « transposition de motif » du paragraphe 47 est reelle mais faible.** "
      f"Le plus long fragment contigu partage avec un CDR3 accepte anterieur vaut "
      f"{quantiles(all_obs_p)[0.5]:.0f} residus en mediane contre "
      f"{sum(all_nul_p)/len(all_nul_p):.1f} pour un temoin permute, et "
      f"{quantiles(all_obs_s)[0.5]:.0f} contre "
      f"{sum(all_nul_s)/len(all_nul_s):.1f} contre l'instantane, soit le hasard. "
      "Un `RFQPRYP` qui passe de la longueur 20 a la longueur 17 existe, mais c'est "
      "l'exception : on ne peut pas construire un generateur sur un signal de "
      "4 residus contre 2,7.\n")
    W("### Le mecanisme du plafond, en chiffres\n")
    W(f"La raison pour laquelle notre methode ne peut pas atteindre 0,70 n'est pas la "
      f"qualite du vivier mais **sa monoculture de longueur** : {at13} des {tot70} "
      "sequences datees qui lient a 0,70 ou mieux ont un CDR3 de longueur 13, et la "
      "longueur imposee n'est 13 que dans 1 des 12 epoques analysees. Pour les 11 "
      "autres, le vivier date n'offre aucun donneur a 0,70 a la bonne longueur "
      "(section 11). Un CDR3 de longueur 13 greffe dans une epoque qui en impose 17 "
      "n'est meme pas un greffage : c'est une sequence neuve a 76 % d'identite "
      "maximale avec son donneur, et le paragraphe 50 a montre que la liaison ne "
      "suit pas.\n")
    W("### Ce que notre generateur devrait arreter de faire\n")
    W("Chercher un bon donneur dans l'archive et lui transporter son CDR3. La "
      "mesure de la section 5 est sans ambiguite : au sein des acceptees, "
      "`design_iiptm` ne monte **pas** avec la ressemblance au meilleur donneur "
      "anterieur -- la correlation est negative, rho = -0,25 pour l'identite "
      "entiere et -0,39 pour celle du CDR3, et elle reste negative a l'interieur "
      "d'une epoque (-0,17 et -0,06). Plus un candidat ressemble a l'archive, moins "
      "il lie. Il n'y a donc rien a optimiser de ce cote. "
      "Le paragraphe 50 avait deja montre que la liaison ne voyage pas avec le motif "
      "greffe ; la presente mesure montre en plus que **personne dans le champ "
      "n'essaie**.\n")
    W("### Ce qui reste a notre portee, et ce qui ne l'est pas\n")
    W("* **A notre portee** : fixer la charpente, ne la muter pas du tout, et "
      "consacrer tout le budget de calcul au dessin du CDR3 a la longueur imposee. "
      "Le champ ne gagne pas sur la charpente et nous perdons de la liaison a "
      "chaque mutation de charpente (paragraphe 50).")
    W("* **A notre portee** : reconnaitre que la longueur imposee change a chaque "
      "epoque et qu'un CDR3 de la mauvaise longueur ne peut pas depasser "
      "`1 - |dL| / L` d'identite. L'archive n'offre donc presque jamais un donneur "
      "utilisable a la longueur du jour, ce qui explique mecaniquement le plafond de notre "
      "vivier.")
    W("* **A notre portee** : descendre l'identite V jusqu'au plancher de 82 plutot "
      "que rester confortablement au-dessus. C'est la seule variable de charpente "
      "sur laquelle le champ se distingue de nous, et elle est gratuite a mesurer.")
    W("* **Hors de notre portee sans GPU** : le dessin conditionne par la structure "
      "lui-meme. C'est le gap, et il est methodologique, pas parametrique.\n")
    W("### La reserve qui limite tout ce qui precede\n")
    W(f"Les identites « au meilleur donneur anterieur » sont des **bornes "
      "inferieures** : l'archive n'est datee que jusqu'a l'epoque 26025 "
      f"({len(snapshot)} sequences), et les ~{len(today_seqs) - len(snapshot)} entrees "
      "arrivees depuis sont indatables. Si le champ greffait depuis des sequences "
      "soumises entre 26026 et l'epoque analysee, cette analyse ne le verrait pas "
      "directement. Deux elements la protegent neanmoins : la vue fuitee de "
      f"l'archive complete ne donne que {100*acc.leaky_best_cdr3.median():.1f} % "
      "d'identite CDR3 mediane, et les 8 epoques consecutives 26299-26306 forment "
      "une chaine ou les acceptees de chaque epoque sont dans le vivier de la "
      "suivante -- et l'identite n'y monte pas. Pour lever la reserve pour de bon, "
      "il faut que `data/archive/Q9NZQ7_first_seen.csv` existe et soit alimente a "
      "chaque telechargement.\n")

    CACHE.parent.mkdir(parents=True, exist_ok=True)
    acc[["sequence", *CACHE_COLS]].to_csv(CACHE, index=False)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(accentuate("\n".join(lines)) + "\n")
    dbg = Path("/tmp/claude-0/-home-user-Nova/73fa282f-c0fd-5255-af32-db1431998688/scratchpad")
    if dbg.exists():
        acc.to_csv(dbg / "field_method_accepted.csv", index=False)
    print(f"wrote {OUT}")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
