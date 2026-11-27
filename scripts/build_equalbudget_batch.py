"""Le lot a budget egal : permuter un donneur connu, ou tirer des candidats neufs ?

HYPOTHESE. Le §90 a trouve que permuter l'interieur du CDR3 deplace l'iiptm de
0,060 d'ecart-type pour un bruit de mesure de 0,006, sans direction (moyenne
signee -0,004). C'est donc une source de variance gratuite : aucune requete GPU
n'est depensee pour *produire* le candidat. Sur un objectif ou nos 311 sequences
interrogees tiennent entre 0,118 et 0,446, une recherche par meilleur-de-k
consomme de la variance. D'ou l'hypothese : **a budget GPU egal, k permutations
d'un donneur deja note battent k candidats neufs du generateur habituel.**

Le §90 ne pouvait pas en decider : il n'avait pas de bras temoin « candidats
neufs ». C'est la seule question restee ouverte sur l'axe generateur (ETAT §9).

DESIGN. Apparie, et c'est le point : l'epoque impose un gene V et une longueur
de CDR3, donc la comparaison operationnelle est « a V et L donnes, je permute
mon meilleur donneur ou je genere ? ». Huit creneaux (V, L), chacun tire d'un
donneur deja note sur la cible vivante, et dans chaque creneau les deux bras
recoivent **le meme nombre de requetes** :

  calib   (14) — des sequences de l'archive P05231 dont le validateur a publie
                 le score : l'etalonnage de l'instrument (lecture 0)
  donneur (8)  — le donneur lui-meme, re-note : c'est le temoin inter-runs,
                 et il rend la lecture 1 immune au biais de selection
  permut  (24) — 3 permutations de l'interieur de son CDR3 (meme V, meme L)
  neuf    (24) — 3 candidats neufs de `scripts.propose_nanobody` pour ce
                 meme (V, L), generateur `pssm`, `skip_surrogate=True`

70 sequences, soit ~140 min sur 2x T4 au tarif mesure de 2,0 min/sequence.

L'ordre d'execution n'est pas celui du CSV : BoltzGen nomme chaque entree par
l'empreinte SHA-256 de la sequence et parcourt le repertoire, donc un run
interrompu rend un sous-echantillon aleatoire — equilibre entre bras en
esperance, ce qui vaut mieux qu'un entrelacement.

Choix des donneurs, et pourquoi la conformite d'abord. Le record de la cible
(0,4457, §90) porte une identite V de 79,6 pour un plancher a 82 : il n'est pas
soumettable, et permuter son CDR3 ne le rend pas soumettable. Un bras dont
aucune sequence ne peut etre soumise ne mesure rien d'operationnel. Les huit
donneurs sont donc les mieux notes qui, a la fois : passent la couche 1 pour
leur propre (V, L), portent un CDR3 de >= 12 residus (sous 12, le §90 a montre
que la permutation frole le plafond de Jaccard 0,90 du validateur), et sont
distincts entre eux (Jaccard 5-mer < 0,80). Si `propose` ne peut pas servir un
creneau (archive trop pauvre en ce gene ou en cette longueur), le creneau passe
au donneur suivant de la liste.

Pas de `archive_cutoff` : on ne rejoue pas une epoque, on demande ou depenser la
prochaine requete aujourd'hui. Generateur `pssm` et `skip_surrogate=True` parce
que c'est le chemin du mineur (`scripts/mine.py`), et parce que le substitut
vaut 0,145 en prospectif (§83) : il ne doit pas choisir a la place du GPU.

LE BRAS D'ETALONNAGE, ET POURQUOI IL EST DEVENU LA LECTURE PRINCIPALE.

En construisant ce lot, le vivier du generateur a rendu un fait que le depot
n'avait jamais regarde : l'archive publique de P05231 porte 46 118 sequences
notees par les validateurs, de mediane 0,654, dont 92,9 % au-dessus de 0,45 et
43,9 % au-dessus de 0,70. Nos 327 requetes sur la meme cible ont une mediane de
0,181 et un maximum de 0,446 — zero au-dessus de 0,45.

Sur l'ANCIENNE cible, notre reproduction du score etait exacte : 84 sequences
d'archive re-notees chez nous donnent rho = +0,984 et un decalage median de
+0,0012. Mais sur la cible vivante, **aucun run n'a jamais note une sequence
dont on connait le score du validateur** : les six `ref_intact` du §77 n'y sont
pas, ni aucune des 327. L'instrument n'a donc jamais ete etalonne depuis la
bascule du 20 novembre.

Les indices libres ne tranchent pas, et il faut le dire :

  - nos mesures tombent SUR la courbe (iiptm, pae) de l'archive — a iiptm egal,
    notre pae colle a 0,3-0,6 pres dans chaque bande. Mais une cible mal posee
    produit de vrais mauvais complexes, qui se rangent eux aussi au mauvais bout
    de cette courbe : l'accord est donc compatible avec les deux hypotheses ;
  - l'amont est identique au notre, caractere par caractere : meme moteur
    (external_tools/boltzgen inchange entre db9bf19 et HEAD), meme MSA, meme
    yaml, memes constantes de config (verifie contre le HEAD amont) ;
  - mais le §77 a mesure 61H7 — le liant cristallographique de 4O9H, la
    structure que la config elle-meme nomme, 8/8 de recouvrement d'epitope — a
    **0,178**, soit la mediane de tout ce que nous mesurons, et sous 92,9 % des
    soumissions ordinaires de l'archive. Un liant natif qui s'amarre a son
    propre epitope ne note pas sous le tout-venant si le montage est correct.

Deux hypotheses restent, et elles ont la meme consequence pratique :
**(a)** notre montage P05231 deprime tout d'environ 0,45, ou **(b)** les scores
de l'archive P05231 ont ete produits sous une autre configuration que celle
d'aujourd'hui — l'archive ne porte aucune date, le depot le sait depuis le §48.
Dans les deux cas, toute comparaison de nos valeurs absolues a celles du champ
est nulle tant que ce n'est pas tranche, et 14 requetes le tranchent.

(Note pour qui cherchera le bug si (a) gagne : `parse_range` de BoltzGen est
positionnel et 1-indexe dans la chaine parsee (`c_start + start`), pas en
`label_seq_id`. La chaine A de 4O9H compte 279 residus et son premier porte
`label_seq_id` 21, donc `res_index: 21..186` ne selectionne pas les residus 21
a 186 de la sequence canonique. L'ancienne cible disait `1..115`, ou les deux
conventions coincident — c'est pourquoi le probleme ne pouvait pas se voir
avant la bascule. Ce n'est PAS une explication de notre ecart avec les
validateurs, qui font tourner le meme code sur les memes fichiers ; c'est une
piste pour (b), et un fait a connaitre dans tous les cas.)

LECTURES PRE-ENREGISTREES (ecrites avant le run, seuils compris).

  0. PRINCIPALE — l'instrument. 14 sequences de l'archive P05231, deux par
     bande de 0,15-0,30 / 0,30-0,45 / 0,45-0,55 / 0,55-0,65 / 0,65-0,72 /
     0,72-0,78 / 0,78-0,90, notees par nous dans ce run.
     rho(notre iiptm, iiptm de l'archive) et decalage median.
     - rho >= 0,80 et |decalage median| <= 0,05 : l'instrument reproduit le
       validateur sur cette cible. Alors nos candidats sont vraiment dans le
       dernier decile du champ, et le §85 est bien plus severe qu'ecrit.
     - rho >= 0,80 et decalage median <= -0,10 : monotone mais deprime. Les
       rangs restent utilisables, mais tout seuil absolu et toute comparaison
       au champ sont a refaire.
     - rho < 0,80 : l'instrument ne reproduit pas le validateur sur cette
       cible. Les lectures absolues des §77, §87, §89 et §90 doivent etre
       relues, et les lectures 1 a 3 ci-dessous ne valent plus que comme
       comparaisons internes au run.
     Repere : 84 sequences sur l'ancienne cible donnaient rho +0,984 et un
     decalage de +0,0012. Avec 14 sequences etalees sur 0,15-0,85, rho >= 0,80
     est une barre basse qu'un instrument fidele franchit sans effort.

  1. meilleur-de-k, appariee. Pour chacun des 8 creneaux,
     delta = max(3 permutations) - max(3 neufs). On rapporte la moyenne des
     delta, leur mediane, et le nombre de creneaux gagnes par la permutation.
     Verdict : la permutation gagne si moyenne(delta) > +0,014 ET au moins 5
     creneaux sur 8 ; les candidats neufs gagnent si moyenne(delta) < -0,014 ET
     au plus 3 creneaux sur 8 ; sinon aucune difference a ce budget — et
     l'avis operationnel est alors de prefer les neufs, qui n'exigent pas
     de donneur deja paye. Le seuil 0,014 est l'ecart typique entre rivaux
     dans une epoque (§49) : sous lui, rien ne change de rang.
  2. Distribution complete. Mann-Whitney bilateral sur les 24 permutations
     contre les 24 neufs, plus l'ecart des medianes. Difference appelee a
     p < 0,05.
  3. La queue haute, qui est ce qui decide vraiment. Nombre de requetes de
     chaque bras au-dessus de 0,3869 — le meilleur iiptm **soumettable**
     jamais mesure sur cette cible — et le maximum de chaque bras.
  4. TEMOIN — bruit entre runs. Les 8 donneurs re-notes contre leur valeur
     precedente. Passe si moyenne|delta| < 0,027 (2 sigma, §78) et
     |moyenne signee| < 0,014. En cas d'echec, le run n'est pas comparable aux
     precedents et la lecture 4 seule tombe : les lectures 1 a 3 sont
     internes au run et survivent (§3).
  5. Soumettabilite et diversite, gratuites. Part de chaque bras qui passe la
     couche 1, et Jaccard 5-mer maximal de chaque sequence contre toutes les
     autres du lot. Attendu : les neufs a 100 % (propose filtre), les
     permutations a 100 % puisque leurs donneurs sont choisis conformes.
     **Si une permutation casse la conformite, c'est un cout de la strategie
     et il doit etre rapporte** : l'interieur du CDR3 entre dans l'alignement
     qui calcule l'identite V.

Usage :
    python -m scripts.build_equalbudget_batch --out inputs/kaggle_equalbudget.csv
    python -m scripts.build_equalbudget_batch --analyse results/scores_nanobody_equalbudget.csv \
        --batch inputs/kaggle_equalbudget.csv
"""
import argparse

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

from nova_lab import DATA_DIR, REPO_ROOT
from nova_lab.nanobodies import check_sequence, germline_identity, regions
from nova_lab.provenance import CURRENT, assert_current_target, runs_on_current_target
from scripts.build_scramble_batch import jaccard_5mer, permute_cdr3, score_path

N_SLOTS = 8
K = 3                      # requetes par bras et par creneau
N_CALIB_PER_BAND = 2
CALIB_BANDS = [(0.15, 0.30), (0.30, 0.45), (0.45, 0.55), (0.55, 0.65),
               (0.65, 0.72), (0.72, 0.78), (0.78, 0.90)]
MIN_CDR3 = 12              # §90 : sous 12, la permutation frole le plafond 0,90
DONOR_DISTINCT = 0.80      # Jaccard 5-mer maximal entre deux donneurs
GENERATOR = "pssm"         # le defaut du mineur
SEED = 20261127

# Seuils pre-enregistres
RIVAL_GAP = 0.014          # §49, ecart typique entre rivaux
SIGMA_2 = 0.027            # §78, 2x le bruit entre runs
WINS_FOR_PERMUT = 5        # creneaux sur 8
WINS_FOR_FRESH = 3
INCUMBENT = 0.3869         # meilleur iiptm soumettable mesure sur P05231
CALIB_RHO = 0.80           # lecture 0 : l'instrument reproduit-il le validateur
CALIB_OFFSET_OK = 0.05
CALIB_OFFSET_DEPRESSED = -0.10
JACCARD_CEILING = 0.90     # la regle du validateur
MWU_ALPHA = 0.05


# ---------------------------------------------------------------------------
# le vivier : tout ce qui est note sur la cible vivante
# ---------------------------------------------------------------------------
def live_pool(log=print) -> pd.DataFrame:
    """Les sequences notees sur la cible vivante, meilleur iiptm par sequence."""
    names = runs_on_current_target()
    paths = [score_path(n) for n in names]
    assert_current_target(paths, "le vivier de donneurs")
    frames = []
    for name, path in zip(names, paths):
        df = pd.read_csv(path, usecols=lambda c: c in ("sequence", "design_iiptm"))
        df["source_run"] = name
        frames.append(df)
    pool = (pd.concat(frames, ignore_index=True)
            .dropna(subset=["design_iiptm"])
            .sort_values("design_iiptm", ascending=False)
            .drop_duplicates("sequence")
            .reset_index(drop=True))
    log(f"vivier cible vivante : {len(pool)} sequences uniques, "
        f"iiptm {pool['design_iiptm'].min():.3f}-{pool['design_iiptm'].max():.3f}")
    return pool


def annotate(pool: pd.DataFrame, log=print) -> pd.DataFrame:
    """Ajouter CDR3, gene V, identites et verdict de couche 1 a chaque sequence."""
    rows = []
    for r in pool.itertuples():
        reg = regions(r.sequence)
        g = germline_identity(r.sequence)
        cdr3 = reg["CDR3"] if reg else ""
        v_gene = g["V"]["genes"][0] if g else None
        if not cdr3 or not v_gene:
            rows.append({"cdr3": "", "cdr3_length": 0, "v_gene": None,
                         "v_identity": np.nan, "j_identity": np.nan,
                         "layer1": "numerotation indisponible"})
            continue
        # La couche 1 pour son propre creneau : c'est ainsi que le creneau est
        # defini, et l'archive est verifiee plus tard sur les seuls retenus.
        problems = check_sequence(r.sequence, v_gene, len(cdr3))
        rows.append({"cdr3": cdr3, "cdr3_length": len(cdr3), "v_gene": v_gene,
                     "v_identity": g["V"]["identity_pct"],
                     "j_identity": g["J"]["identity_pct"],
                     "layer1": "; ".join(problems) if problems else "ok"})
    out = pd.concat([pool.reset_index(drop=True), pd.DataFrame(rows)], axis=1)
    ok = (out["layer1"] == "ok") & (out["cdr3_length"] >= MIN_CDR3)
    log(f"  conformes couche 1 et CDR3 >= {MIN_CDR3} : {int(ok.sum())} / {len(out)}")
    return out


def select_donors(pool: pd.DataFrame, n: int = N_SLOTS, log=print) -> pd.DataFrame:
    """Les mieux notes qui passent la couche 1, portent un CDR3 >= 12, distincts."""
    eligible = pool[(pool["layer1"] == "ok") & (pool["cdr3_length"] >= MIN_CDR3)]
    picked: list[pd.Series] = []
    for _, cand in eligible.iterrows():
        if any(jaccard_5mer(cand["sequence"], p["sequence"]) >= DONOR_DISTINCT
               for p in picked):
            continue
        picked.append(cand)
        if len(picked) >= n * 2:   # marge : un creneau peut etre inservable
            break
    donors = pd.DataFrame(picked).reset_index(drop=True)
    log(f"  {len(donors)} donneurs candidats retenus (distincts a Jaccard < {DONOR_DISTINCT})")
    for r in donors.itertuples():
        log(f"    iiptm {r.design_iiptm:.4f}  {r.v_gene}  CDR3 {r.cdr3_length}  "
            f"V {r.v_identity:.1f}  J {r.j_identity:.1f}  ({r.source_run})")
    return donors


# ---------------------------------------------------------------------------
# le bras d'etalonnage : des sequences dont le validateur a publie le score
# ---------------------------------------------------------------------------
def calibration_arm(per_band: int = N_CALIB_PER_BAND, seed: int = SEED,
                    log=print) -> pd.DataFrame:
    """Des sequences de l'archive, etalees sur les bandes de score du validateur.

    Elles n'ont jamais ete notees chez nous (aucune des 327 n'est dans
    l'archive), donc chacune est une comparaison directe entre notre instrument
    et le leur, sur la cible vivante, dans ce run.
    """
    path = DATA_DIR / "archive" / f"{CURRENT}_nanobodies.csv"
    if not path.exists():
        raise RuntimeError(
            f"archive absente : {path}. La telecharger avec "
            f"NanobodyArchive('{CURRENT}', with_regions=False).")
    arch = pd.read_csv(path)
    arch = arch[arch["sequence"].str.len().between(90, 150)]
    rng = np.random.default_rng(seed)
    picked: list[dict] = []
    for lo, hi in CALIB_BANDS:
        band = arch[(arch["design_iiptm"] >= lo) & (arch["design_iiptm"] < hi)]
        if band.empty:
            raise RuntimeError(f"bande d'etalonnage {lo}-{hi} vide dans l'archive")
        taken = 0
        for idx in rng.permutation(len(band)):
            row = band.iloc[int(idx)]
            if any(jaccard_5mer(row["sequence"], p["sequence"]) >= DONOR_DISTINCT
                   for p in picked):
                continue
            picked.append({"sequence": row["sequence"], "arm": "calib",
                           "slot": -1, "label": f"calib#{lo:.2f}-{hi:.2f}.{taken}",
                           "group": "calibration",
                           "previous_iiptm": float(row["design_iiptm"]),
                           "previous_run": "archive_validateurs",
                           "archive_iiptm": float(row["design_iiptm"]),
                           "archive_pae": float(row["interaction_pae"]),
                           "band_low": lo, "band_high": hi})
            taken += 1
            if taken >= per_band:
                break
        if taken < per_band:
            raise RuntimeError(f"bande {lo}-{hi} : {taken} sequences distinctes sur {per_band}")
    out = pd.DataFrame(picked)
    log(f"  etalonnage : {len(out)} sequences d'archive, iiptm validateur "
        f"{out['archive_iiptm'].min():.3f}-{out['archive_iiptm'].max():.3f}")
    return out


# ---------------------------------------------------------------------------
# le bras neuf : le generateur du mineur, pour le meme creneau
# ---------------------------------------------------------------------------
def fresh_for_slot(v_gene: str, cdr3_length: int, k: int = K, seed: int = SEED,
                   candidates: int = 1200, log=print) -> list[str]:
    """k candidats neufs pour ce (V, L), par le chemin du mineur.

    Rend une liste vide si l'archive ne peut pas servir le creneau — c'est le
    cas prevu par le design, et le creneau passe alors au donneur suivant.
    """
    from scripts.propose_nanobody import propose
    params = {"epoch": None, "allowed_v": v_gene, "allowed_cdr3_length": cdr3_length}
    try:
        out = propose(params, candidates=candidates, shortlist=k, seed=seed,
                      generator_name=GENERATOR, skip_surrogate=True,
                      archive_cutoff=None, log=lambda *a, **kw: None)
    except Exception as exc:                                   # noqa: BLE001
        log(f"    propose a echoue pour {v_gene}/{cdr3_length} : "
            f"{type(exc).__name__}: {exc}")
        return []
    shortlist = out.get("shortlist")
    if shortlist is None or len(shortlist) < k:
        got = 0 if shortlist is None else len(shortlist)
        log(f"    {v_gene}/{cdr3_length} : {got} candidat(s) neuf(s) sur {k} "
            f"demandes (status {out.get('status')}) — creneau inservable")
        return []
    return list(shortlist["sequence"].head(k))


# ---------------------------------------------------------------------------
# construction
# ---------------------------------------------------------------------------
def build(n_slots: int = N_SLOTS, k: int = K, seed: int = SEED,
          log=print) -> pd.DataFrame:
    calib = calibration_arm(seed=seed, log=log)
    pool = annotate(live_pool(log=log), log=log)
    donors = select_donors(pool, n=n_slots, log=log)

    rows: list[dict] = list(calib.to_dict("records"))
    slot = 0
    for d in donors.itertuples():
        if slot >= n_slots:
            break
        log(f"  creneau {slot} : {d.v_gene} / CDR3 {d.cdr3_length}, "
            f"donneur a {d.design_iiptm:.4f}")
        fresh = fresh_for_slot(d.v_gene, d.cdr3_length, k=k, seed=seed + slot, log=log)
        if not fresh:
            continue
        try:
            # Leve si l'interieur du CDR3 ne porte pas k permutations distinctes.
            perms = permute_cdr3(d.sequence, k, seed=seed + 1000 + slot,
                                 log=lambda *a, **kw: None)
        except ValueError as exc:
            log(f"    permutation impossible ({exc}) — creneau abandonne")
            continue
        rows.append({"sequence": d.sequence, "arm": "donneur", "slot": slot,
                     "label": f"donneur#{slot:02d}", "group": "candidate",
                     "previous_iiptm": d.design_iiptm, "previous_run": d.source_run,
                     "v_gene": d.v_gene, "cdr3_length": d.cdr3_length,
                     "cdr3_original": d.cdr3})
        for i, s in enumerate(perms):
            rows.append({"sequence": s, "arm": "permut", "slot": slot,
                         "label": f"permut#{slot:02d}.{i}", "group": "candidate",
                         "previous_iiptm": np.nan, "previous_run": "",
                         "v_gene": d.v_gene, "cdr3_length": d.cdr3_length,
                         "cdr3_original": d.cdr3})
        for i, s in enumerate(fresh):
            rows.append({"sequence": s, "arm": "neuf", "slot": slot,
                         "label": f"neuf#{slot:02d}.{i}", "group": "candidate",
                         "previous_iiptm": np.nan, "previous_run": "",
                         "v_gene": d.v_gene, "cdr3_length": d.cdr3_length,
                         "cdr3_original": d.cdr3})
        slot += 1

    if slot < n_slots:
        raise RuntimeError(f"seulement {slot} creneaux servables sur {n_slots} demandes")
    batch = pd.DataFrame(rows)

    # Verite sur chaque sequence du lot : couche 1 reellement verifiee (hmmscan
    # present), et la diversite contre tout le reste du lot.
    extra = []
    for r in batch.itertuples():
        reg = regions(r.sequence)
        g = germline_identity(r.sequence)
        # Le bras d'etalonnage n'appartient a aucun creneau : sa couche 1 est
        # lue sans contrainte de gene ni de longueur, qui sont celles de
        # l'epoque ou il a ete soumis et que nous ne connaissons pas.
        v_gene = None if r.arm == "calib" else r.v_gene
        cdr3_len = None if r.arm == "calib" else r.cdr3_length
        problems = check_sequence(r.sequence, v_gene, cdr3_len)
        others = [o for o in batch["sequence"] if o != r.sequence]
        extra.append({
            "cdr3_now": reg["CDR3"] if reg else "",
            "v_identity": g["V"]["identity_pct"] if g else np.nan,
            "j_identity": g["J"]["identity_pct"] if g else np.nan,
            "layer1": "; ".join(problems) if problems else "ok",
            "jaccard_max_batch": max((jaccard_5mer(r.sequence, o) for o in others),
                                     default=0.0),
        })
    batch = pd.concat([batch, pd.DataFrame(extra)], axis=1)

    assert batch["sequence"].is_unique, "une sequence est repetee : BoltzGen la paierait deux fois"
    expected = len(calib) + n_slots * (1 + 2 * k)
    assert len(batch) == expected, f"{len(batch)} sequences, attendu {expected}"

    # Entrelacer : un run interrompu garde les deux bras equilibres.
    order = {"donneur": 0, "permut": 1, "neuf": 2}
    batch["_r"] = batch.groupby(["slot", "arm"]).cumcount()
    batch = (batch.sort_values(["_r", "slot", "arm"], key=lambda c: c.map(order) if c.name == "arm" else c)
             .drop(columns="_r").reset_index(drop=True))

    log(f"\n{len(batch)} sequences : {n_slots} creneaux x (1 donneur + {k} permut + {k} neuf)")
    for arm in ("calib", "donneur", "permut", "neuf"):
        sub = batch[batch["arm"] == arm]
        log(f"  {arm:8s} n={len(sub):2d}  couche 1 ok : {(sub['layer1'] == 'ok').sum()}/{len(sub)}"
            f"  Jaccard max {sub['jaccard_max_batch'].max():.3f}")
    return batch


# ---------------------------------------------------------------------------
# lecture
# ---------------------------------------------------------------------------
def analyse(scores_path: str, batch_path: str, log=print) -> dict:
    """Les 5 lectures pre-enregistrees du docstring."""
    scores = pd.read_csv(scores_path)
    batch = pd.read_csv(batch_path)
    cols = ["sequence"] + [c for c in ("design_iiptm", "design_ptm") if c in scores.columns]
    df = batch.merge(scores[cols], on="sequence", how="left")
    scored = df.dropna(subset=["design_iiptm"])
    log(f"{len(scored)} des {len(df)} sequences notees")
    res: dict = {"n_scored": int(len(scored))}

    permut, neuf = scored[scored["arm"] == "permut"], scored[scored["arm"] == "neuf"]

    # --- 0. l'instrument : nous contre le validateur, sur la cible vivante ---
    calib = scored[scored["arm"] == "calib"]
    if len(calib) >= 6 and "archive_iiptm" in calib.columns:
        rho = float(spearmanr(calib["design_iiptm"], calib["archive_iiptm"]).statistic)
        off = calib["design_iiptm"] - calib["archive_iiptm"]
        res.update(reading0_rho=rho, reading0_offset_median=float(off.median()),
                   reading0_offset_mean=float(off.mean()), reading0_n=int(len(calib)))
        if rho >= CALIB_RHO and abs(off.median()) <= CALIB_OFFSET_OK:
            res["reading0_state"] = "fidele"
        elif rho >= CALIB_RHO and off.median() <= CALIB_OFFSET_DEPRESSED:
            res["reading0_state"] = "monotone mais deprime"
        elif rho >= CALIB_RHO:
            res["reading0_state"] = "monotone, decalage intermediaire"
        else:
            res["reading0_state"] = "ne reproduit pas"
        log(f"\n0. L'INSTRUMENT sur {len(calib)} sequences d'archive : "
            f"rho = {rho:+.3f} (seuil {CALIB_RHO}), decalage median "
            f"{off.median():+.4f} -> {res['reading0_state']}")
        for r in calib.sort_values("archive_iiptm").itertuples():
            log(f"   validateur {r.archive_iiptm:.4f} -> nous {r.design_iiptm:.4f} "
                f"({r.design_iiptm - r.archive_iiptm:+.4f})")
    else:
        res["reading0_state"] = None

    # --- 1. meilleur-de-k, appariee par creneau ---
    deltas = []
    for slot in sorted(scored["slot"].unique()):
        p = permut[permut["slot"] == slot]["design_iiptm"]
        n = neuf[neuf["slot"] == slot]["design_iiptm"]
        if p.empty or n.empty:
            continue
        deltas.append({"slot": int(slot), "best_permut": float(p.max()),
                       "best_neuf": float(n.max()), "delta": float(p.max() - n.max())})
    res["reading1"] = deltas
    if deltas:
        d = np.array([x["delta"] for x in deltas])
        wins = int((d > 0).sum())
        res.update(reading1_mean=float(d.mean()), reading1_median=float(np.median(d)),
                   reading1_wins=wins, reading1_n=len(d))
        log(f"\n1. MEILLEUR-DE-K APPARIE sur {len(d)} creneaux")
        for x in deltas:
            log(f"   creneau {x['slot']} : permut {x['best_permut']:.4f}  "
                f"neuf {x['best_neuf']:.4f}  delta {x['delta']:+.4f} "
                f"{'permut' if x['delta'] > 0 else 'neuf'}")
        log(f"   moyenne {d.mean():+.4f}   mediane {np.median(d):+.4f}   "
            f"creneaux gagnes par permut {wins}/{len(d)}   (seuil {RIVAL_GAP})")

    # --- 2. distribution complete ---
    if len(permut) >= 5 and len(neuf) >= 5:
        u = mannwhitneyu(permut["design_iiptm"], neuf["design_iiptm"], alternative="two-sided")
        res.update(reading2_p=float(u.pvalue),
                   reading2_median_diff=float(permut["design_iiptm"].median()
                                              - neuf["design_iiptm"].median()))
        log(f"\n2. DISTRIBUTION : mediane permut {permut['design_iiptm'].median():.4f} "
            f"contre neuf {neuf['design_iiptm'].median():.4f} "
            f"(ecart {res['reading2_median_diff']:+.4f}), Mann-Whitney p = {u.pvalue:.3f}")

    # --- 3. la queue haute ---
    if len(permut) and len(neuf):
        res.update(reading3_permut_above=int((permut["design_iiptm"] >= INCUMBENT).sum()),
                   reading3_neuf_above=int((neuf["design_iiptm"] >= INCUMBENT).sum()),
                   reading3_permut_max=float(permut["design_iiptm"].max()),
                   reading3_neuf_max=float(neuf["design_iiptm"].max()))
        log(f"\n3. QUEUE HAUTE (seuil {INCUMBENT}, meilleur soumettable connu) : "
            f"permut {res['reading3_permut_above']}/{len(permut)} au-dessus (max "
            f"{res['reading3_permut_max']:.4f}), neuf {res['reading3_neuf_above']}/{len(neuf)} "
            f"(max {res['reading3_neuf_max']:.4f})")

    # --- 4. temoin : bruit entre runs ---
    don = scored[(scored["arm"] == "donneur") & scored["previous_iiptm"].notna()]
    if len(don) >= 3:
        delta = don["design_iiptm"] - don["previous_iiptm"]
        res.update(reading4_mean_abs=float(delta.abs().mean()),
                   reading4_mean_signed=float(delta.mean()),
                   reading4_pass=bool(delta.abs().mean() < SIGMA_2
                                      and abs(delta.mean()) < RIVAL_GAP))
        log(f"\n4. TEMOIN inter-runs sur {len(don)} donneurs : moyenne|delta| "
            f"{delta.abs().mean():.4f} (seuil {SIGMA_2}), moyenne signee "
            f"{delta.mean():+.4f} (seuil {RIVAL_GAP}) -> "
            f"{'PASSE' if res['reading4_pass'] else 'ECHEC'}")
    else:
        res["reading4_pass"] = None

    # --- 5. soumettabilite et diversite ---
    res["reading5"] = {}
    for arm in ("donneur", "permut", "neuf"):
        sub = scored[scored["arm"] == arm]
        if sub.empty:
            continue
        res["reading5"][arm] = {
            "layer1_ok": int((sub["layer1"] == "ok").sum()), "n": int(len(sub)),
            "jaccard_max": float(sub["jaccard_max_batch"].max()),
        }
    res["reading5_diversity_pass"] = bool(scored["jaccard_max_batch"].max() < JACCARD_CEILING)
    log("\n5. SOUMETTABILITE ET DIVERSITE")
    for arm, v in res["reading5"].items():
        log(f"   {arm:8s} couche 1 {v['layer1_ok']}/{v['n']}   "
            f"Jaccard max dans le lot {v['jaccard_max']:.3f}")

    log("\n" + "=" * 60)
    for line in verdict(res):
        log(line)
    return res


def verdict(res: dict) -> list[str]:
    """Le verdict, derive des seuils pre-enregistres et de rien d'autre."""
    out = []
    state = res.get("reading0_state")
    if state == "fidele":
        out.append(f"INSTRUMENT FIDELE sur la cible vivante : rho {res['reading0_rho']:+.3f}, "
                   f"decalage {res['reading0_offset_median']:+.4f}. Nos candidats sont donc "
                   f"vraiment dans le dernier decile du champ, et le §85 est plus severe "
                   f"qu'ecrit : la ou nous plafonnons a 0,446, 92,9 % de l'archive est "
                   f"au-dessus de 0,45.")
    elif state == "monotone mais deprime":
        out.append(f"INSTRUMENT MONOTONE MAIS DEPRIME : rho {res['reading0_rho']:+.3f}, "
                   f"decalage {res['reading0_offset_median']:+.4f}. Les rangs tiennent, "
                   f"mais tout seuil absolu et toute comparaison au champ sont a refaire — "
                   f"y compris le 0,178 de 61H7 au §77.")
    elif state == "monotone, decalage intermediaire":
        out.append(f"INSTRUMENT MONOTONE, decalage intermediaire "
                   f"{res['reading0_offset_median']:+.4f} (entre {CALIB_OFFSET_DEPRESSED} "
                   f"et {CALIB_OFFSET_OK}) : ni fidele ni franchement deprime. A trancher "
                   f"avec plus de sequences avant d'en tirer quoi que ce soit d'absolu.")
    elif state == "ne reproduit pas":
        out.append(f"L'INSTRUMENT NE REPRODUIT PAS LE VALIDATEUR sur cette cible : "
                   f"rho {res['reading0_rho']:+.3f} < {CALIB_RHO}. Les lectures absolues "
                   f"des §77, §87, §89 et §90 doivent etre relues, et les lectures 1 a 3 "
                   f"ci-dessous ne valent plus que comme comparaisons internes au run.")
    mean, wins, n = (res.get("reading1_mean"), res.get("reading1_wins"),
                     res.get("reading1_n"))
    if mean is None:
        out.append("LECTURE 1 IMPOSSIBLE : creneaux incomplets.")
    elif mean > RIVAL_GAP and wins >= WINS_FOR_PERMUT:
        out.append(f"LA PERMUTATION GAGNE : +{mean:.4f} par creneau ({wins}/{n} creneaux). "
                   f"Permuter un donneur note vaut mieux que generer, a budget egal.")
    elif mean < -RIVAL_GAP and wins <= WINS_FOR_FRESH:
        out.append(f"LES CANDIDATS NEUFS GAGNENT : {mean:+.4f} par creneau ({wins}/{n}). "
                   f"La variance de la permutation ne paie pas.")
    else:
        out.append(f"AUCUNE DIFFERENCE A CE BUDGET : {mean:+.4f} par creneau "
                   f"({wins}/{n} creneaux), sous l'ecart entre rivaux ({RIVAL_GAP}). "
                   f"Preferer les candidats neufs, qui n'exigent pas un donneur deja paye.")
    p = res.get("reading2_p")
    if p is not None:
        out.append(f"Distribution : {'difference' if p < MWU_ALPHA else 'pas de difference'} "
                   f"(p = {p:.3f}, ecart des medianes {res['reading2_median_diff']:+.4f}).")
    if res.get("reading3_permut_max") is not None:
        best_arm = ("permut" if res["reading3_permut_max"] >= res["reading3_neuf_max"]
                    else "neuf")
        best = max(res["reading3_permut_max"], res["reading3_neuf_max"])
        out.append(f"Queue haute : le meilleur du lot est {best:.4f} ({best_arm}) ; "
                   f"{'nouveau record soumettable' if best > INCUMBENT else 'sous le record soumettable'} "
                   f"({INCUMBENT}).")
    if res.get("reading4_pass") is False:
        out.append("TEMOIN EN ECHEC : ce run n'est pas comparable aux precedents. "
                   "Les lectures 1 a 3 tiennent (internes au run), la 4 non.")
    if res.get("reading5_diversity_pass") is False:
        out.append("DIVERSITE EN ECHEC : une sequence depasse le plafond de Jaccard 0,90.")
    broken = [a for a, v in res.get("reading5", {}).items() if v["layer1_ok"] < v["n"]]
    if broken:
        out.append(f"COUCHE 1 CASSEE dans : {', '.join(broken)} — "
                   f"ces sequences ne sont pas soumettables.")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="inputs/kaggle_equalbudget.csv")
    ap.add_argument("--analyse")
    ap.add_argument("--batch", default="inputs/kaggle_equalbudget.csv")
    ap.add_argument("--slots", type=int, default=N_SLOTS)
    ap.add_argument("--k", type=int, default=K)
    ap.add_argument("--seed", type=int, default=SEED)
    args = ap.parse_args()

    if args.analyse:
        analyse(args.analyse, args.batch)
        return
    batch = build(n_slots=args.slots, k=args.k, seed=args.seed)
    path = REPO_ROOT / args.out
    batch.to_csv(path, index=False)
    print(f"\necrit {path} ({len(batch)} sequences)")


if __name__ == "__main__":
    main()
