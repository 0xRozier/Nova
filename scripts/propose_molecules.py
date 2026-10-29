"""Propose a 20-molecule submission for an epoch, on CPU.

1. Sample random reagent combinations for the epoch's reaction and score the
   products with the surrogate.
2. Refine for a few rounds: resample around the reagents that produced the best
   predicted products (a greedy cousin of Thompson sampling).
3. Keep only candidates that pass every validator check, including the
   novelty rule against the public archive (Tanimoto < 0.5).
4. Pick the 20 maximizing the predicted total under the entropy minimum.

The surrogate is only a pre-filter: before going live, the shortlist should
be re-scored with Boltz-2 on a GPU (see --shortlist).

Usage:
    python -m scripts.propose_molecules --reaction rxn:41
    python -m scripts.propose_molecules --reaction rxn:41 --samples 50000 --rounds 4
"""
import argparse
import random
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from rdkit import Chem

from nova_lab import DATA_DIR
from nova_lab.chem import ATOM_PAIR, MORGAN, HistoricalIndex, check_molecule, check_submission, fingerprint_entropy
from nova_lab.config import MAX_SIMILARITY_TO_HISTORICAL, MIN_ENTROPY, NUM_MOLECULES, SMALL_MOLECULE_TARGET
from nova_lab.features import featurize
from nova_lab.filters import first_failure
from nova_lab.surrogate import Surrogate
from nova_lab.reactions import make_name, molecules_for_role, products, reaction_info


MAX_K = 4  # alternative products per combination (name suffix :k<k>)


def _make(combo):
    """Every filter-passing product of one reagent combination."""
    rxn_id, mol_ids = combo
    out = []
    for k, smi in enumerate(products(rxn_id, list(mol_ids))[:MAX_K]):
        mol = Chem.MolFromSmiles(smi)
        if mol is not None and not first_failure(mol):
            out.append((make_name(rxn_id, list(mol_ids), k), smi, featurize(smi), tuple(mol_ids)))
    return out


def evaluate(pool, model, rxn_id, combos, seen):
    combos = [c for c in combos if c not in seen]
    seen.update(combos)
    out = [p for batch in pool.map(_make, [(rxn_id, c) for c in combos], chunksize=200) for p in batch]
    if not out:
        return pd.DataFrame(columns=["name", "smiles", "pred", "combo"])
    X = np.vstack([x for _, _, x, _ in out])
    return pd.DataFrame({
        "name": [n for n, _, _, _ in out],
        "smiles": [s for _, s, _, _ in out],
        "pred": model.predict(X),
        "combo": [c for _, _, _, c in out],
    })


def refine_combos(scored: pd.DataFrame, roles: list[list[int]], n: int, rng: random.Random, top_frac=0.02):
    """New combos mixing reagents seen in the best products with random ones."""
    top = scored.nlargest(max(20, int(len(scored) * top_frac)), "pred")
    elite = [sorted({c[i] for c in top["combo"]}) for i in range(len(roles))]
    combos = set()
    while len(combos) < n:
        # Keep most positions from the elite pool, occasionally explore one position.
        explore = rng.randrange(len(roles)) if rng.random() < 0.5 else None
        combos.add(tuple(rng.choice(roles[i]) if i == explore else rng.choice(elite[i]) for i in range(len(roles))))
    return list(combos)


def select(candidates: pd.DataFrame, history: HistoricalIndex, allowed: str, shortlist: int,
           max_reuse: int = 2) -> pd.DataFrame:
    """Walk candidates by predicted score, keep the ones the validator would accept.

    Each reagent may appear in at most `max_reuse` kept candidates: without
    the cap the shortlist collapses onto a few elite reagents and no 20 of
    them reach the entropy minimum.
    """
    kept, inchikeys, uses = [], set(), {}
    for row in candidates.sort_values("pred", ascending=False).itertuples(index=False):
        reagents = list(enumerate(row.combo))
        if any(uses.get(r, 0) >= max_reuse for r in reagents):
            continue
        key = Chem.MolToInchiKey(Chem.MolFromSmiles(row.smiles))
        if key in inchikeys:
            continue
        ok, _ = check_molecule(row.name, allowed, history, row.smiles)
        if not ok:
            continue
        kept.append(row)
        inchikeys.add(key)
        for r in reagents:
            uses[r] = uses.get(r, 0) + 1
        if len(kept) >= shortlist:
            break
    return pd.DataFrame(kept)


def _entropy(bit_counts: np.ndarray, n: int) -> float:
    p = bit_counts / n
    with np.errstate(divide="ignore", invalid="ignore"):
        h = np.where((p > 0) & (p < 1), -p * np.log2(p) - (1 - p) * np.log2(1 - p), 0.0)
    return float(h.mean())


def pick_submission(shortlist: pd.DataFrame, margin: float = 0.01) -> pd.DataFrame | None:
    """Maximize the predicted total subject to the validator's entropy minimum.

    Near-analogues share most fingerprint bits and sink the entropy, so the
    greedy pick trades predicted score against entropy gain with weight `lam`;
    the best total among the weights that clear MIN_ENTROPY + margin wins.
    """
    fps = np.array([ATOM_PAIR.GetFingerprintAsNumPy(Chem.MolFromSmiles(s)) for s in shortlist["smiles"]], float)
    pred = shortlist["pred"].to_numpy()
    best, best_total = None, -np.inf
    for lam in np.concatenate([[0.0], np.geomspace(0.05, 50, 40)]):
        chosen, counts = [], np.zeros(fps.shape[1])
        for step in range(NUM_MOLECULES):
            free = np.setdiff1d(np.arange(len(pred)), chosen)
            if step == 0:
                gain = np.zeros(len(free))
            else:
                base = _entropy(counts, step)
                gain = np.array([_entropy(counts + fps[i], step + 1) for i in free]) - base
            pick = free[np.argmax(pred[free] + lam * gain)]
            chosen.append(pick)
            counts += fps[pick]
        if _entropy(counts, NUM_MOLECULES) >= MIN_ENTROPY + margin and pred[chosen].sum() > best_total:
            best, best_total = chosen, pred[chosen].sum()
    if best is None:
        return None
    return shortlist.iloc[_swap_improve(best, pred, fps, MIN_ENTROPY + margin)]


def _swap_improve(chosen: list, pred: np.ndarray, fps: np.ndarray, threshold: float) -> list:
    """Swap weak members for stronger outsiders while the entropy stays above threshold."""
    chosen = list(chosen)
    counts = fps[chosen].sum(axis=0)
    improved = True
    while improved:
        improved = False
        outsiders = [j for j in np.argsort(-pred) if j not in chosen]
        for j in outsiders:
            for i in sorted(chosen, key=lambda k: pred[k]):
                if pred[j] <= pred[i]:
                    break
                trial = counts - fps[i] + fps[j]
                if _entropy(trial, len(chosen)) >= threshold:
                    chosen[chosen.index(i)] = j
                    counts = trial
                    improved = True
                    break
    return chosen


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reaction", required=True, help="e.g. rxn:41 (see scripts.epoch_info)")
    parser.add_argument("--target", default=SMALL_MOLECULE_TARGET[0])
    parser.add_argument("--samples", type=int, default=30000, help="combos per round")
    parser.add_argument("--rounds", type=int, default=3)
    parser.add_argument("--shortlist", type=int, default=300)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    rxn_id = int(args.reaction.split(":")[1])
    name, _, role_a, role_b, role_c = reaction_info(rxn_id)
    roles = [[m for m, _ in molecules_for_role(r)] for r in (role_a, role_b, role_c) if r]
    space = np.prod([len(r) for r in roles], dtype=float)
    print(f"{args.reaction} ({name}): reagents per role {[len(r) for r in roles]}, ~{space:.2e} products")

    model = Surrogate.load(args.target)
    history = HistoricalIndex(args.target, refresh=True)
    print(f"archive: {len(history.fps)} molecules already submitted for {args.target}")

    seen, scored = set(), []
    with Pool(args.workers) as pool:
        combos = list({tuple(rng.choice(r) for r in roles) for _ in range(args.samples)})
        for round_ in range(args.rounds):
            t = time.time()
            batch = evaluate(pool, model, rxn_id, combos, seen)
            scored.append(batch)
            allc = pd.concat(scored, ignore_index=True)
            print(f"round {round_}: {len(batch)} valid products in {time.time() - t:.0f}s, "
                  f"best pred {allc['pred'].max():.4f}, top-20 mean {allc['pred'].nlargest(20).mean():.4f}")
            combos = refine_combos(allc, roles, args.samples, rng)

    shortlist = select(allc, history, args.reaction, args.shortlist)
    if len(shortlist) < NUM_MOLECULES:
        print(f"only {len(shortlist)} valid novel candidates; increase --samples/--rounds")
        return
    chosen = pick_submission(shortlist)
    if chosen is None:
        print(f"could not reach entropy {MIN_ENTROPY} with {len(shortlist)} candidates; raise --shortlist/--samples")
        return
    problems = check_submission(list(chosen["name"]), args.reaction, history)
    if problems:
        print("submission would be rejected:", problems)
        return

    out_dir = DATA_DIR / "proposals"
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{args.target}_{args.reaction.replace(':', '')}"
    shortlist.drop(columns=["combo"]).to_csv(out_dir / f"{stem}_shortlist.csv", index=False)
    (out_dir / f"{stem}_submission.txt").write_text(",".join(chosen["name"]) + "|~\n")

    sims = [history.max_similarity(MORGAN.GetFingerprint(Chem.MolFromSmiles(s))) for s in chosen["smiles"]]
    print(f"\nproposed {len(chosen)} molecules: predicted total {chosen['pred'].sum():.3f}, "
          f"entropy {fingerprint_entropy(list(chosen['smiles'])):.3f}, "
          f"max archive similarity {max(sims):.2f} (limit {MAX_SIMILARITY_TO_HISTORICAL})")
    print(chosen[["name", "pred"]].to_string(index=False))
    print(f"\nwritten to {out_dir}/{stem}_submission.txt (shortlist of {len(shortlist)} for GPU re-scoring)")


if __name__ == "__main__":
    main()
