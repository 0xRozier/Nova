"""Does a protein language model predict the score better than 3-mers?

FINDINGS 83 left one untried axis on the prediction side. Our features are hand
descriptors and 3-mer counts, and they reach rho 0.260 on the live target --
barely above zero with a fold spread of 0.138. The question here is whether the
limit is the FEATURES or the problem: an ESM-2 embedding knows about protein
structure in a way a 3-mer count cannot, and if the oracle's behaviour has any
structural basis, that is where it would show.

Honest expectation, stated before running: low. The live target has 273 queries
and an embedding has hundreds of dimensions, so the model has more features than
samples and ridge regularisation is doing most of the work. And FINDINGS 77
showed the score does not read binding, so there may be no structural signal to
find. But it costs no GPU and it is the last idea on this axis.

Three representations, because they answer different questions:

  whole   mean-pooled over the entire sequence: does the model see anything?
  cdr3    mean-pooled over the CDR3 only. FINDINGS 77 measured that the CDR3's
          identity moves the score by 6.9x what permuting its interior does, so
          if any region carries the signal it is this one -- and pooling over
          120 residues to predict something driven by 11 of them dilutes it.
  both    concatenated.

Compared against the same grouped cross-validation as
scripts.train_score_surrogate (Jaccard 5-mer >= 0.80, whole clusters held out),
so the numbers are directly comparable to its rho 0.260.

Usage:
    python -m scripts.esm_surrogate                      # the live target
    python -m scripts.esm_surrogate --target 8AOK --limit 400
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from nova_lab import DATA_DIR, REPO_ROOT
from nova_lab.config import NANOBODY_STRUCTURE, PREVIOUS_NANOBODY
from nova_lab.nanobodies import regions
from scripts.train_score_surrogate import (build_features, cluster,
                                           grouped_spearman, load,
                                           noise_ceiling)

CACHE = DATA_DIR / "esm"
# esm2_t12_35M is 150 MB and runs on CPU in seconds; t30_150M is 600 MB and
# about four times slower. With 273 training points the bigger model's extra
# dimensions are not obviously an advantage, so the small one is the default
# and the comparison is available rather than assumed.
MODELS = {"t12": "esm2_t12_35M_UR50D", "t30": "esm2_t30_150M_UR50D"}
SEED = 86


def embed(sequences: list[str], model_name: str, batch: int = 8) -> dict[str, np.ndarray]:
    """{sequence: per-residue representation}. Cached: the download is 150 MB."""
    import torch
    import esm

    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{model_name}.npz"
    have: dict[str, np.ndarray] = {}
    if path.exists():
        with np.load(path, allow_pickle=True) as z:
            have = {k: z[k] for k in z.files}
    missing = [s for s in sequences if s not in have]
    if not missing:
        return have

    torch.set_num_threads(max(1, (__import__("os").cpu_count() or 2)))
    model, alphabet = getattr(esm.pretrained, model_name)()
    model.eval()
    converter = alphabet.get_batch_converter()
    layer = model.num_layers
    with torch.no_grad():
        for start in range(0, len(missing), batch):
            chunk = missing[start:start + batch]
            _, _, tokens = converter([(str(i), s) for i, s in enumerate(chunk)])
            out = model(tokens, repr_layers=[layer])["representations"][layer]
            for i, s in enumerate(chunk):
                # drop the BOS token and the padding: keep exactly len(s) rows
                have[s] = out[i, 1:len(s) + 1].numpy().astype(np.float32)
            print(f"  {min(start + batch, len(missing))}/{len(missing)} encodees",
                  end="\r", flush=True)
    print()
    np.savez_compressed(path, **have)
    return have


def pools(sequences: list[str], reps: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Mean over the whole sequence, and mean over the CDR3 alone."""
    whole, cdr3 = [], []
    width = next(iter(reps.values())).shape[1]
    for s in sequences:
        r = reps[s]
        whole.append(r.mean(axis=0))
        region = (regions(s) or {}).get("CDR3", "")
        at = s.find(region) if region else -1
        # A sequence whose CDR3 cannot be located gets zeros rather than the
        # whole-sequence mean: a silent fallback would make the cdr3 column a
        # copy of the whole column for those rows and flatter the comparison.
        cdr3.append(r[at:at + len(region)].mean(axis=0) if at >= 0 else np.zeros(width))
    return {"whole": np.array(whole), "cdr3": np.array(cdr3),
            "both": np.hstack([np.array(whole), np.array(cdr3)])}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default=NANOBODY_STRUCTURE,
                        choices=[NANOBODY_STRUCTURE, PREVIOUS_NANOBODY["structure"]])
    parser.add_argument("--model", default="t12", choices=sorted(MODELS))
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--limit", type=int, help="subsample, for the big target")
    args = parser.parse_args()

    frame = load(args.target)
    if args.limit and len(frame) > args.limit:
        frame = frame.sample(args.limit, random_state=SEED).reset_index(drop=True)
    sigma, ceiling = noise_ceiling(frame)
    print(f"cible {args.target} : {len(frame)} sequences, plafond de Spearman {ceiling:.3f}")

    print(f"encodage par {MODELS[args.model]} (CPU)...")
    reps = embed(list(frame.sequence), MODELS[args.model])
    represented = pools(list(frame.sequence), reps)
    y = frame.y.to_numpy()
    groups = cluster(list(frame.sequence))
    print(f"{len(np.unique(groups))} grappes a Jaccard 5-mer >= 0,80")

    from sklearn.ensemble import HistGradientBoostingRegressor
    from sklearn.decomposition import PCA
    from sklearn.linear_model import RidgeCV
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    ridge = lambda: make_pipeline(StandardScaler(),
                                  RidgeCV(alphas=np.logspace(-1, 6, 30)))
    # PCA first for the tree: a 480-dimensional input with 273 rows is mostly
    # noise to a histogram splitter.
    reduced = lambda k: make_pipeline(StandardScaler(),
                                      PCA(n_components=k, random_state=SEED),
                                      HistGradientBoostingRegressor(
                                          max_iter=400, learning_rate=0.05,
                                          random_state=SEED))

    desc, _, kmer, _ = build_features(frame)
    baseline = np.hstack([desc, kmer])
    trials = [("reference : descripteurs + 3-mers",
               HistGradientBoostingRegressor(max_iter=500, learning_rate=0.05,
                                             random_state=77), baseline)]
    for name, X in represented.items():
        trials.append((f"ridge sur ESM {name}", ridge(), X))
        trials.append((f"GBM sur ESM {name} (PCA 30)", reduced(min(30, X.shape[1])), X))
    trials.append(("ridge sur ESM both + 3-mers", ridge(),
                   np.hstack([represented["both"], kmer])))

    print(f"\n{'representation':<38} {'rho':>7} {'sd plis':>9} {'dim':>6}")
    results = []
    for label, model, X in trials:
        rho, sd = grouped_spearman(X, y, groups, model, args.folds)
        results.append((label, rho, sd, X.shape[1]))
        print(f"{label:<38} {rho:+7.3f} {sd:9.3f} {X.shape[1]:>6}")

    best = max(results, key=lambda r: r[1])
    reference = results[0][1]
    print(f"\nmeilleur : {best[0]} a rho {best[1]:+.3f}")
    print(f"reference (descripteurs + 3-mers) : {reference:+.3f}")
    verdict = ("l'ESM apporte quelque chose" if best[1] > reference + 2 * best[2]
               else "l'ESM n'apporte rien au-dela de l'ecart entre plis")
    print(f"-> {verdict}")
    out = REPO_ROOT / "results" / f"esm_surrogate_{args.target}_{args.model}.csv"
    pd.DataFrame(results, columns=["representation", "rho", "sd_plis", "dim"]).to_csv(
        out, index=False)
    print(f"-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
