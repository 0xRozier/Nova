"""The NOVA nanobody miner: generate, screen, score, encrypt, publish, commit.

The loop per epoch, in order:

  1. resolve the epoch's imposed V gene and CDR3 length from the chain;
  2. refuse to submit if the subnet's own config has moved under us
     (scripts.check_config) — a changed target makes every candidate worthless;
  3. generate and keep what the offline checker accepts
     (scripts.propose_nanobody);
  4. measure PSH on a LARGE pool and keep the least hydrophobic, then run the
     rest of the developability gate on the survivors. The validator refuses
     PSH > 89 outright, which rejected all 11 production submissions, and PSH
     cannot be predicted from the sequence: the best proxies reach Spearman
     0.30 on the length-corrected residual, swinging from -0.22 to +0.58
     depending on CDR3 length. So it has to be measured. It is cheap — 3.6 s
     per candidate on 4 cores — and only 21% of our candidates come in under
     87.4, so the pool has to be wide for the screen to yield anything;
  5. rank the survivors with BoltzGen, not the surrogate. Measured over 21
     epochs: picking by BoltzGen wins 48% of them, picking by the surrogate
     14% — worse than picking at random, which wins 23%. The surrogate is a
     validity pre-filter, never a ranker;
  6. timelock-encrypt to the round that unlocks at the submission deadline,
     publish the file to GitHub, and commit the pointer on chain.

Nothing is sent anywhere without `--live`. The default run does every step,
writes the file it would publish, and stops before the chain and GitHub.

Environment (same names the reference miner uses):
    GITHUB_REPO_OWNER, GITHUB_REPO_NAME, GITHUB_REPO_BRANCH, GITHUB_REPO_PATH,
    GITHUB_TOKEN

Usage:
    python -m scripts.mine --wallet my --hotkey default            # dry run
    python -m scripts.mine --wallet my --hotkey default --live     # for real
    python -m scripts.mine --wallet my --hotkey default --once     # one epoch
"""
import argparse
import base64
import json
import os
import subprocess
import sys
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

from nova_lab import DATA_DIR
from nova_lab.cdr3_generators import GENERATORS
from nova_lab.config import (EPOCH_LENGTH, NANOBODY_TARGET, NO_SUBMISSION_BLOCKS,
                             PSH_STRICT_REJECTION_THRESHOLD as PSH_STRICT)
from nova_lab import timelock
from nova_lab.chain import current_block, missing_api, require_api
from nova_lab import psh_calibration
from nova_lab.psh_calibration import (SIGMA as PSH_SIGMA, THRESHOLD as PSH_THRESHOLD,
                                      refusal_probability)

JOURNAL = DATA_DIR / "mining" / "journal.jsonl"
BLOCK_SECONDS = 12
# Leave room for encryption, the GitHub PUT and the extrinsic. Measured: BoltzGen
# runs at ~1.6 min/sequence on 2x T4, so the candidate budget follows from this.
SAFETY_BLOCKS = 12
# Share of FRESH generator candidates reading at or under each PSH, measured on
# the install the validator runs: 100 candidates per length, one fold each,
# results/psh_fresh_candidates_by_length.csv.
#
#   CDR3   median PSH   <= 87   <= 89   <= 96
#      8        85.0     62 %    75 %    95 %
#     11        89.0     40 %    50 %    81 %
#     18       100.6     10 %    15 %    31 %
#
# The table this replaces was measured before mkdssp 3.1.4 and the filter's
# process_pdb patch were in place, and read 51 / 28 / 8.9 % for the <= 89
# column. The shape survives -- long CDR3 is where the generator struggles --
# but every rate was understated, and the 96 column it sized the pool against
# does not correspond to anything the validator does.
#
# What has not changed: at CDR3 18 the median candidate reads 100.6 while the
# accepted field reads 83.5. That is a generator problem and no pool size
# fixes it; it only buys the tail.
PSH_PASS_RATE = {8: 0.75, 11: 0.50, 18: 0.15}
# Candidates wanted under the validator's threshold before BoltzGen: the
# shortlist of 8, plus slack for the developability screen, which rejected
# 1 of 8 on epoch 26238.
WANT_UNDER_THRESHOLD = 10


def psh_pass_rate(cdr3_length: int, table: dict) -> float:
    """Interpolate the measured rate, clamped at the ends."""
    lengths = sorted(table)
    if cdr3_length <= lengths[0]:
        return table[lengths[0]]
    if cdr3_length >= lengths[-1]:
        return table[lengths[-1]]
    for lo, hi in zip(lengths, lengths[1:]):
        if lo <= cdr3_length <= hi:
            f = (cdr3_length - lo) / (hi - lo)
            return table[lo] + f * (table[hi] - table[lo])
    return table[lengths[-1]]


PSH_SECONDS_PER_CANDIDATE = 3.5  # 34 sequences in 2 min on 4 cores, measured


def psh_pool_for(cdr3_length: int, override: int | None,
                 ceiling: int) -> tuple[int, float]:
    """Pool size that should yield WANT_UNDER_THRESHOLD candidates under 89.

    Sized against the threshold the validator applies, not against a ceiling of
    our own: at CDR3 18 only 15 % of fresh candidates clear 89, so ten of them
    want ~67 measured, which is 4 minutes of the epoch's 66. At CDR3 8 the rate
    is 75 % and thirteen candidates suffice. Going wider than that buys the
    lower tail, which is worth something here -- unlike under the old banded
    curve, where everything below 96 was indistinguishable -- but it buys it
    slowly, and the epoch's budget is the GPU's, not TNP's.
    """
    rate = psh_pass_rate(cdr3_length, PSH_PASS_RATE)
    if override:
        return min(override, ceiling), rate
    need = int(WANT_UNDER_THRESHOLD / max(rate, 0.05))
    return max(WANT_UNDER_THRESHOLD, min(need, ceiling)), rate


def log(msg: str) -> None:
    print(f"[{datetime.now(timezone.utc):%H:%M:%S}] {msg}", flush=True)


@contextmanager
def timed(label: str, timings: dict):
    """Log and record elapsed seconds for a pipeline stage."""
    t0 = time.monotonic()
    log(f"[{label}] start")
    try:
        yield
    finally:
        elapsed = time.monotonic() - t0
        timings[label] = round(elapsed, 1)
        log(f"[{label}] done in {elapsed:.0f}s")


def retry(fn, *args, max_retries: int = 5, base_delay: float = 30, **kwargs):
    """Call *fn* with retries on transient network errors (429, WebSocket, timeout).

    Exponential backoff: base_delay * 2^attempt seconds between attempts.
    """
    for attempt in range(max_retries + 1):
        try:
            return fn(*args, **kwargs)
        except Exception as exc:
            msg = f"{type(exc).__name__}: {exc}"
            transient = any(tok in msg for tok in ("429", "ConnectionClosed",
                                                   "TimeoutError", "WebSocket",
                                                   "Timeout", "timed out"))
            if not transient or attempt == max_retries:
                raise
            delay = base_delay * (2 ** attempt)
            log(f"transient error (attempt {attempt + 1}/{max_retries + 1}): "
                f"{msg}; retrying in {delay:.0f}s")
            time.sleep(delay)


def github_path() -> str:
    """The `owner/repo/branch[/path]` prefix the commitment points at."""
    owner = os.environ.get("GITHUB_REPO_OWNER")
    name = os.environ.get("GITHUB_REPO_NAME")
    branch = os.environ.get("GITHUB_REPO_BRANCH")
    inner = os.environ.get("GITHUB_REPO_PATH", "")
    if not (owner and name and branch):
        sys.exit("set GITHUB_REPO_OWNER, GITHUB_REPO_NAME and GITHUB_REPO_BRANCH")
    path = f"{owner}/{name}/{branch}"
    if inner:
        path = f"{path}/{inner.strip('/')}"
    if len(path) > 100:
        sys.exit(f"github path is {len(path)} characters; the chain caps the commitment at 100")
    return path


def upload(filename: str, content: str) -> None:
    """PUT the submission file through the GitHub contents API."""
    owner = os.environ["GITHUB_REPO_OWNER"]
    name = os.environ["GITHUB_REPO_NAME"]
    branch = os.environ["GITHUB_REPO_BRANCH"]
    inner = os.environ.get("GITHUB_REPO_PATH", "").strip("/")
    token = os.environ.get("GITHUB_TOKEN")
    if not token:
        sys.exit("GITHUB_TOKEN is not set")
    target = f"{inner}/{filename}" if inner else filename
    url = f"https://api.github.com/repos/{owner}/{name}/contents/{target}"
    headers = {"Authorization": f"Bearer {token}", "Accept": "application/vnd.github+json"}
    existing = requests.get(url, headers=headers, params={"ref": branch}, timeout=15)
    payload = {"message": f"Encrypted submission {filename}",
               "content": base64.b64encode(content.encode()).decode(), "branch": branch}
    if existing.status_code == 200:
        payload["sha"] = existing.json().get("sha")
    resp = requests.put(url, headers=headers, json=payload, timeout=30)
    if resp.status_code not in (200, 201):
        raise RuntimeError(f"GitHub upload failed: {resp.status_code} {resp.text[:200]}")


def config_unchanged() -> bool:
    """Has the subnet moved the rules? A changed target voids every candidate."""
    proc = subprocess.run([sys.executable, "-m", "scripts.check_config"],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        log("the subnet config has moved — refusing to submit:")
        for line in proc.stdout.splitlines():
            if "DIFF" in line or "- " in line:
                log(f"  {line.strip()}")
    return proc.returncode == 0


def measure_psh(sequences: list[str], work: Path) -> dict[str, float]:
    """One PSH measurement per sequence. A fresh work dir forces a new fold.

    Sequences whose model TNP could not annotate with DSSP are dropped rather
    than returned: their PSH is not the validator's quantity, and a reading
    that is not comparable is worse than no reading, because it looks like one.
    """
    from scripts.tnp_developability import parse, run
    work.mkdir(parents=True, exist_ok=True)
    names = {f"cand{i}": s for i, s in enumerate(sequences)}
    (work / "candidates.fasta").write_text(
        "".join(f">{n}\n{s}\n" for n, s in names.items()))
    df = parse(run(work / "candidates.fasta", 4, work), work)
    if "dssp" in df.columns:
        blind = int((~df["dssp"]).sum())
        if blind:
            log(f"  {blind} sequences have a PSH measured without DSSP; "
                f"discarding those readings (run scripts.doctor)")
            df = df[df["dssp"]]
    return {names[n]: v for n, v in zip(df["name"], df["PSH"]) if n in names}


# psh_select was removed here. It returned the `keep` LOWEST PSH, which made
# PSH the selection criterion; psh_admissible replaced it with a ceiling, which
# threw away strong candidates one draw over the line. Neither is here now:
# psh_screen measures and psh_calibration prices, and the two multiply. Leaving
# either in place as dead code would invite it back.


def surrogate_rank(sequences: list[str], target: str) -> dict[str, float] | None:
    """Predicted design_iiptm per candidate, best first (an ordered dict).

    This is the quality half of the selection, and it used to be the whole of
    it in one direction or PSH in the other. Over 84 grafted candidates scored
    by BoltzGen the best read design_iiptm 0.692, above the field's median of
    0.664, while the eight chosen on PSH for epoch 26238 topped out at 0.375:
    the good candidates were discarded before anything looked at them. So
    quality is predicted here and the risk of refusal is priced separately,
    and the two are multiplied (nova_lab.psh_calibration.rank) rather than one
    being used as a gate on the other.

    Predicting design_iiptm costs milliseconds, so the whole valid pool is
    ranked rather than a sample of it. Trained on the field and validated by
    ranking our own candidates: rho = 0.569 (n = 124), picking 8 that read a
    real median of 0.572 against 0.332 at random. See
    scripts.train_iiptm_surrogate for the numbers and their caveats.

    Returns None when the model is absent, so the caller can fall back rather
    than crash in the middle of an epoch.
    """
    import joblib
    from nova_lab.nanobodies import regions_batch
    from nova_lab.nanobody_features import featurize
    path = DATA_DIR / "dataset" / f"{target}_iiptm_surrogate.joblib"
    if not path.exists():
        log(f"no design_iiptm surrogate at {path}; falling back to the PSH order. "
            f"Train it with scripts.train_iiptm_surrogate -- selecting on PSH "
            f"measured a median of 0.335 against 0.572 for the surrogate.")
        return None
    bundle = joblib.load(path)
    regions = regions_batch(sequences)
    usable, rows = [], []
    for sequence, region in zip(sequences, regions):
        if not region:
            continue
        vector = featurize(sequence, region)
        if vector is None:
            continue
        usable.append(sequence)
        rows.append(vector)
    if not rows:
        log("no candidate could be featurized; falling back to the PSH order")
        return None
    import numpy as np
    predicted = bundle["model"].predict(np.vstack(rows))
    order = np.argsort(-predicted)
    log(f"surrogate ranked {len(usable)} candidates on predicted "
        f"{bundle['metric']}: best {predicted[order[0]]:.3f}, "
        f"median {float(np.median(predicted)):.3f} "
        f"(trained on {bundle['n_train']} scored sequences)")
    return {usable[i]: float(predicted[i]) for i in order}


# A first pass is a single draw, so a candidate is only dropped before
# refolding when even a favourable draw leaves it hopeless: one measurement
# deviates by PSH_SIGMA, and this is three of them above the threshold.
REFOLD_IF_UNDER = PSH_THRESHOLD + 3 * PSH_SIGMA


def select_by_expected_value(sequences: list[str], quality: dict[str, float] | None,
                             psh_pool: int, work: Path, replicates: int = 3
                             ) -> tuple[list[str], dict[str, float], int]:
    """The miner's choice of what to submit: quality times P(not refused).

    Predicting design_iiptm costs milliseconds against 3.5 s to fold one
    candidate, so quality decides who gets measured; PSH then decides the order
    among those. Neither half is a gate -- a ceiling discards a strong candidate
    that is one draw over the line, and ignoring PSH is what put eleven
    submissions in the bin.

    This lives in one function because `scripts.observe` has to replay the same
    choice. When the replay selects differently from the miner, the win rate it
    reports is some other miner's, and that is precisely how the 14-29 % figure
    came to be measured on candidates chosen by a gate that was wrong.
    """
    if quality:
        pool = list(quality)[:psh_pool]
    else:
        # No model: every candidate is equally promising as far as we know, so
        # the ordering below is the risk of refusal alone.
        pool = sequences[:psh_pool]
        quality = {s: 1.0 for s in pool}
    psh, folds = psh_screen(pool, work, replicates)
    ordered = psh_calibration.rank({s: quality.get(s, 1.0) for s in pool},
                                   psh, folds)
    return ordered, psh, folds


def psh_screen(sequences: list[str], work: Path,
               replicates: int = 3) -> tuple[dict[str, float], int]:
    """Our PSH mean per candidate, and how many folds it averages.

    Two stages, because folding costs 3.5 s a candidate and most of the pool
    does not deserve three: one pass over everything, then k-1 more on what is
    still in reach of the threshold. Nothing is refused here. The reading is a
    cost, priced by nova_lab.psh_calibration, and the caller ranks on it --
    a ceiling was what this used to be, and the ceiling it inherited (96) came
    from the install whose DSSP step was failing.
    """
    from scripts.tnp_developability import TNP
    if not TNP.exists():
        log(f"TNP not installed at {TNP}; cannot screen PSH at all. The "
            f"validator refuses above {PSH_STRICT} and that rejected every "
            f"production submission, so this is not a screen to skip. "
            f"Run scripts.doctor.")
        return {}, 0
    first = measure_psh(sequences, work / "pass1")
    if not first:
        log("PSH pass 1 returned nothing; the toolchain is not measuring "
            "(scripts.doctor)")
        return {}, 0
    worth_refolding = [s for s in sequences if first.get(s, 999.0) < REFOLD_IF_UNDER]
    log(f"PSH pass 1 on {len(first)}: best {min(first.values()):.1f}, median "
        f"{sorted(first.values())[len(first) // 2]:.1f}, "
        f"{sum(v <= PSH_THRESHOLD for v in first.values())} at or under "
        f"{PSH_THRESHOLD:.0f}; refolding {len(worth_refolding)}")
    samples = {s: [v] for s, v in first.items()}
    for i in range(2, replicates + 1):
        for s, v in measure_psh(worth_refolding, work / f"pass{i}").items():
            samples.setdefault(s, []).append(v)
    mean = {s: sum(v) / len(v) for s, v in samples.items()}
    folds = min(len(samples[s]) for s in worth_refolding) if worth_refolding else 1
    best = min(mean.values())
    log(f"PSH over {folds} fold(s): best {best:.1f} -> "
        f"{refusal_probability(best, folds):.0%} chance the validator refuses it "
        f"(the accepted field reads a median of {psh_calibration.ACCEPTED_MEDIAN})")
    return mean, folds


def archive_screen(sequences: list[str], target: str, allowed_v: str,
                   cdr3_len: int) -> list[str]:
    """The archive-duplicate rule, run on the few candidates PSH left.

    5.6 s per candidate, measured: a difflib alignment against each of the
    14424 archived sequences. On the handful that reach here that is under a
    minute; on the pool propose used to screen it was 28 of the epoch's 67.
    """
    from nova_lab.nanobodies import NanobodyArchive, check_sequence
    archive = NanobodyArchive(target)
    kept = []
    for seq in sequences:
        problems = check_sequence(seq, allowed_v, cdr3_len, archive)
        if problems:
            log(f"  dropped: {problems[0]}")
        else:
            kept.append(seq)
    log(f"archive rule: {len(kept)} of {len(sequences)} are new")
    return kept


def developability_screen(sequences: list[str], work: Path) -> dict[str, bool]:
    """True per sequence if both halves of the validator's gate pass.

    TNP builds the structure; Aggrescan3D then scores that same structure, so
    the two run in one pass. Measured on 84 candidates: TNP rejects 2.4% and
    A3D 39%, against 0% and 4% of the real submissions the validators accepted.
    """
    from scripts.aggrescan_gate import screen as a3d_screen
    from scripts.tnp_developability import TNP, parse, run
    if not TNP.exists():
        log(f"TNP not installed at {TNP}; skipping the developability screen. "
            f"Measured, it and Aggrescan3D reject 39% of our candidates, so "
            f"submitting unscreened throws away roughly two epochs in five.")
        return {}
    work.mkdir(parents=True, exist_ok=True)
    fasta = work / "candidates.fasta"
    names = {f"cand{i}": s for i, s in enumerate(sequences)}
    fasta.write_text("".join(f">{n}\n{s}\n" for n, s in names.items()))
    df = parse(run(fasta, min(4, len(sequences)), work), work)
    if "dssp" in df.columns:
        df = df[df["dssp"]]
    # The validator does not use TNP's own "(RED flag)" labels: it takes the
    # values and compares them against its own config, whose numbers are TNP's
    # rounded to two decimals in both directions. And this is a single fold, so
    # a reading next to a bound is noise -- on the calibration set that refused
    # an accepted submission whose PPC folds read 1.2858, 1.1078, 1.1125.
    # nova_lab.developability holds the bounds and the bands.
    from nova_lab import developability
    from scripts.tnp_developability import METRICS as TNP_METRICS
    tnp, flagged = {}, {}
    for row in df.to_dict("records"):
        sequence = names.get(row["name"])
        if sequence is None:
            continue
        values = {m: row[m] for m in TNP_METRICS if m in row}
        breaches = developability.breaches(values, folds=1)
        tnp[sequence] = not breaches
        if breaches:
            flagged[sequence] = "; ".join(breaches)
    for sequence, why in flagged.items():
        log(f"  TNP refuses {sequence[:28]}...: {why}")

    # A3D needs the structure TNP just built, so this pass yields one more PSH
    # fold for free on exactly the candidates that still matter. It is reported,
    # not enforced: psh_screen already priced the risk over three folds, and
    # refusing here on a single draw is how a candidate one draw over the line
    # gets thrown away. What it is good for is the journal -- this is the last
    # thing measured before the submission leaves.
    psh = dict(zip([names[n] for n in df["name"] if n in names],
                   [v for n, v in zip(df["name"], df["PSH"]) if n in names]))
    over = [s for s, v in psh.items() if v > PSH_STRICT]
    if over:
        log(f"{len(over)} of {len(psh)} candidates read above the validator's "
            f"threshold of {PSH_STRICT} on this fold (median "
            f"{sorted(psh.values())[len(psh) // 2]:.0f}); the accepted field "
            f"reads {psh_calibration.ACCEPTED_MEDIAN}")

    a3d = a3d_screen(names, work / "TNP_output" / "Final_Models", work / "a3d")
    if not a3d:
        log("Aggrescan3D could not run; only the TNP half of the gate applied. "
            "Build it with install/install-aggrescan3d.sh in the nanobody-filter repo.")
    log(f"developability: TNP rejected {sum(1 for ok in tnp.values() if not ok)}, "
        f"A3D rejected {sum(1 for ok in a3d.values() if not ok)}")
    verdict = {s: tnp.get(s, True) and a3d.get(s, True) for s in sequences}
    # Drop over-limit candidates only while some stay: never submit nothing over
    # a gate we cannot yet reproduce.
    survivors = [s for s in sequences if verdict.get(s) and psh.get(s, 0) <= PSH_STRICT]
    if survivors:
        verdict = {s: verdict.get(s, True) and psh.get(s, 0) <= PSH_STRICT for s in sequences}
    return verdict


def boltzgen_rank(sequences: list[str], run_name: str, nova_dir: Path,
                  work: Path) -> pd.DataFrame | None:
    """Rank with the validators' own BoltzGen. None when it cannot run.

    The scores are read back from a path derived from `run_name`, and BoltzGen
    reuses an existing output directory instead of recomputing -- which is why
    the README insists on a fresh RUN per batch. So a leftover scores.csv from
    an earlier attempt satisfies `exists()` and would be ranked as if it were
    this epoch's, submitting a sequence chosen from someone else's numbers. Two
    cheap checks make that impossible: the file has to be newer than the batch
    we just wrote, and it has to cover the sequences we asked about.
    """
    work.mkdir(parents=True, exist_ok=True)
    batch = work / f"{run_name}.csv"
    pd.DataFrame({"sequence": sequences}).to_csv(batch, index=False)
    written = batch.stat().st_mtime
    cmd = [sys.executable, "-m", "scripts.boltzgen_score", "--file", str(batch),
           "--run", run_name, "--nova-dir", str(nova_dir)]
    log(f"BoltzGen on {len(sequences)} candidates...")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    scores = DATA_DIR / "boltzgen" / run_name / "scores.csv"
    if proc.returncode != 0 or not scores.exists():
        log(f"BoltzGen failed (exit {proc.returncode}); last output:")
        for line in (proc.stdout + proc.stderr).strip().splitlines()[-5:]:
            log(f"  {line}")
        return None
    if scores.stat().st_mtime < written:
        log(f"{scores} predates this batch: BoltzGen reused an old output "
            f"directory instead of recomputing. Use a fresh --run. Not ranking "
            f"on another batch's numbers.")
        return None

    from scripts.boltzgen_score import rank_sum
    df = pd.read_csv(scores)
    if "sequence" not in df.columns:
        log(f"{scores} has no sequence column; cannot tell what it scored")
        return None
    scored = set(df["sequence"])
    missing = [s for s in sequences if s not in scored]
    if missing:
        log(f"BoltzGen scored {len(scored & set(sequences))} of {len(sequences)} "
            f"candidates; {len(missing)} came back with no metrics. Ranking the "
            f"rest, which is honest, but a sequence with no score is not a "
            f"sequence that scored badly.")
        df = df[df["sequence"].isin(sequences)]
    if df.empty:
        log("none of this batch's candidates were scored; abandoning the epoch")
        return None
    df = df.copy()
    df["rank_sum"] = rank_sum(df)
    return df.sort_values("rank_sum")


def donor_pool_lag(target: str, epoch: int) -> int | None:
    """How many epochs behind the generator's donor pool is, or None if unknown.

    The pool is the training dataset, not the live archive: it only refreshes
    when scripts.build_nanobody_dataset runs. A stale pool is invisible at
    every other checkpoint, which is how production ran a week behind.
    """
    import pandas as pd
    path = DATA_DIR / "dataset" / f"{target}_nanobodies.parquet"
    if not path.exists():
        return None
    dated = pd.read_parquet(path, columns=["epoch"])["epoch"].dropna()
    return None if dated.empty else int(epoch - dated.max())


def mine_epoch(st, args, params: dict) -> dict:
    """One epoch, start to finish. Returns the journal entry."""
    from scripts.propose_nanobody import propose

    epoch = params["epoch"]
    end_block = (epoch + 1) * EPOCH_LENGTH
    deadline = end_block - NO_SUBMISSION_BLOCKS - SAFETY_BLOCKS
    block = current_block(st)
    left = (deadline - block) * BLOCK_SECONDS
    t0 = time.monotonic()
    timings = {}
    entry = {"epoch": epoch, "started_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
             "allowed_v": params["allowed_v"], "allowed_cdr3_length": params["allowed_cdr3_length"],
             "live": bool(args.live), "submitted": False}

    def _stamp(entry):
        """Always persist timing before returning the journal entry."""
        entry["elapsed_seconds"] = round(time.monotonic() - t0, 1)
        entry["timings"] = timings
        remaining = left - (time.monotonic() - t0)
        log(f"total {entry['elapsed_seconds']:.0f}s, {remaining / 60:.0f} min remaining in budget")
        return entry

    log(f"epoch {epoch}: V {params['allowed_v']}, CDR3 {params['allowed_cdr3_length']}; "
        f"block {block}, deadline {deadline} ({left / 60:.0f} min)")
    if left <= 0:
        entry["status"] = "too_late"
        log("past the submission deadline for this epoch")
        return _stamp(entry)

    if not config_unchanged():
        entry["status"] = "config_drift"
        return _stamp(entry)

    stale = donor_pool_lag(args.target, epoch)
    entry["donor_lag_epochs"] = stale
    if stale is not None and stale > args.max_donor_lag:
        # Measured in production: at epoch 26174 the pool stopped at 26025, and
        # none of the 55 scored submissions of 26173-26174 were in it. The
        # generator then grafts from a generation the competition has moved past.
        log(f"donor pool is {stale} epochs behind ({stale * 72 / 60 / 24:.1f} days); "
            f"refusing to submit. Rebuild it in three steps, in this order:\n"
            f"    python -m scripts.fetch_submissions --last 60   # decrypt what was submitted\n"
            f"    python -m scripts.analyze_nanobodies           # join the published metrics\n"
            f"    python -m scripts.build_nanobody_dataset       # re-download the archive, refeaturize\n"
            f"  The last step alone is not enough: without the first two the new "
            f"sequences carry no epoch and the lag does not move.")
        entry["status"] = "donor_pool_stale"
        return _stamp(entry)
    if stale:
        log(f"donor pool is {stale} epochs behind ({stale * 72 / 60 / 24:.1f} days)")

    cdr3_len = params["allowed_cdr3_length"]
    # Size the PSH pool to the length's measured pass rate rather than to a
    # fixed number. Long CDR3 is the hard case -- 15% under 89 at length 18
    # against 75% at length 8 -- and it is also where we rank best, six wins in
    # nine epochs replayed, so it gets the wider pool.
    psh_pool, rate = psh_pool_for(cdr3_len, args.psh_pool, args.candidates)
    entry["psh_pass_rate"] = round(rate, 3)
    entry["psh_pool"] = psh_pool
    log(f"CDR3 length {cdr3_len}: {rate:.0%} of candidates expected at or under "
        f"{PSH_THRESHOLD:.0f}, so PSH pool {psh_pool} to find about "
        f"{WANT_UNDER_THRESHOLD}")
    cost = psh_pool * PSH_SECONDS_PER_CANDIDATE / 60
    if cost > 30:
        log(f"  that pool costs {cost:.0f} min of CPU; the epoch budget is "
            f"about 66 min, so the run may not finish")

    with timed("propose", timings):
        result = propose(params, args.target, candidates=args.candidates,
                         n_donors=args.donors, n_cdr3=args.cdr3_pool,
                         shortlist=max(psh_pool, args.shortlist),
                         seed=args.seed + epoch,
                         generator_name=args.generator, skip_surrogate=True,
                         skip_archive=True, log=log)
    shortlist = result["shortlist"]
    entry.update({"generated": result["generated"], "valid": result["valid"],
                  "status": result["status"]})
    if shortlist.empty:
        log(f"nothing to submit ({result['status']})")
        return _stamp(entry)

    work = DATA_DIR / "mining" / f"epoch{epoch}"
    sequences = list(shortlist["sequence"])

    with timed("surrogate", timings):
        quality = surrogate_rank(sequences, args.target)
    entry["surrogate_ranked"] = bool(quality)
    with timed("psh_screen", timings):
        ordered, psh, folds = select_by_expected_value(
            sequences, quality, psh_pool, work / "psh", args.psh_replicates)
    sequences = ordered[:args.shortlist]
    if psh:
        best_sequence = ordered[0]
        best = psh.get(best_sequence, float("nan"))
        risk = refusal_probability(best, folds)
        entry["psh_folds"] = folds
        entry["psh_best"] = None if best != best else round(best, 2)
        entry["psh_refusal_estimate"] = None if risk != risk else round(risk, 3)
        measured = [psh[s] for s in sequences if s in psh]
        entry["psh_selected_max"] = round(max(measured), 2) if measured else None
        entry["psh_unmeasured_selected"] = sum(1 for s in sequences if s not in psh)
        log(f"chosen on expected value: {psh_calibration.verdict(best, folds)}")
        if risk == risk and risk > 0.5:
            # Still submitted: a refused submission and no submission both score
            # zero, so any chance of passing beats abstaining. It is logged as
            # the generator's problem, which is what it is -- the accepted field
            # reads a median of 83.5 and our candidates 88.7.
            log(f"  no candidate is safe: the best still has a {risk:.0%} chance "
                f"of refusal. Submitting anyway, since not submitting is 0%, but "
                f"this epoch is a coin toss and the fix is a less hydrophobic "
                f"generator, not a wider pool.")
        if entry["psh_unmeasured_selected"]:
            log(f"  {entry['psh_unmeasured_selected']} of the {len(sequences)} "
                f"selected were never measured; they rank last by construction")

    # The archive rule, now that the pool is a handful rather than hundreds.
    # propose skipped it: at 5.6 s per candidate it costs 28 min on a pool of
    # 300, and PSH discards 97% of them. A duplicate is a certain rejection, so
    # this one cannot be downgraded to a warning -- unlike the developability
    # screen below, there is no "submit the best anyway".
    with timed("archive", timings):
        kept = archive_screen(sequences, args.target, params["allowed_v"], cdr3_len)
    entry["archive_rejected"] = len(sequences) - len(kept)
    if not kept:
        log("every surviving candidate duplicates the archive; nothing to submit")
        entry["status"] = "all_duplicates"
        return _stamp(entry)
    sequences = kept

    with timed("developability", timings):
        verdict = developability_screen(sequences, work / "tnp")
    if verdict:
        kept = [s for s in sequences if verdict.get(s, True)]
        entry["developability_rejected"] = len(sequences) - len(kept)
        log(f"{len(kept)} of {len(sequences)} clear the developability gate")
        sequences = kept or sequences  # never submit nothing over a screen
        if not kept:
            log("every candidate is red; submitting the best anyway rather than skipping")

    with timed("boltzgen", timings):
        ranked = boltzgen_rank(sequences, f"epoch{epoch}", args.nova_dir, work / "boltzgen")
    if ranked is None or ranked.empty:
        log("BoltzGen unavailable or returned nothing. The surrogate wins 14% of "
            "epochs vs 48% for BoltzGen vs 23% random — falling back is worse than "
            "skipping. Abandoning this epoch.")
        entry["status"] = "no_boltzgen"
        return _stamp(entry)
    best = ranked["sequence"].iat[0]
    entry["ranked_by"] = "boltzgen"
    # Persist scores for post-mortem analysis (distinguish "rejected by gate" from "scored low").
    score_cols = [c for c in ranked.columns if c not in ("sequence",)]
    entry["boltzgen_scores"] = ranked[["sequence"] + score_cols].head(20).to_dict("records")
    entry["boltzgen_best_rank_sum"] = float(ranked["rank_sum"].iat[0])
    entry["sequence"] = best
    entry["cdr3"] = shortlist.set_index("sequence")["cdr3"].get(best, "")

    uid = st.get_uid_for_hotkey_on_subnet(args.hotkey_ss58, args.netuid) if args.hotkey_ss58 else args.uid
    if uid is None:
        entry["status"] = "no_uid"
        log("this hotkey is not registered on the subnet")
        return _stamp(entry)
    block = retry(current_block, st)
    reveal = timelock.reveal_round(block, EPOCH_LENGTH, NO_SUBMISSION_BLOCKS)
    filename, content = timelock.submission_file(uid, reveal, f"~|{best}")
    pointer = f"{github_path()}/{filename}"
    entry.update({"uid": uid, "reveal_round": reveal, "pointer": pointer})
    (work / filename).parent.mkdir(parents=True, exist_ok=True)
    (work / filename).write_text(content)
    log(f"UID {uid}, reveal round {reveal}, pointer {pointer}")
    log(f"submission file written to {work / filename}")

    if not args.live:
        entry["status"] = "dry_run"
        log("dry run: nothing committed on chain, nothing pushed to GitHub. Pass --live.")
        return _stamp(entry)

    # Upload FIRST, then commit. A file without a pointer is harmless (orphaned
    # on GitHub); a pointer without a file is fatal (dangling, epoch lost).
    # Epoch 26172 committed on-chain then GitHub returned 404 → dangling pointer.
    log("uploading the encrypted file to GitHub...")
    try:
        upload(filename, content)
    except Exception as error:
        entry.update({"status": "upload_failed", "error": f"{type(error).__name__}: {error}"})
        log(f"upload failed: {error}")
        log("not committing, since the file is not reachable")
        return _stamp(entry)

    log("committing the pointer on chain...")
    try:
        resp = st.set_commitment(wallet=args.wallet_obj, netuid=args.netuid, data=pointer)
    except Exception as error:
        # The chain rate-limits commitments, so committing again too soon raises
        # rather than returning. A miner that must run for weeks cannot die here.
        entry.update({"status": "commit_failed", "error": f"{type(error).__name__}: {error}"})
        log(f"commitment failed: {type(error).__name__}: {error}")
        log("the file is uploaded but nothing points at it (harmless orphan)")
        return _stamp(entry)
    entry["commitment"] = str(resp)
    log(f"commitment: {resp}")
    entry.update({"status": "submitted", "submitted": True})
    log("submitted")
    return _stamp(entry)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--live", action="store_true",
                        help="actually commit on chain and push to GitHub")
    parser.add_argument("--once", action="store_true", help="one epoch, then stop")
    parser.add_argument("--wallet", help="bittensor wallet name")
    parser.add_argument("--hotkey", default="default", help="bittensor hotkey name")
    parser.add_argument("--uid", type=int, help="skip the metagraph lookup and use this UID")
    parser.add_argument("--netuid", type=int, default=68)
    parser.add_argument("--target", default=NANOBODY_TARGET[0])
    parser.add_argument("--nova-dir", type=Path, default=Path("/home/user/nova-upstream"))
    parser.add_argument("--candidates", type=int, default=5000)
    parser.add_argument("--donors", type=int, default=40)
    parser.add_argument("--cdr3-pool", type=int, default=60)
    parser.add_argument("--shortlist", type=int, default=8,
                        help="candidates sent to BoltzGen after the PSH cut; "
                             "~1.6 min each on 2x T4")
    parser.add_argument("--psh-pool", type=int,
                        help="override the measured pool sizing (candidates to "
                             "profile with TNP before selecting, 3.5 s each on 4 "
                             "cores). Left out, the pool is sized from the "
                             "length's measured rate under the validator's "
                             "threshold: ~13 at CDR3 8, ~67 at 18. It used to "
                             "default to 500 AS A FLOOR, which silently "
                             "overrode that sizing")
    parser.add_argument("--psh-replicates", type=int, default=3,
                        help="folds averaged for the PSH of the best candidates; one "
                             "measurement deviates by 1.74, k folds by 1.74/sqrt(k)")
    # --psh-margin is gone: it required 4.6 points of headroom under the
    # validator's 89, so 84.4, on the assumption that our reading and theirs
    # are the same measurement. They are not. Sequences the validators accepted
    # read a median of 92.0 at CDR3 17+ in our pipeline, and the acceptance
    # probability only collapses past 96. The veto is now the measured curve in
    # nova_lab.psh_calibration.
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max-donor-lag", type=int, default=40,
                        help="refuse to submit if the donor pool is older than this "
                             "many epochs (40 ~= 2 days)")
    parser.add_argument("--generator", choices=sorted(GENERATORS), default="pssm",
                        help="CDR3 generation strategy (default: pssm)")
    args = parser.parse_args()

    import bittensor as bt
    from scripts.check_nanobody import resolve_epoch

    args.wallet_obj = None
    args.hotkey_ss58 = None
    if args.wallet:
        args.wallet_obj = bt.wallet(name=args.wallet, hotkey=args.hotkey)
        args.hotkey_ss58 = args.wallet_obj.hotkey.ss58_address
    elif args.live:
        sys.exit("--live needs --wallet")
    if args.uid is None and not args.hotkey_ss58:
        sys.exit("give --wallet (to look the UID up) or --uid")

    st = retry(bt.Subtensor, "finney")
    # An upgrade inside the container (11.3.0 was found installed against a
    # 10.4.1 pin) removes every call the submission path makes. Say so here,
    # not forty minutes into an epoch with the CPU already spent.
    if args.live or args.hotkey_ss58:
        require_api(st)
    elif missing_api(st):
        log(f"warning: this bittensor cannot submit (no {', '.join(missing_api(st))}); "
            f"dry runs only, see nova_lab.chain")
    JOURNAL.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if JOURNAL.exists():
        done = {json.loads(l)["epoch"] for l in JOURNAL.open() if l.strip()}

    while True:
        params = retry(resolve_epoch, None)
        if params["epoch"] in done:
            if args.once:
                # Don't sit through the polling sleep to say nothing: --once on
                # an epoch already in the journal is a no-op, and silence here
                # reads as a run in progress.
                log(f"epoch {params['epoch']} is already in the journal, nothing to do")
                return
            time.sleep(60)
            continue
        try:
            entry = mine_epoch(st, args, params)
        except Exception as error:
            # Journal the failure and carry on: one bad epoch must not end the run.
            entry = {"epoch": params["epoch"], "status": "crashed", "submitted": False,
                     "error": f"{type(error).__name__}: {error}",
                     "started_at": datetime.now(timezone.utc).isoformat(timespec="seconds")}
            log(f"epoch {params['epoch']} crashed: {type(error).__name__}: {error}")
        with JOURNAL.open("a") as f:
            f.write(json.dumps(entry) + "\n")
        done.add(params["epoch"])
        if args.once:
            return


if __name__ == "__main__":
    # bittensor leaves a non-daemon websocket thread behind, so returning from
    # main() does not end the process. Measured: --once journalled its entry in
    # 0.8 s and the process was still alive two minutes later with nothing left
    # to do, so a supervisor waiting on the exit (systemd Restart=, cron, a
    # shell loop) would wait on a miner that had already finished. Exit hard,
    # after flushing what we wrote ourselves.
    code = 0
    try:
        main()
    except SystemExit as exit_request:  # sys.exit() inside main(), message or code
        if isinstance(exit_request.code, str):
            print(exit_request.code, file=sys.stderr)
            code = 1
        elif isinstance(exit_request.code, int):
            code = exit_request.code
    finally:
        sys.stdout.flush()
        sys.stderr.flush()
        os._exit(code)
