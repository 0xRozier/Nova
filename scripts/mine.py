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

JOURNAL = DATA_DIR / "mining" / "journal.jsonl"
BLOCK_SECONDS = 12
# Leave room for encryption, the GitHub PUT and the extrinsic. Measured: BoltzGen
# runs at ~1.6 min/sequence on 2x T4, so the candidate budget follows from this.
SAFETY_BLOCKS = 12
# Share of candidates clearing PSH 89, measured per CDR3 length on 540
# candidates (results/psh_by_generator.csv). These are per-CANDIDATE rates, not
# per-epoch: 0% of the depth-4 epochs at length >= 17 had a passing candidate,
# but with only four candidates each and a 9% rate, finding none is the expected
# outcome (0.91^4 = 69%). Long CDR3 is hard, not hopeless — and it is where we
# rank best, 6 wins of 9 epochs — so the pool is widened instead of skipping.
PSH_PASS_RATE = {8: 0.51, 11: 0.28, 18: 0.09}
PSH_MARGIN_RATE = {8: 0.20, 11: 0.11, 18: 0.017}
WANT_WITH_MARGIN = 5  # candidates we want under the margin before BoltzGen


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


def psh_pool_for(cdr3_length: int, floor: int, ceiling: int) -> tuple[int, float]:
    """Pool size that should yield WANT_WITH_MARGIN candidates under the margin."""
    rate = psh_pass_rate(cdr3_length, PSH_MARGIN_RATE)
    need = int(WANT_WITH_MARGIN / max(rate, 0.002))
    return max(floor, min(need, ceiling)), rate


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


PSH_SIGMA = 4.01  # deviation of ONE PSH measurement, from 10 sequences measured twice


def measure_psh(sequences: list[str], work: Path) -> dict[str, float]:
    """One PSH measurement per sequence. A fresh work dir forces a new fold."""
    from scripts.tnp_developability import parse, run
    work.mkdir(parents=True, exist_ok=True)
    names = {f"cand{i}": s for i, s in enumerate(sequences)}
    (work / "candidates.fasta").write_text(
        "".join(f">{n}\n{s}\n" for n, s in names.items()))
    df = parse(run(work / "candidates.fasta", 4, work))
    return {names[n]: v for n, v in zip(df["name"], df["PSH"]) if n in names}


def psh_select(sequences: list[str], keep: int, work: Path,
               replicates: int = 3) -> tuple[list[str], dict]:
    """Return the `keep` least hydrophobic candidates, PSH averaged over folds.

    PSH is computed on a predicted structure and NanoBodyBuilder2 is not
    deterministic, so one measurement is noisy: ten sequences measured twice
    give a per-measurement deviation of 4.01 against a threshold of 89. Our
    local value therefore cannot be calibrated to the validator's — both are
    draws from the same distribution — and a single reading near the limit is a
    coin toss.

    Averaging k folds divides that deviation by root k, which is what buys back
    a usable margin. Two stages keep it affordable: one cheap pass over the
    whole pool (3.6 s per candidate on 4 cores), then k-1 more passes over the
    best few. For a 300-candidate pool keeping 40, that is about 23 minutes
    against a 66-minute usable epoch, where folding everything k times would
    cost 54.
    """
    from scripts.tnp_developability import TNP
    if not TNP.exists():
        log(f"TNP not installed at {TNP}; cannot measure PSH. The validator "
            f"refuses PSH > {PSH_STRICT} and that rejected every production "
            f"submission, so this is not a screen to skip.")
        return sequences[:keep], {}

    first = measure_psh(sequences, work / "pass1")
    if not first:
        return sequences[:keep], {}
    ranked = sorted(first, key=lambda s: first[s])
    log(f"PSH pass 1 on {len(first)} candidates: best {first[ranked[0]]:.1f}, "
        f"median {first[ranked[len(ranked) // 2]]:.1f}")

    # Re-fold only the promising head, with headroom for the ranking to change.
    head = ranked[:max(keep * 3, keep)]
    samples = {s: [first[s]] for s in head}
    for i in range(2, replicates + 1):
        more = measure_psh(head, work / f"pass{i}")
        for s, v in more.items():
            samples.setdefault(s, []).append(v)
    mean = {s: sum(v) / len(v) for s, v in samples.items()}
    ranked = sorted(mean, key=lambda s: mean[s])
    n = min(len(v) for v in samples.values()) or 1
    log(f"PSH over {replicates} folds on the best {len(head)}: "
        f"best {mean[ranked[0]]:.1f}, deviation of the mean "
        f"{PSH_SIGMA / n ** 0.5:.1f}; keeping the {keep} lowest")
    return ranked[:keep], mean


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
    df = parse(run(fasta, min(4, len(sequences)), work))
    tnp = {names[n]: bool(ok) for n, ok in zip(df["name"], ~df["rejected"]) if n in names}

    # The validator refuses PSH > 89 outright (config.yaml
    # psh_strict_rejection_threshold), far below TNP's own red 155.47. Our PSH is
    # not calibrated to theirs, so this ranks candidates rather than rejecting
    # them: prefer the least hydrophobic, which is the one scale-free signal we
    # have. Our candidates sit at median 98 against 90 for accepted submissions.
    psh = dict(zip([names[n] for n in df["name"] if n in names],
                   [v for n, v in zip(df["name"], df["PSH"]) if n in names]))
    over = [s for s, v in psh.items() if v > PSH_STRICT]
    if over:
        log(f"{len(over)} of {len(psh)} candidates exceed the validator's strict "
            f"PSH limit of {PSH_STRICT} (ours median "
            f"{sorted(psh.values())[len(psh) // 2]:.0f}); keeping the least hydrophobic")

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
    """Rank with the validators' own BoltzGen. None when it cannot run."""
    work.mkdir(parents=True, exist_ok=True)
    batch = work / f"{run_name}.csv"
    pd.DataFrame({"sequence": sequences}).to_csv(batch, index=False)
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
    from scripts.boltzgen_score import rank_sum
    df = pd.read_csv(scores)
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
    # Size the PSH pool to the length's measured pass rate rather than skipping:
    # at length 18 only 1.7% of candidates clear the margin, but 1.7% of 300 is
    # still five, and these are the epochs we rank best in.
    psh_pool, rate = psh_pool_for(cdr3_len, args.psh_pool, args.candidates)
    entry["psh_margin_rate"] = round(rate, 3)
    log(f"CDR3 length {cdr3_len}: {rate:.1%} of candidates expected under the "
        f"margin, so PSH pool {psh_pool} to find about {WANT_WITH_MARGIN}")
    if psh_pool * 3.6 / 60 > 30:
        log(f"  that pool costs {psh_pool * 3.6 / 60:.0f} min of CPU; the epoch "
            f"budget is about 66 min, so the run may not finish")

    with timed("propose", timings):
        result = propose(params, args.target, candidates=args.candidates,
                         n_donors=args.donors, n_cdr3=args.cdr3_pool,
                         shortlist=max(psh_pool, args.shortlist),
                         seed=args.seed + epoch,
                         generator_name=args.generator, skip_surrogate=True, log=log)
    shortlist = result["shortlist"]
    entry.update({"generated": result["generated"], "valid": result["valid"],
                  "status": result["status"]})
    if shortlist.empty:
        log(f"nothing to submit ({result['status']})")
        return _stamp(entry)

    work = DATA_DIR / "mining" / f"epoch{epoch}"
    sequences = list(shortlist["sequence"])

    # PSH first, on the whole pool: it is the gate that rejected every
    # production submission, and selecting on it needs a wide funnel rather
    # than a filter at the end of a narrow one.
    with timed("psh_select", timings):
        sequences, psh = psh_select(sequences, args.shortlist, work / "psh",
                                    args.psh_replicates)
    if psh:
        best = min(psh[s] for s in sequences)
        entry["psh_best"] = round(best, 2)
        entry["psh_selected_max"] = round(max(psh[s] for s in sequences), 2)
        if best > PSH_STRICT - args.psh_margin:
            log(f"best PSH is {best:.1f}; the validator refuses above "
                f"{PSH_STRICT} and our measurement's deviation is 4.0, so this "
                f"epoch has no candidate with margin. Skipping rather than "
                f"spending GPU on a submission that will be refused.")
            entry["status"] = "psh_no_margin"
            return _stamp(entry)

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
    parser.add_argument("--psh-pool", type=int, default=500,
                        help="candidates to measure PSH on before selecting; "
                             "3.6 s each on 4 cores, and only ~21% come in under "
                             "87.4, so a narrow pool yields nothing")
    parser.add_argument("--psh-replicates", type=int, default=3,
                        help="folds averaged for the PSH of the best candidates; "
                             "one measurement deviates by 4.0, k folds by 4.0/sqrt(k)")
    parser.add_argument("--psh-margin", type=float, default=4.6,
                        help="required headroom under the validator's PSH limit; "
                             "4.6 is two deviations of a 3-fold mean (8.0 for one fold)")
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
