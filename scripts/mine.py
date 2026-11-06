"""The NOVA nanobody miner: generate, screen, score, encrypt, publish, commit.

The loop per epoch, in order:

  1. resolve the epoch's imposed V gene and CDR3 length from the chain;
  2. refuse to submit if the subnet's own config has moved under us
     (scripts.check_config) — a changed target makes every candidate worthless;
  3. generate and keep what the offline checker accepts
     (scripts.propose_nanobody);
  4. screen the shortlist through both halves of the validator's developability
     stage: TNP (measured: rejects ~2.4% of our output) and Aggrescan3D, which
     requires a3d_average <= -0.680 and a3d_maximum <= 1.580 and rejects ~39%
     of it, almost always on the peak. Skipping this screen is what turns a
     generated winner into a silent rejection;
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
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

from nova_lab import DATA_DIR
from nova_lab.config import EPOCH_LENGTH, NANOBODY_TARGET, NO_SUBMISSION_BLOCKS
from nova_lab import timelock

JOURNAL = DATA_DIR / "mining" / "journal.jsonl"
BLOCK_SECONDS = 12
# Leave room for encryption, the GitHub PUT and the extrinsic. Measured: BoltzGen
# runs at ~1.6 min/sequence on 2x T4, so the candidate budget follows from this.
SAFETY_BLOCKS = 12


def log(msg: str) -> None:
    print(f"[{datetime.now(timezone.utc):%H:%M:%S}] {msg}", flush=True)


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

    a3d = a3d_screen(names, work / "TNP_output" / "Final_Models", work / "a3d")
    if not a3d:
        log("Aggrescan3D could not run; only the TNP half of the gate applied. "
            "Build it with install/install-aggrescan3d.sh in the nanobody-filter repo.")
    log(f"developability: TNP rejected {sum(1 for ok in tnp.values() if not ok)}, "
        f"A3D rejected {sum(1 for ok in a3d.values() if not ok)}")
    return {s: tnp.get(s, True) and a3d.get(s, True) for s in sequences}


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


def mine_epoch(st, args, params: dict) -> dict:
    """One epoch, start to finish. Returns the journal entry."""
    from scripts.propose_nanobody import propose

    epoch = params["epoch"]
    end_block = (epoch + 1) * EPOCH_LENGTH
    deadline = end_block - NO_SUBMISSION_BLOCKS - SAFETY_BLOCKS
    block = st.get_current_block()
    left = (deadline - block) * BLOCK_SECONDS
    entry = {"epoch": epoch, "started_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
             "allowed_v": params["allowed_v"], "allowed_cdr3_length": params["allowed_cdr3_length"],
             "live": bool(args.live), "submitted": False}
    log(f"epoch {epoch}: V {params['allowed_v']}, CDR3 {params['allowed_cdr3_length']}; "
        f"block {block}, deadline {deadline} ({left / 60:.0f} min)")
    if left <= 0:
        entry["status"] = "too_late"
        log("past the submission deadline for this epoch")
        return entry

    if not config_unchanged():
        entry["status"] = "config_drift"
        return entry

    result = propose(params, args.target, candidates=args.candidates,
                     n_donors=args.donors, n_cdr3=args.cdr3_pool,
                     shortlist=args.shortlist, seed=args.seed + epoch, log=log)
    shortlist = result["shortlist"]
    entry.update({"generated": result["generated"], "valid": result["valid"],
                  "status": result["status"]})
    if shortlist.empty:
        log(f"nothing to submit ({result['status']})")
        return entry

    work = DATA_DIR / "mining" / f"epoch{epoch}"
    sequences = list(shortlist["sequence"])
    verdict = developability_screen(sequences, work / "tnp")
    if verdict:
        kept = [s for s in sequences if verdict.get(s, True)]
        entry["developability_rejected"] = len(sequences) - len(kept)
        log(f"{len(kept)} of {len(sequences)} clear the developability gate")
        sequences = kept or sequences  # never submit nothing over a screen
        if not kept:
            log("every candidate is red; submitting the best anyway rather than skipping")

    ranked = boltzgen_rank(sequences, f"epoch{epoch}", args.nova_dir, work / "boltzgen") \
        if args.scorer == "boltzgen" else None
    if ranked is not None and not ranked.empty:
        best = ranked["sequence"].iat[0]
        entry["ranked_by"] = "boltzgen"
    else:
        best = sequences[0]
        entry["ranked_by"] = "surrogate"
        log("WARNING: ranking by the surrogate. Measured over 21 epochs it wins 14% "
            "of them against 48% for BoltzGen, and 23% for a random pick — it is "
            "worse than chance as a selector.")
    entry["sequence"] = best
    entry["cdr3"] = shortlist.set_index("sequence")["cdr3"].get(best, "")

    uid = st.get_uid_for_hotkey_on_subnet(args.hotkey_ss58, args.netuid) if args.hotkey_ss58 else args.uid
    if uid is None:
        entry["status"] = "no_uid"
        log("this hotkey is not registered on the subnet")
        return entry
    block = st.get_current_block()
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
        return entry

    # The reference miner commits first and uploads after, so a failed upload
    # leaves a dangling pointer rather than an unreferenced file.
    log("committing the pointer on chain...")
    try:
        resp = st.set_commitment(wallet=args.wallet_obj, netuid=args.netuid, data=pointer)
    except Exception as error:
        # The chain rate-limits commitments, so committing again too soon raises
        # rather than returning. A miner that must run for weeks cannot die here.
        entry.update({"status": "commit_failed", "error": f"{type(error).__name__}: {error}"})
        log(f"commitment failed: {type(error).__name__}: {error}")
        log("not uploading, since nothing points at the file")
        return entry
    entry["commitment"] = str(resp)
    log(f"commitment: {resp}")

    log("uploading the encrypted file to GitHub...")
    try:
        upload(filename, content)
    except Exception as error:
        entry.update({"status": "upload_failed", "error": f"{type(error).__name__}: {error}"})
        log(f"upload failed: {error}")
        log("the pointer is committed but the file is missing: this epoch is lost, "
            "and the validator will read a dangling pointer")
        return entry
    entry.update({"status": "submitted", "submitted": True})
    log("submitted")
    return entry


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
    parser.add_argument("--scorer", choices=["boltzgen", "surrogate"], default="boltzgen")
    parser.add_argument("--nova-dir", type=Path, default=Path("/home/user/nova-upstream"))
    parser.add_argument("--candidates", type=int, default=3000)
    parser.add_argument("--donors", type=int, default=40)
    parser.add_argument("--cdr3-pool", type=int, default=60)
    parser.add_argument("--shortlist", type=int, default=8,
                        help="candidates sent to BoltzGen; ~1.6 min each on 2x T4")
    parser.add_argument("--seed", type=int, default=0)
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

    st = bt.Subtensor("finney")
    JOURNAL.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if JOURNAL.exists():
        done = {json.loads(l)["epoch"] for l in JOURNAL.open() if l.strip()}

    while True:
        params = resolve_epoch(None)
        if params["epoch"] in done:
            time.sleep(60)
            if args.once:
                return
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
    main()
