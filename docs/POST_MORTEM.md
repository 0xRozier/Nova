# Post-mortem: SN68 Production Run (Epochs 26172-26185)

UID 22 | 2x RTX 4090 (Lium) | 2026-11-06 17:04 to 2026-11-07 08:10 UTC

## Result

12 epochs attempted, 11 submitted on-chain, **0 scored by validators**.
Zero incentive earned. $13.50 spent on GPU rental.

## Timeline

| Epoch | V gene | CDR3 len | Ranked by | PSH filter | Status | Root cause |
|-------|--------|----------|-----------|------------|--------|------------|
| 26172 | IGHV3-23 | 16 | surrogate | none | **upload_failed** | GitHub 404 |
| 26173 | IGHV3-7 | 15 | surrogate | none | submitted | PSH not filtered |
| 26174 | IGHV3-74 | 8 | BoltzGen | none | submitted | PSH not filtered |
| 26175 | IGHV3-74 | 11 | BoltzGen | none | submitted | PSH not filtered |
| 26176 | IGHV3-74 | 18 | BoltzGen | none | submitted | PSH not filtered |
| 26177 | IGHV3-74 | 18 | BoltzGen | **31/40 over** | submitted | PSH fix deployed, but miscalibrated |
| 26178 | IGHV3-74 | 6 | BoltzGen | 4/40 over | submitted | PSH fix deployed, but miscalibrated |
| 26179 | IGHV3-23 | 15 | BoltzGen | 26/30 over | submitted | PSH fix deployed, but miscalibrated |
| 26180 | IGHV3-23 | 13 | BoltzGen | 22/33 over | submitted | PSH fix deployed, but miscalibrated |
| 26181 | IGHV3-64 | 11 | BoltzGen | 25/40 over | submitted | PSH fix deployed, but miscalibrated |
| *26182-26184* | — | — | — | — | **skipped** | HTTP 429 crash, 5h downtime |
| 26185 | IGHV3-74 | 5 | BoltzGen | 1/40 over | submitted | PSH fix deployed, but miscalibrated |

## Cause of failure: PSH miscalibration

**Every single submission was rejected by the validator's PSH strict rejection
gate (threshold: PSH <= 89).** The nanobody scoring pipeline was active
throughout — other miners received scores in every epoch we participated in.

### Evidence

The validator's W&B logs (run `sjd835bx`, epoch 26209) show our UID being
explicitly rejected:

```
UID 22: contains sequences that are rejected by PSH strict rejection:
[('ENQLVESGGGVVKPGGSLTLSCAASGFTFSSYAMSWIRQAPGKERVWVG...', 89.5651)]
```

The score-sharing API confirms: our sequence hashes return 404 for every epoch
(26173-26185), while other miners' hashes return 200 with full BoltzGen metrics.

### The calibration problem

Our miner and the validator both call TNP (Therapeutic Nanobody Profiler) to
compute PSH. But the values differ between our local installation and the
validator's:

| What we observed | Our local PSH | Validator's PSH |
|------------------|---------------|-----------------|
| Our epoch 26209 submission | ~85 (passed our local filter) | 89.57 (rejected) |
| Accepted submissions median | 90.0 | <89 (by definition) |
| Our candidates median | 95-98 | unknown |

Our "PSH fix" (deployed at epoch 26177) filtered candidates to local PSH <= 89.
But because the validator's PSH computation produces **different values** than
ours, sequences that pass our local filter still fail the validator's gate.

We confirmed this asymmetry empirically: 62% of validator-accepted submissions
score above 89 on our local TNP. The hydrophobicity scale or TNP version is
different. Five Kyte-Doolittle variants were tested locally; none matched.

### What we should have done

Set our local PSH threshold much lower (e.g., 75-80) to create a safety margin
that accounts for the scale mismatch. The filter as implemented (local PSH <= 89)
was necessary but insufficient — it was operating on the wrong scale.

## Secondary issues

### 1. GitHub upload failure (epoch 26172)

The first epoch committed a pointer on-chain, then the GitHub upload failed with
HTTP 404. The submissions repo had just been initialized minutes earlier
(`cleanup` and `test` commits at 17:15). The pointer became a dangling reference.

**Impact**: 1 epoch lost.
**Fix needed**: Verify repo accessibility before committing the pointer.

### 2. HTTP 429 crash (epochs 26182-26184)

After epoch 26181 (03:05 UTC), the miner crashed when connecting to
`archive.chain.opentensor.ai`:

```
websockets.exceptions.InvalidStatus: server rejected WebSocket connection: HTTP 429
```

The process died and was only manually restarted at 07:51 UTC, losing ~5 hours
and 3 epochs.

**Impact**: 3 epochs lost.
**Fix needed**: Retry logic with exponential backoff for subtensor connections.

### 3. Surrogate selector is worse than random

Epochs 26172-26173 were ranked by the surrogate model (no BoltzGen available
yet). The surrogate's Spearman correlation against BoltzGen on our candidates
is **-0.20** — worse than random selection (+0.23 expected). The miner warns:
"it wins 14% of epochs against 48% for BoltzGen, and 23% for a random pick".

**Impact**: 2 epochs with near-random candidate selection.
**Fix needed**: Never submit surrogate-ranked candidates in production.

### 4. Candidate diversity too low

Even when BoltzGen was available, only 3-25 candidates reached BoltzGen scoring
per epoch (after validity, liability, developability, and PSH filters). For
comparison, competitive miners likely score hundreds.

Bottlenecks in the pipeline:

| Stage | Epoch 26179 example | Survival |
|-------|---------------------|----------|
| Generated | 5000 | 100% |
| Validity (germline, CDR3, etc.) | 2353 | 47% |
| Liability pre-filter (count <= 8) | 31 | 0.6% |
| Surrogate shortlist | 30 | 0.6% |
| PSH filter (local <= 89) | 4 | 0.08% |
| Developability (TNP + A3D) | 4 | 0.08% |
| BoltzGen scored | 4 | 0.08% |

The liability pre-filter eliminates 98.7% of valid candidates. Combined with
the PSH filter, this leaves almost nothing for BoltzGen to score.

### 5. BoltzGen scores not persisted

The miner stores sequences in the journal but not their BoltzGen scores. We
cannot retroactively evaluate whether our submitted candidates were competitive.
This makes it impossible to distinguish "rejected by validator gate" from
"scored but ranked last".

### 6. Weight override: all incentive goes to UID 42

`config.yaml` contains `payout.override_uid: 42`. All on-chain weight goes to
UID 42 by design (compound payout mechanism). Real winners are paid off-chain
via `emission-transfer-api.metanova-labs.ai`. Our 0 incentive on-chain is
expected behavior — but we also received 0 off-chain, since we were never scored.

## Validator pipeline: what we pass vs. what we fail

We tested every documented gate against our submitted sequences:

| Gate | Our status | Notes |
|------|-----------|-------|
| Content hash (file integrity) | PASS | SHA256 matches on-chain and GitHub |
| Drand timelock decryption | PASS | Correct reveal round, valid ciphertext |
| Payload format (`~\|{sequence}`) | PASS | Correct format verified |
| Sequence length (90-150) | PASS | All sequences 109-125 AAs |
| Valid amino acids | PASS | Standard 20 AAs only |
| No homopolymers (<=6) | PASS | |
| No signal peptide | PASS | |
| V germline identity (>=82%) | PASS | 82.65% to 84.69% |
| J germline identity (>=86%) | PASS | 92.86% |
| V gene epoch match | PASS | All match allowed_v |
| VHH classification | PASS | H49=E/Q, H50=R, canonical cysteines |
| CDR3 length match | PASS | Exact match to allowed_cdr3_length |
| Archive uniqueness (<95% identity) | PASS | |
| CDR3 novelty (<85% identity) | PASS | |
| HuggingFace archive check | PASS | |
| TNP structure quality | PASS | |
| Aggrescan3D (avg <= -0.680, max <= 1.580) | PASS | avg -0.715, max 1.239 |
| **PSH strict rejection (<=89)** | **FAIL** | Validator sees ~89.5, we see ~85 |

## Competitive landscape (from score-sharing API)

The subnet is active with 26 nanobodies scored per epoch. Top miners:

| Miner | Accept rate | Wins (61 epochs) |
|-------|-------------|-------------------|
| frankji-groundcontrol | 83% | 30 |
| QTTECH2023 | 73% | 13 |
| thomasvangurp | 78% | 5 |
| saon-miner | 59% | 5 |
| bitty-labs | 68% | 2 |
| cczulaika | 40% | 2 |
| uchiha0818 | 26% | 1 |
| aichorx | 5% | 0 |
| **us (UID 22)** | **0%** | **0** |

The accept rate is the primary determinant of success. Being admitted to scoring
is the hard part; ranking within scored submissions is secondary.

## What to fix before the next run

### Critical (blocking)

1. **Calibrate PSH to the validator's scale.** Either:
   - Match the exact TNP version and hydrophobicity scale the validator uses
   - Set a much lower local threshold (PSH <= 75) to create margin
   - Reverse-engineer the mapping by submitting known sequences and comparing

### Important

2. **Persist BoltzGen scores to disk** — without this, we cannot evaluate
   performance retroactively.
3. **Add subtensor retry logic** — the 429 crash lost 3 epochs / 5 hours.
4. **Verify GitHub repo before committing** — the epoch 26172 dangling pointer
   was avoidable.

### Nice to have

5. **Widen the candidate funnel** — relax or restructure the liability pre-filter
   to let more candidates through to BoltzGen.
6. **Drop the surrogate fallback** — if BoltzGen is unavailable, skip the epoch
   rather than submitting a random pick.
7. **Log epoch-level timing** — track how much time budget remains after each
   pipeline stage.

## Financial summary

| Item | Cost |
|------|------|
| Lium GPU rental (2x RTX 4090, 15h) | ~$13.50 |
| Subtensor registration (UID 22) | ~0.1 TAO |
| Revenue earned | 0 TAO |
| **Net loss** | **~$13.50 + 0.1 TAO** |

## Files

All production data is preserved in this repository:

- `journal_final.jsonl` — mining journal (12 entries)
- `mine_live_final.log` — complete miner log (541 lines)
- `scores_*.csv` — BoltzGen observation scores (6 files)
- `boltzgen_*.log` — scoring run logs (5 files)
- `docs/OBSERVATIONS.md` — real-time observations and analysis
