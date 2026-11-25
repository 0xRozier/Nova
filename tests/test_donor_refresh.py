"""Refreshing the donor pool inside an epoch, which the miner used to refuse.

Pool freshness is the only bottleneck this repo has both measured and cured:
+0.075 of design_iiptm on GPU (FINDINGS 49), for about 22 s of decryption per
epoch. Before this, the miner detected a stale pool and spent the epoch telling
a human to run three commands -- so a measured, free win was lost every time it
was needed.

What has to hold:

  - holes are filled newest first, not just the tail. Fetching only the last
    few epochs moves donor_pool_lag to zero while leaving a permanent hole, so
    the lag would report a fresh pool that is missing most of the week;
  - the cost is bounded, and a refresh that does not fit the epoch is declined
    rather than started;
  - a failing step leaves the pool alone and says why.

    python -m pytest tests/test_donor_refresh.py -q
"""
import subprocess

import pytest

import scripts.mine as mine


def test_it_fills_the_newest_holes_not_only_the_tail():
    present = {26300, 26301, 26305}          # 26302-26304 are holes
    todo = mine.epochs_to_decrypt(present, epoch=26310, limit=8)
    assert todo[0] == 26308, "newest first"
    assert 26304 in todo and 26303 in todo and 26302 in todo, "holes included"
    assert 26305 not in todo and 26301 not in todo, "what we have is not refetched"


def test_the_epoch_in_progress_and_the_one_before_are_left_out():
    """Neither has anything decryptable or reliably scored yet."""
    todo = mine.epochs_to_decrypt(set(), epoch=26310, limit=50)
    assert max(todo) == 26308


def test_it_never_looks_further_back_than_the_window():
    todo = mine.epochs_to_decrypt(set(), epoch=26310, limit=1000, window=10)
    assert min(todo) == 26299 and len(todo) == 10


def test_the_limit_bounds_what_one_epoch_pays():
    todo = mine.epochs_to_decrypt(set(), epoch=26310, limit=3)
    assert todo == [26308, 26307, 26306]


def test_a_refresh_that_does_not_fit_is_declined_before_it_starts(monkeypatch):
    called = []
    monkeypatch.setattr(subprocess, "run",
                        lambda *a, **k: called.append(a) or None)
    did, why = mine.refresh_donor_pool("Q9NZQ7", list(range(26300, 26308)),
                                       seconds_left=60, log=lambda *a: None)
    assert did is False and "only 1 min" in why
    assert not called, "nothing may run when the budget says no"


def test_nothing_to_decrypt_is_not_a_failure(monkeypatch):
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: pytest.fail("ran anyway"))
    did, why = mine.refresh_donor_pool("Q9NZQ7", [], 10_000, log=lambda *a: None)
    assert did is False and "nothing" in why


def test_the_three_steps_run_in_the_order_that_moves_the_lag(monkeypatch):
    """Without the first two the new sequences carry no epoch."""
    seen = []

    class Ok:
        returncode = 0
        stdout = stderr = ""

    def fake(cmd, **kwargs):
        seen.append(cmd[2])
        return Ok()

    monkeypatch.setattr(subprocess, "run", fake)
    did, why = mine.refresh_donor_pool("Q9NZQ7", [26300, 26301], 10_000,
                                       log=lambda *a: None)
    assert did is True, why
    assert seen == ["scripts.fetch_submissions", "scripts.analyze_nanobodies",
                    "scripts.build_nanobody_dataset"]


def test_a_failing_step_stops_the_refresh_and_says_why(monkeypatch):
    seen = []

    class Fail:
        returncode = 1
        stdout = ""
        stderr = "the archive node would not answer"

    class Ok:
        returncode = 0
        stdout = stderr = ""

    def fake(cmd, **kwargs):
        seen.append(cmd[2])
        return Fail() if cmd[2] == "scripts.analyze_nanobodies" else Ok()

    monkeypatch.setattr(subprocess, "run", fake)
    did, why = mine.refresh_donor_pool("Q9NZQ7", [26300], 10_000, log=lambda *a: None)
    assert did is False
    assert "analyze_nanobodies failed" in why and "would not answer" in why
    assert "scripts.build_nanobody_dataset" not in seen, "must not refeaturize"


def test_it_decrypts_the_span_the_holes_cover(monkeypatch):
    spans = []

    class Ok:
        returncode = 0
        stdout = stderr = ""

    def fake(cmd, **kwargs):
        if cmd[2] == "scripts.fetch_submissions":
            spans.append((cmd[cmd.index("--from-epoch") + 1],
                          cmd[cmd.index("--to-epoch") + 1]))
        return Ok()

    monkeypatch.setattr(subprocess, "run", fake)
    mine.refresh_donor_pool("Q9NZQ7", [26308, 26304, 26303], 10_000,
                            log=lambda *a: None)
    assert spans == [("26303", "26308")]


# --- FINDINGS 83/84: the shortlist read off the clock ----------------------

def test_the_shortlist_uses_the_room_the_cpu_stages_leave():
    """Epoch 26312's journal: 755 s for everything but BoltzGen, 68 min window."""
    from scripts.mine import shortlist_for_budget
    assert shortlist_for_budget(68 * 60, 755) == 26
    # the old fixed default left two thirds of the dominant lever unused
    assert shortlist_for_budget(68 * 60, 755) > 3 * 8


def test_the_shortlist_never_reaches_zero():
    """A badly scored candidate still beats abstaining, which scores zero."""
    from scripts.mine import shortlist_for_budget
    for seconds_left, elapsed in ((60, 0), (600, 900), (0, 0), (-100, 50)):
        assert shortlist_for_budget(seconds_left, elapsed) == 1


def test_the_shortlist_is_capped():
    from scripts.mine import SHORTLIST_CAP, shortlist_for_budget
    assert shortlist_for_budget(100 * 3600, 0) == SHORTLIST_CAP


def test_the_shortlist_shrinks_as_the_epoch_runs_out():
    from scripts.mine import shortlist_for_budget
    sizes = [shortlist_for_budget(68 * 60, e) for e in (0, 600, 1800, 3000)]
    assert sizes == sorted(sizes, reverse=True)


def test_the_throughput_matches_what_the_runs_measured():
    """211 sequences, three runs: 167.0 + 121.1 + 143.6 min."""
    from scripts.mine import BOLTZGEN_SECONDS_PER_SEQUENCE
    measured = (167.0 + 121.1 + 143.6) * 60 / (85 + 56 + 70)
    assert abs(BOLTZGEN_SECONDS_PER_SEQUENCE - measured) < 2.0
