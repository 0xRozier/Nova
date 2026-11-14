"""The gate that asks the field whether our derivation is still the validator's.

Everything in an epoch is built on two numbers derived from one block hash. If
that derivation drifts from upstream's, every candidate is for the wrong
challenge and every submission is refused -- which looks exactly like a lost
epoch, so nothing says so.

`scripts.check_config` cannot cover this: it compares our pin to a ref, so it
is blind to values still on a branch, and PR 102's were announced for epoch
26331 while living on `rxn-group-filtering`. Worse, the announcement and the
branch's code disagree about whether the nanobody parameters move at all.

So the check does not read code. An accepted nanobody submission must carry the
epoch's V gene and CDR3 length, so the accepted submissions of an epoch are a
witness to its real parameters. These tests pin the decision table of
`confirmed()`, which is what the miner refuses on -- and the case that matters
is a disagreement at or after the switchover, because that is the one that
costs epochs going forward.

    python -m pytest tests/test_verify_challenge.py -q
"""
import pytest

from nova_lab.config import RXN_GROUP_FILTERING_FROM_EPOCH as SWITCH
import scripts.verify_challenge as vc

BEFORE, AFTER = SWITCH - 5, SWITCH + 5


@pytest.fixture
def field(monkeypatch):
    """Let a test say what the field submitted and what we derived."""
    state = {"witnesses": {}, "agrees": True, "raises": None}

    def accepted_by_epoch(target):
        if not state["witnesses"]:
            raise SystemExit("no decrypted submissions")
        return state["witnesses"]

    def agrees(epoch, sequences):
        if state["raises"] is not None:
            raise state["raises"]
        return state["agrees"]

    monkeypatch.setattr(vc, "accepted_by_epoch", accepted_by_epoch)
    monkeypatch.setattr(vc, "agrees", agrees)
    return state


def test_agreement_lets_the_miner_through(field):
    field["witnesses"] = {BEFORE: ["SEQ"], AFTER: ["SEQ", "SEQ2"]}
    ok, why = vc.confirmed("Q9NZQ7", log=lambda *a: None)
    assert ok
    assert f"epoch {AFTER}" in why and "2 accepted" in why


def test_a_disagreement_after_the_switchover_stops_it(field):
    """The expensive case: from here on every candidate is for another challenge."""
    field["witnesses"] = {AFTER: ["SEQ"]}
    field["agrees"] = False
    ok, why = vc.confirmed("Q9NZQ7", log=lambda *a: None)
    assert not ok
    assert "does NOT match" in why and str(SWITCH) in why


def test_a_disagreement_before_the_switchover_is_history(field):
    """It says something about an old epoch, not about the one being mined."""
    field["witnesses"] = {BEFORE: ["SEQ"]}
    field["agrees"] = False
    ok, why = vc.confirmed("Q9NZQ7", log=lambda *a: None)
    assert ok and "history" in why


def test_only_the_newest_epoch_decides(field):
    """An old disagreement must not block once a newer epoch agrees.

    The newest verifiable epoch is the only one that says anything about now,
    and this check already trails reality by one epoch; trailing it by the
    whole archive would make it useless.
    """
    field["witnesses"] = {BEFORE: ["OLD"], AFTER: ["NEW"]}
    seen = []
    field["agrees"] = True
    original = vc.agrees

    def record(epoch, sequences):
        seen.append(epoch)
        return original(epoch, sequences)

    # The fixture already replaced vc.agrees; wrap that replacement so the
    # monkeypatch still unwinds at the end of the test.
    monkeypatch = pytest.MonkeyPatch()
    monkeypatch.setattr(vc, "agrees", record)
    try:
        ok, _ = vc.confirmed("Q9NZQ7", log=lambda *a: None)
    finally:
        monkeypatch.undo()
    assert ok and seen == [AFTER]


def test_no_data_does_not_block(field):
    """A machine that has not decrypted anything yet is not a machine in drift.

    Blocking here would make the miner refuse to start on a fresh install,
    which is a worse failure than the one being guarded against.
    """
    field["witnesses"] = {}
    ok, why = vc.confirmed("Q9NZQ7", log=lambda *a: None)
    assert ok and "check against" in why


def test_a_chain_failure_does_not_block(field):
    """The archive node answers 429 under load; that is not a verdict.

    It is also exactly what happened in production for five hours, and a gate
    that reads an outage as drift would have turned that into a week of
    skipped epochs.
    """
    field["witnesses"] = {AFTER: ["SEQ"]}
    field["raises"] = RuntimeError("HTTP 429")
    ok, why = vc.confirmed("Q9NZQ7", log=lambda *a: None)
    assert ok and "could not reach the chain" in why


def test_the_switchover_epoch_itself_counts_as_after(field):
    field["witnesses"] = {SWITCH: ["SEQ"]}
    field["agrees"] = False
    ok, _ = vc.confirmed("Q9NZQ7", log=lambda *a: None)
    assert not ok
