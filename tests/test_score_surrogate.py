"""The surrogate's honesty guarantees, not its accuracy.

A surrogate is easy to flatter: split at random and near-twin sequences land on
both sides, so the model recites a framework it has already seen. Our batches
are built from arms of 32 sequences differing by one loop, which is exactly the
shape that breaks a random split. These tests pin the two things that make the
reported number mean something -- whole clusters held out, and a
decision-relevant metric that is 0 for a useless model and 1 for a perfect one.
"""
import numpy as np
import pandas as pd
import pytest

from scripts.train_score_surrogate import (cluster, decision_value, descriptors,
                                           grouped_spearman, noise_ceiling)

FRAME = ("EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYRMYWVRQPPGKEREWVSAISAGGGSTYYGDSVKGRFTIS"
         "RDNAKNTVYLQMNSLKPEDTAVYYCA{}WGQGTLVTVSS")


def test_near_twins_share_a_cluster_and_strangers_do_not():
    """One loop apart is the same cluster; a different framework is not.

    Measured, and it is why the default link is 0.80 and not the validator's
    0.90: two shuffles of one CDR3 sit at 5-mer Jaccard 0.854, so under the
    subnet's own diversity rule they are distinct sequences while for a
    cross-validation split they are twins.
    """
    from scripts.train_score_surrogate import VALIDATOR_LINK, kmers
    twins = [FRAME.format(c) for c in ("NRAGWGMGDY", "NRGGMAWGDY", "NRAWGGMGDY")]
    stranger = ("QITLKESGPALVKPTQTLTLTCTFSGFSLSTSGMGVGWIRQPPGKALEWLAHIWWDDDKYYNPSLKSR"
                "LTISKDTSKNQVVLTMTNMDPVDTATYYCARRANYGTSYDYGMDYWGQGTLVTVSS")
    a, b = kmers(twins[0], 5), kmers(twins[1], 5)
    jaccard = len(a & b) / len(a | b)
    assert 0.80 <= jaccard < VALIDATOR_LINK, jaccard
    groups = cluster(twins + [stranger])
    assert len(set(groups[:3])) == 1, "le defaut doit regrouper les jumelles"
    assert groups[3] not in set(groups[:3])


def test_grouped_split_never_puts_a_cluster_on_both_sides():
    """The guarantee the reported rho rests on, checked by construction."""
    from sklearn.dummy import DummyRegressor
    seen: list[set] = []

    class Spy(DummyRegressor):
        def fit(self, X, y, **kw):
            seen.append(set(np.asarray(X)[:, 0].tolist()))
            return super().fit(X, y, **kw)

    groups = np.repeat(np.arange(20), 5)            # 20 clusters of 5
    X = groups.reshape(-1, 1).astype(float)         # feature IS the cluster id
    y = np.random.default_rng(0).normal(size=len(groups))
    grouped_spearman(X, y, groups, Spy(), folds=4)
    assert seen, "aucun pli n'a ete entraine"
    for train_ids in seen:
        held_out = set(np.unique(groups).tolist()) - train_ids
        assert held_out, "un pli a vu toutes les grappes : la validation est vide"
        assert not (train_ids & held_out)


def test_a_useless_model_captures_nothing_and_a_perfect_one_captures_all():
    rng = np.random.default_rng(3)
    y = rng.normal(size=400)
    bands = {"tout": np.ones(len(y), bool)}
    perfect = decision_value(y, y.copy(), bands, pool=10, draws=400)[0]
    noise = decision_value(y, rng.normal(size=len(y)), bands, pool=10, draws=400)[0]
    inverted = decision_value(y, -y, bands, pool=10, draws=400)[0]
    assert perfect["part_captee"] == pytest.approx(1.0, abs=1e-9)
    assert abs(noise["part_captee"]) < 0.15
    assert inverted["part_captee"] < -0.5


def test_the_noise_ceiling_falls_when_the_label_is_noisier():
    base = pd.DataFrame({"y": np.linspace(0, 1, 100), "runs": 1,
                         "spread": [np.nan] * 100})
    base.loc[:49, ["runs", "spread"]] = [2, 0.01]
    quiet = noise_ceiling(base)[1]
    loud = base.copy()
    loud.loc[:49, "spread"] = 0.30
    assert noise_ceiling(loud)[1] < quiet


def test_descriptors_are_finite_for_a_sequence_that_cannot_be_numbered():
    """regions() returns None on junk; the features must not become NaN."""
    d = descriptors("AAAAAAAAAAAAAAAAAAAAAAAA")
    assert all(np.isfinite(v) for v in d.values())
    assert d["cdr3_length"] == 0
