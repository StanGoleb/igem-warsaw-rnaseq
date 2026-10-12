import numpy as np
import pandas as pd
import pytest

from igem_warsaw_rnaseq.config import AnchorConfig, ClassConfig
from igem_warsaw_rnaseq.normalize import anchors, combat, rank, run

CLASSES = ClassConfig(
    window=0.10,
    min_cons_frac=0.95,
    silent_min_median_pct=0.75,
    inducible_min_amplitude=0.30,
    top_pct=0.10,
)


def test_rank_pct_ties_and_absent() -> None:
    tpm = pd.DataFrame(
        {"s1": [10.0, 5.0, 0.0, 0.0], "s2": [1.0, np.nan, 3.0, 2.0]},
        index=["a", "b", "c", "d"],
    )
    pct = rank.rank_pct(tpm)
    assert pct.loc["a", "s1"] == 0.25 and pct.loc["b", "s1"] == 0.5
    assert pct.loc["c", "s1"] == 1.0 and pct.loc["d", "s1"] == 1.0
    assert pct.loc["b", "s2"] == 1.0
    assert pct.loc["c", "s2"] == pytest.approx(1 / 3)


def test_cons_frac() -> None:
    pct = pd.DataFrame([[0.10, 0.12, 0.15, 0.90], [0.1, 0.3, 0.5, 0.7]])
    out = rank.cons_frac(pct, 0.10)
    assert out.iloc[0] == 0.75
    assert out.iloc[1] == 0.25


def test_classify_four_classes() -> None:
    n = 40
    rng = np.random.default_rng(1)
    constitutive = rng.uniform(0.02, 0.06, n)
    silent = rng.uniform(0.90, 0.95, n)
    inducible = np.concatenate([rng.uniform(0.7, 0.9, n - 3), [0.01, 0.02, 0.03]])
    variable = rng.uniform(0.0, 1.0, n)
    tpm = pd.DataFrame(
        1.0 - np.vstack([constitutive, silent, inducible, variable]),
        index=["c", "s", "i", "v"],
        columns=[f"x{i}" for i in range(n)],
    )
    tpm = tpm.rank(axis=0, ascending=False)
    pct = pd.DataFrame(
        np.vstack([constitutive, silent, inducible, variable]),
        index=tpm.index,
        columns=tpm.columns,
    )
    m = pd.DataFrame(
        {
            "median_pct": pct.median(axis=1),
            "p05_pct": pct.quantile(0.05, axis=1),
            "cons_frac": rank.cons_frac(pct, CLASSES.window),
            "induction_amplitude": pct.median(axis=1) - pct.quantile(0.05, axis=1),
        }
    )
    labelled = rank.classify(m, CLASSES)
    assert (
        labelled.loc["c", "class"] == "constitutive" and labelled.loc["c", "is_strong"]
    )
    assert labelled.loc["s", "class"] == "silent"
    assert labelled.loc["i", "class"] == "inducible"
    assert labelled.loc["v", "class"] == "variable"


def test_genorm_prefers_stable_gene() -> None:
    rng = np.random.default_rng(0)
    depth = rng.normal(0, 1, 30)
    panel = pd.DataFrame(
        {
            "stable": 8 + depth,
            "stable2": 9 + depth + rng.normal(0, 0.05, 30),
            "noisy": 8 + depth + rng.normal(0, 2, 30),
        }
    ).T
    m = anchors.genorm_m(panel)
    assert m["noisy"] > m["stable"]
    cfg = AnchorConfig(n_auto=2, min_mean_log2tpm=0.0, candidate_top_n=10)
    chosen, table = anchors.select(panel, [], cfg)
    assert chosen == ["stable", "stable2"] or chosen == ["stable2", "stable"]
    assert set(table["source"]) == {"genorm"}
    chosen, table = anchors.select(panel, ["noisy", "missing"], cfg)
    assert chosen == ["noisy"] and list(table["source"]) == ["config"]


def test_anchor_normalise() -> None:
    log2 = pd.DataFrame(
        {"s1": [5.0, 3.0, 1.0], "s2": [6.0, 4.0, 2.0]}, index=["a", "b", "g"]
    )
    out = anchors.normalise(log2, ["a", "b"])
    assert out.loc["g", "s1"] == pytest.approx(1.0 - 4.0)
    assert out.loc["g", "s2"] == pytest.approx(2.0 - 5.0)
    assert anchors.normalise(log2, []).isna().all().all()


def test_combat_skips_without_two_batches() -> None:
    log2 = pd.DataFrame(
        np.random.default_rng(2).normal(size=(5, 4)), columns=list("abcd")
    )
    batch = pd.Series(["P1", "P1", "P2", "P3"], index=list("abcd"))
    out, entry = combat.run(log2, batch, None, min_batch_size=2)
    assert entry["applied"] is False and out.isna().all().all()
    assert entry["n_batches"] == 1


def test_combat_corrects_batch_offset() -> None:
    rng = np.random.default_rng(3)
    base = rng.normal(5, 1, size=(50, 8))
    data = base.copy()
    data[:, 4:] += 3.0
    log2 = pd.DataFrame(data, columns=[f"s{i}" for i in range(8)])
    batch = pd.Series(["A"] * 4 + ["B"] * 4, index=log2.columns)
    out, entry = combat.run(log2, batch, None, min_batch_size=2)
    assert entry["applied"] is True
    gap = out.iloc[:, 4:].mean().mean() - out.iloc[:, :4].mean().mean()
    assert abs(gap) < 0.5


def test_long_table() -> None:
    layers = {
        "tpm": pd.DataFrame({"s1": [1.0, 2.0]}, index=["g1", "g2"]),
        "log2tpm": pd.DataFrame({"s1": [1.0, 1.58]}, index=["g1", "g2"]),
    }
    table = run.long_table(layers, "gene_id")
    assert list(table.columns) == ["gene_id", "sample_id", "tpm", "log2tpm"]
    assert len(table) == 2
