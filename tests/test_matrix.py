import pandas as pd
import pytest

from igem_warsaw_rnaseq.matrix import build, representative, samples


def _genes(**kw: tuple[float, float]) -> pd.DataFrame:
    return pd.DataFrame(
        {"counts": [v[0] for v in kw.values()], "tpm": [v[1] for v in kw.values()]},
        index=pd.Index(list(kw), name="gene_id"),
    )


def test_merge_runs_sums_counts_and_weights_tpm() -> None:
    a = _genes(g1=(90, 10.0), g2=(10, 90.0))
    b = _genes(g1=(300, 50.0), g2=(0, 50.0))
    merged = samples.merge_runs({"a": a, "b": b})
    assert merged.loc["g1", "counts"] == 390
    assert merged.loc["g1", "tpm"] == pytest.approx((10 * 100 + 50 * 300) / 400)
    assert samples.merge_runs({"a": a}).equals(a[["counts", "tpm"]])


def test_sample_table_groups_lanes_and_checks_reference() -> None:
    qc = pd.DataFrame(
        {"run_id": ["r1", "r2", "r3", "r4"], "decision": ["keep"] * 3 + ["drop"]}
    )
    routing = pd.DataFrame({"run_id": ["r1", "r2", "r3", "r4"], "reference": ["A"] * 4})
    meta = pd.DataFrame(
        {
            "run_accession": ["r1", "r2", "r3", "r4"],
            "experiment_accession": ["e1", "e1", "", "e3"],
            "study_accession": ["P1", "P1", "P2", "P3"],
        }
    )
    table = samples.sample_table(qc, routing, meta).set_index("sample_id")
    assert table.loc["e1", "runs"] == "r1;r2"
    assert table.loc["e1", "n_runs"] == 2
    assert table.loc["r3", "n_runs"] == 1
    assert "e3" not in table.index
    routing.loc[1, "reference"] = "B"
    with pytest.raises(ValueError, match="several references"):
        samples.sample_table(qc, routing, meta)


def test_leaf_owner_longest_prefix() -> None:
    refs = ["sp", "sp_strain"]
    assert representative.leaf_owner("sp_strain_g1", refs) == ("sp_strain", "g1")
    assert representative.leaf_owner("sp_g1", refs) == ("sp", "g1")
    assert representative.leaf_owner("other_g1", refs) is None


TREE = "((A_x:0.1,B_y:0.1)n1:0.5,(A_z:0.1,C_w:0.1)n2:0.5)n0;"


def test_choose_by_tree_picks_orthologous_copy() -> None:
    assert representative.choose_by_tree(TREE, "A", ["x", "z"], ["A", "B", "C"]) in {
        "x",
        "z",
    }
    tree = "((A_x:0.1,B_y:0.1)n1:0.1,(A_z:2.0,C_w:0.1)n2:0.1)n0;"
    assert representative.choose_by_tree(tree, "A", ["x", "z"], ["A", "B", "C"]) == "x"
    assert representative.choose_by_tree(tree, "A", ["x"], ["A", "B", "C"]) is None


def test_choose_fallback_chain() -> None:
    refs = ["A", "B"]
    assert representative.choose("A", ["g"], None, None, {}, refs) == (
        "g",
        "single_copy",
    )
    assert representative.choose(
        "A", ["g1", "g2"], None, {"g1": 1.0, "g2": 5.0}, {}, refs
    ) == (
        "g2",
        "max_mean_tpm",
    )
    assert representative.choose(
        "A", ["g1", "g2"], None, None, {"g1": 10, "g2": 3}, refs
    ) == (
        "g1",
        "longest_protein",
    )


def test_orthogroup_matrix_uses_representatives_and_na() -> None:
    reps = pd.DataFrame(
        {
            "orthogroup": ["OG1", "OG2", "OG1"],
            "reference": ["A", "A", "B"],
            "gene_id": ["a1", "a2", "b1"],
            "criterion": ["single_copy"] * 3,
            "n_paralogs_collapsed": [0, 0, 0],
        }
    )
    per_ref = {
        "A": pd.DataFrame(
            {"s1": [1.0, 2.0]}, index=pd.Index(["a1", "a2"], name="gene_id")
        ),
        "B": pd.DataFrame({"s2": [3.0]}, index=pd.Index(["b1"], name="gene_id")),
    }
    m = build.orthogroup_matrix(reps, per_ref, ["OG1", "OG2", "OG3"])
    assert m.loc["OG1", "s1"] == 1.0 and m.loc["OG1", "s2"] == 3.0
    assert m.loc["OG2", "s1"] == 2.0 and pd.isna(m.loc["OG2", "s2"])
    assert m.loc["OG3"].isna().all()
