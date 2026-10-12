import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from igem_warsaw_rnaseq.config import ConfigError, ReferenceSpec, load_config
from igem_warsaw_rnaseq.quant import collect, routing, salmon
from igem_warsaw_rnaseq.quant.runner import QuantResult


def _spec(name: str, tax_id: int | None) -> ReferenceSpec:
    return ReferenceSpec(name, "", None, None, None, tax_id, ())


def test_most_specific_reference_prefers_deepest_taxon() -> None:
    refs = {5543: "genus", 63577: "species", 452589: "strain"}
    lineage = "1;2759;5543;63577;452589"
    assert routing.most_specific_reference(lineage, refs) == "strain"
    assert routing.most_specific_reference("1;2759;5543;5547", refs) == "genus"
    assert routing.most_specific_reference("1;2759;4751", refs) is None
    assert routing.most_specific_reference("", refs) is None


def test_route_overrides_and_taxonomy() -> None:
    kept = pd.DataFrame(
        {
            "run_id": ["A", "B", "C"],
            "tax_lineage": ["1;5543;63577", "1;5543;63577", "1;9606"],
            "reference": ["", "", ""],
        }
    )
    refs = {"atro": _spec("atro", 63577), "other": _spec("other", None)}
    table = routing.route(kept, refs, {"B": "other"})
    assert table.set_index("run_id").loc["A", "reference"] == "atro"
    assert table.set_index("run_id").loc["A", "routing_source"] == "taxonomy"
    assert table.set_index("run_id").loc["B", "reference"] == "other"
    assert table.set_index("run_id").loc["B", "routing_source"] == "override"
    assert table.set_index("run_id").loc["C", "status"] == routing.NO_REFERENCE
    assert list(routing.routable(table)["run_id"]) == ["A", "B"]


def test_load_overrides_rejects_unknown_reference(tmp_path: Path) -> None:
    f = tmp_path / "o.tsv"
    f.write_text("run_id\treference\nA\tnope\n")
    with pytest.raises(ConfigError, match="unknown references"):
        routing.load_overrides(f, {"atro": _spec("atro", 1)})


def test_duplicate_tax_id_rejected(tmp_path: Path) -> None:
    cfg = tmp_path / "c.yaml"
    cfg.write_text(
        "references:\n  a: {ncbi_accession: X, tax_id: 5}\n  b: {ncbi_accession: Y, tax_id: 5}\n"
    )
    with pytest.raises(ConfigError, match="share tax_id"):
        load_config(cfg)


def test_quant_config_defaults(tmp_path: Path) -> None:
    cfg = tmp_path / "c.yaml"
    cfg.write_text("quant: {}\n")
    q = load_config(cfg).quant
    assert q.bootstraps == 30
    assert q.min_pct_mapped == 5.0
    assert q.crossmap_min_median_tpm is None


def _write_salmon_outputs(out: Path, names: list[str], boots: np.ndarray) -> None:
    (out / "aux_info" / "bootstrap").mkdir(parents=True)
    with gzip.open(out / "aux_info" / "bootstrap" / "names.tsv.gz", "wt") as fh:
        fh.write("\t".join(names) + "\n")
    with gzip.open(out / "aux_info" / "bootstrap" / "bootstraps.gz", "wb") as fh:
        fh.write(boots.astype(np.float64).tobytes())
    (out / "aux_info" / "meta_info.json").write_text(
        '{"num_processed": 1000, "num_mapped": 400, "percent_mapped": 40.0,'
        ' "num_decoy_fragments": 100, "num_bootstraps": 3}'
    )
    (out / "lib_format_counts.json").write_text('{"expected_format": "ISR"}')
    (out / "quant.sf").write_text(
        "Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
        "t1\t100\t80\t10\t5\nt2\t100\t80\t30\t15\nt3\t100\t80\t0\t0\n"
    )


def test_stats_bootstraps_and_gene_aggregation(tmp_path: Path) -> None:
    boots = np.array([[4, 14, 0], [6, 16, 0], [5, 15, 0]])
    _write_salmon_outputs(tmp_path, ["t1", "t2", "t3"], boots)
    stats = salmon.read_stats(tmp_path)
    assert stats.pct_mapped == 40.0
    assert stats.pct_decoy == 10.0
    assert stats.detected_libtype == "ISR"

    b = salmon.load_bootstraps(tmp_path)
    assert b.shape == (3, 3)
    genes = salmon.aggregate_genes(
        tmp_path / "quant.sf", {"t1": "g1", "t2": "g1", "t3": "g2"}, b
    ).set_index("gene_id")
    assert genes.loc["g1", "counts"] == 20
    assert genes.loc["g1", "tpm"] == 40
    assert genes.loc["g1", "boot_cv"] == pytest.approx(
        np.std([18, 22, 20], ddof=1) / 20
    )
    assert np.isnan(genes.loc["g2", "boot_cv"])


def test_read_stats_without_lib_format_counts(tmp_path: Path) -> None:
    (tmp_path / "aux_info").mkdir()
    (tmp_path / "aux_info" / "meta_info.json").write_text(
        '{"num_processed": 10, "num_mapped": 0, "percent_mapped": 0.0}'
    )
    stats = salmon.read_stats(tmp_path)
    assert stats.pct_mapped == 0.0
    assert stats.detected_libtype == ""
    assert stats.pct_decoy == 0.0
    assert salmon.load_bootstraps(tmp_path) is None


def test_build_table_decisions() -> None:
    rt = pd.DataFrame({"run_id": ["ok", "low", "fail", "todo"], "reference": ["r"] * 4})
    results = {
        "ok": QuantResult("ok", "r", 100, 80, 80.0, 1.0, "ISR", 5000),
        "low": QuantResult("low", "r", 100, 2, 2.0, 30.0, "IU", 10),
        "fail": QuantResult("fail", "r", drop_reason="salmon_failed"),
    }
    table = collect.build_table(rt, results, {"ok": "PRJ1"}, 5.0).set_index("run_id")
    assert table.loc["ok", "decision"] == "keep"
    assert table.loc["ok", "bioproject"] == "PRJ1"
    assert table.loc["low", "drop_reason"] == collect.LOW_MAPPING_RATE
    assert table.loc["fail", "drop_reason"] == "salmon_failed"
    assert table.loc["todo", "decision"] == "pending"
    assert (
        collect.build_table(rt, results, {}, None)
        .set_index("run_id")
        .loc["low", "decision"]
        == "keep"
    )


def test_embed_pca_and_spearman() -> None:
    rng = np.random.default_rng(0)
    base = rng.random(50)
    matrix = pd.DataFrame(
        {"a": base, "b": base + 0.01 * rng.random(50), "c": rng.random(50)}
    )
    out = collect.embed(matrix)
    assert (
        out.loc["a", "median_spearman_to_others"]
        > out.loc["c", "median_spearman_to_others"]
    )
    assert abs(out.loc["a", "pca1"] - out.loc["b", "pca1"]) < abs(
        out.loc["a", "pca1"] - out.loc["c", "pca1"]
    )
    single = collect.embed(matrix[["a"]])
    assert single.isna().all().all()
