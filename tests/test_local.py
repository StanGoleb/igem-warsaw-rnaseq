import gzip
from pathlib import Path

import pandas as pd
import pytest

from igem_warsaw_rnaseq.config import ConfigError, ReferenceSpec, load_config
from igem_warsaw_rnaseq.qc import decision, local
from igem_warsaw_rnaseq.quant import routing

REFS = {"a": ReferenceSpec("a", "", None, None, None, 1, ())}


def _fastq(path: Path) -> Path:
    with gzip.open(path, "wt") as fh:
        fh.write("@r1\nACGTACGTAC\n+\nFFFFFFFFFF\n")
    return path


def test_load_local_runs_schema_and_layout(tmp_path: Path) -> None:
    _fastq(tmp_path / "s_1.fq.gz")
    _fastq(tmp_path / "s_2.fq.gz")
    _fastq(tmp_path / "t.fq.gz")
    table_path = tmp_path / "runs_local.tsv"
    table_path.write_text(
        "run_id\tfastq_1\tfastq_2\treference\tbatch\tsample\n"
        "L1\ts_1.fq.gz\ts_2.fq.gz\ta\tB1\tS1\n"
        "L2\tt.fq.gz\t\ta\tB1\t\n"
        "L3\tmissing.fq.gz\t\ta\tB2\t\n"
    )
    meta = local.load_local_runs(table_path, REFS).set_index("run_accession")
    assert meta.loc["L1", "library_layout"] == "PAIRED"
    assert meta.loc["L2", "library_layout"] == "SINGLE"
    assert meta.loc["L1", "experiment_accession"] == "S1"
    assert meta.loc["L2", "experiment_accession"] == "L2"
    assert meta.loc["L1", "study_accession"] == "B1"
    assert meta.loc["L1", "fastq_ftp"].split(";")[0] == str(tmp_path / "s_1.fq.gz")
    assert (meta["source"] == "local").all() and not meta["in_ena"].any()

    rows = meta.reset_index().to_dict("records")
    reasons = {
        r["run_accession"]: decision.metadata_drop_reason(r, {"RNA-Seq"}, set())
        for r in rows
    }
    assert reasons == {"L1": None, "L2": None, "L3": decision.MISSING_FILE}


def test_load_local_runs_validation(tmp_path: Path) -> None:
    p = tmp_path / "bad.tsv"
    p.write_text("run_id\tfastq_1\treference\tbatch\nL1\tx\tnope\tB\n")
    with pytest.raises(ConfigError, match="unknown references"):
        local.load_local_runs(p, REFS)
    p.write_text("run_id\tfastq_1\treference\tbatch\nL1\tx\ta\tB\nL1\ty\ta\tB\n")
    with pytest.raises(ConfigError, match="duplicated"):
        local.load_local_runs(p, REFS)
    p.write_text("run_id\tfastq_1\n")
    with pytest.raises(ConfigError, match="missing columns"):
        local.load_local_runs(p, REFS)


def test_route_prefers_declared_reference() -> None:
    kept = pd.DataFrame(
        {
            "run_id": ["E", "L"],
            "tax_lineage": ["1;5", ""],
            "reference": ["", "a"],
        }
    )
    refs = {"a": ReferenceSpec("a", "", None, None, None, 5, ())}
    table = routing.route(kept, refs, {}).set_index("run_id")
    assert table.loc["E", "routing_source"] == routing.SOURCE_TAXONOMY
    assert table.loc["L", "routing_source"] == routing.SOURCE_LOCAL
    assert table.loc["L", "reference"] == "a"


def test_config_accepts_local_with_or_without_ena(tmp_path: Path) -> None:
    (tmp_path / "runs.txt").write_text("SRR1\n")
    (tmp_path / "local.tsv").write_text("run_id\tfastq_1\treference\tbatch\n")
    cfg = load_config(_cfg(tmp_path, "qc:\n  runs: {local_file: local.tsv}\n"))
    assert cfg.qc.runs.local_file == tmp_path / "local.tsv"
    cfg = load_config(
        _cfg(
            tmp_path,
            "qc:\n  runs: {accessions_file: runs.txt, local_file: local.tsv}\n",
        )
    )
    assert (
        cfg.qc.runs.accessions_file is not None and cfg.qc.runs.local_file is not None
    )
    with pytest.raises(ConfigError, match="alternatives"):
        load_config(
            _cfg(tmp_path, "qc:\n  runs: {accessions_file: runs.txt, taxon_id: 1}\n")
        )
    with pytest.raises(ConfigError, match="set accessions_file"):
        load_config(_cfg(tmp_path, "qc:\n  runs: {}\n"))


def _cfg(tmp_path: Path, body: str) -> Path:
    p = tmp_path / "c.yaml"
    p.write_text(body)
    return p
