import json
from pathlib import Path

import pandas as pd
import pytest

from igem_warsaw_rnaseq.download import select_fastq
from igem_warsaw_rnaseq.final import build as final_build
from igem_warsaw_rnaseq.matrix import samples
from igem_warsaw_rnaseq.qc import decision
from igem_warsaw_rnaseq.qc.metadata import load_metadata
from igem_warsaw_rnaseq.qc.paths import QcPaths
from igem_warsaw_rnaseq.qc.runner import RunResult
from igem_warsaw_rnaseq.quant.runner import QuantResult
from igem_warsaw_rnaseq.ref.promoters import reverse_complement


def test_gene_orthogroups_covers_every_member() -> None:
    og = pd.DataFrame(
        {"A": ["a1, a2", "a3", ""], "B": ["b1", "", "b2"]},
        index=pd.Index(["OG0", "OG1", "OG2"], name="Orthogroup"),
    )
    m = final_build.gene_orthogroups(og)
    assert m[("A", "a2")] == "OG0"
    assert m[("B", "b2")] == "OG2"
    assert ("A", "a4") not in m
    assert len(m) == 5


def test_numeric_gene_ids_stay_strings(tmp_path: Path) -> None:
    t = tmp_path / "r.tsv"
    t.write_text("gene_id\tcounts\ttpm\tboot_cv\n123\t5\t1.0\t0.1\n456\t0\t0.0\t\n")
    table = pd.read_csv(t, sep="\t", index_col="gene_id", dtype={"gene_id": str})
    assert list(table.index) == ["123", "456"]


def test_quant_cache_respects_reference(tmp_path: Path) -> None:
    r = QuantResult("run", "old_ref", pct_mapped=50.0)
    r.save(tmp_path / "run.json")
    cached = QuantResult.load(tmp_path / "run.json")
    assert cached.reference != "new_ref"


def test_qc_transient_failures_are_retried(tmp_path: Path) -> None:
    from igem_warsaw_rnaseq.qc.runner import RETRIED

    assert decision.DOWNLOAD_FAILED in RETRIED
    assert decision.MD5_MISMATCH in RETRIED
    assert decision.FASTP_FAILED in RETRIED
    assert decision.NO_READS_AFTER_FASTP not in RETRIED
    r = RunResult("x", drop_reason=decision.DOWNLOAD_FAILED)
    r.save(tmp_path / "x.json")
    assert json.loads((tmp_path / "x.json").read_text())["drop_reason"] in RETRIED


def test_load_metadata_fills_missing_columns(tmp_path: Path) -> None:
    (tmp_path / "runs_metadata.tsv").write_text(
        "run_accession\tin_ena\tdrop_reason\nR1\tTrue\t\n"
    )
    table = load_metadata(QcPaths(tmp_path))
    assert list(table["source"]) == ["ena"]
    assert list(table["reference"]) == [""]
    assert table["in_ena"].tolist() == [True]


def test_reverse_complement_iupac() -> None:
    assert reverse_complement("ACGTN") == "NACGT"
    assert reverse_complement("RYKMBDHV") == "BDHVKMRY"
    assert reverse_complement("acgtry") == "ryacgt"


def test_select_fastq_mismatched_md5_count() -> None:
    files, paired = select_fastq(["a_1.fastq.gz", "a_2.fastq.gz"], ["m1"], "PAIRED")
    assert files == [] and paired is False


def test_sample_table_rejects_unknown_run() -> None:
    qc = pd.DataFrame({"run_id": ["r1"], "decision": ["keep"]})
    routing = pd.DataFrame({"run_id": ["r1"], "reference": ["A"]})
    meta = pd.DataFrame(
        {
            "run_accession": ["other"],
            "experiment_accession": ["e"],
            "study_accession": ["P"],
        }
    )
    with pytest.raises(ValueError, match="missing from runs_metadata"):
        samples.sample_table(qc, routing, meta)
