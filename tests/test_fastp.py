import json
from pathlib import Path

from igem_warsaw_rnaseq.qc.fastp import parse_report


def _report(tmp_path: Path, body: dict) -> Path:
    path = tmp_path / "run.json"
    path.write_text(json.dumps(body))
    return path


def test_parse_report_paired_counts_both_mates(tmp_path: Path) -> None:
    body = {
        "summary": {
            "before_filtering": {"total_reads": 140552136, "q30_rate": 0.900488},
            "after_filtering": {
                "total_reads": 139565014,
                "q30_rate": 0.934599,
                "read1_mean_length": 140,
                "read2_mean_length": 130,
            },
        },
        "adapter_cutting": {"adapter_trimmed_reads": 71213403},
    }
    m = parse_report(_report(tmp_path, body))
    assert m.reads_raw == 140552136
    assert m.reads_after_fastp == 139565014
    assert m.pct_q30 == 90.0488
    assert round(m.pct_adapter, 2) == 50.67
    assert m.read_length_after_fastp == 135.0


def test_parse_report_without_adapter_section(tmp_path: Path) -> None:
    body = {
        "summary": {
            "before_filtering": {"total_reads": 1000, "q30_rate": 0.5},
            "after_filtering": {"total_reads": 0, "q30_rate": 0.0},
        }
    }
    m = parse_report(_report(tmp_path, body))
    assert m.pct_adapter == 0.0
    assert m.reads_after_fastp == 0
    assert m.read_length_after_fastp == 0.0


def test_parse_report_zero_reads(tmp_path: Path) -> None:
    body = {
        "summary": {
            "before_filtering": {"total_reads": 0, "q30_rate": 0.0},
            "after_filtering": {"total_reads": 0, "q30_rate": 0.0},
        },
        "adapter_cutting": {"adapter_trimmed_reads": 0},
    }
    assert parse_report(_report(tmp_path, body)).pct_adapter == 0.0
