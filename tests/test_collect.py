import pandas as pd

from igem_warsaw_rnaseq.config import QcGates
from igem_warsaw_rnaseq.qc import decision
from igem_warsaw_rnaseq.qc.collect import COLUMNS, build_table
from igem_warsaw_rnaseq.qc.runner import RunResult

GATES_OFF = QcGates(min_pct_q30=None, max_pct_adapter=None)


def _metadata(*rows: dict) -> pd.DataFrame:
    base = {
        "study_accession": "PRJNA1",
        "scientific_name": "Trichoderma atroviride",
        "strain": "P1",
        "library_layout": "PAIRED",
        "library_strategy": "RNA-Seq",
        "library_source": "TRANSCRIPTOMIC",
        "library_selection": "PolyA",
        "instrument_platform": "ILLUMINA",
        "drop_reason": "",
    }
    return pd.DataFrame([{**base, **r} for r in rows])


def test_build_table_states() -> None:
    metadata = _metadata(
        {"run_accession": "KEEP"},
        {"run_accession": "META", "drop_reason": decision.UNSUPPORTED_PLATFORM},
        {"run_accession": "BAD"},
        {"run_accession": "TODO"},
    )
    results = {
        "KEEP": RunResult("KEEP", "SINGLE", True, 100, 90, 95.0, 20.0, 2.5, None),
        "BAD": RunResult("BAD", "PAIRED", False, drop_reason=decision.MD5_MISMATCH),
    }
    table = build_table(metadata, results, GATES_OFF).set_index("run_id")

    assert list(table.reset_index().columns) == list(COLUMNS)
    assert table.loc["KEEP", "decision"] == decision.KEEP
    assert table.loc["KEEP", "layout"] == "SINGLE"
    assert table.loc["KEEP", "pct_rrna"] == 2.5
    assert table.loc["META", "decision"] == decision.DROP
    assert table.loc["META", "drop_reason"] == decision.UNSUPPORTED_PLATFORM
    assert pd.isna(table.loc["META", "reads_raw"])
    assert table.loc["BAD", "drop_reason"] == decision.MD5_MISMATCH
    assert not table.loc["BAD", "md5_ok"]
    assert table.loc["TODO", "decision"] == decision.PENDING
    assert table.loc["TODO", "layout"] == "PAIRED"


def test_gates_apply_only_to_processed_runs() -> None:
    metadata = _metadata({"run_accession": "LOWQ"})
    results = {"LOWQ": RunResult("LOWQ", "SINGLE", True, 100, 90, 70.0, 1.0)}
    gates = QcGates(min_pct_q30=80.0, max_pct_adapter=None)
    table = build_table(metadata, results, gates).set_index("run_id")
    assert table.loc["LOWQ", "decision"] == decision.DROP
    assert table.loc["LOWQ", "drop_reason"] == decision.LOW_Q30
    assert build_table(metadata, results, GATES_OFF).loc[0, "decision"] == decision.KEEP
