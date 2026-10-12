"""Assemble runs_qc.tsv from metadata and per-run results; MultiQC report."""

from __future__ import annotations

import logging
import subprocess

import pandas as pd

from igem_warsaw_rnaseq.config import QcConfig, QcGates
from igem_warsaw_rnaseq.qc import decision
from igem_warsaw_rnaseq.qc.metadata import load_metadata
from igem_warsaw_rnaseq.qc.paths import QcPaths
from igem_warsaw_rnaseq.qc.runner import RunResult

log = logging.getLogger(__name__)

COLUMNS: tuple[str, ...] = (
    "run_id",
    "bioproject",
    "species",
    "strain",
    "layout",
    "library_strategy",
    "library_source",
    "library_selection",
    "instrument_platform",
    "md5_ok",
    "reads_raw",
    "reads_after_fastp",
    "read_length_after_fastp",
    "pct_adapter",
    "pct_q30",
    "pct_rrna",
    "decision",
    "drop_reason",
)


def build_table(
    metadata: pd.DataFrame, results: dict[str, RunResult], gates: QcGates
) -> pd.DataFrame:
    rows = []
    for meta in metadata.to_dict("records"):
        run = meta["run_accession"]
        result = results.get(run)
        rows.append(
            {
                "run_id": run,
                "bioproject": meta["study_accession"],
                "species": meta["scientific_name"],
                "strain": meta["strain"],
                "layout": result.layout if result else meta["library_layout"],
                "library_strategy": meta["library_strategy"],
                "library_source": meta["library_source"],
                "library_selection": meta["library_selection"],
                "instrument_platform": meta["instrument_platform"],
                "md5_ok": result.md5_ok if result else None,
                "reads_raw": result.reads_raw if result else None,
                "reads_after_fastp": result.reads_after_fastp if result else None,
                "read_length_after_fastp": (
                    result.read_length_after_fastp if result else None
                ),
                "pct_adapter": result.pct_adapter if result else None,
                "pct_q30": result.pct_q30 if result else None,
                "pct_rrna": result.pct_rrna if result else None,
                **_decide(meta["drop_reason"], result, gates),
            }
        )
    table = pd.DataFrame(rows, columns=list(COLUMNS))
    return table.astype(
        {"md5_ok": "boolean", "reads_raw": "Int64", "reads_after_fastp": "Int64"}
    )


def _decide(
    metadata_reason: str, result: RunResult | None, gates: QcGates
) -> dict[str, str | None]:
    if metadata_reason:
        return {"decision": decision.DROP, "drop_reason": metadata_reason}
    if result is None:
        return {"decision": decision.PENDING, "drop_reason": None}
    reason = result.drop_reason or decision.gate_drop_reason(
        result.pct_q30, result.pct_adapter, gates
    )
    if reason:
        return {"decision": decision.DROP, "drop_reason": reason}
    return {"decision": decision.KEEP, "drop_reason": None}


def collect(cfg: QcConfig) -> pd.DataFrame:
    paths = QcPaths(cfg.output_dir)
    metadata = load_metadata(paths)
    results = (
        {p.stem: RunResult.load(p) for p in sorted(paths.results_dir.glob("*.json"))}
        if paths.results_dir.is_dir()
        else {}
    )
    table = build_table(metadata, results, cfg.gates)
    table.to_csv(paths.qc_table, sep="\t", index=False)
    counts = table["decision"].value_counts().to_dict()
    log.info("runs_qc.tsv: %s", ", ".join(f"{k}={v}" for k, v in counts.items()))
    return table


def run_multiqc(cfg: QcConfig) -> None:
    paths = QcPaths(cfg.output_dir)
    inputs = [d for d in (paths.fastp_dir, paths.rrna_dir) if d.is_dir()]
    if not inputs:
        log.warning("no fastp/rrna outputs yet; skipping MultiQC")
        return
    subprocess.run(
        [
            "multiqc",
            "--force",
            "--no-data-dir",
            "--outdir",
            str(paths.root),
            "--filename",
            paths.multiqc_report.name,
            *map(str, inputs),
        ],
        check=True,
    )
