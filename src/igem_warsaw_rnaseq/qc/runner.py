"""Per-run processing: download, verify, trim, screen.  Idempotent per run."""

from __future__ import annotations

import json
import logging
import shutil
import subprocess
from collections.abc import Mapping
from dataclasses import asdict, dataclass
from pathlib import Path

from igem_warsaw_rnaseq.config import QcConfig
from igem_warsaw_rnaseq.download import download, ena_url, md5sum, select_fastq
from igem_warsaw_rnaseq.qc import decision, fastp, rrna
from igem_warsaw_rnaseq.qc.paths import QcPaths

log = logging.getLogger(__name__)

RETRIED = frozenset(
    {decision.DOWNLOAD_FAILED, decision.MD5_MISMATCH, decision.FASTP_FAILED}
)


@dataclass
class RunResult:
    run_id: str
    layout: str | None = None
    md5_ok: bool | None = None
    reads_raw: int | None = None
    reads_after_fastp: int | None = None
    pct_q30: float | None = None
    pct_adapter: float | None = None
    pct_rrna: float | None = None
    read_length_after_fastp: float | None = None
    drop_reason: str | None = None

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(asdict(self), indent=2))

    @classmethod
    def load(cls, path: Path) -> RunResult:
        return cls(**json.loads(path.read_text()))


def process_run(
    row: Mapping[str, str], cfg: QcConfig, force: bool = False
) -> RunResult:
    paths = QcPaths(cfg.output_dir)
    run = row["run_accession"]
    result_path = paths.result(run)
    if result_path.exists() and not force:
        cached = RunResult.load(result_path)
        # Transient failures (download, checksum, fastp) are retried on rerun.
        if cached.drop_reason not in RETRIED:
            log.info("%s: cached", run)
            return cached
    result = _process(row, cfg, paths)
    result.save(result_path)
    log.info("%s: %s", run, result.drop_reason or "ok")
    return result


def _process(row: Mapping[str, str], cfg: QcConfig, paths: QcPaths) -> RunResult:
    run = row["run_accession"]
    local = row.get("source", "ena") == "local"
    if local:
        raw = [Path(p) for p in row["fastq_ftp"].split(";") if p]
        paired = len(raw) == 2
    else:
        files, paired = select_fastq(
            row["fastq_ftp"].split(";"),
            row["fastq_md5"].split(";"),
            row["library_layout"],
        )
        if not files:
            return RunResult(run_id=run, drop_reason=decision.NO_FASTQ)
    result = RunResult(run_id=run, layout="PAIRED" if paired else "SINGLE")

    work = paths.work(run)
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    try:
        if not local:
            raw = []
            for i, (url, md5) in enumerate(files, start=1):
                dst = work / f"{run}_{i}.fastq.gz"
                if not download(ena_url(url), dst):
                    result.drop_reason = decision.DOWNLOAD_FAILED
                    return result
                if md5sum(dst) != md5:
                    result.md5_ok = False
                    result.drop_reason = decision.MD5_MISMATCH
                    return result
                raw.append(dst)
            result.md5_ok = True

        trimmed = paths.trimmed(run, paired)
        paths.trimmed_dir.mkdir(parents=True, exist_ok=True)
        paths.fastp_dir.mkdir(parents=True, exist_ok=True)
        try:
            fastp.run_fastp(
                raw,
                trimmed,
                paths.fastp_json(run),
                work / "fastp.html",
                cfg.fastp.threads,
                cfg.fastp.extra_args,
            )
        except subprocess.CalledProcessError:
            _unlink(trimmed)
            result.drop_reason = decision.FASTP_FAILED
            return result
    finally:
        shutil.rmtree(work, ignore_errors=True)

    metrics = fastp.parse_report(paths.fastp_json(run))
    result.reads_raw = metrics.reads_raw
    result.reads_after_fastp = metrics.reads_after_fastp
    result.pct_q30 = metrics.pct_q30
    result.pct_adapter = metrics.pct_adapter
    result.read_length_after_fastp = metrics.read_length_after_fastp
    if metrics.reads_after_fastp == 0:
        _unlink(trimmed)
        result.drop_reason = decision.NO_READS_AFTER_FASTP
        return result

    if cfg.rrna.enabled:
        log_path = rrna.run_sortmerna(
            trimmed,
            cfg.rrna.references,
            paths.rrna_index_dir,
            paths.rrna_log_prefix(run),
            paths.work(f"{run}_rrna"),
            cfg.rrna.threads,
        )
        result.pct_rrna = rrna.parse_pct_rrna(log_path)
    return result


def _unlink(files: list[Path]) -> None:
    for f in files:
        f.unlink(missing_ok=True)
