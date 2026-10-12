"""Output layout of the QC stage."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class QcPaths:
    root: Path

    @property
    def metadata(self) -> Path:
        return self.root / "runs_metadata.tsv"

    @property
    def qc_table(self) -> Path:
        return self.root / "runs_qc.tsv"

    @property
    def multiqc_report(self) -> Path:
        return self.root / "multiqc_report.html"

    @property
    def trimmed_dir(self) -> Path:
        return self.root / "fastq_trimmed"

    @property
    def fastp_dir(self) -> Path:
        return self.root / "fastp"

    @property
    def rrna_dir(self) -> Path:
        return self.root / "rrna"

    @property
    def rrna_index_dir(self) -> Path:
        return self.root / "rrna_index"

    @property
    def results_dir(self) -> Path:
        return self.root / "runs"

    @property
    def work_dir(self) -> Path:
        return self.root / "_work"

    def trimmed(self, run: str, paired: bool) -> list[Path]:
        if paired:
            return [self.trimmed_dir / f"{run}_{i}.fastq.gz" for i in (1, 2)]
        return [self.trimmed_dir / f"{run}.fastq.gz"]

    def fastp_json(self, run: str) -> Path:
        return self.fastp_dir / f"{run}.json"

    def rrna_log_prefix(self, run: str) -> Path:
        return self.rrna_dir / run

    def result(self, run: str) -> Path:
        return self.results_dir / f"{run}.json"

    def work(self, run: str) -> Path:
        return self.work_dir / run
