"""Output layout of the quantification stage."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class QuantPaths:
    root: Path

    @property
    def routing(self) -> Path:
        return self.root / "routing.tsv"

    @property
    def salmon_dir(self) -> Path:
        return self.root / "salmon"

    @property
    def gene_dir(self) -> Path:
        return self.root / "quant_gene"

    @property
    def results_dir(self) -> Path:
        return self.root / "runs"

    @property
    def qc_table(self) -> Path:
        return self.root / "quant_qc.tsv"

    @property
    def gene_flags(self) -> Path:
        return self.root / "gene_flags.tsv"

    @property
    def plots_dir(self) -> Path:
        return self.root / "qc_plots"

    def salmon(self, run: str) -> Path:
        return self.salmon_dir / run

    def gene_table(self, run: str) -> Path:
        return self.gene_dir / f"{run}.tsv"

    def result(self, run: str) -> Path:
        return self.results_dir / f"{run}.json"
