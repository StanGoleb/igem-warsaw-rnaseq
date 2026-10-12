"""Output layout of the normalisation stage."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class NormalizePaths:
    root: Path

    @property
    def expr(self) -> Path:
        return self.root / "expr_normalized.parquet"

    @property
    def expr_orthogroup(self) -> Path:
        return self.root / "expr_normalized_orthogroup.parquet"

    @property
    def gene_summary(self) -> Path:
        return self.root / "gene_summary.parquet"

    @property
    def orthogroup_summary(self) -> Path:
        return self.root / "orthogroup_summary.parquet"

    @property
    def anchors(self) -> Path:
        return self.root / "anchors.tsv"

    @property
    def combat_log(self) -> Path:
        return self.root / "combat_log.json"

    @property
    def plots_dir(self) -> Path:
        return self.root / "plots"

    def gene_map(self, space: str, conserved: bool = False) -> Path:
        suffix = "_conserved" if conserved else ""
        return self.plots_dir / f"map_{space}{suffix}.png"

    def gene_plot(self, gene: str) -> Path:
        return self.plots_dir / f"gene_{gene}.png"
