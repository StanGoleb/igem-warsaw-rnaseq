"""Output layout of the matrix stage."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class MatrixPaths:
    root: Path

    @property
    def samples(self) -> Path:
        return self.root / "samples.tsv"

    def counts(self, reference: str) -> Path:
        return self.root / f"matrix_counts_raw_{reference}.parquet"

    def tpm(self, reference: str) -> Path:
        return self.root / f"matrix_tpm_{reference}.parquet"

    @property
    def orthogroup_counts(self) -> Path:
        return self.root / "matrix_orthogroup_counts.parquet"

    @property
    def orthogroup_tpm(self) -> Path:
        return self.root / "matrix_orthogroup_tpm.parquet"

    @property
    def representative(self) -> Path:
        return self.root / "orthogroup_representative.tsv"
