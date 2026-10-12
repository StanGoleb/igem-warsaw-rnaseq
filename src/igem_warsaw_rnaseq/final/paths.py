"""Output layout of the final dataset (mirrors the layout previously delivered)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

LAYERS: tuple[str, ...] = ("rank", "anchor", "combat")


@dataclass(frozen=True)
class FinalPaths:
    root: Path

    @property
    def all_genes_dir(self) -> Path:
        return self.root / "all_genes_scope"

    @property
    def og_dir(self) -> Path:
        return self.root / "og_scope"

    def matrix(self, scope: str, layer: str) -> Path:
        return self.root / scope / f"{layer}.tsv.gz"

    @property
    def samples_metadata(self) -> Path:
        return self.root / "samples_metadata.tsv"

    @property
    def genes_metadata(self) -> Path:
        return self.root / "genes_metadata.tsv"

    @property
    def gene_index(self) -> Path:
        return self.root / "orthology" / "gene_index.tsv"

    @property
    def representative(self) -> Path:
        return self.root / "orthology" / "orthogroup_representative.tsv"

    @property
    def sequences(self) -> Path:
        return self.root / "sequences" / "sequences.tsv.gz"

    @property
    def promoters_fasta(self) -> Path:
        return self.root / "sequences" / "promoters.fasta"

    @property
    def promoters_bed(self) -> Path:
        return self.root / "sequences" / "promoters.bed"

    @property
    def split(self) -> Path:
        return self.root / "split_assignment.tsv"

    @property
    def schema_validation(self) -> Path:
        return self.root / "schema_validation.json"
