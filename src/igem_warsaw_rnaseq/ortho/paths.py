"""Output layout of the orthology stage."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class OrthoPaths:
    root: Path

    @property
    def orthofinder_dir(self) -> Path:
        return self.root / "orthofinder"

    @property
    def orthogroups(self) -> Path:
        return self.orthofinder_dir / "Orthogroups.tsv"

    @property
    def gene_counts(self) -> Path:
        return self.orthofinder_dir / "Orthogroups.GeneCount.tsv"

    @property
    def single_copy(self) -> Path:
        return self.orthofinder_dir / "Orthogroups_SingleCopyOrthologues.txt"

    @property
    def gene_trees(self) -> Path:
        return self.orthofinder_dir / "gene_trees.tsv"

    @property
    def orthofinder_species_tree(self) -> Path:
        return self.orthofinder_dir / "SpeciesTree_rooted.nwk"

    @property
    def ortho_qc(self) -> Path:
        return self.root / "ortho_qc.tsv"

    @property
    def tree_dir(self) -> Path:
        return self.root / "species_tree"

    @property
    def species_tree(self) -> Path:
        return self.root / "species_tree.nwk"

    @property
    def patristic(self) -> Path:
        return self.root / "patristic_distances.tsv"

    @property
    def promoter_conservation(self) -> Path:
        return self.root / "promoter_conservation.tsv"
