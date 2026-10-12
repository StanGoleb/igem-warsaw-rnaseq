"""Output layout of the reference stage."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class RefPaths:
    root: Path
    name: str

    @property
    def dir(self) -> Path:
        return self.root / self.name

    @property
    def source_genome(self) -> Path:
        return self.dir / "source" / "genome.fna.gz"

    @property
    def source_annotation(self) -> Path:
        return self.dir / "source" / "annotation.gff3.gz"

    @property
    def transcripts(self) -> Path:
        return self.dir / "transcripts.fa"

    @property
    def proteins(self) -> Path:
        return self.dir / "proteins.fa"

    @property
    def t2g(self) -> Path:
        return self.dir / "t2g.tsv"

    @property
    def decoys(self) -> Path:
        return self.dir / "decoys.txt"

    @property
    def gentrome(self) -> Path:
        return self.dir / "gentrome.fa"

    @property
    def salmon_index(self) -> Path:
        return self.dir / "salmon_index"

    @property
    def busco_dir(self) -> Path:
        return self.dir / "busco"

    @property
    def busco_downloads(self) -> Path:
        return self.root / "busco_downloads"

    @property
    def ref_qc(self) -> Path:
        return self.dir / "ref_qc.json"

    def promoters(self, upstream_bp: int) -> Path:
        return self.dir / f"promoters_{upstream_bp}bp.fa"

    @property
    def work(self) -> Path:
        return self.dir / "_work"
