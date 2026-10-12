"""Four tiled windows per gene: flank_up | promoter | gene | flank_down, 5'->3'."""

from __future__ import annotations

from dataclasses import dataclass

from igem_warsaw_rnaseq.ref.promoters import CdsSpan, reverse_complement

GROUPS: tuple[str, ...] = ("FlankLeft", "Promotor", "Gene", "FlankRight")


@dataclass(frozen=True)
class Windows:
    flank_up: str
    promoter: str
    gene: str
    flank_down: str
    promoter_bed: tuple[int, int]  # 0-based half-open on the forward strand

    def by_group(self) -> dict[str, str]:
        return dict(
            zip(GROUPS, (self.flank_up, self.promoter, self.gene, self.flank_down))
        )


def windows(contig: str, span: CdsSpan, promoter_bp: int, flank_bp: int) -> Windows:
    """Clip at contig ends; minus-strand genes are mirrored and reverse-complemented."""
    n = len(contig)
    if span.strand == "-":
        atg = span.end  # forward coordinate just past the CDS = start-codon side
        gene_end = span.start - 1
        promoter = (atg, min(n, atg + promoter_bp))
        flank_up = (promoter[1], min(n, promoter[1] + flank_bp))
        gene = (gene_end, atg)
        flank_down = (max(0, gene_end - flank_bp), gene_end)
        cut = [
            reverse_complement(contig[a:b])
            for a, b in (flank_up, promoter, gene, flank_down)
        ]
        return Windows(*cut, promoter)
    atg = span.start - 1
    gene_end = span.end
    promoter = (max(0, atg - promoter_bp), atg)
    flank_up = (max(0, promoter[0] - flank_bp), promoter[0])
    gene = (atg, gene_end)
    flank_down = (gene_end, min(n, gene_end + flank_bp))
    cut = [contig[a:b] for a, b in (flank_up, promoter, gene, flank_down)]
    return Windows(*cut, promoter)
