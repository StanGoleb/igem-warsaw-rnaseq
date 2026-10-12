"""Promoters: the bases immediately 5' of each gene's start codon.

Anchored on the ATG because it is the one landmark every annotation source
provides for every gene; transcript starts (UTRs) are annotated inconsistently.
The transcript is the one chosen for proteins.fa, so promoters, proteins and
orthogroups share gene keys 1:1.
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path

from igem_warsaw_rnaseq.ref.sequences import read_fasta

_COMPLEMENT = str.maketrans(
    "ACGTURYKMSWBDHVNacgturykmswbdhvn", "TGCAAYRMKSWVHDBNtgcaayrmkswvhdbn"
)


@dataclass(frozen=True)
class CdsSpan:
    contig: str
    strand: str
    start: int  # 1-based, first CDS base
    end: int  # 1-based, last CDS base


def parse_cds_table(lines: list[str]) -> dict[str, CdsSpan]:
    """``gffread --table @id,@chr,@strand,@cds`` rows -> transcript -> CDS extent."""
    spans: dict[str, CdsSpan] = {}
    for line in lines:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 4 or fields[3] in ("", "."):
            continue
        segments = [seg.split("-") for seg in fields[3].split(",")]
        starts = [int(s[0]) for s in segments]
        ends = [int(s[1]) for s in segments]
        spans[fields[0]] = CdsSpan(fields[1], fields[2], min(starts), max(ends))
    return spans


def cds_spans(gff: Path) -> dict[str, CdsSpan]:
    table = subprocess.run(
        ["gffread", str(gff), "--table", "@id,@chr,@strand,@cds"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    return parse_cds_table(table.splitlines())


def representative_transcripts(proteins_fa: Path) -> dict[str, str]:
    """gene -> transcript from ``>gene transcript=<id>`` headers of proteins.fa."""
    mapping = {}
    with proteins_fa.open() as fh:
        for line in fh:
            if line.startswith(">"):
                gene, _, rest = line[1:].partition(" ")
                mapping[gene] = rest.strip().removeprefix("transcript=")
    return mapping


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def promoter_sequence(contig_seq: str, span: CdsSpan, upstream_bp: int) -> str:
    """Up to ``upstream_bp`` bases 5' of the start codon, on the gene's strand."""
    if span.strand == "-":
        return reverse_complement(contig_seq[span.end : span.end + upstream_bp])
    return contig_seq[max(0, span.start - 1 - upstream_bp) : span.start - 1]


def write_promoters(
    genome: Path, gff: Path, proteins_fa: Path, out: Path, upstream_bp: int
) -> int:
    """Write ``>gene transcript=<id>`` promoter records; returns the record count."""
    spans = cds_spans(gff)
    contigs = dict(read_fasta(genome))
    n = 0
    with out.open("w") as fh:
        for gene, transcript in representative_transcripts(proteins_fa).items():
            span = spans.get(transcript)
            if span is None or span.contig not in contigs:
                continue
            seq = promoter_sequence(contigs[span.contig], span, upstream_bp).upper()
            if not seq:
                continue
            fh.write(f">{gene} transcript={transcript}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")
            n += 1
    return n
