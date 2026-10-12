"""Transcriptome, transcript-to-gene map and one-protein-per-gene proteome (gffread)."""

from __future__ import annotations

import gzip
import shutil
import subprocess
from collections.abc import Iterable, Iterator
from pathlib import Path


def plain_copy(src: Path, dst: Path) -> Path:
    """Write ``src`` to ``dst`` decompressed (gffread cannot read gzip)."""
    opener = gzip.open if src.name.endswith(".gz") else open
    with opener(src, "rb") as fi, dst.open("wb") as fo:
        shutil.copyfileobj(fi, fo)
    return dst


def plain_gff(src: Path, dst: Path) -> Path:
    """Decompress a GFF3 and drop features of unknown strand (``?``).

    NCBI marks trans-spliced genes (e.g. plant mitochondrial nad1) with strand
    ``?``, which gffread rejects outright; such features have no single genomic
    locus and are skipped.
    """
    opener = gzip.open if src.name.endswith(".gz") else open
    with opener(src, "rt") as fi, dst.open("w") as fo:
        for line in fi:
            fields = line.split("\t")
            if len(fields) >= 7 and fields[6] == "?":
                continue
            fo.write(line)
    return dst


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Yield (first header token, sequence) for each record; gzip transparent."""
    name: str | None = None
    chunks: list[str] = []
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def gene_key(locus_tag: str, protein_id: str, gene_id: str) -> str:
    """Gene identifier used across the pipeline.

    ``locus_tag`` (NCBI) is preferred, then JGI ``proteinId``, then the GFF gene
    id stripped of the ``gene-``/``gene_`` prefixes NCBI adds.
    """
    if locus_tag and locus_tag != ".":
        return locus_tag
    if protein_id and protein_id != ".":
        return protein_id
    return gene_id.removeprefix("gene-").removeprefix("gene_")


def parse_t2g_table(lines: Iterable[str]) -> list[tuple[str, str]]:
    """(transcript, gene) pairs from ``gffread --table @id,locus_tag,proteinId,@geneid``."""
    pairs = []
    for line in lines:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 4:
            continue
        pairs.append((fields[0], gene_key(fields[1], fields[2], fields[3])))
    return pairs


def longest_per_gene(
    proteins: Iterable[tuple[str, str]], t2g: dict[str, str]
) -> list[tuple[str, str, str]]:
    """One (gene, transcript, sequence) per gene: the longest protein, first on ties."""
    best: dict[str, tuple[str, str]] = {}
    for transcript, seq in proteins:
        gene = t2g.get(transcript)
        if gene is None:
            continue
        # gffread writes '.' for untranslatable codons and '*' for stops;
        # OrthoFinder rejects '.', so internal ones become X (unknown residue).
        seq = seq.rstrip("*.").replace(".", "X").replace("*", "X")
        if gene not in best or len(seq) > len(best[gene][1]):
            best[gene] = (transcript, seq)
    return [(gene, tid, seq) for gene, (tid, seq) in best.items()]


def build(
    genome: Path, gff: Path, transcripts_out: Path, t2g_out: Path, proteins_out: Path
) -> None:
    """Run gffread on plain-text genome + GFF and write the three outputs."""
    subprocess.run(
        ["gffread", str(gff), "-g", str(genome), "-w", str(transcripts_out)],
        check=True,
    )
    table = subprocess.run(
        ["gffread", str(gff), "--table", "@id,locus_tag,proteinId,@geneid"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    pairs = parse_t2g_table(table.splitlines())
    with t2g_out.open("w") as fh:
        fh.write("transcript\tgene\n")
        for transcript, gene in pairs:
            fh.write(f"{transcript}\t{gene}\n")

    all_proteins = proteins_out.with_suffix(".all.tmp")
    subprocess.run(
        ["gffread", str(gff), "-g", str(genome), "-y", str(all_proteins)],
        check=True,
    )
    chosen = longest_per_gene(read_fasta(all_proteins), dict(pairs))
    all_proteins.unlink()
    with proteins_out.open("w") as fh:
        for gene, transcript, seq in chosen:
            fh.write(f">{gene} transcript={transcript}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")
