"""ref_qc.json: annotation/genome consistency and completeness summary."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

from igem_warsaw_rnaseq.ref.busco import BuscoResult
from igem_warsaw_rnaseq.ref.sequences import read_fasta


def gff_seqids(gff: Path) -> set[str]:
    ids = set()
    with gff.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            ids.add(line.split("\t", 1)[0])
    return ids


def contig_check(genome: Path, gff: Path) -> tuple[int, list[str]]:
    """(number of genome records, GFF seqids absent from the genome)."""
    fasta_ids = {name for name, _ in read_fasta(genome)}
    missing = sorted(gff_seqids(gff) - fasta_ids)
    return len(fasta_ids), missing


def count_t2g(t2g: Path) -> tuple[int, int]:
    """(n_transcripts, n_genes) from t2g.tsv."""
    genes = set()
    n = 0
    with t2g.open() as fh:
        next(fh)
        for line in fh:
            n += 1
            genes.add(line.rstrip("\n").split("\t")[1])
    return n, len(genes)


def tool_version(cmd: list[str]) -> str:
    out = subprocess.run(cmd, capture_output=True, text=True, check=False)
    text = (out.stdout or out.stderr).strip()
    return text.splitlines()[0] if text else ""


def build(
    name: str,
    annotation_source: str,
    n_contigs: int,
    missing_seqids: list[str],
    n_transcripts: int,
    n_genes: int,
    n_decoys: int | None,
    busco: BuscoResult | None,
) -> dict[str, object]:
    return {
        "reference": name,
        "annotation_source": annotation_source,
        "n_contigs": n_contigs,
        "n_genes": n_genes,
        "n_transcripts": n_transcripts,
        "contig_names_match": not missing_seqids,
        "gff_seqids_missing_from_genome": missing_seqids,
        "salmon_decoys": n_decoys,
        "busco": busco.to_dict() if busco else None,
        "tool_versions": {
            "gffread": tool_version(["gffread", "--version"]),
            "salmon": tool_version(["salmon", "--version"]),
            "busco": tool_version(["busco", "--version"]) if busco else None,
        },
    }


def write(qc: dict[str, object], path: Path) -> None:
    path.write_text(json.dumps(qc, indent=2) + "\n")
