"""NCBI Datasets: download an assembly's genome FASTA and GFF by accession."""

from __future__ import annotations

import gzip
import shutil
import subprocess
import zipfile
from pathlib import Path

DOWNLOAD_URL = (
    "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download"
    "?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF"
)


def select_members(names: list[str]) -> tuple[str, str]:
    """Pick the genome FASTA and GFF members of a datasets zip."""
    genome = next(
        (
            n
            for n in names
            if n.endswith("_genomic.fna") and "cds_from_genomic" not in n
        ),
        None,
    )
    gff = next((n for n in names if n.endswith("genomic.gff")), None)
    if genome is None or gff is None:
        raise ValueError(f"datasets zip lacks genome FASTA or GFF: {names}")
    return genome, gff


def fetch_assembly(accession: str, genome_out: Path, gff_out: Path) -> None:
    """Download ``accession`` and write gzipped genome and GFF to the given paths."""
    genome_out.parent.mkdir(parents=True, exist_ok=True)
    zip_path = genome_out.parent / f"{accession}.zip"
    # curl, not urllib: the datasets API streams chunked responses that urllib
    # intermittently truncates (IncompleteRead).
    subprocess.run(
        [
            "curl",
            "-sS",
            "-L",
            "--fail",
            "--retry",
            "4",
            "--retry-delay",
            "3",
            "-o",
            str(zip_path),
            DOWNLOAD_URL.format(acc=accession),
        ],
        check=True,
    )
    try:
        with zipfile.ZipFile(zip_path) as zf:
            genome, gff = select_members(zf.namelist())
            for member, out in ((genome, genome_out), (gff, gff_out)):
                # Write to a temporary name so an interrupted extraction never
                # leaves a truncated file that a rerun would accept as complete.
                tmp = out.with_suffix(out.suffix + ".part")
                with zf.open(member) as src, gzip.open(tmp, "wb") as dst:
                    shutil.copyfileobj(src, dst)
                tmp.replace(out)
    finally:
        zip_path.unlink(missing_ok=True)
