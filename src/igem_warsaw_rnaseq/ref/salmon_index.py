"""Decoy-aware salmon index: transcriptome plus the genome flagged as decoy."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from igem_warsaw_rnaseq.ref.sequences import read_fasta


def write_decoys(genome: Path, decoys_out: Path) -> int:
    """Decoy names are the genome record ids; returns their count."""
    names = [name for name, _ in read_fasta(genome)]
    decoys_out.write_text("".join(f"{n}\n" for n in names))
    return len(names)


def write_gentrome(transcripts: Path, genome: Path, gentrome_out: Path) -> None:
    """Transcripts followed by the genome: salmon requires decoys last."""
    with gentrome_out.open("wb") as out:
        for part in (transcripts, genome):
            with part.open("rb") as src:
                shutil.copyfileobj(src, out)


def build(
    transcripts: Path,
    index_dir: Path,
    kmer: int,
    threads: int,
    decoys: Path | None = None,
    keep_duplicates: bool = True,
) -> None:
    """Identical transcripts are kept (not collapsed) so every gene stays quantifiable."""
    if index_dir.exists():
        shutil.rmtree(index_dir)
    cmd = [
        "salmon",
        "index",
        "-t",
        str(transcripts),
        "-i",
        str(index_dir),
        "-k",
        str(kmer),
        "-p",
        str(threads),
    ]
    if decoys is not None:
        cmd += ["-d", str(decoys)]
    if keep_duplicates:
        cmd.append("--keepDuplicates")
    subprocess.run(cmd, check=True)
