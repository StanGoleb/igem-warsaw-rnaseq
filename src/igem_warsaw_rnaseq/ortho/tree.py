"""Species tree from a supermatrix of single-copy orthologues (MAFFT, trimAl, IQ-TREE)."""

from __future__ import annotations

import logging
import shutil
import subprocess
from collections.abc import Mapping
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from igem_warsaw_rnaseq.ortho import alignment, newick
from igem_warsaw_rnaseq.ortho.orthofinder import (
    load_orthogroups,
    members,
    single_copy_ids,
)
from igem_warsaw_rnaseq.ref.sequences import read_fasta

log = logging.getLogger(__name__)


def single_copy_sequences(
    orthogroups_tsv: Path, single_copy_txt: Path, proteomes: Mapping[str, Path]
) -> dict[str, dict[str, str]]:
    """orthogroup -> {reference: protein sequence} for single-copy orthogroups."""
    table = load_orthogroups(orthogroups_tsv)
    wanted = single_copy_ids(single_copy_txt)
    sequences = {name: dict(read_fasta(path)) for name, path in proteomes.items()}
    out: dict[str, dict[str, str]] = {}
    for og in wanted:
        row = table.loc[og]
        out[og] = {
            ref: sequences[ref][members(row[ref])[0]] for ref in proteomes if row[ref]
        }
    return out


def align_trim(args: tuple[str, dict[str, str], Path]) -> tuple[str, dict[str, str]]:
    og, records, work = args
    work = work / og
    aligned = alignment.mafft(records, work)
    alignment.write_fasta(aligned, work / "aln.fa")
    subprocess.run(
        [
            "trimal",
            "-in",
            str(work / "aln.fa"),
            "-out",
            str(work / "trim.fa"),
            "-automated1",
        ],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    trimmed = alignment.read_aligned(work / "trim.fa")
    shutil.rmtree(work, ignore_errors=True)
    return og, trimmed


def concatenate(
    alignments: Mapping[str, Mapping[str, str]], references: list[str]
) -> tuple[dict[str, str], list[tuple[str, int, int]]]:
    """Supermatrix per reference plus (orthogroup, start, end) partitions (1-based)."""
    matrix = {ref: "" for ref in references}
    partitions = []
    position = 1
    for og in sorted(alignments):
        block = alignments[og]
        width = len(next(iter(block.values())))
        if width == 0:
            continue
        for ref in references:
            matrix[ref] += block.get(ref, "-" * width)
        partitions.append((og, position, position + width - 1))
        position += width
    return matrix, partitions


def iqtree(
    supermatrix: Path,
    prefix: Path,
    model: str,
    bootstrap: int,
    threads: int,
    n_sequences: int,
) -> Path:
    cmd = [
        "iqtree",
        "-s",
        str(supermatrix),
        "-m",
        model,
        "-T",
        str(threads),
        "--prefix",
        str(prefix),
        "-quiet",
        "-redo",
    ]
    if bootstrap and n_sequences >= 4:
        cmd += ["-B", str(max(bootstrap, 1000))]
    subprocess.run(cmd, check=True)
    return prefix.with_suffix(".treefile")


def two_leaf_tree(matrix: Mapping[str, str]) -> str:
    """IQ-TREE needs three sequences; two references get the p-distance split evenly."""
    (a, seq_a), (b, seq_b) = matrix.items()
    identity = alignment.pair_identity(seq_a, seq_b)
    distance = 1.0 - identity / 100.0 if identity is not None else 0.0
    return f"({a}:{distance / 2:.6g},{b}:{distance / 2:.6g});"


def build(
    orthogroups_tsv: Path,
    single_copy_txt: Path,
    proteomes: Mapping[str, Path],
    work: Path,
    tree_out: Path,
    distances_out: Path,
    model: str,
    bootstrap: int,
    threads: int,
) -> int:
    """Full pipeline; returns the number of orthogroups in the supermatrix."""
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    sequences = single_copy_sequences(orthogroups_tsv, single_copy_txt, proteomes)
    log.info("aligning %d single-copy orthogroups", len(sequences))
    with ProcessPoolExecutor(max_workers=threads) as pool:
        trimmed = dict(
            pool.map(align_trim, [(og, recs, work) for og, recs in sequences.items()])
        )
    matrix, partitions = concatenate(trimmed, list(proteomes))
    supermatrix = work / "supermatrix.fa"
    alignment.write_fasta(matrix, supermatrix)
    with (work / "partitions.tsv").open("w") as fh:
        fh.write("orthogroup\tstart\tend\n")
        for og, start, end in partitions:
            fh.write(f"{og}\t{start}\t{end}\n")
    log.info(
        "supermatrix: %d references x %d columns",
        len(matrix),
        len(next(iter(matrix.values()))),
    )
    if len(matrix) < 2:
        raise ValueError("a species tree needs at least two references")
    if len(matrix) == 2:
        tree_out.write_text(two_leaf_tree(matrix) + "\n")
    else:
        treefile = iqtree(
            supermatrix, work / "iqtree", model, bootstrap, threads, len(matrix)
        )
        shutil.copyfile(treefile, tree_out)
    tree = newick.parse(tree_out.read_text())
    newick.patristic_distances(tree).to_csv(
        distances_out, sep="\t", float_format="%.6g"
    )
    return len(partitions)
