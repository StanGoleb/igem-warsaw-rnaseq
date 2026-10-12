"""OrthoFinder over all references; keep the small tables, drop the bulk."""

from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
from collections.abc import Mapping
from pathlib import Path

import pandas as pd

from igem_warsaw_rnaseq.ortho.paths import OrthoPaths

log = logging.getLogger(__name__)

QC_COLUMNS: tuple[str, ...] = (
    "orthogroup",
    "n_species",
    "n_genes",
    "is_single_copy",
    "max_copies_in_species",
    "species_specific",
)


def run(
    proteomes: Mapping[str, Path],
    paths: OrthoPaths,
    threads: int,
    gene_tree_method: str,
) -> Path:
    """Run OrthoFinder on ``{reference: proteins.fa}`` in a temporary directory.

    OrthoFinder writes tens of thousands of working files; only the orthogroup
    tables, gene trees and statistics are kept under ``paths.orthofinder_dir``.
    """
    with tempfile.TemporaryDirectory(prefix="igem_warsaw_rnaseq_orthofinder_") as tmp:
        proteome_dir, run_dir = Path(tmp) / "proteomes", Path(tmp) / "run"
        proteome_dir.mkdir()
        for name, fasta in proteomes.items():
            shutil.copyfile(fasta, proteome_dir / f"{name}.fa")
        subprocess.run(
            [
                "orthofinder",
                "-f",
                str(proteome_dir),
                "-o",
                str(run_dir),
                "-t",
                str(threads),
                "-a",
                str(max(1, min(threads, 4))),
                "-M",
                gene_tree_method,
            ],
            check=True,
        )
        results = sorted(run_dir.glob("Results_*"))
        if not results:
            raise FileNotFoundError(f"{run_dir}: no Results_* directory")
        finalize(results[-1], paths)
    return paths.orthofinder_dir


def finalize(results: Path, paths: OrthoPaths) -> None:
    """Copy orthogroup tables, statistics and trees out of a Results_* directory."""
    out = paths.orthofinder_dir
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    for name in (
        "Orthogroups.tsv",
        "Orthogroups.GeneCount.tsv",
        "Orthogroups_SingleCopyOrthologues.txt",
        "Orthogroups_UnassignedGenes.tsv",
    ):
        shutil.copyfile(results / "Orthogroups" / name, out / name)
    species_tree = results / "Species_Tree" / "SpeciesTree_rooted.txt"
    if species_tree.exists():
        shutil.copyfile(species_tree, paths.orthofinder_species_tree)
    stats = results / "Comparative_Genomics_Statistics"
    if stats.is_dir():
        shutil.copytree(stats, out / "Comparative_Genomics_Statistics")
    resolved = results / "Resolved_Gene_Trees" / "Resolved_Gene_Trees.txt"
    if resolved.exists():
        write_gene_trees(parse_resolved_trees(resolved.read_text()), paths.gene_trees)


def parse_resolved_trees(text: str) -> list[tuple[str, str]]:
    """``OG0000001: (newick);`` lines -> (orthogroup, newick) pairs."""
    pairs = []
    for line in text.splitlines():
        og, sep, newick = line.partition(": ")
        if sep and newick.strip():
            pairs.append((og.strip(), newick.strip()))
    return pairs


def write_gene_trees(pairs: list[tuple[str, str]], out: Path) -> None:
    with out.open("w") as fh:
        fh.write("orthogroup\tnewick\n")
        for og, newick in pairs:
            fh.write(f"{og}\t{newick}\n")


def load_orthogroups(path: Path) -> pd.DataFrame:
    """Orthogroups.tsv as orthogroup x species, cells = comma-separated gene ids."""
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    return table.set_index("Orthogroup")


def members(cell: str) -> list[str]:
    return [g.strip() for g in cell.split(",") if g.strip()]


def ortho_qc(gene_counts: pd.DataFrame) -> pd.DataFrame:
    """Per-orthogroup summary from Orthogroups.GeneCount.tsv."""
    counts = gene_counts.drop(columns=["Total"]).astype(int)
    return pd.DataFrame(
        {
            "orthogroup": counts.index,
            "n_species": (counts > 0).sum(axis=1).to_numpy(),
            "n_genes": counts.sum(axis=1).to_numpy(),
            "is_single_copy": (counts == 1).all(axis=1).to_numpy(),
            "max_copies_in_species": counts.max(axis=1).to_numpy(),
            "species_specific": ((counts > 0).sum(axis=1) == 1).to_numpy(),
        },
        columns=list(QC_COLUMNS),
    )


def load_gene_counts(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col="Orthogroup")


def single_copy_ids(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]
