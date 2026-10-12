"""Build per-reference and orthogroup-space matrices."""

from __future__ import annotations

import logging
from collections.abc import Mapping

import pandas as pd

from igem_warsaw_rnaseq.matrix import representative, samples
from igem_warsaw_rnaseq.matrix.paths import MatrixPaths
from igem_warsaw_rnaseq.ortho.orthofinder import load_orthogroups, members
from igem_warsaw_rnaseq.ortho.paths import OrthoPaths
from igem_warsaw_rnaseq.qc.metadata import load_metadata
from igem_warsaw_rnaseq.qc.paths import QcPaths
from igem_warsaw_rnaseq.quant.paths import QuantPaths
from igem_warsaw_rnaseq.quant.routing import load_routing
from igem_warsaw_rnaseq.ref.paths import RefPaths
from igem_warsaw_rnaseq.ref.sequences import read_fasta

log = logging.getLogger(__name__)


def reference_matrices(
    sample_rows: pd.DataFrame, quant: QuantPaths
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(counts, tpm) genes x samples for one reference's samples."""
    counts, tpm = {}, {}
    for sample, runs in zip(sample_rows["sample_id"], sample_rows["runs"], strict=True):
        tables = {
            run: pd.read_csv(
                quant.gene_table(run),
                sep="\t",
                index_col="gene_id",
                dtype={"gene_id": str},
            )
            for run in runs.split(";")
        }
        merged = samples.merge_runs(tables)
        counts[sample], tpm[sample] = merged["counts"], merged["tpm"]
    return (
        pd.DataFrame(counts).fillna(0.0).rename_axis("gene_id"),
        pd.DataFrame(tpm).fillna(0.0).rename_axis("gene_id"),
    )


def representatives(
    orthogroups: pd.DataFrame,
    gene_trees: Mapping[str, str],
    mean_tpm: Mapping[str, pd.Series],
    protein_lengths: Mapping[str, Mapping[str, int]],
) -> pd.DataFrame:
    references = list(orthogroups.columns)
    rows = []
    for og, row in orthogroups.iterrows():
        for ref in references:
            genes = members(row[ref])
            if not genes:
                continue
            gene, criterion = representative.choose(
                ref,
                genes,
                gene_trees.get(og),
                mean_tpm.get(ref),
                protein_lengths[ref],
                references,
            )
            rows.append(
                {
                    "orthogroup": og,
                    "reference": ref,
                    "gene_id": gene,
                    "criterion": criterion,
                    "n_paralogs_collapsed": len(genes) - 1,
                }
            )
    return representative.table(rows)


def orthogroup_matrix(
    reps: pd.DataFrame,
    per_reference: Mapping[str, pd.DataFrame],
    all_orthogroups: list[str],
) -> pd.DataFrame:
    """Orthogroups x samples; NA where the sample's reference has no gene in the OG."""
    blocks = []
    for ref, matrix in per_reference.items():
        rep = reps[reps["reference"] == ref].set_index("orthogroup")["gene_id"]
        block = matrix.reindex(rep.to_numpy())
        block.index = rep.index
        blocks.append(block)
    if not blocks:
        return pd.DataFrame(index=pd.Index(all_orthogroups, name="orthogroup"))
    return pd.concat(blocks, axis=1).reindex(all_orthogroups).rename_axis("orthogroup")


def build(
    paths: MatrixPaths,
    qc: QcPaths,
    quant: QuantPaths,
    ortho: OrthoPaths,
    ref_root,
) -> pd.DataFrame:
    paths.root.mkdir(parents=True, exist_ok=True)
    quant_qc = pd.read_csv(quant.qc_table, sep="\t", dtype=str, keep_default_na=False)
    table = samples.sample_table(
        quant_qc, load_routing(quant.routing), load_metadata(qc)
    )
    if table.empty:
        raise ValueError("no samples: every run was dropped in QC or quantification")
    table.to_csv(paths.samples, sep="\t", index=False)

    orthogroups = load_orthogroups(ortho.orthogroups)
    trees = pd.read_csv(ortho.gene_trees, sep="\t", dtype=str)
    gene_trees = dict(zip(trees["orthogroup"], trees["newick"], strict=True))

    per_counts, per_tpm, mean_tpm = {}, {}, {}
    for ref, rows in table.groupby("reference"):
        counts, tpm = reference_matrices(rows, quant)
        counts.to_parquet(paths.counts(str(ref)))
        tpm.to_parquet(paths.tpm(str(ref)))
        per_counts[str(ref)], per_tpm[str(ref)] = counts, tpm
        mean_tpm[str(ref)] = tpm.mean(axis=1)
        log.info("%s: %d genes x %d samples", ref, len(tpm), tpm.shape[1])

    protein_lengths = {
        ref: {
            gene: len(seq) for gene, seq in read_fasta(RefPaths(ref_root, ref).proteins)
        }
        for ref in orthogroups.columns
    }
    reps = representatives(orthogroups, gene_trees, mean_tpm, protein_lengths)
    reps.to_csv(paths.representative, sep="\t", index=False)

    ogs = list(orthogroups.index)
    orthogroup_matrix(reps, per_counts, ogs).to_parquet(paths.orthogroup_counts)
    og_tpm = orthogroup_matrix(reps, per_tpm, ogs)
    og_tpm.to_parquet(paths.orthogroup_tpm)
    log.info(
        "orthogroup space: %d orthogroups x %d samples; criteria: %s",
        len(ogs),
        og_tpm.shape[1],
        reps["criterion"].value_counts().to_dict(),
    )
    return table
