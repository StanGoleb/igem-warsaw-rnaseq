"""Assemble the deliverable from stages 2-6."""

from __future__ import annotations

import json
import logging
import shutil
import tempfile
from collections.abc import Mapping
from pathlib import Path

import pandas as pd

from igem_warsaw_rnaseq.config import FinalConfig, RefConfig, ReferenceSpec
from igem_warsaw_rnaseq.final import sequences, splits, validate
from igem_warsaw_rnaseq.final.paths import LAYERS, FinalPaths
from igem_warsaw_rnaseq.matrix.paths import MatrixPaths
from igem_warsaw_rnaseq.normalize.paths import NormalizePaths
from igem_warsaw_rnaseq.normalize.run import load_conditions
from igem_warsaw_rnaseq.ortho.orthofinder import load_orthogroups, members
from igem_warsaw_rnaseq.ortho.paths import OrthoPaths
from igem_warsaw_rnaseq.quant.paths import QuantPaths
from igem_warsaw_rnaseq.ref.paths import RefPaths
from igem_warsaw_rnaseq.ref.prepare import source_files
from igem_warsaw_rnaseq.ref.promoters import cds_spans, representative_transcripts
from igem_warsaw_rnaseq.ref.sequences import plain_gff, read_fasta

log = logging.getLogger(__name__)

LAYER_SOURCE = {"rank": "rank_pct", "anchor": "anchor", "combat": "combat"}


def gene_orthogroups(orthogroups: pd.DataFrame) -> pd.Series:
    """(reference, gene_id) -> orthogroup for every member gene, not just representatives."""
    pairs = [
        ((ref, gene), og)
        for og, row in orthogroups.iterrows()
        for ref in orthogroups.columns
        for gene in members(row[ref])
    ]
    index = pd.MultiIndex.from_tuples(
        [p for p, _ in pairs], names=["reference", "gene_id"]
    )
    return pd.Series([og for _, og in pairs], index=index, dtype=str)


def gene_index(references: Mapping[str, RefPaths]) -> pd.DataFrame:
    """Deterministic dataset-wide key ``<reference>_G<n>`` for every gene of every reference."""
    rows = []
    for ref in sorted(references):
        t2g = pd.read_csv(references[ref].t2g, sep="\t", dtype=str)
        for i, gene in enumerate(sorted(set(t2g["gene"])), start=1):
            rows.append({"ID": f"{ref}_G{i:06d}", "reference": ref, "gene_id": gene})
    return pd.DataFrame(rows, columns=["ID", "reference", "gene_id"])


def wide_layers(long: pd.DataFrame, key: list[str]) -> dict[str, pd.DataFrame]:
    out = {}
    for layer, source in LAYER_SOURCE.items():
        # unstack keeps all-NA rows (pivot_table would drop them), so every
        # layer has the same row set.
        out[layer] = long.set_index([*key, "sample_id"])[source].unstack("sample_id")
        out[layer].columns.name = None
    return out


def build(
    paths: FinalPaths,
    cfg: FinalConfig,
    ref_cfg: RefConfig,
    specs: Mapping[str, ReferenceSpec],
    matrix: MatrixPaths,
    normalized: NormalizePaths,
    ortho: OrthoPaths,
    quant: QuantPaths,
    conditions_file: Path | None,
) -> dict[str, object]:
    for d in (
        paths.all_genes_dir,
        paths.og_dir,
        paths.gene_index.parent,
        paths.sequences.parent,
    ):
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)
    references = {name: RefPaths(ref_cfg.output_dir, name) for name in specs}
    index = gene_index(references)
    index.to_csv(paths.gene_index, sep="\t", index=False)
    key = index.set_index(["reference", "gene_id"])["ID"]

    reps = pd.read_csv(matrix.representative, sep="\t", dtype=str)
    reps.to_csv(paths.representative, sep="\t", index=False)
    gene_og = gene_orthogroups(load_orthogroups(ortho.orthogroups))

    # --- all-genes scope: wide matrices from the long normalised table
    expr = pd.read_parquet(normalized.expr)
    expr["ID"] = key.reindex(
        pd.MultiIndex.from_frame(expr[["reference", "gene_id"]])
    ).to_numpy()
    expr["orthogroup"] = (
        gene_og.reindex(pd.MultiIndex.from_frame(expr[["reference", "gene_id"]]))
        .fillna("")
        .to_numpy()
    )
    row_key = ["ID", "reference", "gene_id", "orthogroup"]
    layers = wide_layers(expr, row_key)
    samples_table = pd.read_csv(matrix.samples, sep="\t", dtype=str)
    if samples_table.empty:
        raise ValueError("no samples: every run was dropped in QC or quantification")
    og_layers = wide_layers(pd.read_parquet(normalized.expr_orthogroup), ["orthogroup"])
    report = []
    for scope, scope_layers in (("all_genes_scope", layers), ("og_scope", og_layers)):
        for layer in LAYERS:
            table = scope_layers[layer]
            table.to_csv(paths.matrix(scope, layer), sep="\t", float_format="%.6g")
            bounds = (0.0, 1.0) if layer == "rank" else (None, None)
            report.append(validate.check_matrix(f"{scope}/{layer}", table, *bounds))

    # --- samples metadata
    samples = samples_table.copy()
    quant_qc = pd.read_csv(quant.qc_table, sep="\t", dtype={"run_id": str})
    per_run = quant_qc.set_index("run_id")
    samples["mapping_rate"] = [
        per_run.loc[r.split(";"), "mapping_rate"].mean() for r in samples["runs"]
    ]
    samples["pct_decoy"] = [
        per_run.loc[r.split(";"), "pct_decoy"].mean() for r in samples["runs"]
    ]
    samples["n_genes_detected"] = [
        int(per_run.loc[r.split(";"), "n_genes_detected"].max())
        for r in samples["runs"]
    ]
    corrected = layers["combat"].notna().any(axis=0)
    samples["batch_corrected"] = (
        samples["sample_id"].map(corrected).fillna(False).astype(bool)
    )
    conditions = load_conditions(conditions_file)
    if conditions is not None:
        samples = samples.merge(conditions.reset_index(), on="sample_id", how="left")
    else:
        samples["condition"] = ""
        samples["description"] = ""
    samples["n_runs"] = samples["n_runs"].astype(int)
    samples.to_csv(paths.samples_metadata, sep="\t", index=False, float_format="%.6g")
    report.append(validate.check("samples_metadata", validate.SAMPLES, samples))

    # --- genes metadata + splits
    summary = pd.read_parquet(normalized.gene_summary)
    genes = index.merge(summary, on=["reference", "gene_id"], how="left")
    genes["orthogroup"] = gene_og.reindex(
        pd.MultiIndex.from_frame(genes[["reference", "gene_id"]])
    ).to_numpy()
    crit = reps.set_index(["reference", "gene_id"])["criterion"]
    genes["representative_criterion"] = crit.reindex(
        pd.MultiIndex.from_frame(genes[["reference", "gene_id"]])
    ).to_numpy()
    if ortho.promoter_conservation.exists():
        cons = pd.read_csv(ortho.promoter_conservation, sep="\t")
        genes = genes.merge(
            cons[
                [
                    "orthogroup",
                    "prot_mean_pid",
                    "prom_mean_pid",
                    "prom_background_pid",
                    "prom_above_background",
                ]
            ],
            on="orthogroup",
            how="left",
        )
    if quant.gene_flags.exists():
        flags = pd.read_csv(quant.gene_flags, sep="\t")
        genes = genes.merge(
            flags[["reference", "gene_id", "median_boot_cv", "crossmap_flag"]],
            on=["reference", "gene_id"],
            how="left",
        )
    groups = sorted(set(genes["orthogroup"].dropna()))
    singletons = [
        f"{r}:{g}"
        for r, g, og in zip(genes["reference"], genes["gene_id"], genes["orthogroup"])
        if pd.isna(og)
    ]
    split_table = splits.assign(groups + singletons, cfg.split)
    split_table.to_csv(paths.split, sep="\t", index=False)
    split_of = split_table.set_index("orthogroup")["split"]
    genes["split"] = [
        split_of[og] if not pd.isna(og) else split_of[f"{r}:{g}"]
        for r, g, og in zip(genes["reference"], genes["gene_id"], genes["orthogroup"])
    ]
    genes.to_csv(paths.genes_metadata, sep="\t", index=False, float_format="%.6g")
    report.append(
        validate.check(
            "genes_metadata", validate.GENES, genes.fillna({"orthogroup": pd.NA})
        )
    )

    # --- sequences
    seq_table = write_sequences(
        references, specs, key, gene_og, ref_cfg.promoter_bp, cfg.flank_bp, paths
    )
    report.append(validate.check("sequences", validate.SEQUENCES, seq_table))

    result = {"tables": report, "passed": all(r["passed"] for r in report)}
    paths.schema_validation.write_text(json.dumps(result, indent=2) + "\n")
    log.info(
        "final dataset: %d samples, %d genes, validation %s",
        len(samples),
        len(genes),
        "passed" if result["passed"] else "FAILED",
    )
    return result


def write_sequences(
    references: Mapping[str, RefPaths],
    specs: Mapping[str, ReferenceSpec],
    key: pd.Series,
    gene_og: pd.Series,
    promoter_bp: int,
    flank_bp: int,
    paths: FinalPaths,
) -> pd.DataFrame:
    rows = []
    with paths.promoters_fasta.open("w") as fasta, paths.promoters_bed.open("w") as bed:
        for ref in sorted(references):
            genome_src, gff_src = source_files(specs[ref], references[ref])
            with tempfile.TemporaryDirectory(dir=paths.root) as tmp:
                spans = cds_spans(plain_gff(gff_src, Path(tmp) / "annotation.gff"))
            contigs = dict(read_fasta(genome_src))
            for gene, transcript in representative_transcripts(
                references[ref].proteins
            ).items():
                span = spans.get(transcript)
                if span is None or span.contig not in contigs:
                    continue
                w = sequences.windows(contigs[span.contig], span, promoter_bp, flank_bp)
                gid = key.get((ref, gene))
                og = gene_og.get((ref, gene), "")
                for group, seq in w.by_group().items():
                    rows.append(
                        (
                            ref,
                            seq.upper(),
                            group,
                            gid,
                            og if isinstance(og, str) else "",
                            gene,
                        )
                    )
                if w.promoter:
                    fasta.write(f">{gid} {ref} {gene}\n{w.promoter.upper()}\n")
                    bed.write(
                        f"{span.contig}\t{w.promoter_bed[0]}\t{w.promoter_bed[1]}\t{gid}\t0\t{span.strand}\n"
                    )
    table = pd.DataFrame(
        rows, columns=["SPECIES", "SEQUENCE", "GROUP", "ID", "OG_ID", "GENE_ID"]
    )
    table.to_csv(paths.sequences, sep="\t", index=False)
    return table
