"""Prepare one reference end to end.  Idempotent: skipped when ref_qc.json exists."""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

from igem_warsaw_rnaseq import ncbi
from igem_warsaw_rnaseq.config import RefConfig, ReferenceSpec
from igem_warsaw_rnaseq.ref import busco, qc, salmon_index, sequences
from igem_warsaw_rnaseq.ref.paths import RefPaths

log = logging.getLogger(__name__)


def source_files(spec: ReferenceSpec, paths: RefPaths) -> tuple[Path, Path]:
    """(genome, annotation) a reference was or will be built from."""
    if spec.ncbi_accession:
        return paths.source_genome, paths.source_annotation
    return spec.genome, spec.annotation


def prepare(spec: ReferenceSpec, cfg: RefConfig, force: bool = False) -> RefPaths:
    paths = RefPaths(cfg.output_dir, spec.name)
    if paths.ref_qc.exists() and not force:
        log.info("%s: prepared, skipping", spec.name)
        return paths
    paths.dir.mkdir(parents=True, exist_ok=True)
    for stale in paths.dir.glob("promoters_*bp.fa"):
        stale.unlink()

    if spec.ncbi_accession and not (
        paths.source_genome.exists() and paths.source_annotation.exists()
    ):
        log.info("%s: downloading %s from NCBI", spec.name, spec.ncbi_accession)
        ncbi.fetch_assembly(
            spec.ncbi_accession, paths.source_genome, paths.source_annotation
        )
    genome_src, gff_src = source_files(spec, paths)

    if paths.work.exists():
        shutil.rmtree(paths.work)
    paths.work.mkdir()
    try:
        genome = sequences.plain_copy(genome_src, paths.work / "genome.fa")
        gff = sequences.plain_gff(gff_src, paths.work / "annotation.gff")

        n_contigs, missing = qc.contig_check(genome, gff)
        if missing:
            log.warning(
                "%s: %d GFF seqids absent from the genome (e.g. %s)",
                spec.name,
                len(missing),
                missing[:3],
            )

        log.info("%s: gffread", spec.name)
        sequences.build(genome, gff, paths.transcripts, paths.t2g, paths.proteins)
        n_transcripts, n_genes = qc.count_t2g(paths.t2g)

        n_decoys = None
        index_input = paths.transcripts
        decoys = None
        if cfg.salmon.decoy:
            n_decoys = salmon_index.write_decoys(genome, paths.decoys)
            salmon_index.write_gentrome(paths.transcripts, genome, paths.gentrome)
            index_input, decoys = paths.gentrome, paths.decoys
        log.info("%s: salmon index", spec.name)
        salmon_index.build(
            index_input,
            paths.salmon_index,
            cfg.salmon.kmer,
            cfg.salmon.threads,
            decoys,
            cfg.salmon.keep_duplicates,
        )
    finally:
        shutil.rmtree(paths.work, ignore_errors=True)

    busco_result = None
    if cfg.busco.enabled:
        log.info("%s: busco (%s)", spec.name, cfg.busco.lineage)
        summary = busco.run(
            paths.proteins,
            cfg.busco.lineage,
            paths.busco_dir,
            paths.busco_downloads,
            cfg.busco.threads,
        )
        busco_result = busco.parse_summary(summary, cfg.busco.lineage)

    qc.write(
        qc.build(
            spec.name,
            spec.annotation_source,
            n_contigs,
            missing,
            n_transcripts,
            n_genes,
            n_decoys,
            busco_result,
        ),
        paths.ref_qc,
    )
    log.info("%s: done (%d genes, %d transcripts)", spec.name, n_genes, n_transcripts)
    return paths
