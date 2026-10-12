"""Command-line entry point: ``igem-warsaw-rnaseq <stage> <command> --config FILE``."""

from __future__ import annotations

import argparse
import logging
import shutil
import sys
import tempfile
from pathlib import Path

import pandas as pd

from igem_warsaw_rnaseq import conditions, tools
from igem_warsaw_rnaseq.config import (
    Config,
    ConfigError,
    MatrixConfig,
    NormalizeConfig,
    OrthoConfig,
    QcConfig,
    QuantConfig,
    RefConfig,
    load_config,
)
from igem_warsaw_rnaseq.final import build as final_build
from igem_warsaw_rnaseq.final.paths import FinalPaths
from igem_warsaw_rnaseq.matrix import build as matrix_build
from igem_warsaw_rnaseq.matrix.paths import MatrixPaths
from igem_warsaw_rnaseq.normalize import plots as normalize_plots
from igem_warsaw_rnaseq.normalize import run as normalize_run
from igem_warsaw_rnaseq.normalize.paths import NormalizePaths
from igem_warsaw_rnaseq.ortho import conservation as ortho_conservation
from igem_warsaw_rnaseq.ortho import orthofinder
from igem_warsaw_rnaseq.ortho import tree as ortho_tree
from igem_warsaw_rnaseq.ortho.orthofinder import load_orthogroups
from igem_warsaw_rnaseq.ortho.paths import OrthoPaths
from igem_warsaw_rnaseq.qc import collect as qc_collect
from igem_warsaw_rnaseq.qc import decision as qc_decision
from igem_warsaw_rnaseq.qc import metadata as qc_metadata
from igem_warsaw_rnaseq.qc import rrna
from igem_warsaw_rnaseq.qc.paths import QcPaths
from igem_warsaw_rnaseq.qc.runner import process_run
from igem_warsaw_rnaseq.quant import collect as quant_collect
from igem_warsaw_rnaseq.quant import routing
from igem_warsaw_rnaseq.quant.paths import QuantPaths
from igem_warsaw_rnaseq.quant.runner import process_run as quant_process_run
from igem_warsaw_rnaseq.ref import promoters as ref_promoters
from igem_warsaw_rnaseq.ref.paths import RefPaths
from igem_warsaw_rnaseq.ref.prepare import prepare as ref_prepare
from igem_warsaw_rnaseq.ref.prepare import source_files
from igem_warsaw_rnaseq.ref.sequences import plain_gff


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s"
    )
    try:
        cfg = load_config(Path(args.config))
        return args.func(cfg, args)
    except (ConfigError, ValueError) as exc:
        parser.error(str(exc))
    except FileNotFoundError as exc:
        parser.error(f"{exc.filename}: not found (run the preceding stage first?)")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="igem-warsaw-rnaseq")
    stages = parser.add_subparsers(dest="stage", required=True)

    qc = stages.add_parser("qc", help="stage 1: read QC")
    commands = qc.add_subparsers(dest="command", required=True)

    p = commands.add_parser("metadata", help="fetch ENA metadata, apply metadata gates")
    _add_config(p)
    p.set_defaults(func=cmd_qc_metadata)

    p = commands.add_parser("rrna-index", help="build the SortMeRNA index once")
    _add_config(p)
    p.set_defaults(func=cmd_qc_rrna_index)

    p = commands.add_parser("run", help="process runs: download, verify, fastp, rRNA")
    _add_config(p)
    which = p.add_mutually_exclusive_group(required=True)
    which.add_argument("--run", metavar="ACCESSION", help="one run")
    which.add_argument(
        "--task", type=int, metavar="N", help="N-th eligible run (1-based)"
    )
    which.add_argument(
        "--all", action="store_true", help="every eligible run, in order"
    )
    which.add_argument("--list", action="store_true", help="print eligible accessions")
    p.add_argument("--force", action="store_true", help="reprocess cached runs")
    p.set_defaults(func=cmd_qc_run)

    p = commands.add_parser("collect", help="write runs_qc.tsv and the MultiQC report")
    _add_config(p)
    p.add_argument("--no-multiqc", action="store_true")
    p.set_defaults(func=cmd_qc_collect)

    ref = stages.add_parser("ref", help="stage 2: reference preparation")
    commands = ref.add_subparsers(dest="command", required=True)

    p = commands.add_parser(
        "prepare", help="transcriptome, proteome, salmon index, ref_qc.json"
    )
    _add_config(p)
    which = p.add_mutually_exclusive_group(required=True)
    which.add_argument("--name", metavar="REFERENCE", help="one reference")
    which.add_argument("--all", action="store_true", help="every configured reference")
    p.add_argument("--force", action="store_true", help="rebuild prepared references")
    p.set_defaults(func=cmd_ref_prepare)

    p = commands.add_parser("promoters", help="promoter FASTA per prepared reference")
    _add_config(p)
    p.add_argument("--force", action="store_true", help="rewrite existing files")
    p.set_defaults(func=cmd_ref_promoters)

    quant = stages.add_parser("quant", help="stage 3: quantification")
    commands = quant.add_subparsers(dest="command", required=True)

    p = commands.add_parser("route", help="assign QC-kept runs to references")
    _add_config(p)
    p.set_defaults(func=cmd_quant_route)

    p = commands.add_parser("run", help="salmon quant + gene tables")
    _add_config(p)
    which = p.add_mutually_exclusive_group(required=True)
    which.add_argument("--run", metavar="ACCESSION", help="one run")
    which.add_argument(
        "--task", type=int, metavar="N", help="N-th routed run (1-based)"
    )
    which.add_argument("--all", action="store_true", help="every routed run")
    which.add_argument("--list", action="store_true", help="print routed accessions")
    p.add_argument("--force", action="store_true", help="requantify cached runs")
    p.set_defaults(func=cmd_quant_run)

    p = commands.add_parser("collect", help="quant_qc.tsv, gene_flags.tsv, plots")
    _add_config(p)
    p.add_argument("--no-plots", action="store_true")
    p.set_defaults(func=cmd_quant_collect)

    ortho = stages.add_parser("ortho", help="stage 4: orthology and phylogeny")
    commands = ortho.add_subparsers(dest="command", required=True)

    p = commands.add_parser(
        "run", help="OrthoFinder over all references + ortho_qc.tsv"
    )
    _add_config(p)
    p.add_argument("--force", action="store_true", help="rerun OrthoFinder")
    p.set_defaults(func=cmd_ortho_run)

    p = commands.add_parser(
        "tree", help="species tree from single-copy orthologues + patristic distances"
    )
    _add_config(p)
    p.add_argument("--force", action="store_true", help="rebuild the tree")
    p.set_defaults(func=cmd_ortho_tree)

    p = commands.add_parser(
        "conservation", help="protein vs promoter identity per single-copy orthogroup"
    )
    _add_config(p)
    p.add_argument("--force", action="store_true", help="recompute")
    p.set_defaults(func=cmd_ortho_conservation)

    matrix = stages.add_parser("matrix", help="stage 5: expression matrices")
    commands = matrix.add_subparsers(dest="command", required=True)
    p = commands.add_parser(
        "build", help="per-reference and orthogroup matrices, samples, representatives"
    )
    _add_config(p)
    p.set_defaults(func=cmd_matrix_build)

    normalize = stages.add_parser("normalize", help="stage 6: normalisation")
    commands = normalize.add_subparsers(dest="command", required=True)
    p = commands.add_parser(
        "run", help="ranks, anchors, ComBat and gene classes; long-format tables"
    )
    _add_config(p)
    p.set_defaults(func=cmd_normalize_run)

    p = commands.add_parser("plots", help="gene maps per reference and for orthogroups")
    _add_config(p)
    p.set_defaults(func=cmd_normalize_plots)

    p = commands.add_parser("plot-gene", help="map and rank histograms of one gene")
    _add_config(p)
    p.add_argument("--gene", required=True, metavar="GENE_ID")
    p.add_argument(
        "--reference", metavar="NAME", help="needed when the id occurs in several"
    )
    p.set_defaults(func=cmd_normalize_plot_gene)

    conditions = stages.add_parser("conditions", help="condition layer helpers")
    commands = conditions.add_subparsers(dest="command", required=True)
    p = commands.add_parser(
        "export", help="per-sample ENA/BioSample texts + empty conditions template"
    )
    _add_config(p)
    p.set_defaults(func=cmd_conditions_export)

    final = stages.add_parser("final", help="stage 7: final dataset")
    commands = final.add_subparsers(dest="command", required=True)
    p = commands.add_parser(
        "build", help="wide matrices, metadata, sequences, splits, schema validation"
    )
    _add_config(p)
    p.set_defaults(func=cmd_final_build)
    return parser


def _add_config(p: argparse.ArgumentParser) -> None:
    p.add_argument("--config", required=True, help="pipeline YAML")


def _qc(cfg: Config) -> QcConfig:
    if cfg.qc is None:
        raise ConfigError("config has no 'qc' section")
    return cfg.qc


def _ref(cfg: Config) -> RefConfig:
    if cfg.ref is None:
        raise ConfigError("config has no 'ref' section")
    return cfg.ref


def _quant(cfg: Config) -> QuantConfig:
    if cfg.quant is None:
        raise ConfigError("config has no 'quant' section")
    return cfg.quant


def _select(
    table: pd.DataFrame, args: argparse.Namespace, what: str
) -> pd.DataFrame | None:
    """Rows chosen by --run / --task / --all; None after handling --list."""
    if args.list:
        print("\n".join(table["run_id"]))
        return None
    if args.run:
        selected = table[table["run_id"] == args.run]
        if selected.empty:
            raise ConfigError(f"{args.run}: not {what}")
        return selected
    if args.task is not None:
        if not 1 <= args.task <= len(table):
            raise ConfigError(f"--task {args.task}: out of range 1..{len(table)}")
        return table.iloc[[args.task - 1]]
    return table


def cmd_qc_metadata(cfg: Config, args: argparse.Namespace) -> int:
    qc = _qc(cfg)
    table = qc_metadata.build_metadata(qc, cfg.references)
    n_ok = (table["drop_reason"] == "").sum()
    print(f"{len(table)} runs, {n_ok} eligible -> {QcPaths(qc.output_dir).metadata}")
    return 0


def cmd_qc_rrna_index(cfg: Config, args: argparse.Namespace) -> int:
    qc = _qc(cfg)
    tools.require("qc rrna-index")
    if not qc.rrna.enabled:
        print("qc.rrna.enabled is false; nothing to do", file=sys.stderr)
        return 1
    rrna.build_index(
        qc.rrna.references, QcPaths(qc.output_dir).rrna_index_dir, qc.rrna.threads
    )
    return 0


def cmd_qc_run(cfg: Config, args: argparse.Namespace) -> int:
    qc = _qc(cfg)
    needed = ["fastp"]
    if qc.runs.accessions_file is not None or qc.runs.taxon_id is not None:
        needed.append("curl")
    if qc.rrna.enabled:
        needed.append("sortmerna")
    tools.require("qc run", tuple(needed))
    runs = qc_metadata.eligible(qc_metadata.load_metadata(QcPaths(qc.output_dir)))
    selected = _select(
        runs.rename(columns={"run_accession": "run_id"}), args, "an eligible run"
    )
    if selected is None:
        return 0
    for row in selected.rename(columns={"run_id": "run_accession"}).to_dict("records"):
        process_run(row, qc, force=args.force)
    return 0


def cmd_qc_collect(cfg: Config, args: argparse.Namespace) -> int:
    qc = _qc(cfg)
    if not args.no_multiqc:
        tools.require("qc collect")
    table = qc_collect.collect(qc)
    if not args.no_multiqc:
        qc_collect.run_multiqc(qc)
    print(f"{len(table)} runs -> {QcPaths(qc.output_dir).qc_table}")
    return 0


def cmd_ref_prepare(cfg: Config, args: argparse.Namespace) -> int:
    ref = _ref(cfg)
    tools.require("ref prepare", ("busco",) if ref.busco.enabled else ())
    if args.name:
        if args.name not in cfg.references:
            print(f"{args.name}: not in 'references'", file=sys.stderr)
            return 1
        specs = [cfg.references[args.name]]
    else:
        specs = list(cfg.references.values())
    if not specs:
        print("no references configured", file=sys.stderr)
        return 1
    for spec in specs:
        paths = ref_prepare(spec, ref, force=args.force)
        print(f"{spec.name} -> {paths.dir}")
    return 0


def cmd_quant_route(cfg: Config, args: argparse.Namespace) -> int:
    quant, qc = _quant(cfg), _qc(cfg)
    qc_paths = QcPaths(qc.output_dir)
    qc_table = pd.read_csv(
        qc_paths.qc_table, sep="\t", dtype=str, keep_default_na=False
    )
    kept = qc_table[qc_table["decision"] == qc_decision.KEEP][["run_id"]]
    metadata = qc_metadata.load_metadata(qc_paths)
    kept = kept.merge(
        metadata[["run_accession", "tax_lineage", "reference"]].rename(
            columns={"run_accession": "run_id"}
        ),
        on="run_id",
        how="left",
    ).fillna({"tax_lineage": "", "reference": ""})
    overrides = (
        routing.load_overrides(quant.routing_overrides, cfg.references)
        if quant.routing_overrides
        else {}
    )
    table = routing.route(kept, cfg.references, overrides)
    paths = QuantPaths(quant.output_dir)
    paths.root.mkdir(parents=True, exist_ok=True)
    table.to_csv(paths.routing, sep="\t", index=False)
    n_ok = (table["status"] == routing.OK).sum()
    print(f"{len(table)} kept runs, {n_ok} routed -> {paths.routing}")
    return 0


def cmd_quant_run(cfg: Config, args: argparse.Namespace) -> int:
    tools.require("quant run")
    quant, qc, ref = _quant(cfg), _qc(cfg), _ref(cfg)
    routed = routing.routable(
        routing.load_routing(QuantPaths(quant.output_dir).routing)
    )
    selected = _select(routed, args, "a routed run")
    if selected is None:
        return 0
    for row in selected.to_dict("records"):
        quant_process_run(row, quant, qc.output_dir, ref.output_dir, force=args.force)
    return 0


def cmd_quant_collect(cfg: Config, args: argparse.Namespace) -> int:
    quant, qc = _quant(cfg), _qc(cfg)
    table = quant_collect.collect(
        quant, QcPaths(qc.output_dir), make_plots=not args.no_plots
    )
    print(f"{len(table)} runs -> {QuantPaths(quant.output_dir).qc_table}")
    return 0


def _ortho(cfg: Config) -> OrthoConfig:
    if cfg.ortho is None:
        raise ConfigError("config has no 'ortho' section")
    return cfg.ortho


def _orthofinder_matches(paths: OrthoPaths, prepared: dict[str, RefPaths]) -> None:
    columns = set(load_orthogroups(paths.orthogroups).columns)
    if set(prepared) != columns:
        raise ConfigError(
            "references differ from the OrthoFinder run "
            f"({sorted(set(prepared) ^ columns)}); rerun `ortho run --force`"
        )


def _prepared(cfg: Config) -> dict[str, RefPaths]:
    ref = _ref(cfg)
    prepared = {
        name: RefPaths(ref.output_dir, name)
        for name in cfg.references
        if RefPaths(ref.output_dir, name).ref_qc.exists()
    }
    if not prepared:
        raise ConfigError(
            "no prepared references (run `igem-warsaw-rnaseq ref prepare` first)"
        )
    return prepared


def cmd_ref_promoters(cfg: Config, args: argparse.Namespace) -> int:
    tools.require("ref promoters")
    ref = _ref(cfg)
    for name, paths in _prepared(cfg).items():
        out = paths.promoters(ref.promoter_bp)
        if out.exists() and not args.force:
            print(f"{name}: {out.name} present")
            continue
        genome, gff = source_files(cfg.references[name], paths)
        with tempfile.TemporaryDirectory(dir=paths.dir) as tmp:
            n = ref_promoters.write_promoters(
                genome,
                plain_gff(gff, Path(tmp) / "annotation.gff"),
                paths.proteins,
                out,
                ref.promoter_bp,
            )
        print(f"{name}: {n} promoters -> {out}")
    return 0


def cmd_ortho_run(cfg: Config, args: argparse.Namespace) -> int:
    tools.require("ortho run")
    ortho = _ortho(cfg)
    paths = OrthoPaths(ortho.output_dir)
    prepared = _prepared(cfg)
    if paths.orthogroups.exists() and not args.force:
        print(f"{paths.orthogroups} present; use --force to rerun")
        return 0
    paths.root.mkdir(parents=True, exist_ok=True)
    for stale in (
        paths.species_tree,
        paths.patristic,
        paths.promoter_conservation,
        paths.ortho_qc,
    ):
        stale.unlink(missing_ok=True)
    shutil.rmtree(paths.tree_dir, ignore_errors=True)
    orthofinder.run(
        {name: p.proteins for name, p in prepared.items()},
        paths,
        ortho.threads,
        ortho.gene_tree_method,
    )
    qc = orthofinder.ortho_qc(orthofinder.load_gene_counts(paths.gene_counts))
    qc.to_csv(paths.ortho_qc, sep="\t", index=False)
    print(
        f"{len(qc)} orthogroups, {int(qc['is_single_copy'].sum())} single-copy "
        f"-> {paths.orthofinder_dir}"
    )
    return 0


def cmd_ortho_tree(cfg: Config, args: argparse.Namespace) -> int:
    tools.require("ortho tree")
    ortho = _ortho(cfg)
    paths = OrthoPaths(ortho.output_dir)
    if paths.species_tree.exists() and not args.force:
        print(f"{paths.species_tree} present; use --force to rebuild")
        return 0
    _orthofinder_matches(paths, _prepared(cfg))
    n = ortho_tree.build(
        paths.orthogroups,
        paths.single_copy,
        {name: p.proteins for name, p in _prepared(cfg).items()},
        paths.tree_dir,
        paths.species_tree,
        paths.patristic,
        ortho.tree_model,
        ortho.tree_bootstrap,
        ortho.threads,
    )
    print(f"species tree from {n} orthogroups -> {paths.species_tree}")
    return 0


def cmd_ortho_conservation(cfg: Config, args: argparse.Namespace) -> int:
    tools.require("ortho conservation")
    ortho, ref = _ortho(cfg), _ref(cfg)
    paths = OrthoPaths(ortho.output_dir)
    if paths.promoter_conservation.exists() and not args.force:
        print(f"{paths.promoter_conservation} present; use --force to recompute")
        return 0
    prepared = _prepared(cfg)
    _orthofinder_matches(paths, prepared)
    missing = [
        n for n, p in prepared.items() if not p.promoters(ref.promoter_bp).exists()
    ]
    if missing:
        raise ConfigError(
            f"promoters missing for {missing}; run `igem-warsaw-rnaseq ref promoters`"
        )
    table = ortho_conservation.run(
        paths.orthogroups,
        paths.single_copy,
        {name: p.proteins for name, p in prepared.items()},
        {name: p.promoters(ref.promoter_bp) for name, p in prepared.items()},
        paths.root / "_conservation_work",
        paths.promoter_conservation,
        ortho.conservation_seed,
        ortho.threads,
    )
    above = int(table["prom_above_background"].fillna(False).sum())
    print(
        f"{len(table)} orthogroups, {above} with promoter identity above background "
        f"-> {paths.promoter_conservation}"
    )
    return 0


def _matrix(cfg: Config) -> MatrixConfig:
    if cfg.matrix is None:
        raise ConfigError("config has no 'matrix' section")
    return cfg.matrix


def cmd_matrix_build(cfg: Config, args: argparse.Namespace) -> int:
    matrix, qc, quant, ortho, ref = (
        _matrix(cfg),
        _qc(cfg),
        _quant(cfg),
        _ortho(cfg),
        _ref(cfg),
    )
    paths = MatrixPaths(matrix.output_dir)
    table = matrix_build.build(
        paths,
        QcPaths(qc.output_dir),
        QuantPaths(quant.output_dir),
        OrthoPaths(ortho.output_dir),
        ref.output_dir,
    )
    print(f"{len(table)} samples -> {paths.root}")
    return 0


def cmd_final_build(cfg: Config, args: argparse.Namespace) -> int:
    tools.require("final build")
    if cfg.final is None:
        raise ConfigError("config has no 'final' section")
    normalize = _normalize(cfg)
    paths = FinalPaths(cfg.final.output_dir)
    result = final_build.build(
        paths,
        cfg.final,
        _ref(cfg),
        cfg.references,
        MatrixPaths(_matrix(cfg).output_dir),
        NormalizePaths(normalize.output_dir),
        OrthoPaths(_ortho(cfg).output_dir),
        QuantPaths(_quant(cfg).output_dir),
        normalize.conditions_file,
    )
    print(
        f"final dataset -> {paths.root} (schema validation "
        f"{'passed' if result['passed'] else 'FAILED'})"
    )
    return 0 if result["passed"] else 1


def cmd_conditions_export(cfg: Config, args: argparse.Namespace) -> int:
    qc = QcPaths(_qc(cfg).output_dir)
    table = conditions.export(
        qc_metadata.load_metadata(qc),
        qc.root / "conditions_raw.tsv",
        qc.root / "conditions_template.tsv",
    )
    print(f"{len(table)} samples -> {qc.root / 'conditions_raw.tsv'}")
    return 0


def _normalize(cfg: Config) -> NormalizeConfig:
    if cfg.normalize is None:
        raise ConfigError("config has no 'normalize' section")
    return cfg.normalize


def cmd_normalize_run(cfg: Config, args: argparse.Namespace) -> int:
    normalize, matrix = _normalize(cfg), _matrix(cfg)
    paths = NormalizePaths(normalize.output_dir)
    samples = normalize_run.run(
        paths, MatrixPaths(matrix.output_dir), cfg.references, normalize
    )
    print(f"{len(samples)} samples -> {paths.root}")
    return 0


def _conserved_genes(cfg: Config) -> set[tuple[str, str]] | None:
    """(reference, gene) pairs whose orthogroup promoter identity is above background."""
    if cfg.ortho is None:
        return None
    ortho = OrthoPaths(cfg.ortho.output_dir)
    if not ortho.promoter_conservation.exists():
        return None
    cons = pd.read_csv(ortho.promoter_conservation, sep="\t")
    above = set(cons.loc[cons["prom_above_background"] == True, "orthogroup"])
    gene_og = final_build.gene_orthogroups(load_orthogroups(ortho.orthogroups))
    return {key for key, og in gene_og.items() if og in above}


def cmd_normalize_plots(cfg: Config, args: argparse.Namespace) -> int:
    paths = NormalizePaths(_normalize(cfg).output_dir)
    paths.plots_dir.mkdir(parents=True, exist_ok=True)
    conserved = _conserved_genes(cfg)
    summary = pd.read_parquet(paths.gene_summary)
    samples = pd.read_parquet(paths.expr, columns=["reference", "sample_id"])
    for ref, table in summary.groupby("reference"):
        table = table.set_index("gene_id")
        n = samples.loc[samples["reference"] == ref, "sample_id"].nunique()
        normalize_plots.gene_map(
            table, paths.gene_map(str(ref)), f"{ref}: gene map", "gene", n
        )
        if conserved is not None:
            keep = [g for g in table.index if (str(ref), g) in conserved]
            normalize_plots.gene_map(
                table.loc[keep],
                paths.gene_map(str(ref), conserved=True),
                f"{ref}: genes with conserved promoters",
                "gene",
                n,
            )
    og = pd.read_parquet(paths.orthogroup_summary).set_index("orthogroup")
    n_all = pd.read_parquet(paths.expr_orthogroup, columns=["sample_id"])[
        "sample_id"
    ].nunique()
    normalize_plots.gene_map(
        og, paths.gene_map("orthogroup"), "orthogroups: gene map", "orthogroup", n_all
    )
    if conserved is not None:
        cons = pd.read_csv(
            OrthoPaths(cfg.ortho.output_dir).promoter_conservation, sep="\t"
        )
        keep = og.index.intersection(
            cons.loc[cons["prom_above_background"] == True, "orthogroup"]
        )
        normalize_plots.gene_map(
            og.loc[keep],
            paths.gene_map("orthogroup", conserved=True),
            "orthogroups with conserved promoters",
            "orthogroup",
            n_all,
        )
    print(f"plots -> {paths.plots_dir}")
    return 0


def cmd_normalize_plot_gene(cfg: Config, args: argparse.Namespace) -> int:
    tools.require("ref promoters")
    normalize, ref_cfg = _normalize(cfg), _ref(cfg)
    paths = NormalizePaths(normalize.output_dir)
    summary = pd.read_parquet(paths.gene_summary)
    hits = summary[summary["gene_id"] == args.gene]
    if args.reference:
        hits = hits[hits["reference"] == args.reference]
    if hits.empty:
        raise ConfigError(f"{args.gene}: not in gene_summary.parquet")
    if len(hits) > 1:
        raise ConfigError(
            f"{args.gene}: present in {sorted(hits['reference'])}; pass --reference"
        )
    row = hits.iloc[0]
    reference = str(row["reference"])
    ref_summary = summary[summary["reference"] == reference].set_index("gene_id")
    expr = pd.read_parquet(paths.expr, columns=["gene_id", "reference", "rank_pct"])
    strength = (
        1.0
        - expr.loc[
            (expr["gene_id"] == args.gene) & (expr["reference"] == reference),
            "rank_pct",
        ]
    )

    ref_paths = RefPaths(ref_cfg.output_dir, reference)
    t2g = pd.read_csv(ref_paths.t2g, sep="\t", dtype=str)
    _, gff = source_files(cfg.references[reference], ref_paths)
    with tempfile.TemporaryDirectory() as tmp:
        descriptions = normalize_plots.gene_descriptions(
            plain_gff(gff, Path(tmp) / "annotation.gff"),
            dict(zip(t2g["transcript"], t2g["gene"], strict=True)),
        )

    orthogroup = og_row = og_strength = None
    if cfg.ortho is not None and OrthoPaths(cfg.ortho.output_dir).orthogroups.exists():
        gene_og = final_build.gene_orthogroups(
            load_orthogroups(OrthoPaths(cfg.ortho.output_dir).orthogroups)
        )
        orthogroup = gene_og.get((reference, args.gene))
        if orthogroup is not None:
            og_summary = pd.read_parquet(paths.orthogroup_summary).set_index(
                "orthogroup"
            )
            og_expr = pd.read_parquet(
                paths.expr_orthogroup, columns=["orthogroup", "rank_pct"]
            )
            og_row = og_summary.loc[orthogroup]
            og_strength = (
                1.0 - og_expr.loc[og_expr["orthogroup"] == orthogroup, "rank_pct"]
            )

    paths.plots_dir.mkdir(parents=True, exist_ok=True)
    out = paths.gene_plot(args.gene)
    normalize_plots.gene_figure(
        args.gene,
        reference,
        descriptions.get(args.gene, ""),
        row,
        strength,
        ref_summary,
        orthogroup,
        og_row,
        og_strength,
        out,
    )
    print(f"{args.gene} ({row['class']}) -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
