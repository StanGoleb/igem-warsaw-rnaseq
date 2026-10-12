"""Normalise every reference's matrix and the orthogroup matrix; write long tables."""

from __future__ import annotations

import json
import logging
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from igem_warsaw_rnaseq.config import NormalizeConfig, ReferenceSpec
from igem_warsaw_rnaseq.matrix.paths import MatrixPaths
from igem_warsaw_rnaseq.normalize import anchors, combat, rank
from igem_warsaw_rnaseq.normalize.paths import NormalizePaths

log = logging.getLogger(__name__)

SUMMARY_COLUMNS: tuple[str, ...] = (
    "n_samples",
    "median_pct",
    "p05_pct",
    "p95_pct",
    "band",
    "rank_iqr",
    "rank_max",
    "rank_max_category",
    "cons_frac",
    "induction_amplitude",
    "class",
    "is_strong",
)


def load_conditions(path) -> pd.DataFrame | None:
    if path is None:
        return None
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for col in ("sample_id", "condition"):
        if col not in table.columns:
            raise ValueError(f"{path}: missing column '{col}'")
    return table.drop_duplicates("sample_id").set_index("sample_id")


def normalise_space(
    tpm: pd.DataFrame,
    configured_anchors: Sequence[str],
    batch: pd.Series,
    condition: pd.Series | None,
    cfg: NormalizeConfig,
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame, dict[str, object]]:
    """All normalised layers for one matrix (genes or orthogroups x samples)."""
    log2 = np.log2(tpm + 1.0)
    chosen, anchor_table = anchors.select(log2, configured_anchors, cfg.anchors)
    layers = {
        "tpm": tpm,
        "log2tpm": log2,
        "rank_pct": rank.rank_pct(tpm),
        "anchor": anchors.normalise(log2, chosen),
    }
    if cfg.combat.enabled:
        layers["combat"], entry = combat.run(
            log2, batch, condition, cfg.combat.min_batch_size
        )
    else:
        layers["combat"] = pd.DataFrame(np.nan, index=tpm.index, columns=tpm.columns)
        entry = {"applied": False, "reason_skipped": "disabled in config"}
    return layers, anchor_table, entry


def summarise(
    tpm: pd.DataFrame, condition: pd.Series | None, cfg: NormalizeConfig
) -> pd.DataFrame:
    m = rank.classify(rank.metrics(tpm, cfg.classes), cfg.classes)
    best = rank.best_sample(tpm)
    m["rank_max_category"] = best.map(condition) if condition is not None else pd.NA
    return m[list(SUMMARY_COLUMNS)]


def long_table(layers: Mapping[str, pd.DataFrame], id_name: str) -> pd.DataFrame:
    stacked = {name: layer.stack(future_stack=True) for name, layer in layers.items()}
    table = pd.DataFrame(stacked)
    table.index = table.index.set_names([id_name, "sample_id"])
    return table.reset_index()


def run(
    paths: NormalizePaths,
    matrices: MatrixPaths,
    references: Mapping[str, ReferenceSpec],
    cfg: NormalizeConfig,
) -> pd.DataFrame:
    paths.root.mkdir(parents=True, exist_ok=True)
    samples = pd.read_csv(matrices.samples, sep="\t", dtype=str).set_index("sample_id")
    if samples.empty:
        raise ValueError("no samples: every run was dropped in QC or quantification")
    conditions = load_conditions(cfg.conditions_file)
    condition = (
        conditions["condition"].reindex(samples.index)
        if conditions is not None
        else None
    )

    gene_tables, gene_summaries, anchor_rows, combat_log = [], [], [], {}
    for ref, rows in samples.groupby("reference"):
        tpm = pd.read_parquet(matrices.tpm(str(ref)))[rows.index]
        batch = rows["bioproject"]
        layers, anchor_table, entry = normalise_space(
            tpm,
            references[str(ref)].anchor_genes,
            batch,
            condition.reindex(rows.index) if condition is not None else None,
            cfg,
        )
        table = long_table(layers, "gene_id")
        table.insert(2, "reference", str(ref))
        gene_tables.append(table)
        summary = summarise(tpm, condition, cfg)
        summary.insert(0, "reference", str(ref))
        gene_summaries.append(summary.rename_axis("gene_id").reset_index())
        anchor_rows.append(anchor_table.assign(space=str(ref)))
        combat_log[str(ref)] = entry
        log.info(
            "%s: %d genes x %d samples; anchors %s; combat %s",
            ref,
            len(tpm),
            tpm.shape[1],
            list(anchor_table["gene"]),
            entry["applied"],
        )

    pd.concat(gene_tables, ignore_index=True).to_parquet(paths.expr)
    pd.concat(gene_summaries, ignore_index=True).to_parquet(paths.gene_summary)

    og_tpm = pd.read_parquet(matrices.orthogroup_tpm)[samples.index]
    og_anchors = _orthogroup_anchors(matrices, references)
    layers, anchor_table, entry = normalise_space(
        og_tpm, og_anchors, samples["bioproject"], condition, cfg
    )
    long_table(layers, "orthogroup").to_parquet(paths.expr_orthogroup)
    summarise(og_tpm, condition, cfg).rename_axis(
        "orthogroup"
    ).reset_index().to_parquet(paths.orthogroup_summary)
    anchor_rows.append(anchor_table.assign(space="orthogroup"))
    combat_log["orthogroup"] = entry
    log.info(
        "orthogroup space: %d x %d; anchors %s; combat %s",
        len(og_tpm),
        og_tpm.shape[1],
        list(anchor_table["gene"]),
        entry["applied"],
    )

    pd.concat(anchor_rows, ignore_index=True)[
        ["space", "gene", "source", "genorm_M"]
    ].to_csv(paths.anchors, sep="\t", index=False, float_format="%.4g")
    paths.combat_log.write_text(json.dumps(combat_log, indent=2) + "\n")
    return samples


def _orthogroup_anchors(
    matrices: MatrixPaths, references: Mapping[str, ReferenceSpec]
) -> list[str]:
    """Orthogroups of configured anchor genes (any reference); empty -> automatic."""
    configured = {g for spec in references.values() for g in spec.anchor_genes}
    if not configured:
        return []
    reps = pd.read_csv(matrices.representative, sep="\t", dtype=str)
    hits = reps[reps["gene_id"].isin(configured)]
    return sorted(set(hits["orthogroup"]))
