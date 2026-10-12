"""quant_qc.tsv, gene_flags.tsv and QC plots from per-run results."""

from __future__ import annotations

import logging
from collections.abc import Mapping

import numpy as np
import pandas as pd

from igem_warsaw_rnaseq.config import QuantConfig
from igem_warsaw_rnaseq.qc import decision as qc_decision
from igem_warsaw_rnaseq.qc.metadata import load_metadata
from igem_warsaw_rnaseq.qc.paths import QcPaths
from igem_warsaw_rnaseq.quant import plots
from igem_warsaw_rnaseq.quant.paths import QuantPaths
from igem_warsaw_rnaseq.quant.routing import load_routing, routable
from igem_warsaw_rnaseq.quant.runner import QuantResult

log = logging.getLogger(__name__)

LOW_MAPPING_RATE = "low_mapping_rate"

COLUMNS: tuple[str, ...] = (
    "run_id",
    "reference",
    "bioproject",
    "mapping_rate",
    "pct_decoy",
    "detected_libtype",
    "n_genes_detected",
    "pca1",
    "pca2",
    "median_spearman_to_others",
    "decision",
    "drop_reason",
)

FLAG_COLUMNS: tuple[str, ...] = (
    "reference",
    "gene_id",
    "median_tpm",
    "median_boot_cv",
    "n_runs_detected",
    "crossmap_flag",
)


def build_table(
    routing: pd.DataFrame,
    results: Mapping[str, QuantResult],
    bioprojects: Mapping[str, str],
    min_pct_mapped: float | None,
) -> pd.DataFrame:
    rows = []
    for run, reference in zip(routing["run_id"], routing["reference"], strict=True):
        r = results.get(run)
        rows.append(
            {
                "run_id": run,
                "reference": reference,
                "bioproject": bioprojects.get(run, ""),
                "mapping_rate": r.pct_mapped if r else None,
                "pct_decoy": r.pct_decoy if r else None,
                "detected_libtype": r.detected_libtype if r else None,
                "n_genes_detected": r.n_genes_detected if r else None,
                "pca1": None,
                "pca2": None,
                "median_spearman_to_others": None,
                **_decide(r, min_pct_mapped),
            }
        )
    table = pd.DataFrame(rows, columns=list(COLUMNS))
    return table.astype({"n_genes_detected": "Int64"})


def _decide(
    r: QuantResult | None, min_pct_mapped: float | None
) -> dict[str, str | None]:
    if r is None:
        return {"decision": qc_decision.PENDING, "drop_reason": None}
    if r.drop_reason:
        return {"decision": qc_decision.DROP, "drop_reason": r.drop_reason}
    if min_pct_mapped is not None and r.pct_mapped < min_pct_mapped:
        return {"decision": qc_decision.DROP, "drop_reason": LOW_MAPPING_RATE}
    return {"decision": qc_decision.KEEP, "drop_reason": None}


def log_tpm_matrix(paths: QuantPaths, runs: list[str]) -> pd.DataFrame:
    """Genes x runs matrix of log2(TPM + 1)."""
    columns = {
        run: pd.read_csv(
            paths.gene_table(run), sep="\t", index_col="gene_id", dtype={"gene_id": str}
        )["tpm"]
        for run in runs
    }
    return np.log2(pd.DataFrame(columns).fillna(0.0) + 1.0)


def embed(matrix: pd.DataFrame) -> pd.DataFrame:
    """Per run: PCA scores (PC1, PC2) and median Spearman correlation to other runs."""
    runs = list(matrix.columns)
    out = pd.DataFrame(
        index=runs, columns=["pca1", "pca2", "median_spearman_to_others"], dtype=float
    )
    if len(runs) < 2:
        return out
    x = matrix.loc[matrix.var(axis=1) > 0].T
    centered = x - x.mean(axis=0)
    u, s, _ = np.linalg.svd(centered.to_numpy(), full_matrices=False)
    scores = u * s
    if scores.shape[1] > 0:
        out["pca1"] = scores[:, 0]
    if scores.shape[1] > 1:
        out["pca2"] = scores[:, 1]
    corr = matrix.corr(method="spearman")
    off_diagonal = corr.mask(np.eye(len(corr), dtype=bool))
    out["median_spearman_to_others"] = off_diagonal.median(axis=1)
    return out


def gene_flags(
    paths: QuantPaths, reference: str, runs: list[str], cfg: QuantConfig
) -> pd.DataFrame:
    tables = {
        run: pd.read_csv(
            paths.gene_table(run), sep="\t", index_col="gene_id", dtype={"gene_id": str}
        )
        for run in runs
    }
    tpm = pd.DataFrame({r: t["tpm"] for r, t in tables.items()})
    counts = pd.DataFrame({r: t["counts"] for r, t in tables.items()})
    cv = pd.DataFrame({r: t["boot_cv"] for r, t in tables.items()})
    detected = counts >= cfg.detected_min_reads
    flags = pd.DataFrame(
        {
            "reference": reference,
            "gene_id": tpm.index,
            "median_tpm": tpm.median(axis=1).to_numpy(),
            "median_boot_cv": cv.where(detected).median(axis=1).to_numpy(),
            "n_runs_detected": detected.sum(axis=1).to_numpy(),
        }
    )
    flags["crossmap_flag"] = pd.array([None] * len(flags), dtype="boolean")
    if (
        cfg.crossmap_min_median_tpm is not None
        and cfg.crossmap_min_median_boot_cv is not None
    ):
        flags["crossmap_flag"] = (
            flags["median_tpm"] >= cfg.crossmap_min_median_tpm
        ) & (flags["median_boot_cv"] >= cfg.crossmap_min_median_boot_cv)
    return flags[list(FLAG_COLUMNS)]


def collect(
    cfg: QuantConfig, qc_root: QcPaths, make_plots: bool = True
) -> pd.DataFrame:
    paths = QuantPaths(cfg.output_dir)
    routing = routable(load_routing(paths.routing))
    results = (
        {p.stem: QuantResult.load(p) for p in sorted(paths.results_dir.glob("*.json"))}
        if paths.results_dir.is_dir()
        else {}
    )
    metadata = load_metadata(qc_root)
    bioprojects = dict(zip(metadata["run_accession"], metadata["study_accession"]))
    table = build_table(routing, results, bioprojects, cfg.min_pct_mapped)

    flag_frames = []
    for reference, group in table[table["decision"] == qc_decision.KEEP].groupby(
        "reference"
    ):
        runs = list(group["run_id"])
        matrix = log_tpm_matrix(paths, runs)
        scores = embed(matrix)
        for col in scores.columns:
            table.loc[table["run_id"].isin(runs), col] = table.loc[
                table["run_id"].isin(runs), "run_id"
            ].map(scores[col])
        flag_frames.append(gene_flags(paths, str(reference), runs, cfg))
        if make_plots:
            plots.write(
                paths.plots_dir,
                str(reference),
                table[table["run_id"].isin(runs)],
                matrix,
            )

    table.to_csv(paths.qc_table, sep="\t", index=False, float_format="%.6g")
    flags = (
        pd.concat(flag_frames, ignore_index=True)
        if flag_frames
        else pd.DataFrame(columns=list(FLAG_COLUMNS))
    )
    flags.to_csv(paths.gene_flags, sep="\t", index=False, float_format="%.6g")
    counts = table["decision"].value_counts().to_dict()
    log.info("quant_qc.tsv: %s", ", ".join(f"{k}={v}" for k, v in counts.items()))
    return table
