"""Samples (ENA experiments) from QC-kept runs, and lane merging."""

from __future__ import annotations

from collections.abc import Mapping

import pandas as pd

COLUMNS: tuple[str, ...] = ("sample_id", "reference", "bioproject", "runs", "n_runs")


def sample_table(
    quant_qc: pd.DataFrame, routing: pd.DataFrame, metadata: pd.DataFrame
) -> pd.DataFrame:
    """One row per experiment accession; its kept runs share one reference."""
    kept = quant_qc[quant_qc["decision"] == "keep"][["run_id"]]
    kept = kept.merge(routing[["run_id", "reference"]], on="run_id", how="left")
    meta = metadata.rename(columns={"run_accession": "run_id"})[
        ["run_id", "experiment_accession", "study_accession"]
    ]
    kept = kept.merge(meta, on="run_id", how="left")
    unknown = kept["run_id"][kept["experiment_accession"].isna()].tolist()
    if unknown:
        raise ValueError(f"runs missing from runs_metadata.tsv: {unknown}")
    kept["experiment_accession"] = kept["experiment_accession"].where(
        kept["experiment_accession"].astype(bool), kept["run_id"]
    )
    rows = []
    for sample, group in kept.groupby("experiment_accession", sort=True):
        references = set(group["reference"])
        if len(references) != 1:
            raise ValueError(
                f"{sample}: runs routed to several references {references}"
            )
        rows.append(
            {
                "sample_id": sample,
                "reference": references.pop(),
                "bioproject": group["study_accession"].iloc[0],
                "runs": ";".join(sorted(group["run_id"])),
                "n_runs": len(group),
            }
        )
    return pd.DataFrame(rows, columns=list(COLUMNS))


def merge_runs(tables: Mapping[str, pd.DataFrame]) -> pd.DataFrame:
    """Combine runs of one library: counts summed, TPM averaged by read weight."""
    if len(tables) == 1:
        return next(iter(tables.values()))[["counts", "tpm"]]
    counts = pd.concat([t["counts"] for t in tables.values()], axis=1).fillna(0.0)
    tpm = pd.concat([t["tpm"] for t in tables.values()], axis=1).fillna(0.0)
    weights = counts.sum(axis=0)
    if weights.sum() == 0:
        weights = pd.Series(1.0, index=weights.index)
    merged_tpm = (tpm * weights.to_numpy()).sum(axis=1) / weights.sum()
    return pd.DataFrame({"counts": counts.sum(axis=1), "tpm": merged_tpm})
