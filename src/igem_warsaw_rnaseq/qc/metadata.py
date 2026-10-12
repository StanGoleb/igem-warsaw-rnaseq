"""Run metadata table: ENA lookup and/or local runs, plus metadata-stage drop reasons."""

from __future__ import annotations

from collections.abc import Mapping

import pandas as pd

from igem_warsaw_rnaseq import ena
from igem_warsaw_rnaseq.config import ConfigError, QcConfig, ReferenceSpec
from igem_warsaw_rnaseq.qc.decision import metadata_drop_reason
from igem_warsaw_rnaseq.qc.local import load_local_runs
from igem_warsaw_rnaseq.qc.paths import QcPaths


def build_metadata(
    cfg: QcConfig, references: Mapping[str, ReferenceSpec]
) -> pd.DataFrame:
    frames = []
    if cfg.runs.accessions_file is not None:
        accessions = cfg.runs.accessions_file.read_text().split()
        frames.append(ena.fetch_runs(accessions).assign(source="ena", reference=""))
    elif cfg.runs.taxon_id is not None:
        frames.append(
            ena.fetch_runs_for_taxon(cfg.runs.taxon_id).assign(
                source="ena", reference=""
            )
        )
    if cfg.runs.local_file is not None:
        frames.append(load_local_runs(cfg.runs.local_file, references))
    table = pd.concat(frames, ignore_index=True)
    duplicated = table["run_accession"][table["run_accession"].duplicated()].tolist()
    if duplicated:
        raise ConfigError(f"run ids present more than once: {duplicated}")
    table["drop_reason"] = [
        metadata_drop_reason(
            row, cfg.allowed_library_strategies, cfg.excluded_platforms
        )
        or ""
        for row in table.to_dict("records")
    ]
    paths = QcPaths(cfg.output_dir)
    paths.root.mkdir(parents=True, exist_ok=True)
    table.to_csv(paths.metadata, sep="\t", index=False)
    return table


def load_metadata(paths: QcPaths) -> pd.DataFrame:
    table = pd.read_csv(paths.metadata, sep="\t", dtype=str, keep_default_na=False)
    table["in_ena"] = table["in_ena"].eq("True")
    for column, default in (("source", "ena"), ("reference", "")):
        if column not in table.columns:
            table[column] = default
    return table


def eligible(metadata: pd.DataFrame) -> pd.DataFrame:
    """Runs that passed the metadata gates, in table order."""
    return metadata[metadata["drop_reason"] == ""].reset_index(drop=True)
