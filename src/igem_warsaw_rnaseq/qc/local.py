"""Local FASTQ runs, described by a table and mapped onto the ENA metadata schema."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import pandas as pd

from igem_warsaw_rnaseq.config import ConfigError, ReferenceSpec
from igem_warsaw_rnaseq.ena import RUN_FIELDS

REQUIRED: tuple[str, ...] = ("run_id", "fastq_1", "reference", "batch")
OPTIONAL: tuple[str, ...] = ("fastq_2", "sample")


def load_local_runs(
    path: Path, references: Mapping[str, ReferenceSpec]
) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = [c for c in REQUIRED if c not in table.columns]
    if missing:
        raise ConfigError(f"{path}: missing columns {missing}")
    for col in OPTIONAL:
        if col not in table.columns:
            table[col] = ""
    duplicated = table["run_id"][table["run_id"].duplicated()].tolist()
    if duplicated:
        raise ConfigError(f"{path}: duplicated run_id {duplicated}")
    unknown = sorted(set(table["reference"]) - set(references))
    if unknown:
        raise ConfigError(f"{path}: unknown references {unknown}")

    rows = []
    for r in table.to_dict("records"):
        files = [
            str(_resolve(p, path.parent)) for p in (r["fastq_1"], r["fastq_2"]) if p
        ]
        row = dict.fromkeys(RUN_FIELDS, "")
        row.update(
            run_accession=r["run_id"],
            experiment_accession=r["sample"] or r["run_id"],
            study_accession=r["batch"],
            library_layout="PAIRED" if len(files) == 2 else "SINGLE",
            library_strategy="RNA-Seq",
            library_source="TRANSCRIPTOMIC",
            fastq_ftp=";".join(files),
            in_ena=False,
            source="local",
            reference=r["reference"],
        )
        rows.append(row)
    return pd.DataFrame(rows)


def _resolve(value: str, base: Path) -> Path:
    p = Path(value).expanduser()
    return p if p.is_absolute() else (base / p).resolve()
