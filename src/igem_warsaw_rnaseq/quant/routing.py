"""Assign each QC-kept run to a reference: overrides first, then ENA taxonomy."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import pandas as pd

from igem_warsaw_rnaseq.config import ConfigError, ReferenceSpec

OK = "ok"
NO_REFERENCE = "no_reference"
SOURCE_OVERRIDE = "override"
SOURCE_TAXONOMY = "taxonomy"
SOURCE_LOCAL = "local"

COLUMNS: tuple[str, ...] = ("run_id", "reference", "routing_source", "status")


def most_specific_reference(
    tax_lineage: str, refs_by_taxid: Mapping[int, str]
) -> str | None:
    """The reference whose tax_id is deepest in the run's root-to-leaf lineage."""
    for token in reversed(tax_lineage.split(";")):
        token = token.strip()
        if token.isdigit() and int(token) in refs_by_taxid:
            return refs_by_taxid[int(token)]
    return None


def load_overrides(
    path: Path, references: Mapping[str, ReferenceSpec]
) -> dict[str, str]:
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for col in ("run_id", "reference"):
        if col not in table.columns:
            raise ConfigError(f"{path}: missing column '{col}'")
    unknown = sorted(set(table["reference"]) - set(references))
    if unknown:
        raise ConfigError(f"{path}: unknown references {unknown}")
    return dict(zip(table["run_id"], table["reference"], strict=True))


def route(
    kept: pd.DataFrame,
    references: Mapping[str, ReferenceSpec],
    overrides: Mapping[str, str],
) -> pd.DataFrame:
    """``kept`` needs columns run_id, tax_lineage and reference (QC-kept runs only).

    Precedence: overrides file, the run's own ``reference`` (local runs), ENA taxonomy.
    """
    refs_by_taxid = {
        spec.tax_id: name for name, spec in references.items() if spec.tax_id
    }
    rows = []
    for run, lineage, declared in zip(
        kept["run_id"], kept["tax_lineage"], kept["reference"], strict=True
    ):
        if run in overrides:
            rows.append((run, overrides[run], SOURCE_OVERRIDE, OK))
            continue
        if declared:
            rows.append((run, declared, SOURCE_LOCAL, OK))
            continue
        ref = most_specific_reference(lineage, refs_by_taxid)
        if ref is None:
            rows.append((run, "", "", NO_REFERENCE))
        else:
            rows.append((run, ref, SOURCE_TAXONOMY, OK))
    return pd.DataFrame(rows, columns=list(COLUMNS))


def load_routing(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def routable(routing: pd.DataFrame) -> pd.DataFrame:
    return routing[routing["status"] == OK].reset_index(drop=True)
