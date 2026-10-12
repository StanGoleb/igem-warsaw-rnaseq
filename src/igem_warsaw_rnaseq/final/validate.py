"""pandera schemas for the deliverable; results collected into schema_validation.json."""

from __future__ import annotations

import pandas as pd
import pandera.pandas as pa

from igem_warsaw_rnaseq.final.sequences import GROUPS
from igem_warsaw_rnaseq.final.splits import SPLITS

SAMPLES = pa.DataFrameSchema(
    {
        "sample_id": pa.Column(str, unique=True),
        "reference": pa.Column(str),
        "bioproject": pa.Column(str, nullable=True),
        "n_runs": pa.Column(int, pa.Check.ge(1)),
        "mapping_rate": pa.Column(float, pa.Check.in_range(0, 100), nullable=True),
    },
    strict=False,
)

GENES = pa.DataFrameSchema(
    {
        "ID": pa.Column(str, unique=True),
        "reference": pa.Column(str),
        "gene_id": pa.Column(str),
        "orthogroup": pa.Column(str, nullable=True),
        "class": pa.Column(
            str,
            pa.Check.isin(["constitutive", "silent", "inducible", "variable"]),
            nullable=True,
        ),
        "split": pa.Column(str, pa.Check.isin(list(SPLITS))),
    },
    strict=False,
)

SEQUENCES = pa.DataFrameSchema(
    {
        "SPECIES": pa.Column(str),
        "SEQUENCE": pa.Column(str, pa.Check.str_matches(r"^[ACGTNRYKMSWBDHV]*$")),
        "GROUP": pa.Column(str, pa.Check.isin(list(GROUPS))),
        "ID": pa.Column(str),
        "OG_ID": pa.Column(str, nullable=True),
        "GENE_ID": pa.Column(str),
    },
    strict=False,
)


def check(
    name: str, schema: pa.DataFrameSchema, table: pd.DataFrame
) -> dict[str, object]:
    try:
        schema.validate(table, lazy=True)
        return {"table": name, "passed": True, "n_rows": len(table), "errors": []}
    except pa.errors.SchemaErrors as exc:
        failures = exc.failure_cases
        errors = (
            failures[["column", "check", "failure_case"]]
            .head(20)
            .astype(str)
            .to_dict("records")
        )
        return {
            "table": name,
            "passed": False,
            "n_rows": len(table),
            "errors": errors,
        }


def check_matrix(
    name: str, matrix: pd.DataFrame, lower: float | None, upper: float | None
) -> dict[str, object]:
    values = matrix.select_dtypes("number")
    problems = []
    if not matrix.index.is_unique:
        problems.append("duplicate row keys")
    if lower is not None and (values < lower).any().any():
        problems.append(f"values below {lower}")
    if upper is not None and (values > upper).any().any():
        problems.append(f"values above {upper}")
    return {
        "table": name,
        "passed": not problems,
        "n_rows": len(matrix),
        "n_samples": int(values.shape[1]),
        "errors": problems,
    }
