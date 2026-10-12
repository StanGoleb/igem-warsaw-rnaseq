"""ENA Portal API client: run-level metadata."""

from __future__ import annotations

import io
import urllib.parse
import urllib.request
from collections.abc import Iterable

import pandas as pd

SEARCH_URL = "https://www.ebi.ac.uk/ena/portal/api/search"
TIMEOUT_SECONDS = 120
ACCESSION_CHUNK = 100

RUN_FIELDS: tuple[str, ...] = (
    "run_accession",
    "experiment_accession",
    "study_accession",
    "sample_accession",
    "scientific_name",
    "tax_id",
    "tax_lineage",
    "strain",
    "library_layout",
    "library_strategy",
    "library_source",
    "library_selection",
    "instrument_platform",
    "instrument_model",
    "read_count",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
    "study_title",
    "experiment_title",
    "sample_title",
    "sample_description",
    "sample_alias",
    "library_name",
)


def fetch_runs(accessions: Iterable[str]) -> pd.DataFrame:
    """Metadata for the given run accessions, in input order.

    Accessions unknown to ENA are kept as rows with ``in_ena == False``.
    """
    wanted = list(dict.fromkeys(a.strip() for a in accessions if a.strip()))
    frames = [
        _search({"includeAccessions": ",".join(wanted[i : i + ACCESSION_CHUNK])})
        for i in range(0, len(wanted), ACCESSION_CHUNK)
    ]
    found = pd.concat(frames, ignore_index=True) if frames else _empty()
    found["in_ena"] = True
    table = found.set_index("run_accession").reindex(wanted)
    table["in_ena"] = table["in_ena"].eq(True)
    return table.fillna("").reset_index()


def fetch_runs_for_taxon(taxon_id: int) -> pd.DataFrame:
    """All TRANSCRIPTOMIC runs under an NCBI taxon (ENA ``tax_tree``)."""
    table = _search(
        {"query": f'tax_tree({taxon_id}) AND library_source="TRANSCRIPTOMIC"'}
    )
    table["in_ena"] = True
    return table


def _search(params: dict[str, str]) -> pd.DataFrame:
    query = urllib.parse.urlencode(
        {
            "result": "read_run",
            "fields": ",".join(RUN_FIELDS),
            "format": "tsv",
            "limit": "0",
            **params,
        }
    )
    with urllib.request.urlopen(f"{SEARCH_URL}?{query}", timeout=TIMEOUT_SECONDS) as r:
        text = r.read().decode("utf-8")
    if not text.strip():
        return _empty()
    return pd.read_csv(io.StringIO(text), sep="\t", dtype=str, keep_default_na=False)


def _empty() -> pd.DataFrame:
    return pd.DataFrame({field: pd.Series(dtype=str) for field in RUN_FIELDS})
