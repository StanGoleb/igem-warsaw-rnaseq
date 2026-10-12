"""Raw material for the condition layer: every text ENA/BioSample holds per sample.

The pipeline does not categorise conditions.  It exports what a curator (or a
language model driven by one) needs, and defines the file it reads back:
``conditions.tsv`` with ``sample_id`` (ENA experiment), ``condition`` (category)
and ``description`` (normalised prose); extra columns pass through to the final
dataset.
"""

from __future__ import annotations

import urllib.request
import xml.etree.ElementTree as ET
from collections.abc import Iterable
from pathlib import Path

import pandas as pd

BROWSER_XML_URL = "https://www.ebi.ac.uk/ena/browser/api/xml/{accessions}"
TIMEOUT_SECONDS = 120
CHUNK = 100

RAW_COLUMNS: tuple[str, ...] = (
    "sample_id",
    "bioproject",
    "sample_accession",
    "runs",
    "study_title",
    "experiment_title",
    "sample_title",
    "sample_description",
    "sample_alias",
    "library_name",
    "biosample_title",
    "biosample_attributes",
)
TEMPLATE_COLUMNS: tuple[str, ...] = ("sample_id", "condition", "description")


def parse_sample_xml(text: str) -> dict[str, tuple[str, str]]:
    """accession -> (title, 'tag=value; ...') from ENA browser SAMPLE XML."""
    out: dict[str, tuple[str, str]] = {}
    for sample in ET.fromstring(text).findall("SAMPLE"):
        attrs = [
            f"{a.findtext('TAG', '')}={a.findtext('VALUE', '')}"
            for a in sample.findall(".//SAMPLE_ATTRIBUTE")
            if a.findtext("TAG")
        ]
        out[sample.get("accession", "")] = (
            sample.findtext("TITLE", "") or "",
            "; ".join(attrs),
        )
    return out


def biosample_attributes(accessions: Iterable[str]) -> dict[str, tuple[str, str]]:
    wanted = sorted({a for a in accessions if a})
    out: dict[str, tuple[str, str]] = {}
    for i in range(0, len(wanted), CHUNK):
        url = BROWSER_XML_URL.format(accessions=",".join(wanted[i : i + CHUNK]))
        with urllib.request.urlopen(url, timeout=TIMEOUT_SECONDS) as r:
            out.update(parse_sample_xml(r.read().decode("utf-8")))
    return out


def raw_table(
    metadata: pd.DataFrame, attributes: dict[str, tuple[str, str]]
) -> pd.DataFrame:
    """One row per experiment (sample) with every text field, runs joined."""
    rows = []
    found = metadata[metadata["in_ena"] | (metadata["source"] == "local")].copy()
    found["experiment_accession"] = found["experiment_accession"].where(
        found["experiment_accession"].astype(bool), found["run_accession"]
    )
    for sample, group in found.groupby("experiment_accession", sort=True):
        first = group.iloc[0]
        title, attrs = attributes.get(first["sample_accession"], ("", ""))
        rows.append(
            {
                "sample_id": sample,
                "bioproject": first["study_accession"],
                "sample_accession": first["sample_accession"],
                "runs": ";".join(sorted(group["run_accession"])),
                "study_title": first["study_title"],
                "experiment_title": first["experiment_title"],
                "sample_title": first["sample_title"],
                "sample_description": first["sample_description"],
                "sample_alias": first["sample_alias"],
                "library_name": first["library_name"],
                "biosample_title": title,
                "biosample_attributes": attrs,
            }
        )
    return pd.DataFrame(rows, columns=list(RAW_COLUMNS))


def export(metadata: pd.DataFrame, raw_out: Path, template_out: Path) -> pd.DataFrame:
    table = raw_table(metadata, biosample_attributes(metadata["sample_accession"]))
    table.to_csv(raw_out, sep="\t", index=False)
    template = pd.DataFrame(
        {"sample_id": table["sample_id"], "condition": "", "description": ""}
    )
    template[list(TEMPLATE_COLUMNS)].to_csv(template_out, sep="\t", index=False)
    return table
