"""Keep/drop rules for the QC stage.

A run is dropped only for mechanical failures and metadata mismatches.  Quality
metrics (pct_q30, pct_adapter, pct_rrna, read counts) are always reported; the
two optional gates are off unless set in the configuration.
"""

from __future__ import annotations

from collections.abc import Mapping
from collections.abc import Set as AbstractSet
from pathlib import Path

from igem_warsaw_rnaseq.config import QcGates

KEEP = "keep"
DROP = "drop"
PENDING = "pending"

NOT_IN_ENA = "not_in_ena"
MISSING_FILE = "missing_file"
UNSUPPORTED_PLATFORM = "unsupported_platform"
NOT_RNASEQ = "not_rnaseq"
NO_FASTQ = "no_fastq"
DOWNLOAD_FAILED = "download_failed"
MD5_MISMATCH = "md5_mismatch"
FASTP_FAILED = "fastp_failed"
NO_READS_AFTER_FASTP = "no_reads_after_fastp"
LOW_Q30 = "low_q30"
HIGH_ADAPTER = "high_adapter"


def metadata_drop_reason(
    row: Mapping[str, object],
    allowed_strategies: AbstractSet[str],
    excluded_platforms: AbstractSet[str],
) -> str | None:
    if row.get("source", "ena") == "local":
        files = [p for p in str(row["fastq_ftp"]).split(";") if p]
        if not files or not all(Path(p).is_file() for p in files):
            return MISSING_FILE
        return None
    if not row["in_ena"]:
        return NOT_IN_ENA
    if row["instrument_platform"] in excluded_platforms:
        return UNSUPPORTED_PLATFORM
    if (
        row["library_source"] != "TRANSCRIPTOMIC"
        or row["library_strategy"] not in allowed_strategies
    ):
        return NOT_RNASEQ
    if not row["fastq_ftp"]:
        return NO_FASTQ
    return None


def gate_drop_reason(
    pct_q30: float | None, pct_adapter: float | None, gates: QcGates
) -> str | None:
    if (
        gates.min_pct_q30 is not None
        and pct_q30 is not None
        and pct_q30 < gates.min_pct_q30
    ):
        return LOW_Q30
    if (
        gates.max_pct_adapter is not None
        and pct_adapter is not None
        and pct_adapter > gates.max_pct_adapter
    ):
        return HIGH_ADAPTER
    return None
