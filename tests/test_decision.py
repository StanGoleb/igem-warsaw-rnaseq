import pytest

from igem_warsaw_rnaseq.config import QcGates
from igem_warsaw_rnaseq.qc import decision

ALLOWED = frozenset({"RNA-Seq"})
EXCLUDED = frozenset({"ABI_SOLID"})


def _row(**overrides: object) -> dict[str, object]:
    row: dict[str, object] = {
        "in_ena": True,
        "instrument_platform": "ILLUMINA",
        "library_source": "TRANSCRIPTOMIC",
        "library_strategy": "RNA-Seq",
        "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/x.fastq.gz",
    }
    row.update(overrides)
    return row


@pytest.mark.parametrize(
    ("overrides", "expected"),
    [
        ({}, None),
        ({"in_ena": False}, decision.NOT_IN_ENA),
        ({"instrument_platform": "ABI_SOLID"}, decision.UNSUPPORTED_PLATFORM),
        ({"library_strategy": "miRNA-Seq"}, decision.NOT_RNASEQ),
        ({"library_source": "GENOMIC"}, decision.NOT_RNASEQ),
        ({"fastq_ftp": ""}, decision.NO_FASTQ),
    ],
)
def test_metadata_drop_reason(
    overrides: dict[str, object], expected: str | None
) -> None:
    assert (
        decision.metadata_drop_reason(_row(**overrides), ALLOWED, EXCLUDED) == expected
    )


def test_metadata_gate_order_platform_before_strategy() -> None:
    row = _row(instrument_platform="ABI_SOLID", library_strategy="WGS")
    assert (
        decision.metadata_drop_reason(row, ALLOWED, EXCLUDED)
        == decision.UNSUPPORTED_PLATFORM
    )


def test_gates_off_by_default() -> None:
    gates = QcGates(min_pct_q30=None, max_pct_adapter=None)
    assert decision.gate_drop_reason(10.0, 95.0, gates) is None


def test_gates_when_set() -> None:
    gates = QcGates(min_pct_q30=80.0, max_pct_adapter=10.0)
    assert decision.gate_drop_reason(79.9, 0.0, gates) == decision.LOW_Q30
    assert decision.gate_drop_reason(80.0, 10.0, gates) is None
    assert decision.gate_drop_reason(90.0, 10.1, gates) == decision.HIGH_ADAPTER
    assert decision.gate_drop_reason(None, None, gates) is None
