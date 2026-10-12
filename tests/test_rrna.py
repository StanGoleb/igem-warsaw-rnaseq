from pathlib import Path

import pytest

from igem_warsaw_rnaseq.qc.rrna import parse_pct_rrna

LOG = """
 Results:
    Total reads = 26566730
    Total reads for de novo clustering = 0
    Total reads passing E-value threshold = 1013431 (3.81)
    Total reads failing E-value threshold = 25553299 (96.19)
"""


def test_parse_pct_rrna(tmp_path: Path) -> None:
    log = tmp_path / "run.log"
    log.write_text(LOG)
    assert parse_pct_rrna(log) == 3.81


def test_parse_pct_rrna_missing_line(tmp_path: Path) -> None:
    log = tmp_path / "run.log"
    log.write_text("nothing here\n")
    with pytest.raises(ValueError):
        parse_pct_rrna(log)
