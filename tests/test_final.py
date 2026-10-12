import pandas as pd

from igem_warsaw_rnaseq.config import SplitConfig
from igem_warsaw_rnaseq.final import sequences, splits, validate
from igem_warsaw_rnaseq.ref.promoters import CdsSpan, reverse_complement

CONTIG = "".join("ACGT"[(i * 7) % 4] for i in range(200))


def test_windows_plus_strand_clipped() -> None:
    w = sequences.windows(
        CONTIG, CdsSpan("c", "+", 51, 80), promoter_bp=10, flank_bp=20
    )
    assert w.promoter == CONTIG[40:50]
    assert w.flank_up == CONTIG[20:40]
    assert w.gene == CONTIG[50:80]
    assert w.flank_down == CONTIG[80:100]
    assert w.promoter_bed == (40, 50)
    near_start = sequences.windows(CONTIG, CdsSpan("c", "+", 6, 20), 10, 20)
    assert near_start.promoter == CONTIG[0:5] and near_start.flank_up == ""


def test_windows_minus_strand_mirrors() -> None:
    w = sequences.windows(
        CONTIG, CdsSpan("c", "-", 51, 80), promoter_bp=10, flank_bp=20
    )
    assert w.promoter == reverse_complement(CONTIG[80:90])
    assert w.flank_up == reverse_complement(CONTIG[90:110])
    assert w.gene == reverse_complement(CONTIG[50:80])
    assert w.flank_down == reverse_complement(CONTIG[30:50])
    assert w.promoter_bed == (80, 90)
    assert list(w.by_group()) == list(sequences.GROUPS)


def test_split_fractions_and_determinism() -> None:
    cfg = SplitConfig(train=0.8, val=0.1, test=0.1, seed=7)
    groups = [f"OG{i}" for i in range(100)]
    a = splits.assign(groups, cfg)
    b = splits.assign(groups, cfg)
    assert a.equals(b)
    assert a["split"].value_counts().to_dict() == {"train": 80, "val": 10, "test": 10}
    assert set(a["orthogroup"]) == set(groups)


def test_validate_reports_failures() -> None:
    table = pd.DataFrame(
        {
            "SPECIES": ["a"],
            "SEQUENCE": ["ACGTX"],
            "GROUP": ["Promotor"],
            "ID": ["a_G1"],
            "OG_ID": [None],
            "GENE_ID": ["g"],
        }
    )
    result = validate.check("sequences", validate.SEQUENCES, table)
    assert result["passed"] is False and result["errors"]
    ok = validate.check_matrix(
        "m", pd.DataFrame({"s": [0.0, 1.0]}, index=["a", "b"]), 0.0, 1.0
    )
    assert ok["passed"] is True
    bad = validate.check_matrix(
        "m", pd.DataFrame({"s": [0.0, 2.0]}, index=["a", "a"]), 0.0, 1.0
    )
    assert bad["errors"] == ["duplicate row keys", "values above 1.0"]
