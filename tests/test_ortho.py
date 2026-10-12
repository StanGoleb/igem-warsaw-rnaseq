from pathlib import Path

import pandas as pd
import pytest

from igem_warsaw_rnaseq.ortho import alignment, conservation, newick, orthofinder, tree
from igem_warsaw_rnaseq.ref import promoters


def test_parse_cds_table_and_promoter_both_strands() -> None:
    rows = [
        "t_plus\tchr1\t+\t11-20,31-40",
        "t_minus\tchr1\t-\t51-60,71-80",
        "t_nocds\tchr1\t+\t.",
    ]
    spans = promoters.parse_cds_table(rows)
    assert set(spans) == {"t_plus", "t_minus"}
    assert spans["t_plus"] == promoters.CdsSpan("chr1", "+", 11, 40)
    assert spans["t_minus"] == promoters.CdsSpan("chr1", "-", 51, 80)

    contig = "".join(chr(ord("A") + (i % 4)) for i in range(100))
    contig = contig.translate(str.maketrans("ABCD", "ACGT"))
    plus = promoters.promoter_sequence(contig, spans["t_plus"], 5)
    assert plus == contig[5:10]
    minus = promoters.promoter_sequence(contig, spans["t_minus"], 5)
    assert minus == promoters.reverse_complement(contig[80:85])
    assert promoters.promoter_sequence(contig, spans["t_plus"], 50) == contig[:10]


def test_representative_transcripts(tmp_path: Path) -> None:
    fa = tmp_path / "p.fa"
    fa.write_text(">g1 transcript=t1\nMA\n>g2 transcript=t2\nMK\n")
    assert promoters.representative_transcripts(fa) == {"g1": "t1", "g2": "t2"}


def test_newick_patristic() -> None:
    root = newick.parse("((A:0.1,B:0.2)90:0.3,C:0.4);")
    assert newick.leaves(root) == ["A", "B", "C"]
    d = newick.patristic_distances(root)
    assert d.loc["A", "B"] == pytest.approx(0.3)
    assert d.loc["A", "C"] == pytest.approx(0.8)
    assert d.loc["B", "C"] == pytest.approx(0.9)
    assert d.loc["C", "A"] == d.loc["A", "C"]


def test_pair_identity_gap_excluded() -> None:
    assert alignment.pair_identity("AC-GT", "ACTG-") == pytest.approx(100.0)
    assert alignment.pair_identity("ACGT", "ACGA") == pytest.approx(75.0)
    assert alignment.pair_identity("--", "AA") is None
    assert alignment.mean_pairwise_identity(
        {"a": "AAAA", "b": "AAAT", "c": "AATT"}
    ) == pytest.approx((75 + 50 + 75) / 3)


def test_concatenate_partitions() -> None:
    matrix, parts = tree.concatenate(
        {"OG2": {"x": "AA", "y": "AC"}, "OG1": {"x": "GGG"}}, ["x", "y"]
    )
    assert matrix == {"x": "GGGAA", "y": "---AC"}
    assert parts == [("OG1", 1, 3), ("OG2", 4, 5)]


def test_two_leaf_tree_uses_p_distance() -> None:
    text = tree.two_leaf_tree({"a": "AAAA-A", "b": "AAAT-A"})
    assert text == "(a:0.1,b:0.1);"
    d = newick.patristic_distances(newick.parse(text))
    assert d.loc["a", "b"] == pytest.approx(0.2)


def test_ortho_qc() -> None:
    counts = pd.DataFrame(
        {"a": [1, 2, 0], "b": [1, 1, 0], "c": [1, 0, 3], "Total": [3, 3, 3]},
        index=pd.Index(["OG0", "OG1", "OG2"], name="Orthogroup"),
    )
    qc = orthofinder.ortho_qc(counts).set_index("orthogroup")
    assert list(qc.loc["OG0"]) == [3, 3, True, 1, False]
    assert list(qc.loc["OG1"]) == [2, 3, False, 2, False]
    assert list(qc.loc["OG2"]) == [1, 3, False, 3, True]


def test_parse_resolved_trees() -> None:
    text = "OG0000000: (a_g1:0.1,b_g2:0.2)n0;\n\nOG0000001: (a_g3:0.3,b_g4:0.4)n0;\n"
    assert orthofinder.parse_resolved_trees(text) == [
        ("OG0000000", "(a_g1:0.1,b_g2:0.2)n0;"),
        ("OG0000001", "(a_g3:0.3,b_g4:0.4)n0;"),
    ]


def test_members_and_shuffle_deterministic() -> None:
    assert orthofinder.members("g1, g2,g3 ") == ["g1", "g2", "g3"]
    import random

    assert conservation.shuffled("ACGTACGT", random.Random(3)) == conservation.shuffled(
        "ACGTACGT", random.Random(3)
    )
    assert sorted(conservation.shuffled("ACGTACGT", random.Random(3))) == sorted(
        "ACGTACGT"
    )
