import json
from pathlib import Path

import pytest

from igem_warsaw_rnaseq import ncbi
from igem_warsaw_rnaseq.config import ConfigError, load_config
from igem_warsaw_rnaseq.ref import busco, qc, salmon_index, sequences


def test_gene_key_priority() -> None:
    assert sequences.gene_key("LOCUS_1", "123", "gene-X") == "LOCUS_1"
    assert sequences.gene_key(".", "123", "gene-X") == "123"
    assert sequences.gene_key("", "", "gene-X") == "X"
    assert sequences.gene_key("", ".", "gene_Y") == "Y"
    assert sequences.gene_key("", "", "Z") == "Z"


def test_parse_t2g_table_skips_short_rows() -> None:
    lines = ["rna-1\tLOC1\t.\tgene-1", "bad line", "rna-2\t.\t42\tgene-2"]
    assert sequences.parse_t2g_table(lines) == [("rna-1", "LOC1"), ("rna-2", "42")]


def test_longest_per_gene_first_on_ties_and_strips_stop() -> None:
    proteins = [("t1", "MAAA*"), ("t2", "MAAAAA"), ("t3", "MAAAAA"), ("t9", "MQ")]
    t2g = {"t1": "g1", "t2": "g1", "t3": "g1", "t4": "g2"}
    assert sequences.longest_per_gene(proteins, t2g) == [("g1", "t2", "MAAAAA")]


def test_longest_per_gene_masks_internal_stops_and_dots() -> None:
    proteins = [("t1", "MA.A*AA*."), ("t2", "MAA")]
    t2g = {"t1": "g1", "t2": "g2"}
    assert sequences.longest_per_gene(proteins, t2g) == [
        ("g1", "t1", "MAXAXAA"),
        ("g2", "t2", "MAA"),
    ]


def test_plain_gff_drops_unknown_strand(tmp_path: Path) -> None:
    src = tmp_path / "a.gff"
    src.write_text(
        "##gff-version 3\n"
        "c\tsrc\tgene\t1\t9\t.\t+\t.\tID=g1\n"
        "c\tsrc\tgene\t10\t20\t.\t?\t.\tID=g2;exception=trans-splicing\n"
        "c\tsrc\tgene\t30\t40\t.\t-\t.\tID=g3\n"
    )
    out = sequences.plain_gff(src, tmp_path / "b.gff")
    kept = [line.split("\t")[8].strip() for line in out.read_text().splitlines()[1:]]
    assert kept == ["ID=g1", "ID=g3"]


def test_read_fasta_and_decoys(tmp_path: Path) -> None:
    genome = tmp_path / "g.fa"
    genome.write_text(">chr1 desc\nACGT\nAC\n>chr2\nGG\n")
    assert list(sequences.read_fasta(genome)) == [("chr1", "ACGTAC"), ("chr2", "GG")]
    decoys = tmp_path / "decoys.txt"
    assert salmon_index.write_decoys(genome, decoys) == 2
    assert decoys.read_text() == "chr1\nchr2\n"


def test_gentrome_order(tmp_path: Path) -> None:
    tx = tmp_path / "tx.fa"
    tx.write_text(">t1\nAAA\n")
    genome = tmp_path / "g.fa"
    genome.write_text(">chr1\nCCC\n")
    out = tmp_path / "gentrome.fa"
    salmon_index.write_gentrome(tx, genome, out)
    assert out.read_text() == ">t1\nAAA\n>chr1\nCCC\n"


def test_contig_check(tmp_path: Path) -> None:
    genome = tmp_path / "g.fa"
    genome.write_text(">chr1\nACGT\n>chr2\nACGT\n")
    gff = tmp_path / "a.gff"
    gff.write_text(
        "##gff-version 3\nchr1\tsrc\tgene\t1\t4\t.\t+\t.\tID=g\nchrX\tsrc\tgene\t1\t2\t.\t+\t.\tID=h\n"
    )
    assert qc.contig_check(genome, gff) == (2, ["chrX"])


def test_count_t2g(tmp_path: Path) -> None:
    t2g = tmp_path / "t2g.tsv"
    t2g.write_text("transcript\tgene\nt1\tg1\nt2\tg1\nt3\tg2\n")
    assert qc.count_t2g(t2g) == (3, 2)


def test_busco_parse_summary(tmp_path: Path) -> None:
    summary = tmp_path / "short_summary.specific.fungi_odb12.busco.json"
    summary.write_text(
        json.dumps(
            {
                "results": {
                    "Complete percentage": 98.5,
                    "Single copy percentage": 98.2,
                    "Multi copy percentage": 0.3,
                    "Fragmented percentage": 0.5,
                    "Missing percentage": 1.0,
                    "n_markers": 758,
                }
            }
        )
    )
    result = busco.parse_summary(summary, "fungi_odb12")
    assert result.complete == 98.5
    assert result.n_markers == 758
    assert result.to_dict()["lineage"] == "fungi_odb12"


def test_ncbi_select_members() -> None:
    names = [
        "ncbi_dataset/data/GCF_1/cds_from_genomic.fna",
        "ncbi_dataset/data/GCF_1/GCF_1_ASM_genomic.fna",
        "ncbi_dataset/data/GCF_1/genomic.gff",
    ]
    assert ncbi.select_members(names) == (names[1], names[2])
    with pytest.raises(ValueError):
        ncbi.select_members(names[:1])


def _cfg(tmp_path: Path, body: str) -> Path:
    p = tmp_path / "c.yaml"
    p.write_text(body)
    return p


def test_references_config(tmp_path: Path) -> None:
    (tmp_path / "g.fa").write_text(">c\nA\n")
    (tmp_path / "a.gff").write_text("")
    cfg = load_config(
        _cfg(
            tmp_path,
            "ref: {output_dir: refs}\n"
            "references:\n"
            "  a: {ncbi_accession: GCF_1, annotation_source: NCBI}\n"
            "  b: {genome: g.fa, annotation: a.gff}\n",
        )
    )
    assert cfg.qc is None
    assert cfg.ref.output_dir == tmp_path / "refs"
    assert cfg.references["a"].ncbi_accession == "GCF_1"
    assert cfg.references["b"].genome == tmp_path / "g.fa"
    assert cfg.ref.busco.enabled is False


@pytest.mark.parametrize(
    "entry",
    ["{}", "{ncbi_accession: X, genome: g.fa, annotation: a.gff}", "{genome: g.fa}"],
)
def test_reference_requires_one_source(tmp_path: Path, entry: str) -> None:
    (tmp_path / "g.fa").write_text("")
    (tmp_path / "a.gff").write_text("")
    with pytest.raises(ConfigError, match="references.a"):
        load_config(_cfg(tmp_path, f"references:\n  a: {entry}\n"))


def test_busco_requires_lineage(tmp_path: Path) -> None:
    with pytest.raises(ConfigError, match="lineage"):
        load_config(_cfg(tmp_path, "ref:\n  busco: {enabled: true}\n"))
