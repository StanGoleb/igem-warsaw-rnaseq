import hashlib
from pathlib import Path

from igem_warsaw_rnaseq.download import ena_url, md5sum, select_fastq

BASE = "ftp.sra.ebi.ac.uk/vol1/fastq/SRR1/SRR1"


def test_select_paired_mates_over_orphans() -> None:
    urls = [f"{BASE}.fastq.gz", f"{BASE}_1.fastq.gz", f"{BASE}_2.fastq.gz"]
    files, paired = select_fastq(urls, ["m0", "m1", "m2"], "PAIRED")
    assert paired is True
    assert files == [(urls[1], "m1"), (urls[2], "m2")]


def test_select_paired_with_single_file_falls_back_to_single() -> None:
    files, paired = select_fastq([f"{BASE}.fastq.gz"], ["m0"], "PAIRED")
    assert paired is False
    assert files == [(f"{BASE}.fastq.gz", "m0")]


def test_select_single_prefers_unmated_file() -> None:
    urls = [f"{BASE}_1.fastq.gz", f"{BASE}.fastq.gz"]
    files, paired = select_fastq(urls, ["m1", "m0"], "SINGLE")
    assert paired is False
    assert files == [(urls[1], "m0")]


def test_md5sum(tmp_path: Path) -> None:
    path = tmp_path / "x.bin"
    path.write_bytes(b"abc" * 100000)
    assert md5sum(path) == hashlib.md5(b"abc" * 100000).hexdigest()


def test_ena_url() -> None:
    assert (
        ena_url("ftp.sra.ebi.ac.uk/vol1/x.gz") == "https://ftp.sra.ebi.ac.uk/vol1/x.gz"
    )
