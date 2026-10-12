"""FASTQ retrieval from ENA: resumable download and checksum verification."""

from __future__ import annotations

import gzip
import hashlib
import subprocess
from pathlib import Path

RESUME_ATTEMPTS = 20
HEAD_TIMEOUT_SECONDS = 30


def ena_url(ftp_path: str) -> str:
    """ENA lists FASTQ paths without a scheme (``ftp.sra.ebi.ac.uk/...``)."""
    return f"https://{ftp_path}"


def select_fastq(
    urls: list[str], md5s: list[str], layout: str
) -> tuple[list[tuple[str, str]], bool]:
    """Choose the (url, md5) files to process and whether they form a pair.

    ENA may list three files for a PAIRED run (orphans plus ``_1``/``_2``); the
    mates are used.  A PAIRED run with a single file is treated as SINGLE.
    """
    if len(urls) != len(md5s):
        return [], False
    files = list(zip(urls, md5s, strict=True))
    if layout.upper() == "PAIRED":
        r1 = next((f for f in files if "_1.fastq" in f[0]), None)
        r2 = next((f for f in files if "_2.fastq" in f[0]), None)
        if r1 and r2:
            return [r1, r2], True
        if len(files) >= 2:
            return files[:2], True
        return files[:1], False
    unmated = next(
        (f for f in files if "_1.fastq" not in f[0] and "_2.fastq" not in f[0]), None
    )
    return [unmated or files[0]], False


def download(url: str, dst: Path) -> bool:
    """Download ``url`` to ``dst``; True only when a complete, readable gzip landed.

    ENA closes connections mid-transfer often enough that a plain retry never
    finishes large files, so each attempt resumes (``-C -``) until the local size
    reaches the remote Content-Length.
    """
    target = _content_length(url)
    for _ in range(RESUME_ATTEMPTS):
        subprocess.run(
            [
                "curl",
                "-sSL",
                "-C",
                "-",
                "--retry",
                "3",
                "--retry-all-errors",
                "--speed-limit",
                "50000",
                "--speed-time",
                "30",
                "-o",
                str(dst),
                url,
            ],
            check=False,
        )
        size = dst.stat().st_size if dst.exists() else 0
        if target is not None:
            if size >= target:
                return _gzip_ok(dst)
        elif size > 0 and _gzip_ok(dst):
            return True
    return False


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _content_length(url: str) -> int | None:
    head = subprocess.run(
        ["curl", "-sIL", "--max-time", str(HEAD_TIMEOUT_SECONDS), url],
        capture_output=True,
        text=True,
        check=False,
    ).stdout
    lengths = [
        line.split(":", 1)[1].strip()
        for line in head.splitlines()
        if line.lower().startswith("content-length:")
    ]
    try:
        return int(lengths[-1]) if lengths else None
    except ValueError:
        return None


def _gzip_ok(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        with gzip.open(path, "rb") as fh:
            return bool(fh.read(1024))
    except OSError:
        return False
