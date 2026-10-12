"""SortMeRNA: fraction of reads matching rRNA references."""

from __future__ import annotations

import re
import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path

_PASSING_RE = re.compile(r"Total reads passing E-value threshold = \d+ \(([\d.]+)\)")


def build_index(references: Sequence[Path], idx_dir: Path, threads: int) -> None:
    """Build the reference index once; per-run calls then reuse it."""
    idx_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "sortmerna",
        "--index",
        "1",
        "--idx-dir",
        str(idx_dir),
        "--threads",
        str(threads),
    ]
    for ref in references:
        cmd += ["--ref", str(ref)]
    subprocess.run(cmd, check=True)


def run_sortmerna(
    reads: Sequence[Path],
    references: Sequence[Path],
    idx_dir: Path,
    log_prefix: Path,
    workdir: Path,
    threads: int,
) -> Path:
    """Screen ``reads``; returns the log written as ``<log_prefix>.log``."""
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)
    log_prefix.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "sortmerna",
        "--index",
        "0",
        "--idx-dir",
        str(idx_dir),
        "--kvdb",
        str(workdir / "kvdb"),
        "--workdir",
        str(workdir),
        "--aligned",
        str(log_prefix),
        "--threads",
        str(threads),
    ]
    for ref in references:
        cmd += ["--ref", str(ref)]
    for r in reads:
        cmd += ["--reads", str(r)]
    subprocess.run(cmd, check=True)
    shutil.rmtree(workdir, ignore_errors=True)
    return log_prefix.with_suffix(".log")


def parse_pct_rrna(log_path: Path) -> float:
    match = _PASSING_RE.search(log_path.read_text())
    if match is None:
        raise ValueError(f"{log_path}: no 'Total reads passing E-value threshold' line")
    return float(match.group(1))
