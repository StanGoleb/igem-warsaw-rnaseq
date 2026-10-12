"""MAFFT alignment helpers and gap-excluded pairwise identity."""

from __future__ import annotations

import subprocess
from collections.abc import Mapping
from itertools import combinations
from pathlib import Path


def write_fasta(records: Mapping[str, str], path: Path) -> None:
    with path.open("w") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n{seq}\n")


def read_aligned(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = None
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                records[name] = ""
            elif name is not None:
                records[name] += line.strip()
    return records


def mafft(records: Mapping[str, str], work: Path, threads: int = 1) -> dict[str, str]:
    """Align ``records`` with ``mafft --auto``; returns name -> aligned sequence."""
    work.mkdir(parents=True, exist_ok=True)
    inp, out = work / "in.fa", work / "out.fa"
    write_fasta(records, inp)
    with out.open("w") as fh:
        subprocess.run(
            ["mafft", "--auto", "--quiet", "--thread", str(threads), str(inp)],
            stdout=fh,
            check=True,
        )
    return read_aligned(out)


def pair_identity(a: str, b: str) -> float | None:
    """Percent identity over columns where neither sequence has a gap."""
    same = compared = 0
    for x, y in zip(a, b, strict=True):
        if x == "-" or y == "-":
            continue
        compared += 1
        same += x.upper() == y.upper()
    return 100.0 * same / compared if compared else None


def mean_pairwise_identity(aligned: Mapping[str, str]) -> float | None:
    values = [
        v
        for a, b in combinations(aligned.values(), 2)
        if (v := pair_identity(a, b)) is not None
    ]
    return sum(values) / len(values) if values else None
