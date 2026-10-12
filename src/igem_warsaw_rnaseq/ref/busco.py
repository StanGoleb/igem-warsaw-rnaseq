"""BUSCO completeness of the proteome."""

from __future__ import annotations

import json
import shutil
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path

KEEP = frozenset({"full_table.tsv", "missing_busco_list.tsv"})


@dataclass(frozen=True)
class BuscoResult:
    lineage: str
    complete: float
    single_copy: float
    multi_copy: float
    fragmented: float
    missing: float
    n_markers: int

    def to_dict(self) -> dict[str, float | int | str]:
        return asdict(self)


def run(
    proteins: Path, lineage: str, out_dir: Path, download_path: Path, threads: int
) -> Path:
    """Run BUSCO in protein mode; returns the short-summary JSON."""
    if out_dir.exists():
        shutil.rmtree(out_dir)
    subprocess.run(
        [
            "busco",
            "-i",
            str(proteins),
            "-m",
            "proteins",
            "-l",
            lineage,
            "-o",
            out_dir.name,
            "--out_path",
            str(out_dir.parent),
            "--download_path",
            str(download_path),
            "-c",
            str(threads),
            "-f",
        ],
        check=True,
    )
    _trim(out_dir)
    summaries = sorted(out_dir.glob("short_summary.*.json"))
    if not summaries:
        raise FileNotFoundError(f"{out_dir}: no short_summary.*.json written by BUSCO")
    return summaries[0]


def parse_summary(summary_json: Path, lineage: str) -> BuscoResult:
    results = json.loads(summary_json.read_text())["results"]
    return BuscoResult(
        lineage=lineage,
        complete=float(results["Complete percentage"]),
        single_copy=float(results["Single copy percentage"]),
        multi_copy=float(results["Multi copy percentage"]),
        fragmented=float(results["Fragmented percentage"]),
        missing=float(results["Missing percentage"]),
        n_markers=int(results["n_markers"]),
    )


def _trim(out_dir: Path) -> None:
    """Keep summaries and the full table; drop per-marker outputs (thousands of files)."""
    for child in out_dir.iterdir():
        if child.name.startswith("short_summary."):
            continue
        if child.is_dir():
            for inner in child.iterdir():
                if inner.name in KEEP:
                    shutil.move(str(inner), out_dir / inner.name)
            shutil.rmtree(child, ignore_errors=True)
        elif child.name not in KEEP:
            child.unlink()
