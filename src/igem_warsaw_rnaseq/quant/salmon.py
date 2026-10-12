"""salmon quant: run, trim the output directory, parse stats, aggregate to genes."""

from __future__ import annotations

import gzip
import json
import shutil
import subprocess
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

KEEP_IN_OUT = frozenset({"quant.sf", "aux_info", "lib_format_counts.json"})
KEEP_IN_AUX = frozenset({"meta_info.json", "bootstrap"})


@dataclass(frozen=True)
class MappingStats:
    num_processed: int
    num_mapped: int
    pct_mapped: float
    pct_decoy: float
    detected_libtype: str
    num_bootstraps: int


def run(
    index: Path, reads: Sequence[Path], out_dir: Path, threads: int, bootstraps: int
) -> None:
    if out_dir.exists():
        shutil.rmtree(out_dir)
    cmd = [
        "salmon",
        "quant",
        "-i",
        str(index),
        "-l",
        "A",
        "-p",
        str(threads),
        "--numBootstraps",
        str(bootstraps),
        "-o",
        str(out_dir),
        "-q",
    ]
    if len(reads) == 2:
        cmd += ["-1", str(reads[0]), "-2", str(reads[1])]
    else:
        cmd += ["-r", str(reads[0])]
    completed = subprocess.run(cmd, check=False)
    # With bootstraps enabled salmon exits 1 when no read maps, yet still writes
    # quant.sf and meta_info.json; that is a 0% mapping result, not a failure.
    if completed.returncode != 0 and not _has_outputs(out_dir):
        raise subprocess.CalledProcessError(completed.returncode, cmd)
    _trim(out_dir)


def _has_outputs(out_dir: Path) -> bool:
    return (out_dir / "quant.sf").exists() and (
        out_dir / "aux_info" / "meta_info.json"
    ).exists()


def _trim(out_dir: Path) -> None:
    """Keep quant.sf, mapping stats and bootstraps; drop the rest (file budget)."""
    for child in out_dir.iterdir():
        if child.name not in KEEP_IN_OUT:
            shutil.rmtree(child) if child.is_dir() else child.unlink()
    aux = out_dir / "aux_info"
    for child in aux.iterdir():
        if child.name not in KEEP_IN_AUX:
            shutil.rmtree(child) if child.is_dir() else child.unlink()


def read_stats(out_dir: Path) -> MappingStats:
    meta = json.loads((out_dir / "aux_info" / "meta_info.json").read_text())
    libfmt_path = out_dir / "lib_format_counts.json"
    libtype = (
        str(json.loads(libfmt_path.read_text())["expected_format"])
        if libfmt_path.exists()
        else ""
    )
    processed = int(meta["num_processed"])
    decoy = int(meta.get("num_decoy_fragments", 0))
    return MappingStats(
        num_processed=processed,
        num_mapped=int(meta["num_mapped"]),
        pct_mapped=float(meta["percent_mapped"]),
        pct_decoy=100.0 * decoy / processed if processed else 0.0,
        detected_libtype=libtype,
        num_bootstraps=int(meta.get("num_bootstraps", 0)),
    )


def load_bootstraps(out_dir: Path) -> pd.DataFrame | None:
    """Bootstrap read counts as (replicate x transcript); None when not written."""
    boot = out_dir / "aux_info" / "bootstrap"
    if not (boot / "bootstraps.gz").exists():
        return None
    with gzip.open(boot / "names.tsv.gz", "rt") as fh:
        names = fh.read().rstrip("\n").split("\t")
    with gzip.open(boot / "bootstraps.gz") as fh:
        values = np.frombuffer(fh.read(), dtype=np.float64)
    return pd.DataFrame(values.reshape(-1, len(names)), columns=names)


def aggregate_genes(
    quant_sf: Path, t2g: Mapping[str, str], bootstraps: pd.DataFrame | None
) -> pd.DataFrame:
    """Per-gene counts and TPM (isoforms summed) plus bootstrap CV of the counts."""
    quant = pd.read_csv(
        quant_sf, sep="\t", usecols=["Name", "TPM", "NumReads"], dtype={"Name": str}
    )
    quant["gene_id"] = quant["Name"].map(t2g)
    quant = quant.dropna(subset=["gene_id"])
    genes = (
        quant.groupby("gene_id")[["NumReads", "TPM"]]
        .sum()
        .rename(columns={"NumReads": "counts", "TPM": "tpm"})
    )
    genes["boot_cv"] = np.nan
    if bootstraps is not None:
        per_gene = bootstraps.T.groupby(bootstraps.columns.map(t2g)).sum().T
        mean = per_gene.mean(axis=0)
        std = per_gene.std(axis=0, ddof=1)
        cv = (std / mean).where(mean > 0)
        genes["boot_cv"] = cv.reindex(genes.index)
    return genes.reset_index()


def load_t2g(path: Path) -> dict[str, str]:
    table = pd.read_csv(path, sep="\t", dtype=str)
    return dict(zip(table["transcript"], table["gene"], strict=True))
