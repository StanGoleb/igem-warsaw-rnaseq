"""Per-orthogroup protein vs promoter conservation, with a shuffled-promoter background."""

from __future__ import annotations

import logging
import random
import shutil
from collections.abc import Mapping
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

from igem_warsaw_rnaseq.ortho import alignment
from igem_warsaw_rnaseq.ortho.orthofinder import (
    load_orthogroups,
    members,
    single_copy_ids,
)
from igem_warsaw_rnaseq.ref.sequences import read_fasta

log = logging.getLogger(__name__)

COLUMNS: tuple[str, ...] = (
    "orthogroup",
    "n_aligned",
    "prot_mean_pid",
    "prom_mean_pid",
    "prom_background_pid",
    "prom_above_background",
)


def shuffled(seq: str, rng: random.Random) -> str:
    chars = list(seq)
    rng.shuffle(chars)
    return "".join(chars)


def conserve_one(
    args: tuple[str, dict[str, str], dict[str, str], int, Path],
) -> dict[str, object]:
    og, proteins, promoters, seed, work = args
    work = work / og
    rng = random.Random(seed)
    prot = alignment.mean_pairwise_identity(alignment.mafft(proteins, work / "prot"))
    prom = alignment.mean_pairwise_identity(alignment.mafft(promoters, work / "prom"))
    background = alignment.mean_pairwise_identity(
        alignment.mafft(
            {k: shuffled(v, rng) for k, v in promoters.items()}, work / "bg"
        )
    )
    shutil.rmtree(work, ignore_errors=True)
    above = None
    if prom is not None and background is not None:
        above = prom > background
    return {
        "orthogroup": og,
        "n_aligned": len(promoters),
        "prot_mean_pid": prot,
        "prom_mean_pid": prom,
        "prom_background_pid": background,
        "prom_above_background": above,
    }


def run(
    orthogroups_tsv: Path,
    single_copy_txt: Path,
    proteomes: Mapping[str, Path],
    promoters: Mapping[str, Path],
    work: Path,
    out: Path,
    seed: int,
    threads: int,
) -> pd.DataFrame:
    table = load_orthogroups(orthogroups_tsv)
    prot_seqs = {ref: dict(read_fasta(p)) for ref, p in proteomes.items()}
    prom_seqs = {ref: dict(read_fasta(p)) for ref, p in promoters.items()}
    jobs = []
    for i, og in enumerate(single_copy_ids(single_copy_txt)):
        row = table.loc[og]
        genes = {ref: members(row[ref])[0] for ref in proteomes if row[ref]}
        prots = {ref: prot_seqs[ref][g] for ref, g in genes.items()}
        proms = {
            ref: prom_seqs[ref][g] for ref, g in genes.items() if g in prom_seqs[ref]
        }
        if len(proms) < 2:
            continue
        jobs.append((og, prots, proms, seed + i, work))
    log.info("conservation: %d orthogroups", len(jobs))
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    with ProcessPoolExecutor(max_workers=threads) as pool:
        rows = list(pool.map(conserve_one, jobs, chunksize=8))
    shutil.rmtree(work, ignore_errors=True)
    result = pd.DataFrame(rows, columns=list(COLUMNS))
    result.to_csv(out, sep="\t", index=False, float_format="%.4g")
    return result
