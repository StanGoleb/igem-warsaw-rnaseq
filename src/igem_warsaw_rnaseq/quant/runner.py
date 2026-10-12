"""Per-run quantification against the routed reference.  Idempotent per run."""

from __future__ import annotations

import json
import logging
import subprocess
from collections.abc import Mapping
from dataclasses import asdict, dataclass
from pathlib import Path

from igem_warsaw_rnaseq.config import QuantConfig
from igem_warsaw_rnaseq.qc.paths import QcPaths
from igem_warsaw_rnaseq.qc.runner import RunResult
from igem_warsaw_rnaseq.quant import salmon
from igem_warsaw_rnaseq.quant.paths import QuantPaths
from igem_warsaw_rnaseq.ref.paths import RefPaths

log = logging.getLogger(__name__)

SALMON_FAILED = "salmon_failed"


@dataclass
class QuantResult:
    run_id: str
    reference: str
    num_processed: int | None = None
    num_mapped: int | None = None
    pct_mapped: float | None = None
    pct_decoy: float | None = None
    detected_libtype: str | None = None
    n_genes_detected: int | None = None
    drop_reason: str | None = None

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(asdict(self), indent=2))

    @classmethod
    def load(cls, path: Path) -> QuantResult:
        return cls(**json.loads(path.read_text()))


def process_run(
    route: Mapping[str, str],
    cfg: QuantConfig,
    qc_root: Path,
    ref_root: Path,
    force: bool = False,
) -> QuantResult:
    run, reference = route["run_id"], route["reference"]
    paths = QuantPaths(cfg.output_dir)
    result_path = paths.result(run)
    if result_path.exists() and not force:
        cached = QuantResult.load(result_path)
        if cached.reference == reference and cached.drop_reason is None:
            log.info("%s: cached", run)
            return cached

    qc = QcPaths(qc_root)
    ref = RefPaths(ref_root, reference)
    layout = RunResult.load(qc.result(run)).layout
    reads = qc.trimmed(run, paired=(layout == "PAIRED"))
    result = QuantResult(run_id=run, reference=reference)

    log.info("%s: salmon quant vs %s (%s)", run, reference, layout)
    try:
        salmon.run(
            ref.salmon_index, reads, paths.salmon(run), cfg.threads, cfg.bootstraps
        )
    except subprocess.CalledProcessError:
        result.drop_reason = SALMON_FAILED
        result.save(result_path)
        return result

    stats = salmon.read_stats(paths.salmon(run))
    genes = salmon.aggregate_genes(
        paths.salmon(run) / "quant.sf",
        salmon.load_t2g(ref.t2g),
        salmon.load_bootstraps(paths.salmon(run)),
    )
    paths.gene_dir.mkdir(parents=True, exist_ok=True)
    genes.to_csv(paths.gene_table(run), sep="\t", index=False, float_format="%.6g")

    result.num_processed = stats.num_processed
    result.num_mapped = stats.num_mapped
    result.pct_mapped = stats.pct_mapped
    result.pct_decoy = stats.pct_decoy
    result.detected_libtype = stats.detected_libtype
    result.n_genes_detected = int((genes["counts"] >= cfg.detected_min_reads).sum())
    result.save(result_path)
    log.info(
        "%s: %.1f%% mapped, %d genes detected",
        run,
        stats.pct_mapped,
        result.n_genes_detected,
    )
    return result
