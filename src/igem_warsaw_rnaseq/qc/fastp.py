"""fastp: adapter/quality trimming and report parsing."""

from __future__ import annotations

import json
import subprocess
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class FastpMetrics:
    reads_raw: int
    reads_after_fastp: int
    pct_q30: float
    pct_adapter: float
    read_length_after_fastp: float


def run_fastp(
    inputs: Sequence[Path],
    outputs: Sequence[Path],
    report_json: Path,
    report_html: Path,
    threads: int,
    extra_args: Sequence[str] = (),
) -> None:
    cmd = [
        "fastp",
        "--thread",
        str(threads),
        "--json",
        str(report_json),
        "--html",
        str(report_html),
        "--in1",
        str(inputs[0]),
        "--out1",
        str(outputs[0]),
    ]
    if len(inputs) == 2:
        cmd += ["--in2", str(inputs[1]), "--out2", str(outputs[1])]
    cmd += list(extra_args)
    subprocess.run(cmd, check=True)


def parse_report(report_json: Path) -> FastpMetrics:
    """Metrics from a fastp JSON report.

    Read counts follow fastp: for paired-end data both mates are counted, so
    ``reads_raw`` is twice the number of pairs.
    """
    report = json.loads(report_json.read_text())
    before = report["summary"]["before_filtering"]
    after = report["summary"]["after_filtering"]
    total = int(before["total_reads"])
    trimmed = int(report.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0))
    lengths = [
        float(after[k])
        for k in ("read1_mean_length", "read2_mean_length")
        if k in after
    ]
    return FastpMetrics(
        reads_raw=total,
        reads_after_fastp=int(after["total_reads"]),
        pct_q30=100.0 * float(before["q30_rate"]),
        pct_adapter=100.0 * trimmed / total if total else 0.0,
        read_length_after_fastp=sum(lengths) / len(lengths) if lengths else 0.0,
    )
