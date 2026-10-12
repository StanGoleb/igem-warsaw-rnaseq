"""External tools each command depends on, checked before any work starts."""

from __future__ import annotations

import shutil

from igem_warsaw_rnaseq.config import ConfigError

REQUIRED: dict[str, tuple[str, ...]] = {
    "qc metadata": (),
    "qc rrna-index": ("sortmerna",),
    "qc run": ("curl", "fastp"),
    "qc collect": ("multiqc",),
    "ref prepare": ("curl", "gffread", "salmon"),
    "ref promoters": ("gffread",),
    "quant route": (),
    "quant run": ("salmon",),
    "quant collect": (),
    "ortho run": ("orthofinder",),
    "ortho tree": ("mafft", "trimal", "iqtree"),
    "ortho conservation": ("mafft",),
    "matrix build": (),
    "normalize run": (),
    "final build": ("gffread",),
    "conditions export": (),
}


def require(command: str, extra: tuple[str, ...] = ()) -> None:
    missing = [t for t in (*REQUIRED[command], *extra) if shutil.which(t) is None]
    if missing:
        raise ConfigError(
            f"{command}: required tool(s) not found on PATH: {', '.join(missing)}"
        )
