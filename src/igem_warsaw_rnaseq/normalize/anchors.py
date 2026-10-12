"""Anchor (reference-gene) normalisation: configured genes or geNorm's most stable."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd

from igem_warsaw_rnaseq.config import AnchorConfig

SOURCE_CONFIG = "config"
SOURCE_GENORM = "genorm"


def candidates(log2: pd.DataFrame, cfg: AnchorConfig) -> list[str]:
    """Genes measured in every sample, expressed above the floor, top by mean."""
    complete = log2.loc[log2.notna().all(axis=1)]
    means = complete.mean(axis=1)
    expressed = means[means >= cfg.min_mean_log2tpm]
    return (
        expressed.sort_values(ascending=False).head(cfg.candidate_top_n).index.tolist()
    )


def genorm_m(panel: pd.DataFrame) -> pd.Series:
    """geNorm stability M (Vandesompele 2002): mean SD of log-ratios to other genes."""
    values = panel.to_numpy(dtype=float)
    n = values.shape[0]
    m = np.empty(n)
    for j in range(n):
        v = (values - values[j]).std(axis=1, ddof=1)
        m[j] = v.sum() / (n - 1) if n > 1 else np.nan
    return pd.Series(m, index=panel.index, name="genorm_M")


def select(
    log2: pd.DataFrame, configured: Sequence[str], cfg: AnchorConfig
) -> tuple[list[str], pd.DataFrame]:
    """(anchor genes, table with gene, source, genorm_M)."""
    present = [g for g in configured if g in log2.index]
    if present:
        panel = log2.loc[present]
        m = genorm_m(panel) if len(present) > 1 else pd.Series(np.nan, index=present)
        table = pd.DataFrame(
            {"gene": present, "source": SOURCE_CONFIG, "genorm_M": m.to_numpy()}
        )
        return present, table
    pool = candidates(log2, cfg)
    if len(pool) < 2:
        return [], pd.DataFrame(columns=["gene", "source", "genorm_M"])
    m = genorm_m(log2.loc[pool]).sort_values()
    chosen = m.head(cfg.n_auto)
    table = pd.DataFrame(
        {"gene": chosen.index, "source": SOURCE_GENORM, "genorm_M": chosen.to_numpy()}
    )
    return chosen.index.tolist(), table


def normalise(log2: pd.DataFrame, anchors: Sequence[str]) -> pd.DataFrame:
    """log2(TPM+1) minus the per-sample mean log2(TPM+1) of the anchors."""
    if not anchors:
        return pd.DataFrame(np.nan, index=log2.index, columns=log2.columns)
    offset = log2.loc[list(anchors)].mean(axis=0)
    return log2.sub(offset, axis=1)
