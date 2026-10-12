"""Within-sample rank percentiles, per-gene rank metrics and expression classes."""

from __future__ import annotations

import numpy as np
import pandas as pd

from igem_warsaw_rnaseq.config import ClassConfig

CONSTITUTIVE = "constitutive"
SILENT = "silent"
INDUCIBLE = "inducible"
VARIABLE = "variable"


def rank_pct(tpm: pd.DataFrame, ties: str = "max") -> pd.DataFrame:
    """rank / genes ranked per sample; 0 = strongest.  Absent genes (NA) -> 1.0.

    ``ties="max"`` sends every zero-TPM gene to exactly 1.0 (the dataset convention);
    ``ties="average"`` is what the class metrics use.
    """
    ranks = tpm.rank(ascending=False, method=ties)
    pct = ranks.div(ranks.notna().sum(axis=0), axis=1)
    return pct.fillna(1.0)


def cons_frac(pct: pd.DataFrame, width: float) -> pd.Series:
    """Fraction of a gene's samples inside its densest rank window of ``width``."""
    values = pct.to_numpy()
    out = np.empty(values.shape[0])
    idx = np.arange(values.shape[1])
    for i in range(values.shape[0]):
        row = np.sort(values[i][~np.isnan(values[i])])
        n = row.size
        if n == 0:
            out[i] = np.nan
            continue
        end = np.searchsorted(row, row + width, side="right")
        out[i] = (end - idx[:n]).max() / n
    return pd.Series(out, index=pct.index)


def metrics(tpm: pd.DataFrame, cfg: ClassConfig) -> pd.DataFrame:
    """Per-gene rank metrics on ties-averaged percentiles (absent samples ignored)."""
    ranks = tpm.rank(ascending=False, method="average")
    pct = ranks.div(ranks.notna().sum(axis=0), axis=1)
    p05 = pct.quantile(0.05, axis=1)
    p95 = pct.quantile(0.95, axis=1)
    median = pct.median(axis=1)
    return pd.DataFrame(
        {
            "n_samples": pct.notna().sum(axis=1),
            "median_pct": median,
            "p05_pct": p05,
            "p95_pct": p95,
            "band": p95 - p05,
            "rank_iqr": pct.quantile(0.75, axis=1) - pct.quantile(0.25, axis=1),
            "rank_max": pct.min(axis=1),
            "cons_frac": cons_frac(pct, cfg.window),
            "induction_amplitude": median - p05,
        }
    )


def classify(m: pd.DataFrame, cfg: ClassConfig) -> pd.DataFrame:
    """constitutive / silent / inducible / variable, plus the is_strong filter.

    Constitutive = always in the same rank range (stability), whatever the level;
    level only separates usable (constitutive) from off (silent).  Inducible is
    defined by the swing from baseline to peak, not by a median cutoff.
    """
    always = m["cons_frac"] >= cfg.min_cons_frac
    off = m["median_pct"] > cfg.silent_min_median_pct
    inducible = (
        ~always
        & (m["p05_pct"] <= cfg.top_pct)
        & (m["induction_amplitude"] >= cfg.inducible_min_amplitude)
    )
    label = pd.Series(VARIABLE, index=m.index)
    label[inducible] = INDUCIBLE
    label[always & off] = SILENT
    label[always & ~off] = CONSTITUTIVE
    return m.assign(
        **{
            "class": label,
            "is_strong": (label == CONSTITUTIVE) & (m["median_pct"] <= cfg.top_pct),
        }
    )


def best_sample(tpm: pd.DataFrame) -> pd.Series:
    """Sample in which each gene ranks best (NA when the gene is never measured)."""
    ranks = tpm.rank(ascending=False, method="average")
    measured = ranks.notna().any(axis=1)
    best = pd.Series(pd.NA, index=ranks.index, dtype=object)
    best[measured] = ranks[measured].idxmin(axis=1)
    return best
