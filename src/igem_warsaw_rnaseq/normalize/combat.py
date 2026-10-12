"""ComBat batch correction (inmoose.pycombat_norm) with batch = bioproject."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


def run(
    log2: pd.DataFrame,
    batch: pd.Series,
    condition: pd.Series | None,
    min_batch_size: int,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Corrected matrix (NA where not correctable) and a log entry.

    Samples in batches smaller than ``min_batch_size`` and genes with missing or
    constant values across the fitted samples are left NA.
    """
    sizes = batch.value_counts()
    usable = batch[batch.map(sizes) >= min_batch_size].index
    entry: dict[str, object] = {
        "applied": False,
        "n_samples_fitted": len(usable),
        "n_samples_total": len(batch),
        "n_batches": int(batch.loc[usable].nunique()),
        "n_conditions": None,
        "design_rank": None,
        "reason_skipped": None,
    }
    out = pd.DataFrame(np.nan, index=log2.index, columns=log2.columns)
    if entry["n_batches"] < 2:
        entry["reason_skipped"] = "fewer than two batches of usable size"
        return out, entry

    data = log2[usable]
    rows = data.notna().all(axis=1) & (data.std(axis=1) > 0)
    data = data.loc[rows]
    covar = None
    if condition is not None:
        levels = condition.reindex(usable).fillna("")
        dummies = pd.get_dummies(levels, drop_first=True, dtype=float)
        entry["n_conditions"] = int(levels.nunique())
        if dummies.shape[1]:
            covar = dummies
            entry["design_rank"] = int(np.linalg.matrix_rank(dummies.to_numpy()))

    from inmoose.pycombat import pycombat_norm

    try:
        corrected = pycombat_norm(data, batch.loc[usable].to_numpy(), covar_mod=covar)
    except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
        if covar is None:
            entry["reason_skipped"] = f"pycombat_norm failed: {exc}"
            return out, entry
        log.warning(
            "ComBat with condition covariates failed (%s); refitting without", exc
        )
        entry["design_rank"] = None
        entry["reason_covariates_dropped"] = str(exc)
        corrected = pycombat_norm(data, batch.loc[usable].to_numpy())

    out.loc[data.index, usable] = pd.DataFrame(
        np.asarray(corrected), index=data.index, columns=usable
    )
    entry["applied"] = True
    entry["n_genes_fitted"] = len(data)
    return out, entry
