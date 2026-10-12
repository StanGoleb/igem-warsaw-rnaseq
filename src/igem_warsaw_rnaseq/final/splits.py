"""Train / val / test assignment by orthogroup, so paralogs never straddle splits."""

from __future__ import annotations

import numpy as np
import pandas as pd

from igem_warsaw_rnaseq.config import SplitConfig

SPLITS: tuple[str, ...] = ("train", "val", "test")


def assign(groups: list[str], cfg: SplitConfig) -> pd.DataFrame:
    """Shuffle ``groups`` with the seed and cut them by the configured fractions."""
    order = list(groups)
    np.random.default_rng(cfg.seed).shuffle(order)
    n = len(order)
    n_train = round(n * cfg.train)
    n_val = round(n * cfg.val)
    labels = ["train"] * n_train + ["val"] * n_val + ["test"] * (n - n_train - n_val)
    table = pd.DataFrame({"orthogroup": order, "split": labels, "seed": cfg.seed})
    return table.sort_values("orthogroup").reset_index(drop=True)
