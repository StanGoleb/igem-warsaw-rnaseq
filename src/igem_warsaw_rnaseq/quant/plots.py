"""QC figures per reference: PCA scatter coloured by bioproject, Spearman heatmap."""

from __future__ import annotations

from pathlib import Path

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

CATEGORICAL = (
    "#2a78d6",
    "#eb6834",
    "#1baf7a",
    "#eda100",
    "#e87ba4",
    "#008300",
    "#4a3aa7",
    "#e34948",
)
OTHER = "#898781"
SEQUENTIAL = ("#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b")
SURFACE, INK, INK_2, MUTED, GRID, AXIS = (
    "#fcfcfb",
    "#0b0b0b",
    "#52514e",
    "#898781",
    "#e1e0d9",
    "#c3c2b7",
)


def write(
    plots_dir: Path, reference: str, qc: pd.DataFrame, matrix: pd.DataFrame
) -> None:
    plots_dir.mkdir(parents=True, exist_ok=True)
    _pca(qc, plots_dir / f"pca_{reference}.png", reference)
    _heatmap(
        matrix.corr(method="spearman"),
        plots_dir / f"spearman_{reference}.png",
        reference,
    )


def _style(ax: plt.Axes) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.grid(color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)


def _pca(qc: pd.DataFrame, out: Path, reference: str) -> None:
    groups = list(dict.fromkeys(sorted(qc["bioproject"].fillna(""))))
    named = groups[: len(CATEGORICAL)]
    colour = {g: CATEGORICAL[i] for i, g in enumerate(named)}
    fig, ax = plt.subplots(figsize=(7, 5), facecolor=SURFACE)
    _style(ax)
    for group in groups:
        sub = qc[qc["bioproject"].fillna("") == group]
        label = group if group in colour else "Other"
        ax.scatter(
            sub["pca1"],
            sub["pca2"],
            s=48,
            color=colour.get(group, OTHER),
            edgecolors=SURFACE,
            linewidths=1.5,
            label=label if label != "Other" else None,
        )
    if len(groups) > len(named):
        ax.scatter(
            [], [], s=48, color=OTHER, label=f"Other ({len(groups) - len(named)})"
        )
    ax.set_xlabel("PC1", color=INK_2)
    ax.set_ylabel("PC2", color=INK_2)
    ax.set_title(
        f"{reference}: PCA of log2(TPM+1), {len(qc)} runs",
        color=INK,
        fontsize=11,
        loc="left",
    )
    if len(groups) > 1:
        ax.legend(
            title="bioproject",
            fontsize=8,
            title_fontsize=8,
            frameon=False,
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            labelcolor=INK_2,
        )
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def _heatmap(corr: pd.DataFrame, out: Path, reference: str) -> None:
    n = len(corr)
    cmap = LinearSegmentedColormap.from_list("sequential", SEQUENTIAL)
    fig, ax = plt.subplots(figsize=(6.5, 5.5), facecolor=SURFACE)
    image = ax.imshow(
        corr.to_numpy(), cmap=cmap, vmin=float(corr.min().min()), vmax=1.0
    )
    ax.set_facecolor(SURFACE)
    show_labels = n <= 40
    ax.set_xticks(range(n) if show_labels else [])
    ax.set_yticks(range(n) if show_labels else [])
    if show_labels:
        ax.set_xticklabels(corr.columns, rotation=90, fontsize=6, color=MUTED)
        ax.set_yticklabels(corr.index, fontsize=6, color=MUTED)
    for side in ax.spines.values():
        side.set_visible(False)
    ax.set_title(
        f"{reference}: Spearman correlation of log2(TPM+1), {n} runs",
        color=INK,
        fontsize=11,
        loc="left",
    )
    bar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.02)
    bar.outline.set_visible(False)
    bar.ax.tick_params(colors=MUTED, labelsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
