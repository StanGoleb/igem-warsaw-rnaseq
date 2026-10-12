"""Gene maps and per-gene rank histograms from the classification of stage 6."""

from __future__ import annotations

import subprocess
from collections.abc import Mapping
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from igem_warsaw_rnaseq.normalize import rank
from igem_warsaw_rnaseq.quant.plots import AXIS, GRID, INK, INK_2, MUTED, SURFACE

CONSTITUTIVE_COLOUR = "#2a78d6"
INDUCIBLE_COLOUR = "#e34948"
OTHER_COLOUR = "#c3c2b7"
HIGHLIGHT = "#eda100"

CAPTION = (
    "Each point is one {unit}, ranked by expression within each of {n} samples "
    "(1 = strongest). x: strength at its weakest (1 - p95), y: at its strongest "
    "(1 - p05). Constitutive genes lie on the diagonal top right, inducible genes "
    "above the diagonal on the left."
)


def gene_map(
    summary: pd.DataFrame,
    out: Path,
    title: str,
    unit: str,
    n_samples: int,
    highlight: str | None = None,
) -> None:
    """Scatter of trough strength vs peak strength coloured by class."""
    x = 1.0 - summary["p95_pct"]
    y = 1.0 - summary["p05_pct"]
    cls = summary["class"]
    fig, ax = plt.subplots(figsize=(7, 7), facecolor=SURFACE)
    _style(ax)
    groups = (
        ("other", ~cls.isin([rank.CONSTITUTIVE, rank.INDUCIBLE]), OTHER_COLOUR, 4),
        ("constitutive", cls == rank.CONSTITUTIVE, CONSTITUTIVE_COLOUR, 10),
        ("inducible", cls == rank.INDUCIBLE, INDUCIBLE_COLOUR, 10),
    )
    for label, mask, colour, size in groups:
        ax.scatter(
            x[mask],
            y[mask],
            s=size,
            color=colour,
            alpha=0.9 if label != "other" else 0.5,
            linewidths=0,
            label=f"{label} ({int(mask.sum())})",
        )
    ax.plot([0, 1], [0, 1], color=AXIS, linewidth=0.8, linestyle="--")
    if highlight is not None and highlight in summary.index:
        ax.scatter(
            [x[highlight]],
            [y[highlight]],
            s=160,
            facecolors="none",
            edgecolors=HIGHLIGHT,
            linewidths=2,
            label=highlight,
            zorder=5,
        )
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("trough strength (1 - p95 percentile)", color=INK_2)
    ax.set_ylabel("peak strength (1 - p05 percentile)", color=INK_2)
    ax.set_title(title, color=INK, fontsize=11, loc="left")
    ax.legend(frameon=False, fontsize=8, loc="lower right", labelcolor=INK_2)
    fig.text(
        0.02,
        0.005,
        CAPTION.format(unit=unit, n=n_samples),
        fontsize=7.5,
        color=INK_2,
        wrap=True,
    )
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(out, dpi=150)
    plt.close(fig)


def gene_figure(
    gene: str,
    reference: str,
    description: str,
    gene_row: pd.Series,
    gene_strength: pd.Series,
    reference_summary: pd.DataFrame,
    orthogroup: str | None,
    og_row: pd.Series | None,
    og_strength: pd.Series | None,
    out: Path,
) -> None:
    """Map with the gene highlighted plus its rank-strength histograms."""
    panels = 2 if og_strength is None else 3
    fig = plt.figure(figsize=(5 * panels, 5.4), facecolor=SURFACE)
    axes = fig.subplots(1, panels)
    _map_axis(axes[0], reference_summary, gene, reference)
    _histogram(
        axes[1],
        gene_strength,
        f"{gene} in {reference}",
        _row_caption(gene_row),
    )
    if og_strength is not None and og_row is not None:
        _histogram(
            axes[2],
            og_strength,
            f"{orthogroup} across all samples",
            _row_caption(og_row),
        )
    fig.suptitle(
        f"{gene}: {description}" if description else gene,
        color=INK,
        fontsize=11,
        x=0.02,
        ha="left",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out, dpi=150)
    plt.close(fig)


def _map_axis(ax: plt.Axes, summary: pd.DataFrame, gene: str, reference: str) -> None:
    _style(ax)
    x = 1.0 - summary["p95_pct"]
    y = 1.0 - summary["p05_pct"]
    cls = summary["class"]
    colours = cls.map(
        {rank.CONSTITUTIVE: CONSTITUTIVE_COLOUR, rank.INDUCIBLE: INDUCIBLE_COLOUR}
    ).fillna(OTHER_COLOUR)
    ax.scatter(x, y, s=4, c=colours, alpha=0.6, linewidths=0)
    ax.plot([0, 1], [0, 1], color=AXIS, linewidth=0.8, linestyle="--")
    if gene in summary.index:
        ax.scatter(
            [x[gene]],
            [y[gene]],
            s=160,
            facecolors="none",
            edgecolors=HIGHLIGHT,
            linewidths=2,
            zorder=5,
        )
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("trough strength", color=INK_2)
    ax.set_ylabel("peak strength", color=INK_2)
    ax.set_title(f"gene map, {reference}", color=INK, fontsize=10, loc="left")


def _histogram(ax: plt.Axes, strength: pd.Series, title: str, caption: str) -> None:
    _style(ax)
    cls = caption.split(",")[0]
    colour = {
        rank.CONSTITUTIVE: CONSTITUTIVE_COLOUR,
        rank.INDUCIBLE: INDUCIBLE_COLOUR,
    }.get(cls, MUTED)
    ax.hist(strength.dropna(), bins=np.linspace(0, 1, 26), color=colour)
    ax.set_xlim(0, 1)
    ax.set_xlabel("strength (1 - rank percentile)", color=INK_2)
    ax.set_ylabel("samples", color=INK_2)
    ax.set_title(f"{title}\n{caption}", color=INK, fontsize=9, loc="left")


def _row_caption(row: pd.Series) -> str:
    return (
        f"{row['class']}, cons_frac {row['cons_frac']:.2f}, "
        f"median {1 - row['median_pct']:.2f}, amplitude {row['induction_amplitude']:.2f}, "
        f"n={int(row['n_samples'])}"
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


def gene_descriptions(gff: Path, t2g: Mapping[str, str]) -> dict[str, str]:
    """gene -> 'product (name)' from the annotation, via the transcript table."""
    table = subprocess.run(
        ["gffread", str(gff), "--table", "@id,product,gene"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    out: dict[str, str] = {}
    for line in table.splitlines():
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 3:
            continue
        gene = t2g.get(fields[0])
        if gene is None or gene in out:
            continue
        product = fields[1] if fields[1] != "." else ""
        name = fields[2] if fields[2] != "." else ""
        out[gene] = f"{product} ({name})" if product and name else product or name
    return out
