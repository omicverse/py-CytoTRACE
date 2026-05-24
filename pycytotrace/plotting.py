"""Visualisation — 1:1 ports of CytoTRACE::plotCytoTRACE + plotCytoGenes."""

from __future__ import annotations

import numpy as np
import pandas as pd

from ggplot2_py import (
    aes,
    geom_bar,
    geom_point,
    ggplot,
    labs,
    scale_color_gradientn,
    theme_classic,
)

from .core import CytoTRACE


_RDYLBU_REV = [
    "#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#FFFFBF",
    "#FEE090", "#FDAE61", "#F46D43", "#D73027",
]


def plot_cytotrace(
    result: CytoTRACE,
    embedding: np.ndarray,
    *,
    dim_x: int = 1,
    dim_y: int = 2,
    point_size: float = 1.5,
):
    """Scatter of an external embedding (UMAP / PCA / tSNE) coloured by CytoTRACE score.

    Args:
        result: ``CytoTRACE`` returned by ``cytotrace_run``.
        embedding: (n_cells × ≥2) coords.
        dim_x, dim_y: 1-indexed components.
    """
    i, j = int(dim_x) - 1, int(dim_y) - 1
    arr = np.asarray(embedding, dtype=np.float64)
    df = pd.DataFrame(
        {
            f"Dim{i + 1}": arr[:, i],
            f"Dim{j + 1}": arr[:, j],
            "CytoTRACE": result.cytotrace,
        }
    )
    return (
        ggplot(df, aes(x=df.columns[0], y=df.columns[1], colour="CytoTRACE"))
        + geom_point(size=point_size)
        + scale_color_gradientn(colours=_RDYLBU_REV)
        + theme_classic()
        + labs(x=df.columns[0], y=df.columns[1], colour="CytoTRACE")
    )


def plot_cyto_genes(
    result: CytoTRACE,
    n_top: int = 10,
    n_bottom: int = 10,
):
    """Bar plot of the top N cytotrace-correlated + bottom N anti-correlated genes."""
    cg = result.cyto_genes
    pos = cg.head(n_top)
    neg = cg.tail(n_bottom)
    df = pd.DataFrame(
        {
            "gene": list(pos.index) + list(neg.index),
            "correlation": list(pos.values) + list(neg.values),
        }
    )
    df["gene"] = pd.Categorical(df["gene"], categories=df["gene"], ordered=True)

    return (
        ggplot(df, aes(x="gene", y="correlation"))
        + geom_bar(stat="identity", fill="#4575B4")
        + theme_classic()
        + labs(x="Gene", y="Pearson r with CytoTRACE",
               title=f"Top {n_top} (positive) + bottom {n_bottom} (negative) genes")
    )
