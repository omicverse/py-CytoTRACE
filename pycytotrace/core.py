"""Core CytoTRACE algorithm.

References:
- Gulati, G. S. et al. *Single-cell transcriptional diversity is a hallmark
  of developmental potential.* Science 367, 405–411 (2020).
- Github: https://github.com/gunsagargulati/CytoTRACE (R reference).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.optimize import nnls


@dataclass
class CytoTRACE:
    """Container — mirrors R CytoTRACE() return list.

    Attributes:
        cytotrace: (n,) score in [0,1] — 1 = most stem-like.
        cytotrace_rank: (n,) rank-1...n.
        cyto_genes: (n_genes,) Pearson r with cytotrace, sorted descending.
        gcs: (n,) gene-counts-signature (mean of top-200 gc-correlated genes).
        gcs_genes: (n_genes,) Pearson r with gene counts, sorted descending.
        counts: (n,) number of detectably expressed genes per cell.
        gene_names: list of genes retained.
        cell_names: list of cells retained.
    """
    cytotrace: np.ndarray
    cytotrace_rank: np.ndarray
    cyto_genes: pd.Series
    gcs: np.ndarray
    gcs_genes: pd.Series
    counts: np.ndarray
    gene_names: list[str]
    cell_names: list[str]


def _rescale_01(x: np.ndarray) -> np.ndarray:
    lo, hi = x.min(), x.max()
    if hi - lo < 1e-12:
        return np.zeros_like(x)
    return (x - lo) / (hi - lo)


def _census_normalize(mat_log2: np.ndarray, counts: np.ndarray) -> np.ndarray:
    """Census: rescale each cell to its detected-gene count (Qiu 2017)."""
    xnl = (2.0 ** mat_log2) - 1.0
    cell_sums = xnl.sum(axis=0)
    cell_sums = np.maximum(cell_sums, 1e-12)
    rnorm = xnl * (counts / cell_sums)[None, :]
    return np.log2(rnorm + 1.0)


def _most_variable_genes(matn: np.ndarray, n_top: int = 1000) -> np.ndarray:
    """Top-n_top genes by dispersion (var/mean) among genes expressed in ≥5% cells."""
    n_expr = (matn > 0).sum(axis=1)
    keep = n_expr >= int(0.05 * matn.shape[1])
    if keep.sum() == 0:
        return np.arange(matn.shape[0])
    sub = matn[keep]
    means = sub.mean(axis=1)
    vars_ = sub.var(axis=1)
    disp = vars_ / np.maximum(means, 1e-12)
    order_in_sub = np.argsort(disp)[::-1]
    sub_idx = np.where(keep)[0][order_in_sub[: min(n_top, sub.shape[0])]]
    return np.sort(sub_idx)


def _similarity_matrix(mat: np.ndarray) -> np.ndarray:
    """Pearson cell-cell similarity, zero diag, threshold below mean, row-normalise."""
    # Cells are columns; we want cell-cell correlation
    Xc = mat - mat.mean(axis=0, keepdims=True)
    Xc = Xc / np.maximum(Xc.std(axis=0, keepdims=True), 1e-12)
    C = (Xc.T @ Xc) / mat.shape[0]   # cell × cell
    cutoff = float(C.mean())
    np.fill_diagonal(C, 0.0)
    C[C < 0] = 0.0
    C[C <= cutoff] = 0.0
    row_sums = C.sum(axis=1)
    safe = np.maximum(row_sums, 1e-12)
    D = C / safe[:, None]
    D[row_sums == 0, :] = 0.0
    return D


def _diffuse(S: np.ndarray, score: np.ndarray, alpha: float = 0.9,
             max_iter: int = 10000, tol: float = 1e-6) -> np.ndarray:
    """v_{t+1} = α S v_t + (1-α) v_0."""
    v = score.copy()
    v0 = score.copy()
    for _ in range(max_iter):
        v_new = alpha * S @ v + (1.0 - alpha) * v0
        diff = np.mean(np.abs(v_new - v))
        v = v_new
        if diff < tol:
            break
    return v


def cytotrace_run(
    mat: np.ndarray | pd.DataFrame,
    *,
    gene_names: list[str] | None = None,
    cell_names: list[str] | None = None,
) -> CytoTRACE:
    """End-to-end CytoTRACE pipeline.

    Args:
        mat: (n_genes × n_cells) raw count matrix.
        gene_names / cell_names: optional explicit names.
    """
    if isinstance(mat, pd.DataFrame):
        gene_names = list(mat.index)
        cell_names = list(mat.columns)
        M = mat.to_numpy(dtype=np.float64)
    else:
        M = np.asarray(mat, dtype=np.float64)
        if gene_names is None:
            gene_names = [f"gene_{i}" for i in range(M.shape[0])]
        if cell_names is None:
            cell_names = [f"cell_{i}" for i in range(M.shape[1])]

    # 1. Filter zero-variance genes
    pq = (M.sum(axis=1) == 0) | (M.var(axis=1) == 0)
    M = M[~pq]
    gene_names = [g for g, drop in zip(gene_names, pq) if not drop]

    # 2. Per-cell sequencing depth normalisation → log2
    col_sums = np.maximum(M.sum(axis=0), 1e-12)
    M_cpm = (M / col_sums) * 1e6
    M_log2 = np.log2(M_cpm + 1.0)

    # 3. Gene counts
    counts = (M_log2 > 0).sum(axis=0).astype(np.float64)
    counts = np.maximum(counts, 1.0)

    # 4. Census normalise
    M_norm = _census_normalize(M_log2, counts)

    # 5. Find most-variable genes
    mvg_idx = _most_variable_genes(M_norm, n_top=1000)

    # 6. Cell-cell similarity from MVGs
    S = _similarity_matrix(M_norm[mvg_idx])

    # 7. Gene-counts signature (mean of top 200 genes by Pearson with counts)
    Mn_c = M_norm - M_norm.mean(axis=1, keepdims=True)
    counts_c = counts - counts.mean()
    num = (Mn_c * counts_c[None, :]).sum(axis=1)
    den = np.sqrt((Mn_c ** 2).sum(axis=1) * (counts_c ** 2).sum())
    pearson_gc = num / np.maximum(den, 1e-12)
    pearson_gc = pd.Series(pearson_gc, index=gene_names).sort_values(ascending=False)
    top200 = pearson_gc.head(200).index.tolist()
    top200_idx = [gene_names.index(g) for g in top200]
    gcs = M_norm[top200_idx].mean(axis=0)

    # 8. NNLS regression of GCS onto similarity matrix
    coef, _ = nnls(S, gcs)
    gcs_regressed = S @ coef

    # 9. Diffuse
    gcs_diffused = _diffuse(S, gcs_regressed, alpha=0.9)

    # 10. Rank → [0,1]
    cytotrace_rank = pd.Series(gcs_diffused).rank().to_numpy()
    cytotrace = _rescale_01(cytotrace_rank)

    # 11. Gene-cytotrace correlations
    ct_c = cytotrace - cytotrace.mean()
    num2 = (Mn_c * ct_c[None, :]).sum(axis=1)
    den2 = np.sqrt((Mn_c ** 2).sum(axis=1) * (ct_c ** 2).sum())
    cyto_genes = pd.Series(num2 / np.maximum(den2, 1e-12),
                           index=gene_names).sort_values(ascending=False)

    return CytoTRACE(
        cytotrace=cytotrace,
        cytotrace_rank=cytotrace_rank,
        cyto_genes=cyto_genes,
        gcs=gcs,
        gcs_genes=pearson_gc,
        counts=counts,
        gene_names=gene_names,
        cell_names=cell_names,
    )
