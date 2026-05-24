"""Smoke tests for pycytotrace."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_PORT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PORT))

import pycytotrace


def test_import():
    assert pycytotrace.__version__ == "0.1.0"


def test_pipeline_runs():
    rng = np.random.RandomState(42)
    n_cells = 50
    n_genes = 200
    # Build a gradient: cell i has higher counts the smaller i is (stem-like)
    base = rng.poisson(2.0, (n_genes, n_cells)).astype(float)
    for c in range(n_cells):
        # add (n_cells - c) extra genes detected for the stem cells
        extra_genes = rng.choice(n_genes, size=max(1, n_cells - c), replace=False)
        base[extra_genes, c] += rng.poisson(5.0, size=len(extra_genes))
    res = pycytotrace.cytotrace_run(base)
    assert res.cytotrace.shape == (n_cells,)
    assert res.cytotrace.min() >= 0 and res.cytotrace.max() <= 1
    # CytoTRACE should be positively correlated with gene counts
    from scipy.stats import spearmanr
    rho = spearmanr(res.counts, res.cytotrace)[0]
    assert rho > 0.6, f"cytotrace ~ counts Spearman = {rho:.3f}"
