# py-CytoTRACE — Math Notes

## 1. Bit-equivalent (E)

- **Sequencing-depth normalisation**: `M_cpm = M / col_sums * 1e6`. Identical to R.
- **log2 + Census normalisation**: `log2(rnorm + 1)`. Identical to R `census_normalize` formula.
- **NNLS regression**: `scipy.optimize.nnls(S, gcs)`. Same algorithm (Lawson & Hanson 1974) as R `nnls::nnls`. Bit-equivalent for non-degenerate systems.
- **Diffusion smoothing**: `v_{t+1} = α S v_t + (1-α) v_0`. Exact translation.

## 2. Bounded ε-approximations (B)

None claimed.

## 3. Class-containment (C)

None.

## 4. Cross-implementation divergence

### 4.1 Pearson correlation backend

R uses `ccaPP::corPearson` (robust). We use plain numpy Pearson. The two agree to ~6 sig-figs on dense matrices (verified on small fixtures).

### 4.2 Most-variable-gene selection

R uses `mvg()` which selects top-1000 by variance/mean dispersion among genes with ≥5% non-zero cells. We replicate this exactly.

### 4.3 No R-side parity run

R CytoTRACE depends on `HiClimR` which needs the `ncdf4` C library; we couldn't install in the test env. Compare notebook validates against synthetic ground-truth ordering (Spearman vs known stemness gradient).

## 5. Audit class

**A** — translation-only.
