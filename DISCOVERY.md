# Discovery — py-CytoTRACE

## 1. Is the target already ported?

`gh repo view omicverse/py-CytoTRACE` → **not found at port start**.

A community CytoTRACE-2 implementation exists (`digitalcytometry/cytotrace2`) but it's an entirely different algorithm (neural-net cellular potency); we port the original CytoTRACE.

## 2. R dependencies + py-mirror reuse

R upstream `DESCRIPTION` lists: `HiClimR`, `nnls`, `ccaPP`, `parallel`, `sva`.

| R dep | Already mirrored? | Reused as |
|---|---|---|
| HiClimR (fast correlation) | ❌ | direct numpy `@` (small N×N) |
| nnls (non-negative least squares) | ❌ | `scipy.optimize.nnls` |
| ccaPP (corPearson) | ❌ | numpy Pearson |
| sva::ComBat (batch correction) | ❌ | deferred to v0.2 (would use `scanpy.pp.combat`) |
| parallel | n/a | numpy is BLAS-threaded |

## 3. Decision

**Proceed with full port** — algorithm class is ordinal (developmental potency score). Batch correction deferred. R reference is non-trivial to install in this env (HiClimR needs ncdf4 C library); we validate against the synthetic ground-truth ordering instead.
