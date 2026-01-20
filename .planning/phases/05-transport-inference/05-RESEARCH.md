# Phase 5: Transport + Inference - Research

**Researched:** 2026-01-20
**Domain:** Optimal transport (Sinkhorn), permutation inference, parallel execution
**Confidence:** HIGH

## Summary

Phase 5 verifies two related but distinct components: (1) the Sinkhorn optimal transport algorithm used to map cluster values between subject-specific parcellations and a medoid reference, and (2) permutation-based inference methods (sign-flip maxT) for calibrated p-values.

The implementation already exists in the core R/ directory with both R and C++ backends. The C++ implementation (`src/sinkhorn.cpp`) uses log-domain stabilized Sinkhorn-Knopp iterations with OpenMP parallelization. The R wrapper (`R/dkge-transport.R`) provides caching for warm starts and cost matrix construction. Inference (`R/dkge-inference.R`) implements sign-flip max-T testing with FWER control.

Existing tests (`test-transport-sinkhorn.R`, `test-transport.R`, `test-inference.R`) cover basic functionality but lack: (1) systematic doubly-stochastic verification, (2) convergence failure edge cases, (3) null distribution uniformity checks, and (4) parallel/sequential equivalence tests.

**Primary recommendation:** Create two focused plans - one for Sinkhorn mathematical properties and convergence, one for inference null calibration and parallel equivalence.

## Implementation Locations

### Core Implementation (R/)

| File | Purpose | Exported Functions |
|------|---------|-------------------|
| `R/dkge-transport.R` | Transport to medoid, cost matrix, caching | `dkge_transport_to_medoid_sinkhorn`, `dkge_transport_to_medoid_sinkhorn_cpp`, `dkge_transport_contrasts_to_medoid`, `dkge_prepare_transport`, `dkge_clear_sinkhorn_cache` |
| `R/dkge-inference.R` | Sign-flip max-T, unified inference | `dkge_signflip_maxT`, `dkge_infer`, `dkge_freedman_lane` |
| `R/dkge-mapper.R` | Pluggable mapping strategies | `dkge_mapper_spec`, `fit_mapper`, `predict_mapper`, `apply_mapper` |
| `src/sinkhorn.cpp` | C++ Sinkhorn with OpenMP | `sinkhorn_plan_cpp` (internal) |

### Experimental (future/)

| File | Purpose | Status |
|------|---------|--------|
| `future/dkge-sinkhorn.R` | Pure R Sinkhorn implementation | Not exported, reference only |
| `future/dkge-sinkhorn-cpp.R` | Alternative C++ wrapper | Not exported |

**Key finding:** All transport and inference functions are in **core R/** and are exported. The `future/` directory contains experimental/reference implementations that are NOT part of the package API.

### Existing Tests

| Test File | Coverage |
|-----------|----------|
| `test-transport-sinkhorn.R` | Cross-validation vs T4transport, spatial penalty bias, cache clearing |
| `test-transport.R` | Deterministic operators, R/CPP equivalence, cost matrix penalty |
| `test-inference.R` | Structure validation, transport integration |
| `test-inference-transport.R` | Error handling for mismatched clusters |

## Sinkhorn Algorithm Analysis

### Algorithm Structure

The C++ implementation in `src/sinkhorn.cpp` uses:

```cpp
// Log-domain Sinkhorn-Knopp iterations
// P(i,j) = exp(log_u[i] + logK[i,j] + log_v[j])

for (int it = 0; it < max_iter; ++it) {
  // Update log_u: log_u = log_mu - logsumexp(logK + log_v, axis=1)
  // Update log_v: log_v = log_nu - logsumexp(logK + log_u, axis=0)

  // Check convergence every 5 iterations
  if ((it % 5 == 4) || (it == max_iter - 1)) {
    // Compute plan and marginal errors
    if (max(err1, err2) < tol) break;
  }
}
```

### Parameters

| Parameter | Default | Purpose | Source |
|-----------|---------|---------|--------|
| `epsilon` | 0.05 | Entropic regularization | `dkge_transport_to_medoid_sinkhorn` |
| `max_iter` | 200-300 | Maximum iterations | Varies by function |
| `tol` | 1e-6 | Marginal error tolerance | Per CONTEXT.md decision |
| `lambda_emb` | 1 | Embedding cost weight | Cost matrix construction |
| `lambda_spa` | 0.5 | Spatial cost weight | Cost matrix construction |
| `sigma_mm` | 15 | Spatial rescaling (mm) | Cost matrix construction |

### Convergence Criteria

From `src/sinkhorn.cpp:152-154`:
```cpp
double err1 = arma::abs(row_sums - mu_vec).max();
double err2 = arma::abs(col_sums - nu_vec).max();
if (std::max(err1, err2) < tol) break;
```

**Convergence means:** Maximum absolute marginal error < tolerance (1e-6)

### Warm-Start Caching

The R wrapper (`R/dkge-transport.R:50-108`) implements caching:
- Cache key: dimensions + epsilon + cost matrix fingerprint
- Stores: `log_u`, `log_v`, iterations
- LRU eviction with 64 entry limit

## Mathematical Properties to Test

### Doubly Stochastic (Critical)

Per CONTEXT.md, verify at tolerance 1e-6:

```r
# Row sums = 1
expect_equal(rowSums(plan), mu, tolerance = 1e-6)

# Column sums = 1
expect_equal(colSums(plan), nu, tolerance = 1e-6)

# Non-negative entries
expect_true(all(plan >= 0))
```

### Identity/Deterministic Cases

Per algo.md Algorithm 3:
```
If s = m: y_s <- v_s; continue  # Identity mapping for medoid
```

Test that medoid subject maps to itself exactly.

### Cost-Optimal Transport

For identical source/target distributions with small epsilon, transport should concentrate on diagonal (low cost).

### Marginal Preservation

Transport plan T couples marginals mu and nu:
```
T %*% ones_m = mu  (row marginal)
t(T) %*% ones_n = nu  (column marginal)
```

## Inference Methods Analysis

### Sign-Flip Max-T

From `R/dkge-inference.R:18-54`:

```r
dkge_signflip_maxT <- function(Y, B = 2000, center, tail) {
  # Y: S x Q matrix (subjects x clusters)
  # Observed t-statistics per cluster
  t_obs <- colMeans(Y) / (apply(Y, 2, sd) / sqrt(S) + 1e-12)

  # Permutation null: random sign flips
  for (b in 1:B) {
    Yb <- flips[,b] * Y  # Element-wise sign flip
    t_b <- colMeans(Yb) / (apply(Yb, 2, sd) / sqrt(S) + 1e-12)
    maxnull[b] <- max(abs(t_b))
  }

  # p-values with max-T correction (FWER)
  p <- sapply(abs(t_obs), function(x) (1 + sum(maxnull >= x)) / (B + 1))
}
```

### Null Calibration Requirement

Per algo.md section 9:
> sign-flip max-T on transported subject maps (FWER across clusters)

Per CONTEXT.md:
- 500 permutations for tests (balances precision vs runtime)
- Shuffled labels approach (break true association in structured data)
- Single seed for reproducibility

### Expected Null Behavior

Under the null hypothesis:
- p-values should be uniformly distributed on [0,1]
- Chi-square goodness-of-fit test should pass
- Histogram of p-values should be approximately flat

## Parallel Execution Analysis

### C++ Level (OpenMP)

From `src/sinkhorn.cpp`:
```cpp
#pragma omp parallel for schedule(static)
for (int i = 0; i < n; ++i) {
  // log_u updates
}

#pragma omp parallel for schedule(static)
for (int j = 0; j < m; ++j) {
  // log_v updates
}
```

OpenMP parallelization is automatic when compiled with support.

### R Level (future.apply)

From CONTEXT.md success criteria:
> `parallel=TRUE` produces identical results to sequential execution

The `.dkge_apply()` helper (`R/dkge-utils.R:30-38`) wraps parallel execution:
```r
.dkge_apply <- function(X, FUN, parallel = FALSE, ...) {
  if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE))
      stop("parallel=TRUE requires future.apply package")
    future.apply::future_lapply(X, FUN, ...)
  } else {
    lapply(X, FUN, ...)
  }
}
```

### Parallel Testing Strategy

Per `.planning/research/PITFALLS.md:384-395`:
```r
test_that("parallel and sequential give same results", {
  future::plan(future::sequential)
  result_seq <- dkge_contrast(fit, c, parallel = TRUE)
  future::plan(future::multisession, workers = 2)
  result_par <- dkge_contrast(fit, c, parallel = TRUE)
  expect_equal(result_seq$values, result_par$values)
})
```

## Convergence Edge Cases

### Potential Failure Modes

1. **Very small epsilon**: Near-deterministic transport, numerical underflow
2. **Mismatched mass totals**: Algorithm checks `abs(sum(mu) - sum(nu)) < 1e-6`
3. **Zero/negative weights**: Validation enforces positive weights
4. **Ill-conditioned cost matrix**: Very large cost differences can cause slow convergence

### Non-Convergence Behavior

Per CONTEXT.md:
> Non-convergence behavior: warning + best result (not error)

Current implementation returns last plan without warning. **Gap identified:** Need to verify warning is issued or add test for this behavior.

## Common Pitfalls

### Pitfall 1: Doubly-Stochastic Verification Tolerance

**What goes wrong:** Using different tolerance for verification than algorithm uses internally.
**How to avoid:** Use same 1e-6 tolerance as Sinkhorn convergence criterion.
**Detection:** Test plan marginals explicitly, not just "close enough".

### Pitfall 2: Identity Transport for Medoid

**What goes wrong:** Assuming medoid gets identity plan when it may get scaled version.
**How to avoid:** Verify `plans[[medoid]] == diag(sizes_ref)` exactly.
**Current code:** `plans[[s]] <- diag(sizes_ref, nrow = Q)` - correctly uses identity coupling.

### Pitfall 3: Permutation Test Seed Sensitivity

**What goes wrong:** Test passes with one seed but fails with others.
**How to avoid:** Use single seed per CONTEXT.md decision, but verify algorithm correctness not chance.
**Strategy:** Test uniformity with enough permutations (500) that sampling variance is small.

### Pitfall 4: Parallel Non-Determinism

**What goes wrong:** OpenMP or future.apply gives different results due to floating-point ordering.
**How to avoid:** Test equivalence within tolerance, not exact equality.
**Tolerance:** Per prior phases, use 1e-8 to 1e-10 depending on operation.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Optimal transport | Custom solver | Existing `sinkhorn_plan_cpp` | Log-domain stabilization is subtle |
| Chi-square test | Manual calculation | `chisq.test()` | Standard R function |
| Uniform distribution test | Custom | `stats::ks.test()` or `chisq.test()` | Standard goodness-of-fit |
| Parallel execution | Custom threading | `future.apply` | Handles platform differences |

## Recommended Plan Structure

### Plan 05-01: Sinkhorn Convergence Tests

**Tasks:**
1. Doubly-stochastic verification for standard inputs
2. Identity transport verification for medoid subject
3. Deterministic operator tests (matching embeddings -> diagonal plan)
4. Convergence failure edge cases (small epsilon, max iterations)
5. R/CPP equivalence verification

### Plan 05-02: Null Calibration and Parallel Tests

**Tasks:**
1. Null distribution uniformity test (chi-square on p-value bins)
2. Shuffled-labels null construction
3. Parallel vs sequential equivalence for transport
4. Parallel vs sequential equivalence for inference
5. Integration: transport + inference pipeline

## Open Questions

### Resolved by CONTEXT.md

1. **Convergence tolerance:** 1e-6 (decided)
2. **Permutation count:** 500 for tests (decided)
3. **Null construction:** Shuffled labels (decided)
4. **Non-convergence behavior:** Warning + best result (decided)

### Claude's Discretion

1. **Chi-square significance threshold:** Recommend alpha=0.05 (standard)
2. **Sparsity threshold for near-diagonal:** Recommend >0.9 on diagonal entries
3. **Parallel tolerance:** Recommend 1e-10 (consistent with prior phases)
4. **P-value bin count for chi-square:** Recommend 10 bins (standard choice)

### Unresolved

1. **Warning implementation:** Current C++ code doesn't issue warning on non-convergence. Test should verify expected behavior per CONTEXT.md or document actual behavior.

## Sources

### Primary (HIGH confidence)
- `src/sinkhorn.cpp` - C++ implementation reviewed directly
- `R/dkge-transport.R` - R wrapper reviewed directly
- `R/dkge-inference.R` - Inference implementation reviewed directly
- `data-raw/algo.md` - Algorithm specification (Algorithm 3, section 9)
- `05-CONTEXT.md` - User decisions for this phase

### Secondary (MEDIUM confidence)
- Existing tests (`test-transport*.R`, `test-inference*.R`)
- `.planning/research/PITFALLS.md` - Parallel testing guidance

### Tertiary (LOW confidence)
- None - all findings from authoritative sources

## Metadata

**Confidence breakdown:**
- Implementation locations: HIGH - direct code review
- Algorithm structure: HIGH - direct code review
- Mathematical properties: HIGH - derived from algo.md
- Pitfalls: MEDIUM - based on code patterns and prior phase experience
- Parallel testing: MEDIUM - based on codebase patterns

**Research date:** 2026-01-20
**Valid until:** 30 days (stable codebase)
