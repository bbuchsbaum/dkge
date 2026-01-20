# Phase 2: Fit Layer Correctness - Research

**Researched:** 2026-01-19
**Domain:** DKGE core fitting algorithm - K-orthonormality, pooled design computation, recovery testing
**Confidence:** HIGH

## Summary

Phase 2 focuses on verifying mathematical correctness of `dkge_fit()`, the core algorithm that produces K-orthonormal bases from subject-level GLM coefficients. The codebase already has substantial test infrastructure including `test-fit.R`, `test-fit-helpers.R`, `test-toy-recovery.R`, and `test-dkge-loso.R`. The key mathematical invariant to verify is K-orthonormality: `t(U) %*% K %*% U = I` within tolerance.

The existing toy simulator `dkge_sim_toy()` provides ground-truth generation for recovery tests. The algorithm follows a well-documented specification in `data-raw/algo.md` which defines the exact mathematical steps. Current tests verify helper functions but lack explicit K-orthonormality assertions for the final `dkge_fit()` output.

**Primary recommendation:** Add explicit K-orthonormality property tests to `dkge_fit()` output, create synthetic recovery tests using `dkge_sim_toy()` with known ground truth, and add determinism tests with fixed seeds.

## Standard Stack

The DKGE package uses established patterns for numerical linear algebra testing in R.

### Core Testing Infrastructure
| Component | Source | Purpose | Why Standard |
|-----------|--------|---------|--------------|
| testthat | CRAN | Unit testing framework | De facto R testing standard |
| withr::local_seed() | CRAN | Scoped seed management | Avoids global state pollution |
| dkge_sim_toy() | Internal | Ground truth generation | Produces known K-orthonormal U_true |

### Mathematical Validation Patterns
| Pattern | Implementation | Purpose |
|---------|----------------|---------|
| isSymmetric() | base R | Check matrix symmetry within tolerance |
| eigen(, symmetric=TRUE) | base R | Extract eigenvalues for PSD check |
| crossprod() | base R | Efficient Gram matrix computation |
| max(abs(X - Y)) | base R | Maximum absolute error metric |

### Existing Test Helpers
| Helper | File | Purpose |
|--------|------|---------|
| make_fit_fixture() | test-fit.R | Generate synthetic fitting inputs |
| .make_fit() | test-dkge-loso.R | Build minimal dkge object for LOSO tests |
| rel_err() | test-dkge-loso.R | Relative error metric |
| max_abs() | test-dkge-loso.R | Max absolute deviation |

## Architecture Patterns

### Existing dkge_fit() Architecture

The fit is modularized into four stages in `R/dkge-fit-core.R`:

```
dkge_fit() orchestrates:
  1. .dkge_fit_prepare()   -- Harmonize inputs, compute ruler R, row-standardize betas
  2. .dkge_fit_accumulate() -- Compute Chat in K-metric with subject weights
  3. .dkge_fit_solve()      -- Eigen-decomposition, produce U_hat and U
  4. .dkge_fit_assemble()   -- Build final dkge object with all metadata
```

This modular architecture enables targeted testing of each stage.

### Mathematical Specification (from algo.md)

**Algorithm 1 - DKGE fit (tiny qxq):**
```
1. G_pool <- sum_s X_s^T X_s;  R <- chol(G_pool)
2. For each s: B_tilde_s <- R^T B_s
3. Eigendecompose K -> K^{+-1/2}
4. Accumulate: Chat <- sum_s w_s K^{1/2} (B_tilde_s Omega_s B_tilde_s^T) K^{1/2}
5. Eig(Chat) -> V Lambda V^T
6. U <- K^{-1/2} V[:,1:r]

Key invariant: U^T K U = I_r (K-orthonormality)
```

### Test Pattern: Property-Based Testing

Based on successful patterns in `test-kernel-invariants.R` and `test-dkge-loso.R`:

```r
test_that("dkge_fit output U satisfies K-orthonormality", {
  # Create synthetic data
  fixture <- make_fit_fixture()
  fit <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = 3)

  # K-orthonormality: U^T K U = I_r
  UtKU <- t(fit$U) %*% fit$K %*% fit$U
  expect_equal(UtKU, diag(fit$rank), tolerance = 1e-8)
})
```

### Test Pattern: Recovery Test with Known Ground Truth

Using `dkge_sim_toy()` which produces K-orthonormal `U_true`:

```r
test_that("dkge_fit recovers planted components in high-SNR regime", {
  withr::local_seed(123)
  factors <- list(A = list(L = 3, type = "nominal"))
  sim <- dkge_sim_toy(factors, active_terms = "A", S = 5, P = 50, snr = 20)

  # U_true from sim is already K-orthonormal
  fit <- dkge_fit(sim$B_list, sim$X_list, K = sim$K, rank = 1)

  # Measure subspace agreement via principal angle cosines
  cos_angle <- dkge_cosines_K(sim$U_true[, 1, drop=FALSE],
                              fit$U[, 1, drop=FALSE], sim$K)
  expect_gt(cos_angle[1], 0.95)
})
```

### Anti-Patterns to Avoid

- **Using set.seed() globally:** Pollutes global state; use `withr::local_seed()` instead
- **Hardcoded tolerance without justification:** Use 1e-8 (per Phase 1 decision) with explanation
- **Testing only happy path:** Include edge cases (rank 1, rank = q, near-singular K)
- **Testing final output without intermediate verification:** Test each stage where possible

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| K-orthonormalization | Manual Gram-Schmidt | `dkge_k_orthonormalize()` | Handles K-metric correctly |
| Ground truth generation | Manual U construction | `dkge_sim_toy()` | Guarantees K-orthonormality |
| Subspace comparison | Raw matrix difference | `dkge_cosines_K()` | Handles sign ambiguity, K-metric |
| Symmetric eigen | `eigen()` without `symmetric=TRUE` | `eigen(, symmetric=TRUE)` | Numerical stability |

**Key insight:** K-orthonormality testing must use the K-metric, not Euclidean. The existing `dkge_cosines_K()` function handles this correctly.

## Common Pitfalls

### Pitfall 1: Sign/Rotation Ambiguity in Eigenvectors
**What goes wrong:** Comparing U directly to U_true fails due to arbitrary sign flips
**Why it happens:** Eigenvectors are only defined up to sign (and rotation for repeated eigenvalues)
**How to avoid:** Use `dkge_cosines_K()` which computes principal angles (sign-invariant)
**Warning signs:** Test passes with one seed but fails with another

### Pitfall 2: Forgetting K-Metric in Orthonormality Check
**What goes wrong:** Testing `t(U) %*% U = I` instead of `t(U) %*% K %*% U = I`
**Why it happens:** Euclidean orthonormality is more familiar
**How to avoid:** Always include K in the Gram matrix: `t(U) %*% K %*% U`
**Warning signs:** Orthonormality check passes but LOSO contrasts are wrong

### Pitfall 3: Tolerance Too Tight for Numerical Stability
**What goes wrong:** Tests flaky due to accumulating floating-point errors
**Why it happens:** Multiple matrix multiplications compound errors
**How to avoid:** Use 1e-8 tolerance (established in Phase 1), document choice
**Warning signs:** Tests pass locally but fail in CI or on different platforms

### Pitfall 4: Not Verifying Pooled Design Aggregation
**What goes wrong:** Silently incorrect R when subject designs differ in rank
**Why it happens:** Edge cases in Gram matrix aggregation
**How to avoid:** Test that `t(R) %*% R = sum_s X_s^T X_s` explicitly
**Warning signs:** Weights or Chat values are NaN/Inf

### Pitfall 5: Testing Only Full Rank Cases
**What goes wrong:** Missing edge cases where rank < q or eigen-gaps are small
**Why it happens:** Most synthetic data is full-rank
**How to avoid:** Include tests with rank = 1, rank = q, and near-rank-deficient scenarios
**Warning signs:** Production data with small eigen-gaps produces unexpected results

## Code Examples

Verified patterns from existing codebase:

### K-Orthonormality Check (from test-dkge-loso.R)
```r
# Source: tests/testthat/test-dkge-loso.R lines 119-120
G <- t(out$basis) %*% fit$K %*% out$basis
expect_lt(max_abs(G - diag(r)), 5e-10)
```

### Pooled Gram Verification (from test-fit.R)
```r
# Source: tests/testthat/test-fit.R lines 18-20
ruler <- dkge:::.dkge_compute_shared_ruler(fixture$designs)
expect_equal(ruler$G_pool, Reduce(`+`, lapply(fixture$designs, crossprod)))
expect_equal(t(ruler$R) %*% ruler$R, ruler$G_pool, tolerance = 1e-10)
```

### Recovery Test Pattern (from test-toy-recovery.R)
```r
# Source: tests/testthat/test-toy-recovery.R
sim <- dkge_sim_toy(factors, active_terms = "A", S = 3, P = 12, snr = 10)
fit <- dkge_fit(sim$B_list, designs = sim$X_list, K = sim$K, rank = 1)

# Compare subspaces using K-metric cosines
cos1 <- dkge_cosines_K(sim$U_true[, 1, drop = FALSE],
                       fit$U[, 1, drop = FALSE], sim$K)
expect_gt(cos1[1], 0.98)
```

### U_true K-Orthonormality Verification (from test-toy-generator.R)
```r
# Source: tests/testthat/test-toy-generator.R lines 18-19
UtKU <- t(sim$U_true) %*% sim$K %*% sim$U_true
expect_equal(UtKU, diag(ncol(sim$U_true)), tolerance = 1e-8)
```

### Determinism Test Pattern
```r
# Pattern for verifying deterministic output
test_that("dkge_fit is deterministic with same seed", {
  withr::local_seed(42)
  fixture <- make_fit_fixture(seed = 100)
  fit1 <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = 2)

  withr::local_seed(42)
  fixture <- make_fit_fixture(seed = 100)
  fit2 <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = 2)

  expect_equal(fit1$U, fit2$U)
  expect_equal(fit1$Chat, fit2$Chat)
  expect_equal(fit1$R, fit2$R)
})
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Manual Gram-Schmidt | `dkge_k_orthonormalize()` | Current | Handles K-metric correctly |
| Testing U directly | Using `dkge_cosines_K()` | Current | Sign-ambiguity-aware |
| Global set.seed() | withr::local_seed() | Phase 1 | No side effects |

**Current test coverage gaps for Phase 2:**
- `test-fit.R`: Tests helpers and basic fit, but no explicit K-orthonormality assertion on final U
- `test-fit-helpers.R`: Tests `.dkge_fit_solve` "produces orthonormal basis" - but uses standard orthonormality, not K-orthonormality
- `test-toy-recovery.R`: Tests recovery but only for rank-1 case with single factor

## Open Questions

Things that need validation during implementation:

1. **Edge case: rank = q (full rank)**
   - What we know: Algorithm should work, eigenvalues become full spectrum
   - What's unclear: Whether Chat accumulation is stable for very high rank
   - Recommendation: Add explicit test case with rank = q

2. **Tolerance calibration for high-dimensional K**
   - What we know: 1e-8 works for small q (~10)
   - What's unclear: Appropriate tolerance for q ~ 100
   - Recommendation: Test with q = 50-100 to verify tolerance holds

3. **Determinism across solver methods**
   - What we know: "pooled" solver should be deterministic
   - What's unclear: "jd" solver behavior with random initialization
   - Recommendation: Test determinism for both solvers separately

## Sources

### Primary (HIGH confidence)
- `data-raw/algo.md` - Mathematical specification (read in full)
- `R/dkge-fit.R` - Implementation (read in full, lines 1-313)
- `R/dkge-fit-core.R` - Core implementation stages (read in full, lines 1-497)
- `tests/testthat/test-fit.R` - Existing fit tests (read in full)
- `tests/testthat/test-fit-helpers.R` - Helper tests (read in full)
- `tests/testthat/test-dkge-loso.R` - LOSO tests with K-orthonormality checks (read in full)
- `tests/testthat/test-toy-recovery.R` - Recovery test pattern (read in full)
- `tests/testthat/test-toy-generator.R` - Toy simulator tests (read in full)
- `R/dkge-sim.R` - dkge_sim_toy() implementation (read in full)

### Secondary (MEDIUM confidence)
- `tests/testthat/test-kernel-invariants.R` - Property test patterns from Phase 1
- `tests/testthat/test-data-validation.R` - Edge case test patterns from Phase 1

## Metadata

**Confidence breakdown:**
- K-orthonormality testing: HIGH - Clear mathematical definition, existing code patterns
- Recovery testing: HIGH - `dkge_sim_toy()` well-documented with K-orthonormal U_true
- Determinism testing: HIGH - Standard pattern, no external dependencies
- Tolerance calibration: MEDIUM - 1e-8 established but not validated for large q

**Research date:** 2026-01-19
**Valid until:** 30 days (stable mathematical foundations, unlikely to change)

## Recommended Plan Structure

Based on research findings, Phase 2 should have 3-4 tasks:

### Task 1: Add K-orthonormality property tests to test-fit.R
- Add explicit `t(U) %*% K %*% U = I` assertion for `dkge_fit()` output
- Test across multiple kernel types (identity, ordinal RBF, multi-factor)
- Test rank = 1, rank = 3, rank = q edge cases

### Task 2: Add pooled design verification tests
- Verify `t(R) %*% R = sum_s X_s^T X_s` for various subject counts
- Test with heterogeneous subject designs (different T_s values)
- Verify numerical stability with ill-conditioned designs (add ridge)

### Task 3: Create recovery tests with known ground truth
- Use `dkge_sim_toy()` to generate ground truth
- Test recovery at multiple SNR levels (high: 20, medium: 8, low: 2)
- Verify reconstruction error decreases with SNR

### Task 4: Add determinism tests
- Same seed produces identical U, Chat, R, weights
- Test both "pooled" and "jd" solvers (jd may need special handling)
- Document any non-deterministic components
