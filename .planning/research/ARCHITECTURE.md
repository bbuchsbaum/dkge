# Architecture Patterns: Testing Scientific/Numerical R Packages

**Domain:** Testing structure for scientific R packages (specifically dkge)
**Researched:** 2026-01-19
**Confidence:** HIGH (based on analysis of existing dkge test suite + R package testing best practices)

## Recommended Test Architecture

The dkge package already demonstrates sophisticated test organization. This document synthesizes the existing patterns and recommends a coherent testing architecture for comprehensive coverage.

### Layered Test Hierarchy

Scientific R packages benefit from organizing tests by computational layer, with dependencies flowing upward:

```
Layer 7: Integration/Pipeline Tests
         |
Layer 6: Classification Tests (dkge_classify, CV)
         |
Layer 5: Inference Tests (permutation, sign-flip)
         |
Layer 4: Transport Tests (Sinkhorn, medoid mapping)
         |
Layer 3: Contrast/Cross-Fitting Tests (LOSO, K-fold, analytic)
         |
Layer 2: Fit Layer Tests (dkge_fit, eigendecomposition)
         |
Layer 1: Kernel Layer Tests (design_kernel, kernel_roots)
         |
Layer 0: Data Layer Tests (dkge_subject, dkge_data, alignment)
```

### Component Boundaries

| Component | Responsibility | Tests With |
|-----------|---------------|------------|
| Data constructors | Input validation, effect alignment | Raw matrices, edge cases |
| Kernel functions | K construction, PSD enforcement, roots | Mathematical invariants |
| Fit core | Eigendecomposition, K-orthonormal basis U | Reconstruction, orthonormality |
| Contrast engine | Cross-fitting (LOSO/K-fold), analytic | Equivariance, unbiasedness |
| Transport | Sinkhorn OT, cluster alignment | Doubly stochastic plans |
| Inference | Permutation tests, FWER control | Null calibration, p-value uniformity |
| Classification | CV training, metrics | Label consistency, metric bounds |
| Pipeline | End-to-end orchestration | Integration scenarios |

### Data Flow for Test Fixtures

```
helper-toy.R           helper-fit-fixture.R         helper-weights-null.R
      |                        |                            |
      v                        v                            v
toy_kernel_info()      make_small_fit()            null_betas()
toy_betas()            S, q, P, T params           dkge_test_reliability_prior()
toy_fold_fit()         dkge_fit() wrapper          dkge_test_loso_delta_pvalue()
toy_real_fit()                                     dkge_test_null_uniformity()
      |                        |                            |
      +------------------------+----------------------------+
                               |
                               v
                    Test files use fixtures
                    with known structure
```

## Patterns from Existing dkge Tests

### Pattern 1: Mathematical Invariant Testing

Tests verify mathematical properties that must hold regardless of input data.

**Example from test-design-kernel.R:**
```r
test_that("kernel roots reconstruct original kernel", {
  K <- matrix(c(2, 0.5, 0.5, 1.5), 2, 2)
  roots <- dkge:::.dkge_kernel_roots(K)
  expect_equal(roots$Khalf %*% roots$Khalf, K, tolerance = 1e-10)
  expect_equal(roots$Khalf %*% roots$Kihalf, diag(2), tolerance = 1e-10)
})
```

**Pattern structure:**
1. Construct matrix with known properties
2. Apply function under test
3. Verify mathematical identity (reconstruction, orthogonality, etc.)
4. Use appropriate tolerance (typically 1e-8 to 1e-10)

**Key invariants for dkge:**

| Component | Invariant | Test Form |
|-----------|-----------|-----------|
| Kernel | Symmetry | `K == t(K)` |
| Kernel | PSD | `all(eigen(K)$values >= -tol)` |
| Kernel roots | Reconstruction | `Khalf %*% Khalf == K` |
| Kernel roots | Inverse | `Khalf %*% Kihalf == I` |
| Basis U | K-orthonormality | `t(U) %*% K %*% U == I_r` |
| Contrast | Equivariance | Permuting cluster order permutes output |
| Transport plan | Doubly stochastic | `rowSums(P) == colSums(P) == marginals` |

### Pattern 2: Cross-Fitting Validation

LOSO and K-fold tests verify that held-out subjects are excluded from basis estimation.

**Example from test-dkge-kfold.R:**
```r
test_that(".dkge_contrast_kfold: with k=S equals LOSO and bases are valid", {
  # When k = S, each subject is its own fold -> must match LOSO
  fit <- .make_fit(seed = 42, q = q, r = r, S = S, P = P)
  folds <- dkge_define_folds(fit, type = "subject", k = S, seed = 101)

  kres <- .dkge_contrast_kfold(fit, clist, folds = folds, ...)

  # LOSO baseline
  for (s in seq_len(S)) {
    out <- dkge_loso_contrast(fit, s, cvec, ridge = 0)
    loso_vs[[s]] <- out$v
    loso_Us[[s]] <- out$basis
  }

  # Values must match
  for (s in seq_len(S)) {
    expect_lt(rel_err(kres$values$c1[[s]], loso_vs[[s]]), 1e-12)
  }

  # Bases must be K-orthonormal
  for (s in seq_len(S)) {
    U_fold <- kres$metadata$fold_bases[[fidx]]
    G <- t(U_fold) %*% fit$K %*% U_fold
    expect_lt(max_abs(G - diag(r)), 5e-10)
  }
})
```

**Pattern structure:**
1. Create fit with known parameters
2. Run cross-fitting with k = S (degenerates to LOSO)
3. Compare against explicit LOSO implementation
4. Verify all bases satisfy K-orthonormality

### Pattern 3: Equivariance Testing

Tests that permutations of input produce corresponding permutations of output.

**Example from test-dkge-kfold.R:**
```r
test_that("permuting B_s columns only permutes that subject's v_s; bases unchanged", {
  # Baseline
  base <- .dkge_contrast_kfold(fit, clist, folds = folds, ...)

  # Permute columns of one subject's B matrix
  Pmat <- diag(P)[, sample(P), drop = FALSE]
  fit2 <- fit
  fit2$Btil[[s0]] <- fit$Btil[[s0]] %*% Pmat

  # Rerun
  alt <- .dkge_contrast_kfold(fit2, clist, folds = folds, ...)

  # Bases unchanged (computed on training set)
  expect_lt(max_abs(alt$metadata$fold_bases[[fidx]] -
                    base$metadata$fold_bases[[fidx]]), 1e-12)

  # Values permuted accordingly
  expect_lt(max_abs(v_alt - as.numeric(t(Pmat) %*% v_base)), 1e-12)

  # Other subjects unaffected
  for (s in others) {
    expect_lt(rel_err(alt$values$c1[[s]], base$values$c1[[s]]), 1e-12)
  }
})
```

**Pattern structure:**
1. Run baseline computation
2. Apply known permutation to input
3. Verify output transforms predictably
4. Verify unrelated outputs unchanged

### Pattern 4: Fallback and Edge Case Testing

Tests that verify graceful degradation when analytic methods fail.

**Example from test-analytic.R:**
```r
test_that("dkge_analytic_loso falls back to exact LOSO when tolerance is large", {
  fit <- dkge_fit(toy$betas, toy$designs, K = toy$K, rank = 2)

  exact <- dkge_loso_contrast(fit, s = 2, contrast)
  fb <- dkge_analytic_loso(fit, s = 2, contrast, tol = 10, fallback = TRUE)

  expect_identical(fb$method, "fallback")
  expect_equal(fb$v, exact$v)
  expect_equal(fb$basis, exact$basis)
})

test_that("dkge_contrast analytic metadata records fallback usage", {
  res <- dkge_contrast(fit, contrast, method = "analytic", tol = 10, fallback = TRUE)

  expect_true(all(res$metadata$fallback_rates >= 0))
  expect_equal(res$metadata$fallback_rates[[1]], 1)
  expect_true(all(res$metadata$methods[[1]] == "fallback"))
})
```

**Pattern structure:**
1. Create conditions that trigger fallback
2. Verify fallback produces correct results
3. Verify metadata records fallback usage

### Pattern 5: Null Calibration Testing

Statistical tests verify that p-values are uniform under the null hypothesis.

**Example from helper-weights-null.R:**
```r
dkge_test_null_uniformity <- function(nrep = 30,
                                      nsub = 12,
                                      V = 120,
                                      adapt = c("kenergy_prec", "kenergy", "none"),
                                      n_perm = 199,
                                      seed = 123) {
  pvals <- numeric(nrep)
  for (r in seq_len(nrep)) {
    rng <- rng + 1L
    B_list <- null_betas(nsub = nsub, V = V, seed = rng)
    pvals[r] <- dkge_test_loso_delta_pvalue(B_list, adapt = adapt, ...)
  }
  pvals
}
```

**Pattern structure:**
1. Generate many replicates under the null (no signal)
2. Compute p-values for each replicate
3. Test uniformity (e.g., Kolmogorov-Smirnov test)
4. Mark as `skip_on_cran()` due to computational cost

### Pattern 6: Deterministic Transport Testing

Optimal transport tests with known ground-truth solutions.

**Example from test-transport.R:**
```r
test_that("sinkhorn transport produces deterministic operators", {
  dat <- make_deterministic_sinkhorn_inputs()  # Known permutation structure
  res <- dkge_transport_to_medoid_sinkhorn(
    dat$v_list, dat$A_list, dat$centroids,
    medoid = 1, epsilon = 1e-4, max_iter = 2000, tol = 1e-9
  )

  # Identity for medoid subject
  expect_equal(res$subj_values[1, ], dat$v_list[[1]])
  expect_equal(res$plans[[1]], diag(1, 3))

  # Known permutation for subject 2
  expected_plan_subject2 <- matrix(
    c(0, 0, 1/3,
      1/3, 0, 0,
      0, 1/3, 0), nrow = 3, byrow = TRUE
  )
  expect_equal(res$plans[[2]], expected_plan_subject2, tolerance = 1e-6)
})
```

**Pattern structure:**
1. Construct inputs with known optimal transport solution
2. Run transport algorithm
3. Compare plans against expected matrices
4. Use small epsilon for near-deterministic plans

## Anti-Patterns to Avoid

### Anti-Pattern 1: Testing Implementation, Not Behavior

**Bad:**
```r
test_that("internal variable is set", {
  fit <- dkge_fit(betas, designs, K = K, rank = 2)
  expect_true(exists("temp_var", envir = fit$.internal))  # Testing internals
})
```

**Good:**
```r
test_that("fit produces K-orthonormal basis", {
  fit <- dkge_fit(betas, designs, K = K, rank = 2)
  gram <- t(fit$U) %*% K %*% fit$U
  expect_equal(gram, diag(fit$rank), tolerance = 1e-10)
})
```

### Anti-Pattern 2: Ignoring Numerical Tolerance

**Bad:**
```r
expect_equal(reconstructed, original)  # May fail due to floating point
```

**Good:**
```r
expect_equal(reconstructed, original, tolerance = 1e-10)
```

### Anti-Pattern 3: Overfitting Tests to Specific Seeds

**Bad:**
```r
test_that("random test passes", {
  set.seed(42)  # Only this seed works
  result <- stochastic_function()
  expect_true(result > 0)  # Fails with other seeds
})
```

**Good:**
```r
test_that("invariant holds across random seeds", {
  for (seed in c(1, 42, 123, 999)) {
    set.seed(seed)
    result <- stochastic_function()
    expect_true(invariant_property(result))
  }
})
```

### Anti-Pattern 4: Large Fixtures Without skip_on_cran()

**Bad:**
```r
test_that("large scale test", {
  betas <- make_huge_betas(nsub = 1000, V = 50000)  # Slow
  # ... test code
})
```

**Good:**
```r
test_that("large scale test", {
  skip_on_cran()  # Only run locally or on CI
  betas <- make_huge_betas(nsub = 1000, V = 50000)
  # ... test code
})
```

## Suggested Testing Order

Build tests from the foundation up, ensuring each layer works before testing layers that depend on it.

### Phase 1: Data Layer (Foundation)

**Priority:** HIGH - All other tests depend on correct data handling

Tests to implement:
1. `dkge_subject()` - effect name alignment, omega validation
2. `dkge_data()` - partial effect overlap handling, provenance tracking
3. Effect alignment across subjects with missing effects

**Invariants:**
- Row names of beta match column names of design after alignment
- Provenance accurately records which effects are observed per subject

### Phase 2: Kernel Layer

**Priority:** HIGH - Fit layer requires valid kernels

Tests to implement:
1. `design_kernel()` - PSD enforcement, normalization
2. `kernel_roots()` - Khalf/Kihalf reconstruction
3. Edge cases: rank-deficient kernels, near-zero eigenvalues

**Invariants:**
- K is symmetric
- K is PSD (all eigenvalues >= 0)
- Khalf @ Khalf = K
- Khalf @ Kihalf = I

### Phase 3: Fit Layer

**Priority:** HIGH - Core algorithm

Tests to implement:
1. `dkge_fit()` - eigendecomposition, basis extraction
2. K-orthonormality of U
3. Pooled design Cholesky correctness
4. Subject weight computation

**Invariants:**
- t(U) @ K @ U = I_r
- Reconstruction of Chat from contribs

### Phase 4: Contrast/Cross-Fitting Layer

**Priority:** HIGH - Primary output of package

Tests to implement:
1. LOSO correctness - held-out subject excluded from basis
2. K-fold equivalence when k = S
3. Analytic approximation accuracy and fallback
4. Equivariance under cluster permutation

**Invariants:**
- K-orthonormality of all leave-out bases
- Unbiasedness (no data leakage)

### Phase 5: Transport Layer

**Priority:** MEDIUM - Required for mismatched cluster counts

Tests to implement:
1. Sinkhorn convergence
2. Doubly stochastic constraint satisfaction
3. Deterministic cases with known solutions
4. C++ vs R implementation agreement

**Invariants:**
- Transport plans are doubly stochastic
- Medoid subject maps to itself

### Phase 6: Inference Layer

**Priority:** MEDIUM - Statistical validity

Tests to implement:
1. Sign-flip permutation mechanics
2. Null calibration (p-value uniformity)
3. FWER control

**Invariants:**
- Under null, p-values are uniform
- False positive rate controlled at alpha

### Phase 7: Classification Layer

**Priority:** MEDIUM - Higher-level analysis

Tests to implement:
1. CV fold construction
2. Label consistency
3. Metric bounds (accuracy in [0,1], etc.)
4. Fallback when glmnet unavailable

**Invariants:**
- Predictions only use training data
- Probabilities sum to 1

### Phase 8: Integration/Pipeline Tests

**Priority:** LOW - Run after component tests pass

Tests to implement:
1. `dkge_pipeline()` end-to-end
2. Multi-step workflows
3. Error propagation

## Test File Organization Recommendation

```
tests/testthat/
  |
  +-- helper-toy.R              # Shared fixtures for all tests
  +-- helper-fit-fixture.R      # Fit-specific fixtures
  +-- helper-weights-null.R     # Null calibration fixtures
  +-- helper-transport.R        # Transport fixtures
  |
  +-- test-data.R               # Layer 0: dkge_subject, dkge_data
  +-- test-design-kernel.R      # Layer 1: design_kernel, kernel_roots
  +-- test-fit.R                # Layer 2: dkge_fit core
  +-- test-fit-helpers.R        # Layer 2: Internal fit helpers
  +-- test-contrast.R           # Layer 3: dkge_contrast
  +-- test-dkge-loso.R          # Layer 3: LOSO specifics
  +-- test-dkge-kfold.R         # Layer 3: K-fold specifics
  +-- test-analytic.R           # Layer 3: Analytic approximation
  +-- test-transport.R          # Layer 4: Sinkhorn OT
  +-- test-inference.R          # Layer 5: Permutation tests
  +-- test-classify.R           # Layer 6: Classification
  +-- test-pipeline.R           # Layer 7: Integration
```

## Scalability Considerations

| Test Type | Local | CI | CRAN |
|-----------|-------|-----|------|
| Unit tests (fast) | Always | Always | Always |
| Mathematical invariants | Always | Always | Always |
| Cross-fitting (small) | Always | Always | Always |
| Null calibration | Always | Always | `skip_on_cran()` |
| Large-scale stress | Manual | Nightly | Skip |
| Benchmark comparisons | Manual | Weekly | Skip |

## Build Order Implications for Roadmap

Based on the layered architecture, test development should proceed:

1. **Phase 1: Data + Kernel Tests** - Foundation must be solid
   - Complete test-data.R coverage
   - Complete test-design-kernel.R coverage
   - Estimated: 1-2 days

2. **Phase 2: Fit Layer Tests** - Core algorithm
   - Complete test-fit.R coverage
   - Add mathematical invariant tests
   - Estimated: 2-3 days

3. **Phase 3: Contrast Tests** - Primary functionality
   - Unify LOSO/K-fold/analytic testing
   - Add equivariance tests
   - Estimated: 3-4 days

4. **Phase 4: Transport + Inference** - Can parallelize
   - Transport: deterministic cases
   - Inference: null calibration
   - Estimated: 2-3 days

5. **Phase 5: Classification + Integration** - Final layer
   - Classification metrics
   - Pipeline integration
   - Estimated: 2-3 days

**Total estimated effort:** 10-15 days for comprehensive coverage

## Sources

- [R Packages (2e) - Testing Basics](https://r-pkgs.org/testing-basics.html)
- [R Packages (2e) - Testing Design](https://r-pkgs.org/testing-design.html)
- [R Packages (2e) - Advanced Testing](https://r-pkgs.org/testing-advanced.html)
- [testthat package documentation](https://testthat.r-lib.org/)
- [The 4 Layers of Testing Every R Package Needs](https://www.r-bloggers.com/2025/07/the-4-layers-of-testing-every-r-package-needs/)
- [Organizing tests in R packages - cynkra](https://cynkra.com/blog/2025-03-04-refactoring-test-files/)
- Analysis of existing dkge test suite (60+ test files examined)
