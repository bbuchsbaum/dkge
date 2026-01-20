# Phase 3: Cross-Fitting Validation - Research

**Researched:** 2026-01-19
**Domain:** Leave-one-subject-out (LOSO) and K-fold cross-fitting for DKGE
**Confidence:** HIGH

## Summary

This research investigated the cross-fitting validation infrastructure in the DKGE package. The package implements three cross-fitting strategies for unbiased contrast estimation: iterative LOSO (gold standard), K-fold cross-validation (efficient), and analytic LOSO (fast approximation using eigenvalue perturbation theory).

The codebase has well-structured implementations in `R/dkge-loso.R`, `R/dkge-kfold.R`, and `R/dkge-analytic.R`, with shared fold-building logic in `R/dkge-folds.R`. The core mathematical operation is subtracting a subject's contribution from the pooled covariance (`Chat_minus = Chat - w_s * S_s`) before re-computing the basis via eigen-decomposition.

**Key findings:**
1. LOSO and K-fold implementations share infrastructure via `.dkge_build_fold_bases()` and `.dkge_fold_weight_context()`
2. Analytic LOSO has 5 documented fallback trigger conditions, only 1 currently has direct test coverage (eigengap)
3. K-fold with k=S subjects should mathematically equal LOSO - this equivalence is partially tested but not comprehensively verified for data leakage

**Primary recommendation:** Tests should verify data leakage prevention (held-out subject excluded from basis computation), fallback path coverage for analytic LOSO, and K-fold/LOSO equivalence.

## Standard Stack

The established libraries/tools for this domain:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| base R | any | eigen(), chol(), backsolve() | Core linear algebra |
| testthat | 3.x | Unit testing framework | R standard |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| future.apply | any | Parallel subject processing | Optional parallelism |

### Alternatives Considered
None - testing uses standard R/testthat patterns established in Phase 1-2.

## Architecture Patterns

### Code Organization
```
R/
├── dkge-loso.R          # dkge_loso_contrast() - single subject LOSO
├── dkge-kfold.R         # dkge_define_folds(), .dkge_contrast_kfold()
├── dkge-analytic.R      # dkge_analytic_loso(), .dkge_contrast_analytic_impl()
├── dkge-folds.R         # .dkge_build_fold_bases() - shared infrastructure
├── dkge-contrast.R      # dkge_contrast() - unified entry point
└── dkge-weights.R       # .dkge_fold_weight_context() - Chat computation
```

### Pattern 1: Iterative LOSO (Gold Standard)
**What:** Re-compute basis excluding each subject via full eigen-decomposition
**Implementation:** `dkge_loso_contrast()` in `R/dkge-loso.R`
**Key operation:**
```r
# From algo.md Section 4
Chat_minus <- Chat - w_s * S_s           # Remove subject contribution
eig_minus <- eigen(Chat_minus)           # Full O(q^3) eigen
U_minus <- Kihalf %*% eig_minus$vectors[, 1:r]  # K-orthonormal basis
```
**Used by:** `dkge_contrast(method = "loso")`

### Pattern 2: K-Fold Cross-Fitting
**What:** Split subjects into K folds, recompute basis excluding each fold
**Implementation:** `.dkge_contrast_kfold()` in `R/dkge-kfold.R`
**Key insight:** When k=S, each fold holds out exactly one subject, so K-fold should equal LOSO
**Used by:** `dkge_contrast(method = "kfold")`

### Pattern 3: Analytic LOSO Approximation
**What:** First-order eigenvalue perturbation to approximate held-out basis
**Implementation:** `dkge_analytic_loso()` in `R/dkge-analytic.R`
**Approximation formulas:**
```r
# Eigenvalue shift
delta_lambda_j <- -w_s * H[j, j]  # where H = V^T S_s V

# Eigenvector rotation (first-order perturbation)
coeffs[-j] <- -w_s * H[-j, j] / gaps[-j]  # gaps = lambda_j - lambda_k
delta_v <- V %*% coeffs
```
**Complexity:** O(q^2 * r) vs O(q^3) for full eigen
**Used by:** `dkge_contrast(method = "analytic")`

### Pattern 4: Fallback Trigger Conditions
**What:** Conditions that force analytic LOSO to fall back to iterative LOSO
**Implementation:** `dkge_analytic_loso()` lines 65-118
**Conditions:**
1. `solver_not_pooled` - JD solver used instead of pooled eigen (line 66-76)
2. `nonuniform_voxel_weights` - Voxel weights present (line 78-90)
3. `missing_full_decomposition` - Full eigen vectors/values not stored (line 98-107)
4. `dimension_mismatch` - V or lambda dimensions wrong (line 109-118)
5. `eigengap` - Eigenvalue gap too small (line 146-154)
6. `perturbation_magnitude` - Perturbation coefficients too large (line 165-173)

### Anti-Patterns to Avoid
- **Computing Chat_minus incorrectly:** Must use subject weight `w_s` and contribution `S_s`
- **Forgetting symmetrization:** Chat_minus must be symmetrized before eigen
- **Missing K-orthonormality check:** Held-out basis U_minus must satisfy U^T K U = I

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Fold definition | Custom splitting logic | `dkge_define_folds()` | Handles reproducibility, coverage checking |
| Held-out Chat | Manual subtraction | `.dkge_fold_weight_context()` | Handles weights, missingness, ridge |
| Basis building | Direct eigen calls | `.dkge_build_fold_bases()` | Handles alignment, caching, loaders |

**Key insight:** The cross-fitting infrastructure is well-factored. Tests should use existing helpers where possible but may need to verify internals.

## Common Pitfalls

### Pitfall 1: Data Leakage via Subject Contribution
**What goes wrong:** Held-out subject's data influences the basis used to evaluate them
**Why it happens:** Chat_minus computed incorrectly, or wrong subject excluded
**How to avoid:** Verify that injecting extreme values in held-out subject does not affect U_minus
**Warning signs:** U_minus equal to full U when it should differ

### Pitfall 2: K-fold vs LOSO Mismatch with k=S
**What goes wrong:** K-fold with k=S doesn't produce identical results to LOSO
**Why it happens:** Different code paths, different weighting, or fold assignment issues
**How to avoid:** Test numerical equivalence when k=S
**Warning signs:** Small numerical differences beyond floating-point tolerance

### Pitfall 3: Analytic Fallback Not Triggered
**What goes wrong:** Analytic approximation used in unstable regime
**Why it happens:** Threshold too permissive, condition checks missing
**How to avoid:** Test each fallback condition explicitly
**Warning signs:** Large discrepancy between analytic and iterative results

### Pitfall 4: Analytic vs Iterative Tolerance
**What goes wrong:** Analytic approximation marked as "matching" when it doesn't
**Why it happens:** Tolerance too loose, wrong metric used
**How to avoid:** Define explicit tolerance for "analytic matches iterative"
**Warning signs:** High correlation but large absolute errors

## Code Examples

### LOSO Contrast (from algo.md Algorithm 2)
```r
# Source: dkge-loso.R lines 13-37
dkge_loso_contrast <- function(fit, s, c, ridge = 0) {
  train_ids <- setdiff(seq_len(length(fit$Btil)), s)
  ctx <- .dkge_fold_weight_context(fit, train_ids, ridge = ridge)
  Chat_minus <- ctx$Chat

  eig_minus <- eigen(Chat_minus, symmetric = TRUE)
  r <- ncol(fit$U)
  Uminus <- fit$Kihalf %*% eig_minus$vectors[, seq_len(r), drop = FALSE]

  c_tilde <- backsolve(fit$R, c, transpose = FALSE)
  alpha <- t(Uminus) %*% fit$K %*% c_tilde

  A_s <- t(fit$Btil[[s]]) %*% fit$K %*% Uminus
  v_s <- as.numeric(A_s %*% alpha)

  list(v = v_s, alpha = alpha, basis = Uminus, evals = eig_minus$values)
}
```

### K-fold Equivalence Test Pattern
```r
# Source: test-toy-kfold.R (existing pattern)
test_that("K-fold with k = S matches LOSO", {
  # Setup data
  loso <- dkge_contrast(fit, contrast, method = "loso", align = FALSE)
  kf <- dkge_contrast(fit, contrast, method = "kfold", folds = S, align = FALSE)

  # Compare values per subject
  for (s in seq_len(S)) {
    expect_equal(loso$values[[1]][[s]], kf$values[[1]][[s]],
                 tolerance = 1e-10)
  }
})
```

### Data Leakage Test Pattern
```r
# Recommended pattern for Phase 3 tests
test_that("LOSO basis excludes held-out subject", {
  # Get baseline LOSO basis for subject s
  out1 <- dkge_loso_contrast(fit, s = s, c = contrast)

  # Inject extreme values in subject s
  fit_extreme <- fit
  fit_extreme$Btil[[s]] <- fit$Btil[[s]] * 1000
  fit_extreme$contribs[[s]] <- fit$contribs[[s]] * 1000

  # Recompute LOSO for subject s
  out2 <- dkge_loso_contrast(fit_extreme, s = s, c = contrast)

  # Basis should be IDENTICAL (subject s excluded from both)
  expect_equal(out1$basis, out2$basis, tolerance = 1e-12)
})
```

### Fallback Trigger Test Pattern
```r
# Recommended pattern for testing fallback conditions
test_that("Analytic LOSO falls back when solver is not pooled", {
  fit <- make_fit()
  fit$solver <- "jd"  # Non-pooled solver

  result <- dkge_analytic_loso(fit, s = 1, c = contrast)

  expect_equal(result$method, "fallback")
  expect_equal(result$diagnostic$reason, "solver_not_pooled")
})
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Direct Chat subtraction | `.dkge_fold_weight_context()` | Recent refactor | Centralized fold logic |
| Manual fold building | `.dkge_build_fold_bases()` | Recent refactor | Shared by LOSO and K-fold |

**Current implementation status:**
- LOSO: Fully implemented and tested for basic correctness
- K-fold: Fully implemented, partial equivalence testing with LOSO
- Analytic: Implemented with fallback, eigengap fallback tested, other fallbacks untested

## Open Questions

Things that couldn't be fully resolved:

1. **Tolerance for analytic vs iterative comparison**
   - What we know: test-toy-analytic.R uses cosine similarity > 0.98
   - What's unclear: Is this the canonical tolerance? Should there be absolute error bounds?
   - Recommendation: Use 0.98 cosine for HIGH SNR, allow degradation for LOW SNR; document in tests

2. **Weight handling during cross-fitting**
   - What we know: `.dkge_fold_weight_context()` handles complex weight scenarios
   - What's unclear: How well are edge cases tested (per-subject weights, missingness)?
   - Recommendation: Basic tests first, weight edge cases in later phase if needed

## Success Criteria (from Phase Definition)

Phase 3 tests must verify:

1. **LOSO basis U^{-s} differs from full basis U** - HIGH confidence path exists
2. **Extreme values in held-out subject do not affect U^{-s}** - No existing test, critical for data leakage
3. **K-fold with k=S equals LOSO** - Partial test exists, needs strengthening
4. **All analytic fallback paths have coverage** - Only eigengap tested, 4 untested:
   - `solver_not_pooled`
   - `nonuniform_voxel_weights`
   - `missing_full_decomposition`
   - `perturbation_magnitude`
5. **Analytic approximation matches iterative within tolerance** - Existing test uses 0.98 cosine

## Test Coverage Gaps

### Currently Tested
- LOSO: K-orthonormality, ridge behavior, column permutation invariance
- K-fold: k=S equivalence to LOSO (basic), fold definition
- Analytic: Low-leverage accuracy, eigengap fallback, basis permutation invariance

### Not Tested (Phase 3 Scope)
- **Data leakage prevention:** Held-out subject's extreme values don't affect U_minus
- **Analytic fallback paths:**
  - `solver_not_pooled` - JD solver triggers fallback
  - `nonuniform_voxel_weights` - Voxel weights trigger fallback
  - `missing_full_decomposition` - Missing eigen info triggers fallback
  - `perturbation_magnitude` - Large perturbation triggers fallback
- **LOSO basis differs from full:** Explicit verification that U_minus != U
- **K-fold with k=S numerical equivalence:** Tighter tolerance verification

## Sources

### Primary (HIGH confidence)
- `R/dkge-loso.R` - LOSO implementation
- `R/dkge-kfold.R` - K-fold implementation
- `R/dkge-analytic.R` - Analytic LOSO with fallback conditions
- `R/dkge-folds.R` - Shared fold infrastructure
- `R/dkge-weights.R` - Weight context computation
- `data-raw/algo.md` Section 4 - Algorithm 2 LOSO specification

### Secondary (MEDIUM confidence)
- `tests/testthat/test-dkge-loso.R` - Existing LOSO tests
- `tests/testthat/test-dkge-kfold.R` - Existing K-fold tests
- `tests/testthat/test_dkge_analytic.R` - Existing analytic tests
- `tests/testthat/test-toy-kfold.R` - K-fold=S equivalence test

### Tertiary (LOW confidence)
None - all findings verified against codebase.

## Metadata

**Confidence breakdown:**
- LOSO implementation: HIGH - Source code clearly matches algo.md
- K-fold implementation: HIGH - Source code well-structured
- Analytic implementation: HIGH - Source code with explicit fallback conditions
- Test coverage gaps: HIGH - Verified by grep for fallback reasons

**Research date:** 2026-01-19
**Valid until:** 2026-02-19 (30 days - stable R package code)
