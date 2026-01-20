# Phase 1: Data + Kernel Foundation - Research

**Researched:** 2026-01-19
**Domain:** R package testing, data validation, kernel mathematics
**Confidence:** HIGH

## Summary

This research examines the current state of `dkge_subject()`, `dkge_data()`, `design_kernel()`, and `kernel_roots()` to identify test coverage gaps and prioritize hardening efforts for Phase 1.

The existing test suites provide solid coverage for happy paths but lack systematic edge case testing. Key gaps include: (1) partial effect overlap edge cases in alignment logic, (2) kernel symmetry/PSD verification at construction time, (3) numerical tolerance boundaries for kernel root reconstruction, and (4) input validation error message quality.

**Primary recommendation:** Add property-based tests that verify mathematical invariants (symmetry, PSD, reconstruction identity) across randomized inputs, plus explicit edge case tests for the fragile areas identified in CONCERNS.md.

## Standard Stack

### Core Testing Framework
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| testthat | 3.x | Unit testing | Standard R testing framework, already in Suggests |
| covr | N/A | Coverage reporting | Standard for R package coverage analysis |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| withr | N/A | Local state management | For seed isolation in random tests |

**Installation:**
Already present in DESCRIPTION Suggests; no new dependencies needed.

## Existing Test Coverage Analysis

### test-data.R (126 lines, 11 tests)

**What's covered:**
1. `dkge_subject()` aligns effect names and sets defaults
2. `dkge_subject()` validates omega dimensions
3. `dkge_subject()` list method and defaults work
4. `dkge_subject()` default method errors on unsupported types
5. `dkge_data()` bundles raw matrices and normalizes IDs
6. `dkge_data()` aligns partial effect overlaps and records provenance
7. `dkge_data()` respects provided subject IDs and omega
8. `dkge_data()` errors on mismatched effects
9. `dkge()` high-level stores inputs when requested
10. `dkge()` accepts dkge_data and omega overrides

**What's missing:**
- Empty input lists (betas = list(), designs = list())
- Single-subject case edge handling
- Subjects with zero clusters (P = 0)
- Non-matrix inputs (data.frames, tibbles)
- Unicode/special characters in effect names
- Very long effect names or IDs
- Duplicate effect names across subjects
- Effect ordering invariance (same effects, different orders -> same result)
- Dimension mismatch between beta rows and design columns
- NA values in beta matrices
- Inf values in beta matrices or designs

### test-design-kernel.R (64 lines, 6 tests)

**What's covered:**
1. Cell-basis kernel is PSD and sized correctly
2. Effect-basis map matches block metadata
3. Normalization removes dependence on rho scaling
4. Contrast helpers have expected shapes (sum_contrasts, helmert_contrasts)
5. kernel_roots and alignment behave correctly
6. kernel_roots reports clamped eigenvalues

**What's missing:**
- Symmetry verification for constructed kernels
- Kernel with all zero eigenvalues (numerically degenerate)
- Very small kernels (1x1, 2x2)
- Very large kernels (100x100+) for performance regression
- Circular factor type with L=1 (edge case)
- Continuous factor with single value
- Negative rho values (should error)
- Missing required factor fields (L for discrete, values for continuous)
- Block structure ordering effects
- Jitter interaction with numerical stability
- Reconstruction tolerance: `Khalf %*% Khalf = K` for all kernel types

### test-align-effects.R (136 lines, 3 tests)

**What's covered:**
1. Nystrom completion produces union with PSD aligned kernels
2. Shrinkage and intersection modes behave as expected
3. Fold-aware alignment keeps training/test separation

**What's missing:**
- Zero overlap between subjects (disjoint effects)
- All subjects have identical effects (degenerate alignment)
- Single subject with all effects
- Empty effects list for a subject
- Effect ordering permutation invariance
- Numerical stability with nearly-singular group kernel

## Architecture Patterns

### Recommended Test Structure
```
tests/testthat/
  test-data.R                  # dkge_subject, dkge_data constructors
  test-design-kernel.R         # design_kernel, kernel_roots
  test-align-effects.R         # dkge_align_effects, dkge_align_subjects_to_union
  test-data-validation.R       # NEW: Input validation error paths
  test-kernel-invariants.R     # NEW: Mathematical property tests
  test-alignment-invariants.R  # NEW: Effect alignment invariants
```

### Pattern 1: Property-Based Testing for Mathematical Invariants
**What:** Test mathematical properties that must hold for ALL valid inputs
**When to use:** Kernel symmetry, PSD, reconstruction identity
**Example:**
```r
test_that("kernel_roots reconstruction holds for random PSD matrices", {
  set.seed(42)
  for (q in c(3, 10, 50)) {
    # Generate random PSD matrix
    L <- matrix(rnorm(q * q), q, q)
    K <- crossprod(L)
    K <- (K + t(K)) / 2

    roots <- kernel_roots(K)
    reconstructed <- roots$Khalf %*% roots$Khalf

    expect_equal(reconstructed, K, tolerance = 1e-8,
                 label = paste("q =", q))
  }
})
```

### Pattern 2: Parameterized Edge Cases
**What:** Single test function covering multiple edge case inputs
**When to use:** Input validation, boundary conditions
**Example:**
```r
test_that("dkge_subject rejects invalid inputs with clear messages", {
  design <- matrix(1:10, 5, 2, dimnames = list(NULL, c("A", "B")))

  cases <- list(
    list(beta = "not_a_matrix", msg = "matrix"),
    list(beta = matrix(1:4, 2, 2), msg = "match"),  # wrong rownames
    list(beta = matrix(1:6, 3, 2), msg = "match")   # wrong dimensions
  )

  for (case in cases) {
    expect_error(
      dkge_subject(case$beta, design = design),
      case$msg,
      info = paste("Case:", deparse(case$beta))
    )
  }
})
```

### Pattern 3: Ordering Invariance Tests
**What:** Verify results don't depend on input ordering when they shouldn't
**When to use:** Effect alignment, subject ordering
**Example:**
```r
test_that("effect alignment is invariant to effect ordering", {
  set.seed(123)
  beta1 <- matrix(rnorm(3 * 10), 3, 10, dimnames = list(c("A", "B", "C"), NULL))
  beta2 <- matrix(rnorm(3 * 10), 3, 10, dimnames = list(c("B", "A", "C"), NULL))
  design1 <- matrix(rnorm(20 * 3), 20, 3, dimnames = list(NULL, c("A", "B", "C")))
  design2 <- matrix(rnorm(20 * 3), 20, 3, dimnames = list(NULL, c("B", "A", "C")))

  data1 <- dkge_data(list(beta1), list(design1))
  data2 <- dkge_data(list(beta2[c("A", "B", "C"), ]), list(design2[, c("A", "B", "C")]))

  expect_equal(data1$betas, data2$betas)
})
```

### Anti-Patterns to Avoid
- **Testing implementation, not behavior:** Don't test internal function structure; test observable outputs
- **Fragile numeric comparisons:** Always use `tolerance` parameter for floating-point comparisons
- **Seed-dependent failures:** Use `withr::with_seed()` to isolate random state

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Matrix symmetry checks | Manual comparison | `isSymmetric(K, tol = 1e-8)` | Handles tolerance, faster |
| PSD verification | Manual eigenvalue check | `eigen(K, symmetric = TRUE)$values >= -tol` | Numerically stable |
| Effect name matching | Manual string comparison | `setequal()` + `match()` | Already used in codebase |
| Random test data | Fully custom generators | Fixtures like `make_subject_fixture()` | Consistency, reuse |

**Key insight:** The codebase already has good helper functions (`make_subject_fixture`, `make_dkge_fixture`). Extend these rather than creating parallel infrastructure.

## Common Pitfalls

### Pitfall 1: Floating-Point Comparison Without Tolerance
**What goes wrong:** Tests fail intermittently due to machine epsilon differences
**Why it happens:** Direct `expect_equal` without tolerance on matrix operations
**How to avoid:** Always specify `tolerance = 1e-8` or appropriate level for matrix equality
**Warning signs:** Tests pass locally, fail on CI; tests fail after unrelated changes

### Pitfall 2: Missing Effect Name Validation
**What goes wrong:** Misaligned effects produce silent wrong results
**Why it happens:** Effect matching uses `match()` which returns NA for missing, but downstream code may not check
**How to avoid:** Add explicit validation tests for effect name mismatches
**Warning signs:** Code calls `match()` without `anyNA()` check

### Pitfall 3: Eigenvalue Clamping Without Warning Visibility
**What goes wrong:** Small eigenvalues silently clamped, affecting reconstruction
**Why it happens:** `kernel_roots()` warns when >10% clamped, but threshold may miss edge cases
**How to avoid:** Test reconstruction identity explicitly; verify warning threshold is appropriate
**Warning signs:** Reconstruction tests pass with loose tolerance, fail with tight tolerance

### Pitfall 4: Effect Ordering Assumptions
**What goes wrong:** Code assumes effects arrive in same order across subjects
**Why it happens:** `.align_effects()` reorders, but not all paths go through it
**How to avoid:** Add ordering invariance tests; verify all entry points normalize ordering
**Warning signs:** Results differ based on order effects are listed in design matrix columns

## Code Examples

### Kernel Reconstruction Test
```r
# Source: Based on test-design-kernel.R lines 46-55
test_that("kernel_roots reconstruction satisfies Khalf %*% Khalf = K", {
  factors <- list(A = list(L = 3, type = "nominal"))
  K <- design_kernel(factors, basis = "cell", normalize = "none")$K

  roots <- kernel_roots(K)
  reconstructed <- roots$Khalf %*% roots$Khalf

  expect_equal(reconstructed, K, tolerance = 1e-8)
  # Also verify symmetry
  expect_true(isSymmetric(roots$Khalf, tol = 1e-10))
  expect_true(isSymmetric(roots$Kihalf, tol = 1e-10))
})
```

### Effect Alignment Invariance Test
```r
# Source: Based on R/dkge-data.R lines 313-350
test_that("dkge_data produces identical results regardless of effect order", {
  set.seed(999)
  # Create two subjects with shuffled effect ordering
  effects_a <- c("eff1", "eff2", "eff3")
  effects_b <- c("eff3", "eff1", "eff2")  # permuted

  beta_a <- matrix(rnorm(3 * 10), 3, 10, dimnames = list(effects_a, NULL))
  beta_b <- matrix(rnorm(3 * 10), 3, 10, dimnames = list(effects_b, NULL))

  design_a <- matrix(rnorm(20 * 3), 20, 3, dimnames = list(NULL, effects_a))
  design_b <- matrix(rnorm(20 * 3), 20, 3, dimnames = list(NULL, effects_b))

  data1 <- dkge_data(list(beta_a, beta_b), list(design_a, design_b))

  # All betas should be aligned to same effect ordering
  expect_identical(rownames(data1$betas[[1]]), rownames(data1$betas[[2]]))
  expect_identical(colnames(data1$designs[[1]]), colnames(data1$designs[[2]]))
  expect_identical(data1$effects, rownames(data1$betas[[1]]))
})
```

### Input Validation Error Quality Test
```r
# Source: Best practice for API-01 requirement
test_that("dkge_subject provides clear error messages for invalid inputs", {
  design <- matrix(1:10, 5, 2, dimnames = list(NULL, c("eff1", "eff2")))

  # Wrong type
  expect_error(
    dkge_subject(data.frame(a = 1:2, b = 3:4), design = design),
    class = "error"
  )

  # Dimension mismatch
  beta_wrong_rows <- matrix(1:6, 3, 2, dimnames = list(c("a", "b", "c"), NULL))
  expect_error(
    dkge_subject(beta_wrong_rows, design = design),
    "match"
  )

  # Mismatched effect names
  beta_wrong_names <- matrix(1:4, 2, 2, dimnames = list(c("x", "y"), NULL))
  expect_error(
    dkge_subject(beta_wrong_names, design = design),
    "match"
  )
})
```

### PSD Verification Test
```r
# Source: Mathematical requirement from data-raw/algo.md
test_that("design_kernel produces PSD matrices for all factor types", {
  factor_specs <- list(
    list(type = "nominal", L = 4),
    list(type = "ordinal", L = 4, l = 1.0),
    list(type = "circular", L = 4, l = 1.0),
    list(type = "continuous", values = c(1, 2, 3, 4))
  )

  for (spec in factor_specs) {
    factors <- list(A = spec)
    K <- design_kernel(factors, basis = "cell", normalize = "none")$K

    # Verify symmetry
    expect_true(isSymmetric(K, tol = 1e-10),
                info = paste("Type:", spec$type))

    # Verify PSD (all eigenvalues >= 0)
    eig <- eigen(K, symmetric = TRUE)$values
    expect_true(all(eig >= -1e-10),
                info = paste("Type:", spec$type, "min eig:", min(eig)))
  }
})
```

## Key Invariants to Test

Based on data-raw/algo.md mathematical specification:

| Invariant | Location | Priority | Test Strategy |
|-----------|----------|----------|---------------|
| K = K^T (symmetry) | design_kernel | HIGH | Property test across factor types |
| K is PSD | design_kernel | HIGH | Eigenvalue check >= -epsilon |
| Khalf %*% Khalf = K | kernel_roots | HIGH | Reconstruction within 1e-8 |
| Kihalf %*% K %*% Kihalf = I | kernel_roots | MEDIUM | Identity reconstruction |
| Effects aligned across subjects | dkge_data | HIGH | Ordering invariance test |
| Provenance tracks coverage | dkge_data | MEDIUM | obs_mask, pair_counts accuracy |

## Edge Cases to Cover

### dkge_subject / dkge_data
| Edge Case | Current Coverage | Priority |
|-----------|------------------|----------|
| Empty betas list | None | HIGH |
| Single subject | Partial | MEDIUM |
| Zero clusters (P = 0) | None | HIGH |
| NA/Inf in beta | None | HIGH |
| Duplicate effect names | None | MEDIUM |
| Very long effect names | None | LOW |
| Unicode in effect names | None | LOW |
| Non-matrix input types | Partial | MEDIUM |

### design_kernel
| Edge Case | Current Coverage | Priority |
|-----------|------------------|----------|
| 1x1 kernel | None | MEDIUM |
| 2x2 kernel | None | MEDIUM |
| All-zero rho | None | MEDIUM |
| Negative rho | None | HIGH |
| Missing L field | None | HIGH |
| Missing values for continuous | None | HIGH |
| Single-level factor (L=1) | None | MEDIUM |

### kernel_roots
| Edge Case | Current Coverage | Priority |
|-----------|------------------|----------|
| Already diagonal K | None | MEDIUM |
| Near-singular K | Partial | HIGH |
| Exact zeros on diagonal | None | HIGH |
| Large negative eigenvalue (invalid input) | None | HIGH |

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| testthat 2.x `context()` | testthat 3.x (no context) | testthat 3.0 | Package uses edition 3 |
| Manual tolerance checks | `expect_equal(tolerance=)` | Long ago | Already followed |
| Global seed setting | `withr::local_seed()` | Best practice | Should adopt |

**Deprecated/outdated:**
- `context()` calls: Still present in test-design-kernel.R line 3; should be removed for edition 3

## Recommended Plan Structure

### Task Group 1: Input Validation Hardening
1. Add validation tests for dkge_subject edge cases
2. Add validation tests for dkge_data edge cases
3. Verify error message quality meets API-01

### Task Group 2: Mathematical Invariant Tests
1. Add symmetry tests for design_kernel across all factor types
2. Add PSD tests for design_kernel across all factor types
3. Add reconstruction tests for kernel_roots
4. Add inverse reconstruction tests for kernel_roots

### Task Group 3: Effect Alignment Hardening
1. Add ordering invariance tests
2. Add partial overlap edge case tests
3. Add provenance accuracy tests

### Task Group 4: Coverage Cleanup
1. Remove deprecated context() calls
2. Add missing edge cases to existing tests
3. Run covr and document coverage improvement

## Open Questions

1. **Numerical tolerance thresholds**
   - What we know: kernel_roots uses 1e-10 default jitter; tests use various tolerances
   - What's unclear: Is 1e-8 reconstruction tolerance tight enough? Too tight?
   - Recommendation: Test with range of tolerances, document expected precision

2. **Large kernel performance**
   - What we know: Typical q = 10-100 per algo.md
   - What's unclear: At what q does kernel_roots become slow?
   - Recommendation: Add timing-based regression test for q = 100

## Sources

### Primary (HIGH confidence)
- `/Users/bbuchsbaum/code/dkge/R/dkge-data.R` - Source code analysis
- `/Users/bbuchsbaum/code/dkge/R/design-kernel.R` - Source code analysis
- `/Users/bbuchsbaum/code/dkge/tests/testthat/test-data.R` - Existing coverage
- `/Users/bbuchsbaum/code/dkge/tests/testthat/test-design-kernel.R` - Existing coverage
- `/Users/bbuchsbaum/code/dkge/data-raw/algo.md` - Mathematical specification
- `/Users/bbuchsbaum/code/dkge/.planning/codebase/CONCERNS.md` - Known fragile areas

### Secondary (MEDIUM confidence)
- testthat edition 3 documentation - R package testing standards

## Metadata

**Confidence breakdown:**
- Existing coverage analysis: HIGH - Direct source inspection
- Missing test gaps: HIGH - Systematic comparison to spec
- Edge case identification: MEDIUM - Based on common patterns
- Performance thresholds: LOW - No benchmarks run

**Research date:** 2026-01-19
**Valid until:** 30 days (stable domain, internal codebase)
