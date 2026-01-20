---
phase: 04-numerical-edge-cases
verified: 2026-01-20T04:15:00Z
status: passed
score: 4/4 must-haves verified
re_verification: false
---

# Phase 4: Numerical Edge Cases Verification Report

**Phase Goal:** Package handles degenerate inputs gracefully without silent failures
**Verified:** 2026-01-20T04:15:00Z
**Status:** passed
**Re-verification:** No - initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | Rank-deficient input matrices produce informative errors or graceful degradation | VERIFIED | `.dkge_check_rank()` warns with subject ID; 5 tests pass in test-numerical-robustness.R |
| 2 | Near-singular matrices do not cause NaN propagation | VERIFIED | `.dkge_check_condition()` warns at threshold 1e8; fit completes with finite values; 3 tests pass |
| 3 | Tests pass with multiple different random seeds (not seed-dependent) | VERIFIED | 59 tests in test-multi-seed-robustness.R pass with seeds {1,2,3,5,7,42,123,999,2024} |
| 4 | Effect alignment handles partial overlap between subjects correctly | VERIFIED | sparse subject warnings; fit succeeds with partial overlap; 3 tests pass |

**Score:** 4/4 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `R/dkge-utils.R` | Numerical utility functions (min 50 lines) | VERIFIED | 171 lines; contains `.dkge_check_rank`, `.dkge_check_condition`, `.dkge_voxel_exclusion_mask` |
| `tests/testthat/test-numerical-robustness.R` | Edge case tests (min 150 lines) | VERIFIED | 416 lines; 41 tests covering rank deficiency, ill-conditioning, NaN handling, partial overlap, minimum subjects |
| `tests/testthat/test-multi-seed-robustness.R` | Multi-seed tests (min 100 lines) | VERIFIED | 345 lines; 59 tests covering fit consistency, recovery stability, edge case seed-independence, tolerance verification |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| R/dkge-data.R | R/dkge-utils.R | `.dkge_check_rank()` called in `dkge_subject.matrix()` | WIRED | Line 105: `.dkge_check_rank(design, beta, subject_id = id)` |
| R/dkge-fit.R | R/dkge-utils.R | `.dkge_check_condition()` called before Cholesky | WIRED | Line 26: `.dkge_check_condition(G_pool, threshold = 1e8, name = "pooled Gram matrix")` |
| R/dkge-data.R | - | Minimum 2 subjects check | WIRED | Lines 316-317: `stop("At least 2 subjects required for group analysis.")` |
| R/dkge-align-data.R | - | Sparse subject warning (>50% missing) | WIRED | Lines 82-92: warning for sparse effect coverage |
| R/dkge-fit-core.R | - | effective_rank/rank_reduced metadata | WIRED | Lines 212-265, 317-386, 518-519: tracks and stores rank metadata |

### Requirements Coverage

| Requirement | Status | Notes |
|-------------|--------|-------|
| EDGE-01: Edge case coverage | SATISFIED | 41 edge case tests + 59 multi-seed tests |
| FRAG-01: Fragile areas hardened | SATISFIED | Numerical checks integrated, warnings identify culprits |

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| (none) | - | - | - | No blocking anti-patterns found |

### Human Verification Required

None - all checks can be verified programmatically through the test suite.

### Test Suite Verification

**Numerical robustness tests:**
```
[ FAIL 0 | WARN 2 | SKIP 0 | PASS 41 ]
```
(Warnings are deprecation warnings from multivarious::prep(), not test failures)

**Multi-seed robustness tests:**
```
[ FAIL 0 | WARN 55 | SKIP 0 | PASS 59 ]
```
(Warnings are deprecation warnings from multivarious::prep(), not test failures)

## Summary

Phase 4 goal is fully achieved. The package now:

1. **Detects rank-deficient inputs**: `.dkge_check_rank()` identifies rank deficiency in design and beta matrices, warning with subject IDs for debugging

2. **Handles ill-conditioned matrices**: `.dkge_check_condition()` warns when condition number exceeds 1e8 threshold before Cholesky decomposition, preventing silent NaN propagation

3. **Is seed-independent**: 59 multi-seed tests verify fit determinism, recovery stability, and edge case consistency across 5+ different random seeds

4. **Handles partial effect overlap**: Effect alignment embeds missing effects as zeros, warns about sparse subjects (<50% coverage), and maintains correct pair count provenance

All numerical checks emit informative warnings identifying the culprit (subject ID, matrix name, specific values) while allowing graceful degradation. The `effective_rank` and `rank_reduced` metadata are stored in fit objects for downstream inspection.

---

*Verified: 2026-01-20T04:15:00Z*
*Verifier: Claude (gsd-verifier)*
