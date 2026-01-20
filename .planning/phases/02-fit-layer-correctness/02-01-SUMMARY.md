---
phase: 02-fit-layer-correctness
plan: 01
subsystem: testing
tags: [eigendecomposition, K-orthonormality, Gram-matrix, linear-algebra]

# Dependency graph
requires:
  - phase: 01-data-kernel-foundation
    provides: dkge_fit() core implementation with K-metric eigendecomposition
provides:
  - K-orthonormality property tests (t(U) %*% K %*% U = I_r)
  - Pooled Gram matrix verification tests (R^T R = sum_s X_s^T X_s)
  - Mathematical invariant test coverage for fit layer
affects: [02-02-reconstruction-variance, 03-contrast-inference, 04-pipeline-transport]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "K-orthonormality assertion: max(abs(t(U) %*% K %*% U - diag(r))) < 1e-8"
    - "Gram verification: Reduce(`+`, lapply(designs, crossprod))"
    - "withr::local_seed() for reproducible property tests"

key-files:
  created: []
  modified:
    - tests/testthat/test-fit.R

key-decisions:
  - "Use 1e-8 tolerance for K-orthonormality (consistent with Phase 1 standard)"
  - "Test 9 K-orthonormality configurations: identity, RBF, multi-factor kernels plus rank edge cases"
  - "Use 1e-10 tolerance for Gram matrix verification (tighter for exact linear algebra)"

patterns-established:
  - "K-orthonormality test pattern: UtKU <- t(fit$U) %*% fit$K %*% fit$U; expect_lt(max(abs(UtKU - diag(fit$rank))), 1e-8)"
  - "Pooled Gram test pattern: expect_equal(ruler$G_pool, Reduce(`+`, lapply(designs, crossprod)))"

# Metrics
duration: 2min
completed: 2026-01-20
---

# Phase 2 Plan 1: K-orthonormality and Pooled Design Property Tests Summary

**12 new mathematical property tests verifying U^T K U = I_r across 9 kernel/rank configurations and R^T R = sum_s X_s^T X_s for pooled design aggregation**

## Performance

- **Duration:** 2 min
- **Started:** 2026-01-20T02:43:48Z
- **Completed:** 2026-01-20T02:45:30Z
- **Tasks:** 3
- **Files modified:** 1

## Accomplishments

- Added 9 K-orthonormality property tests covering identity, RBF, and multi-factor kernels
- Added rank edge case tests: rank=1, rank=q, rank=q-1
- Added K-orthonormality tests with ridge regularization and subject weighting schemes (MFA, energy)
- Added 3 pooled Gram verification tests for heterogeneous subjects, many subjects, and row standardization
- Full test suite passes: 898 tests, 0 failures

## Task Commits

Each task was committed atomically:

1. **Task 1: Add K-orthonormality property tests** - `cd738fe` (test)
2. **Task 2: Add pooled design verification tests** - `89bfcd8` (test)
3. **Task 3: Verify full test suite** - verification only, no code changes

## Files Created/Modified

- `tests/testthat/test-fit.R` - Added 12 new test_that() blocks for mathematical property verification

## Decisions Made

- **Tolerance standards:** K-orthonormality uses 1e-8 (established in Phase 1), Gram matrix verification uses 1e-10 (tighter for exact linear algebra operations)
- **Kernel coverage:** Test identity, RBF/ordinal, and multi-factor block diagonal kernels to ensure generality
- **Subject weighting coverage:** Test both "mfa_sigma1" and "energy" methods to verify K-orthonormality holds under all weighting schemes

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None - all tests passed on first implementation.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- K-orthonormality invariant verified for dkge_fit() output
- Pooled design aggregation verified for correct Gram matrix computation
- Ready for 02-02-PLAN.md: reconstruction variance tests
- All 898 tests passing, no regressions introduced

---
*Phase: 02-fit-layer-correctness*
*Completed: 2026-01-20*
