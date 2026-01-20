---
phase: 03-cross-fitting-validation
plan: 01
subsystem: testing
tags: [loso, kfold, cross-fitting, data-leakage, testthat]

# Dependency graph
requires:
  - phase: 02-fit-layer-correctness
    provides: "Verified dkge_fit() correctness and K-orthonormality"
provides:
  - "Data leakage prevention tests for LOSO"
  - "K-fold/LOSO equivalence verification"
  - "Cross-fitting mechanics validation"
affects: [03-02-analytic-fallback-coverage]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Extreme value injection for leakage detection"
    - "Subspace comparison via singular values"

key-files:
  created:
    - tests/testthat/test-cross-fitting-leakage.R
  modified: []

key-decisions:
  - "Chat_minus tolerance 1e-7 due to floating-point accumulation"
  - "Basis comparison via cosine similarity for sign-flip invariance"
  - "Multi-column basis comparison via SVD singular values"

patterns-established:
  - "Data leakage test: inject extreme values in held-out subject, verify basis unchanged"
  - "K-fold=S LOSO equivalence: compare per-subject values with 1e-10 tolerance"

# Metrics
duration: 3min
completed: 2026-01-20
---

# Phase 3 Plan 1: Cross-Fitting Leakage Tests Summary

**LOSO and K-fold cross-fitting verified: held-out subjects properly excluded from basis computation with 6 tests covering data leakage, equivalence, and mechanics**

## Performance

- **Duration:** 3 min
- **Started:** 2026-01-20T03:05:51Z
- **Completed:** 2026-01-20T03:08:43Z
- **Tasks:** 3
- **Files modified:** 1

## Accomplishments

- Verified LOSO basis U_minus differs from full basis U for subjects with non-trivial contribution
- Verified extreme value injection in held-out subject does not affect U_minus (critical data leakage test)
- Verified K-fold with k=S produces numerically identical results to LOSO (tolerance 1e-10)
- Verified held-out bases match between K-fold and LOSO methods
- Full test suite passes (986 tests, 0 failures)

## Task Commits

Each task was committed atomically:

1. **Task 1: LOSO Data Leakage Prevention Tests** - `ea518ee` (test)
2. **Task 2: K-fold LOSO Equivalence Tests** - `4e5d7fe` (test)
3. **Task 3: Run Full Test Suite** - verification only, no commit

## Files Created/Modified

- `tests/testthat/test-cross-fitting-leakage.R` - 6 new tests (347 lines)
  - Test 1: LOSO basis differs from full basis
  - Test 2: Extreme value injection does not affect U_minus
  - Test 3: Different subjects get different held-out bases
  - Test 4: K-fold k=S produces identical values to LOSO
  - Test 5: K-fold and LOSO held-out bases match
  - Test 6: Chat_minus equals Chat minus held-out contribution

## Decisions Made

1. **Chat_minus tolerance 1e-7:** Due to floating-point accumulation when summing multiple subject contributions, the Chat_minus comparison uses 1e-7 instead of 1e-10. The critical tests (basis and eigenvalue equality) still use 1e-10.

2. **Basis comparison via cosine/SVD:** For rank-1 bases, use absolute cosine similarity to account for sign flips. For multi-column bases, use SVD singular values of the Gram matrix to verify subspace agreement.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Cross-fitting leakage prevention verified
- Ready for Plan 2: Analytic fallback path coverage testing
- Foundation established for analytic LOSO approximation validation

---
*Phase: 03-cross-fitting-validation*
*Completed: 2026-01-20*
