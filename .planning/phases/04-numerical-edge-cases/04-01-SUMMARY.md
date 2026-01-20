---
phase: 04-numerical-edge-cases
plan: 01
subsystem: testing
tags: [numerical-stability, rank-deficiency, condition-number, graceful-degradation]

# Dependency graph
requires:
  - phase: 02-fit-layer-correctness
    provides: Core dkge_fit() and fit pipeline functions
provides:
  - Internal numerical utility functions (.dkge_check_rank, .dkge_check_condition, .dkge_voxel_exclusion_mask)
  - Rank deficiency detection with warnings identifying culprit subject
  - Ill-conditioning detection with 1e8 threshold warnings
  - Minimum 2-subject requirement for group analysis
  - Sparse subject warnings for >50% missing effects
  - Effective rank tracking and rank reduction metadata in fit objects
  - 41 edge case tests covering all numerical robustness scenarios
affects: [05-inference-framework, 06-prediction-api]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Graceful degradation: warn and reduce rank rather than error"
    - "Culprit identification: warnings include subject IDs and matrix names"
    - "Metadata tracking: effective_rank and rank_reduced stored in fit objects"

key-files:
  created:
    - tests/testthat/test-numerical-robustness.R
  modified:
    - R/dkge-utils.R
    - R/dkge-data.R
    - R/dkge-fit.R
    - R/dkge-fit-core.R
    - R/dkge-align-data.R
    - tests/testthat/test-data-validation.R
    - tests/testthat/test-data.R

key-decisions:
  - "Minimum 2 subjects required for group analysis (error on single subject)"
  - "Condition number threshold 1e8 for ill-conditioning warnings"
  - "Sparse subject threshold >50% missing effects for warnings"
  - "Rank reduction stores effective_rank and rank_reduced in fit object"

patterns-established:
  - "Numerical checks via internal .dkge_check_* helpers in dkge-utils.R"
  - "Warning messages identify culprit: subject ID, matrix name, specific values"
  - "Graceful degradation: warnings allow continuation, only hard errors for fatal issues"

# Metrics
duration: 6min
completed: 2026-01-20
---

# Phase 04 Plan 01: Numerical Robustness Summary

**Internal numerical utility functions with rank deficiency detection, ill-conditioning warnings, and 41 edge case tests ensuring graceful degradation for degenerate inputs**

## Performance

- **Duration:** 6 min
- **Started:** 2026-01-20T03:51:06Z
- **Completed:** 2026-01-20T03:56:59Z
- **Tasks:** 3
- **Files modified:** 7

## Accomplishments
- Added three internal numerical utility functions for detecting rank deficiency, ill-conditioning, and non-finite values
- Integrated numerical checks into data constructors (dkge_subject, dkge_data) and fit pipeline
- Created comprehensive test suite with 41 tests covering all edge cases from CONTEXT.md
- Established minimum 2-subject requirement for group analysis
- Warnings now identify culprit subject/matrix for easier debugging

## Task Commits

Each task was committed atomically:

1. **Task 1: Add numerical robustness utility functions** - `7ea5bf2` (feat)
2. **Task 2: Integrate checks into data constructors and fit functions** - `9adefa0` (feat)
3. **Task 3: Add numerical edge case tests** - `2c497ff` (test)

## Files Created/Modified
- `R/dkge-utils.R` - Added .dkge_check_rank, .dkge_check_condition, .dkge_voxel_exclusion_mask
- `R/dkge-data.R` - Rank checks in dkge_subject(), 2-subject minimum in dkge_data()
- `R/dkge-fit.R` - Condition number check before Cholesky in shared ruler
- `R/dkge-fit-core.R` - Rank reduction warnings and effective_rank metadata
- `R/dkge-align-data.R` - Sparse subject warnings in effect alignment
- `tests/testthat/test-numerical-robustness.R` - 41 edge case tests (NEW)
- `tests/testthat/test-data-validation.R` - Updated for 2-subject minimum
- `tests/testthat/test-data.R` - Updated for 2-subject minimum

## Decisions Made
- Condition number threshold set to 1e8 (moderate, per CONTEXT.md)
- Sparse subject threshold set to >50% missing effects
- Warnings use call. = FALSE for clean messages without call stack
- effective_rank computed as sum(eigenvalues > 1e-12) consistent with existing tolerance
- rank_reduced flag tracks whether reduction occurred for metadata

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed test-data.R for 2-subject minimum**
- **Found during:** Task 3 (running full test suite)
- **Issue:** Existing test "effect alignment produces correct values after reordering" used single subject
- **Fix:** Added second subject to test data to meet minimum requirement
- **Files modified:** tests/testthat/test-data.R
- **Verification:** Test passes with 2 subjects
- **Committed in:** 2c497ff (Task 3 commit)

**2. [Rule 1 - Bug] Fixed test-data-validation.R for behavior change**
- **Found during:** Task 2 (running data validation tests)
- **Issue:** Test expected single subject to work; now errors
- **Fix:** Updated test to expect error with informative message
- **Files modified:** tests/testthat/test-data-validation.R
- **Verification:** Test validates error message content
- **Committed in:** 9adefa0 (Task 2 commit)

---

**Total deviations:** 2 auto-fixed (2 bugs in existing tests due to new 2-subject requirement)
**Impact on plan:** Necessary fixes for behavioral change. No scope creep.

## Issues Encountered
- Pre-existing deprecation warning from multivarious::prep() appears in test output but does not affect functionality
- R CMD check fails on DESCRIPTION Author/Maintainer fields (pre-existing issue, not related to this plan)

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Numerical robustness infrastructure complete
- All 1036 tests pass (0 failures, 2 skips for optional packages)
- Ready for Plan 02 (multi-seed reproducibility tests) or Phase 05 (inference framework)

---
*Phase: 04-numerical-edge-cases*
*Completed: 2026-01-20*
