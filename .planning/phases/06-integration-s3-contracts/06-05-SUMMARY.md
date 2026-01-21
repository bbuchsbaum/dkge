---
phase: 06-integration-s3-contracts
plan: 05
subsystem: documentation
tags: [roxygen, examples, CRAN, donttest]

# Dependency graph
requires:
  - phase: 06-04
    provides: Example audit identifying 8 \dontrun{} and 7 missing @examples
provides:
  - All \dontrun{} replaced with \donttest{} for slow-but-working examples
  - Key exported functions have runnable @examples
  - 26 R files with @examples (up from 19)
affects: [CRAN-submission]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Use dkge_sim_toy() factorial API for test data in examples"
    - "Wrap slow examples in \\donttest{} instead of \\dontrun{}"

key-files:
  created: []
  modified:
    - R/dkge-targets.R
    - R/dkge-latent-utils.R
    - R/dkge-info-maps.R
    - R/dkge-latent-clf.R
    - R/dkge-analytic.R
    - R/dkge-kfold.R
    - R/dkge-fit-from-kernels.R
    - R/dkge-predict.R
    - R/dkge-transport.R
    - R/dkge-loso.R
    - R/dkge-project.R
    - R/dkge-procrustes.R
    - R/dkge-bootstrap.R
    - R/dkge-cpca.R
    - man/*.Rd (15 files)

key-decisions:
  - "Use dkge_sim_toy factorial API (factors/active_terms/S/P) for examples"
  - "Wrap slow examples in \\donttest{} instead of \\dontrun{}"
  - "Add minimal but runnable @examples showing key API usage"

patterns-established:
  - "Pattern: dkge_sim_toy() for generating test data in examples"

# Metrics
duration: 5min
completed: 2026-01-21
---

# Phase 6 Plan 05: Gap Closure for @examples Summary

**All \dontrun{} replaced with \donttest{}, 7 key functions gained @examples, devtools::run_examples() passes**

## Performance

- **Duration:** 5 min
- **Started:** 2026-01-21T16:37:01Z
- **Completed:** 2026-01-21T16:41:42Z
- **Tasks:** 3
- **Files modified:** 14 R files + 15 man/*.Rd files

## Accomplishments

- Replaced 8 \dontrun{} with \donttest{} in 7 R files
- Added @examples to 7 key exported functions
- All examples verified passing via devtools::run_examples()
- R CMD check examples pass (0 errors)

## Task Commits

Each task was committed atomically:

1. **Task 1: Replace \dontrun{} with \donttest{}** - `ac51e82` (fix)
2. **Task 2: Add @examples to key functions** - `ededc1e` (docs)
3. **Task 3: Regenerate documentation** - `e5899dc` (docs)
4. **Bug fix: Correct dkge_sim_toy API** - `e93296d` (fix)

## Files Created/Modified

**R files with \dontrun -> \donttest (Task 1):**
- `R/dkge-targets.R` - dkge_targets example
- `R/dkge-latent-utils.R` - dkge_project_clusters_to_latent example
- `R/dkge-info-maps.R` - Two examples (dkge_info_map_from_classifier, dkge_info_map_haufe)
- `R/dkge-latent-clf.R` - dkge_cv_train_latent_classifier example
- `R/dkge-analytic.R` - dkge_analytic_loso example
- `R/dkge-kfold.R` - dkge_define_folds example
- `R/dkge-fit-from-kernels.R` - dkge_fit_from_kernels example

**R files with new @examples (Task 2):**
- `R/dkge-predict.R` - dkge_freeze() example
- `R/dkge-transport.R` - dkge_clear_sinkhorn_cache() example
- `R/dkge-loso.R` - dkge_loso_contrast() example (with \donttest{})
- `R/dkge-project.R` - dkge_project_btil() example
- `R/dkge-procrustes.R` - dkge_k_orthonormalize() example
- `R/dkge-bootstrap.R` - dkge_bootstrap_projected() example (with \donttest{})
- `R/dkge-cpca.R` - dkge_projector_K() example

## Decisions Made

1. **Use dkge_sim_toy factorial API** - The function uses `factors`, `active_terms`, `S`, `P` parameters rather than `n_subjects`/`q`. All examples now follow this pattern.

2. **\donttest{} for slow examples** - CRAN requires examples that actually run. \dontrun{} means "never run" while \donttest{} means "skip on CRAN due to time but runs during development."

3. **Minimal but runnable examples** - Added lightweight examples that demonstrate API usage without requiring complex setup (transport, rendering, etc.)

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Incorrect dkge_sim_toy API in examples**
- **Found during:** Task 3 (verification)
- **Issue:** Examples used non-existent parameters `n_subjects`, `q`, `toy$betas`, `toy$designs`
- **Fix:** Updated to use correct API: `factors`, `active_terms`, `S`, `P`, `toy$B_list`, `toy$X_list`
- **Files modified:** R/dkge-predict.R, R/dkge-loso.R, R/dkge-project.R, R/dkge-bootstrap.R
- **Verification:** devtools::run_examples() passes
- **Committed in:** e93296d

---

**Total deviations:** 1 auto-fixed (API mismatch)
**Impact on plan:** Bug fix necessary for correct operation. No scope creep.

## Issues Encountered

None - plan executed with one minor API correction needed.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- All @examples now functional
- R CMD check examples pass (0 errors)
- Example count increased from 19 to 26 files
- Package ready for CRAN submission workflow

---
*Phase: 06-integration-s3-contracts*
*Completed: 2026-01-21*
