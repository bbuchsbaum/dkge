---
phase: 06-integration-s3-contracts
verified: 2026-01-21T16:55:00Z
status: gaps_found
score: 3/5 must-haves verified
re_verification:
  previous_status: gaps_found
  previous_score: 3/5
  gaps_closed:
    - "All \\dontrun{} replaced with \\donttest{} (8 instances fixed)"
  gaps_remaining:
    - "R CMD check has 3 WARNINGs (vignettes + system compiler)"
    - "17 R files with @export still lack @examples"
    - "Test coverage unmeasurable (covr instrumentation broken)"
  regressions:
    - "Vignettes not building (2 new WARNINGs added)"
gaps:
  - truth: "R CMD check passes with 0 errors, 0 warnings, 0 notes"
    status: failed
    reason: "3 WARNINGs: 1 from R system header (compiler), 2 from vignettes not building"
    artifacts:
      - path: "R headers/Boolean.h"
        issue: "System-level warning: -Wfixed-enum-extension (Homebrew clang 20.x)"
      - path: "vignettes/*.Rmd"
        issue: "14 vignettes exist but inst/doc is missing (no built HTML/PDF)"
      - path: "vignettes/dkge-anchors.Rmd"
        issue: "Vignette build fails, preventing package build with vignettes"
    missing:
      - "Build vignettes and commit inst/doc/ directory"
      - "Fix dkge-anchors.Rmd build error"
      - "System warning is outside package control"
  - truth: "All exported functions have working @examples"
    status: partial
    reason: "26/43 files with @export have @examples (60%); 17 files still lack examples"
    artifacts:
      - path: "R/dkge-align-effects.R"
        issue: "Has @export but no @examples"
      - path: "R/dkge-anchor-build.R"
        issue: "Has @export but no @examples"
      - path: "R/dkge-anchor-fit.R"
        issue: "Has @export but no @examples"
      - path: "R/dkge-anchor-targets.R"
        issue: "Has @export but no @examples"
      - path: "R/dkge-components.R"
        issue: "Has @export but no @examples"
      - path: "R/dkge-contrast-validated.R"
        issue: "Has @export but no @examples"
      - path: "R/dkge-cv.R"
        issue: "Has @export but no @examples"
      - path: "R/dkge-input.R"
        issue: "Has @export but no @examples"
      - path: "R/dkge-jd.R"
        issue: "Has @export but no @examples"
      - path: "R/dkge-mapper.R"
        issue: "Has @export but no @examples"
      - path: "R/dkge-neuroim2.R"
        issue: "Has @export but no @examples"
      - path: "R/dkge-plot-suite.R"
        issue: "Has @export but no @examples"
      - path: "R/dkge-plot.R"
        issue: "Has @export but no @examples"
      - path: "R/dkge-render-core.R"
        issue: "Has @export but no @examples"
      - path: "R/dkge-voxel.R"
        issue: "Has @export but no @examples"
      - path: "R/dkge-write.R"
        issue: "Has @export but no @examples"
      - path: "R/hyperdesign-generics.R"
        issue: "Has @export but no @examples"
    missing:
      - "Add @examples to remaining 17 R files with @export"
  - truth: "Test coverage on exported functions exceeds 80%"
    status: failed
    reason: "covr package shows 0% coverage due to Rcpp instrumentation issues"
    artifacts:
      - path: "covr output"
        issue: "All R files show 0.00% coverage due to Rcpp compilation preventing instrumentation"
    missing:
      - "Fix covr instrumentation or use alternative coverage measurement"
      - "Cannot verify 80% threshold with broken coverage tooling"
human_verification:
  - test: "Build vignettes manually"
    expected: "devtools::build_vignettes() completes and inst/doc/ is populated"
    why_human: "Vignette build may require external dependencies or data"
  - test: "Run examples manually with run_examples()"
    expected: "devtools::run_examples() completes without errors"
    why_human: "Already verified passing in automated check"
---

# Phase 6: Integration + S3 Contracts Re-Verification Report

**Phase Goal:** End-to-end workflows and user-facing API behave as documented
**Verified:** 2026-01-21T16:55:00Z
**Status:** gaps_found
**Re-verification:** Yes -- after gap closure attempt (06-05-PLAN)

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | dkge_pipeline() completes successfully on valid multi-subject data | VERIFIED | 1382 tests pass, including 11 pipeline integration tests |
| 2 | All exported S3 methods (print, predict, as.data.frame) dispatch correctly | VERIFIED | 14 print method tests + 19 S3 contract tests pass |
| 3 | R CMD check passes with 0 errors, 0 warnings, 0 notes | FAILED | 0 errors, 3 WARNINGs (1 system, 2 vignettes), 0 notes |
| 4 | All exported functions have working @examples | PARTIAL | 26/43 files with @export have @examples (60%), up from 19 |
| 5 | Test coverage on exported functions exceeds 80% | FAILED | covr shows 0% due to Rcpp instrumentation issues |

**Score:** 3/5 truths verified (2 full passes, 1 partial improvement)

### Gap Closure Analysis

**Previous gaps from 06-VERIFICATION.md:**

1. **\dontrun{} replacement** - CLOSED
   - All 8 `\dontrun{}` instances replaced with `\donttest{}`
   - Verified: `grep -r "\\dontrun" R/*.R` returns empty
   
2. **Missing @examples** - PARTIALLY CLOSED
   - Files with @examples: 19 -> 26 (7 new examples added)
   - Files with @export but no @examples: 24 -> 17 (7 resolved)
   - Remaining: 17 files still lack @examples

3. **R CMD check warnings** - REGRESSED
   - Previous: 1 WARNING (system compiler)
   - Current: 3 WARNINGs (system compiler + 2 vignette issues)
   - New issue: Vignettes exist but inst/doc/ is missing

4. **Test coverage** - UNCHANGED
   - covr instrumentation still broken due to Rcpp

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `tests/testthat/test-integration-pipeline.R` | Pipeline tests | VERIFIED | 345 lines, 11 test blocks |
| `tests/testthat/test-print-methods.R` | Print contracts | VERIFIED | 234 lines, 14 test blocks |
| `tests/testthat/test-s3-contracts.R` | S3 contracts | VERIFIED | 336 lines, 19 test blocks |
| `R/*.R` with @examples | Example coverage | PARTIAL | 26/43 files (60%) |
| `inst/doc/*.html` | Built vignettes | MISSING | Directory does not exist |

### Key Link Verification

| From | To | Via | Status | Details |
|------|-----|-----|--------|---------|
| R/*.R | man/*.Rd | roxygen2 @examples | WIRED | 26 files have examples that appear in .Rd |
| vignettes/*.Rmd | inst/doc/ | knitr build | NOT_WIRED | 14 vignettes not built |
| test-*.R | R/*.R | testthat assertions | WIRED | 1382 tests pass |

### Requirements Coverage

| Requirement | Status | Blocking Issue |
|-------------|--------|----------------|
| API-01 (API contracts verified) | SATISFIED | S3 tests pass |
| CHECK-01 (R CMD check passes) | FAILED | 3 WARNINGs |
| COV-01 (Higher test coverage) | BLOCKED | covr instrumentation broken |

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| vignettes/*.Rmd | - | No built output | WARNING | R CMD check WARNING |
| 17 R files | - | Missing @examples | INFO | Incomplete documentation |

### R CMD Check Results

```
Status: 3 WARNINGs

WARNING 1: Installation warning (system-level)
  /Library/Frameworks/R.framework/Resources/include/R_ext/Boolean.h:62:36:
  warning: unknown warning group '-Wfixed-enum-extension', ignored
  [Homebrew clang 20.x issue, outside package control]

WARNING 2: Files in vignettes/ but no inst/doc
  14 .Rmd files in vignettes/ without corresponding built HTML/PDF

WARNING 3: Package vignettes without single PDF/HTML
  Directory 'inst/doc' does not exist
```

### Human Verification Required

### 1. Build Vignettes
**Test:** Run `devtools::build_vignettes()` manually
**Expected:** Vignettes build successfully, inst/doc/ populated
**Why human:** May require external dependencies, interactive debugging

### 2. Verify Examples Run
**Test:** Run `devtools::run_examples()` and observe
**Expected:** All examples complete without errors
**Why human:** Already verified passing, but good for confidence

## Gaps Summary

Three gaps remain from Phase 6 goals:

1. **R CMD check WARNINGs (3):** One system-level compiler warning (outside package control) and two vignette-related warnings (vignettes exist but not built). The vignette issue is new -- likely caused by a failing vignette (dkge-anchors.Rmd) that prevents the build.

2. **Incomplete @examples (40% missing):** 17 of 43 R files with @export lack @examples. The gap closure plan (06-05) added 7 examples, but the success criterion was "all exported functions" which requires additional work.

3. **Unverifiable test coverage:** The covr package cannot instrument this package due to Rcpp compilation issues. With 1382 tests passing, actual coverage is likely substantial, but the 80% threshold cannot be verified programmatically.

## Progress Since Initial Verification

| Metric | Before 06-05 | After 06-05 | Change |
|--------|--------------|-------------|--------|
| \dontrun{} instances | 8 | 0 | -8 (FIXED) |
| Files with @examples | 19 | 26 | +7 |
| Files missing @examples | 24 | 17 | -7 |
| R CMD check WARNINGs | 1 | 3 | +2 (REGRESSED) |
| Tests passing | 1382 | 1382 | 0 |

---

*Verified: 2026-01-21T16:55:00Z*
*Verifier: Claude (gsd-verifier)*
*Mode: Re-verification after gap closure*
