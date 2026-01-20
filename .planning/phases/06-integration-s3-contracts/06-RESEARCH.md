# Phase 6: Integration + S3 Contracts - Research

**Researched:** 2026-01-20
**Domain:** R package integration testing, S3 contracts, R CMD check compliance, examples
**Confidence:** HIGH

## Summary

This phase ensures the DKGE package works as a complete product with proper R conventions. Research covers three key areas: (1) end-to-end pipeline integration testing, (2) S3 method contract verification, and (3) R CMD check compliance with working examples.

The current state shows:
- Tests pass (0 FAIL, 1177 PASS) but with 288 warnings from deprecated `prep()` calls in multivarious
- R CMD check shows 6 warnings and 3 notes requiring fixes
- No package datasets exist (data/ directory is empty)
- Examples use `\dontrun{}` which CRAN discourages
- Test coverage reporting shows 0% (likely covr configuration issue, tests clearly execute code)

**Primary recommendation:** Focus on fixing R CMD check issues first (documentation mismatches, missing imports), then add integration tests with package datasets, and finally ensure all exported functions have working examples.

## Standard Stack

The established libraries/tools for this domain:

### Core Testing Infrastructure
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| testthat | 3.x | Unit testing framework | Standard R testing, edition 3 used |
| covr | Current | Code coverage measurement | Standard for R package coverage |
| devtools | Current | Package development utilities | R CMD check integration |
| roxygen2 | 7.3.x | Documentation generation | Package already uses this |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| usethis | Current | Package setup utilities | Creating package datasets, GHA |
| pkgdown | Current | Documentation site | Already configured |
| rcmdcheck | Current | Programmatic R CMD check | CI integration |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| testthat | RUnit | testthat is modern standard, RUnit is legacy |
| covr | codecov direct | covr provides local + remote integration |

**Installation:**
All dependencies already available in the package DESCRIPTION (Suggests section).

## Architecture Patterns

### Recommended Test Organization
```
tests/
  testthat/
    helper-*.R        # Shared fixtures (already exists: helper-fit-fixture.R, helper-toy.R)
    test-integration-pipeline.R    # NEW: Full pipeline integration
    test-s3-contracts.R           # NEW: S3 method dispatch verification
    test-examples-work.R          # NEW: Verify examples run without error
```

### Pattern 1: Integration Test with Package Dataset
**What:** Test entire workflow using small bundled dataset
**When to use:** End-to-end pipeline verification
**Example:**
```r
# Source: R Packages book - Testing chapter
test_that("dkge_pipeline completes full workflow", {
  # Use package dataset (to be created)
  data(dkge_example, package = "dkge")

  result <- dkge_pipeline(
    betas = dkge_example$betas,
    designs = dkge_example$designs,
    kernel = dkge_example$kernel,
    contrasts = c(1, -1, 0)
  )

  expect_s3_class(result$fit, "dkge")
  expect_s3_class(result$contrasts, "dkge_contrasts")
  expect_true(!is.null(result$diagnostics))
})
```

### Pattern 2: S3 Method Contract Test
**What:** Verify S3 methods dispatch and return expected types
**When to use:** For each exported S3 method
**Example:**
```r
# Source: testthat S3 testing documentation
test_that("print.dkge_contrasts returns invisibly", {
  fixture <- make_small_fit()
  contrast_obj <- dkge_contrast(fixture$fit, c(1, -1, 0), method = "analytic")

  # Verify print works and returns invisibly
  output <- capture.output(result <- print(contrast_obj))
  expect_s3_class(result, "dkge_contrasts")
  expect_identical(result, contrast_obj)
  expect_true(length(output) > 0)
})

test_that("as.data.frame.dkge_contrasts returns data.frame", {
  fixture <- make_small_fit()
  contrast_obj <- dkge_contrast(fixture$fit, c(1, -1, 0), method = "analytic")

  df <- as.data.frame(contrast_obj)
  expect_s3_class(df, "data.frame")
  expect_true("contrast" %in% names(df))
  expect_true("subject" %in% names(df))
  expect_true("value" %in% names(df))
})
```

### Pattern 3: Minimal Runnable Example
**What:** Small self-contained code demonstrating function
**When to use:** @examples in roxygen documentation
**Example:**
```r
#' @examples
#' # Create minimal test data
#' betas <- replicate(3, matrix(rnorm(5 * 20), 5, 20), simplify = FALSE)
#' designs <- replicate(3, matrix(rnorm(100 * 5), 100, 5,
#'   dimnames = list(NULL, paste0("eff", 1:5))), simplify = FALSE)
#' fit <- dkge(betas, designs = designs, kernel = diag(5), rank = 2)
#' print(fit)
```

### Anti-Patterns to Avoid
- **All examples in `\dontrun{}`:** CRAN discourages, prefer `\donttest{}` for slow operations
- **Tests that depend on external state:** Use `withr::local_*` for isolation
- **Testing R's dispatch mechanism:** Trust S3 dispatch, test output correctness only
- **Ignoring test warnings:** The 288 warnings from `prep()` need upstream fix

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Test data setup | Manual matrix construction in each test | `make_small_fit()` helper | Already exists, consistent fixtures |
| Random reproducibility | Set seed manually everywhere | `withr::local_seed()` | Proper isolation |
| Temp files in tests | Manual file creation | `withr::local_tempfile()` | Auto-cleanup |
| Mock S3 methods | Manual method registration | `testthat::local_mocked_s3_method()` | New in testthat 3.3.0 |
| Coverage tracking | Manual counting | `covr::package_coverage()` | Standard tool |

**Key insight:** The package already has good helper infrastructure (`helper-fit-fixture.R`, `helper-toy.R`). Extend these patterns rather than creating new ones.

## Common Pitfalls

### Pitfall 1: R CMD check Documentation Mismatches
**What goes wrong:** Code and documentation parameters diverge
**Why it happens:** Adding arguments without updating roxygen2 comments
**How to avoid:** Run `devtools::document()` after any function signature change
**Warning signs:** "Codoc mismatches" warning in R CMD check
**Current issue:** `dkge_fit` has undocumented arguments: `solver`, `jd_control`, `jd_mask`, `jd_init`

### Pitfall 2: Missing NAMESPACE Imports
**What goes wrong:** Functions from base packages used without explicit import
**Why it happens:** Base packages auto-loaded in interactive R but not in clean check
**How to avoid:** Add `@importFrom stats p.adjust pt setNames` etc. or use `stats::p.adjust()`
**Warning signs:** "no visible global function definition" NOTE
**Current issue:** `capture.output` needs `@importFrom utils capture.output`

### Pitfall 3: Examples That Don't Run
**What goes wrong:** Examples wrapped in `\dontrun{}` everywhere
**Why it happens:** Developer convenience during initial development
**How to avoid:** Use `\donttest{}` for slow operations, working examples otherwise
**Warning signs:** CRAN policy violation, examples not exercised
**Current issue:** 10+ files have `\dontrun{}` examples

### Pitfall 4: Deprecated Upstream Functions
**What goes wrong:** Package uses deprecated functions from dependencies
**Why it happens:** Dependency updated, package not updated
**How to avoid:** Monitor dependency changelogs, update on major versions
**Warning signs:** "deprecated" warnings in test output
**Current issue:** `multivarious::prep()` deprecated in v0.3.0, should use `fit()`

### Pitfall 5: Hidden Files in Package
**What goes wrong:** Development artifacts included in package tarball
**Why it happens:** Forgot to add to .Rbuildignore
**How to avoid:** Add `.planning` and similar to .Rbuildignore
**Warning signs:** "Found hidden files" NOTE
**Current issue:** `.planning` directory flagged

### Pitfall 6: Macro Errors in Rd Files
**What goes wrong:** LaTeX math macros not recognized
**Why it happens:** Using `\times` without proper escaping in roxygen
**How to avoid:** Use `\\times` or `*` for multiplication, wrap in `\\eqn{}`
**Warning signs:** "unknown macro" warning
**Current issue:** `dkge_input_anchor.Rd` uses `\times` incorrectly

## Code Examples

Verified patterns from official sources and current codebase:

### Package Dataset Creation
```r
# Source: usethis documentation + R Packages book
# Run once in data-raw/dkge_example.R

dkge_example <- list(
  betas = replicate(3, matrix(rnorm(5 * 20), 5, 20,
    dimnames = list(paste0("eff", 1:5), paste0("c", 1:20))), simplify = FALSE),
  designs = replicate(3, {
    X <- matrix(rnorm(100 * 5), 100, 5, dimnames = list(NULL, paste0("eff", 1:5)))
    qr.Q(qr(X))
  }, simplify = FALSE),
  kernel = diag(5)
)

usethis::use_data(dkge_example, overwrite = TRUE)
```

### Fix Documentation Mismatch
```r
# Source: Current dkge_fit.R - Add missing @param entries
#' @param solver Solver for q-space problem: "pooled" or "jd"
#' @param jd_control Control parameters from [dkge_jd_control()]
#' @param jd_mask Optional mask for JD off-diagonal penalty
#' @param jd_init Optional orthogonal initialiser for JD solver
```

### Fix Deprecated Function Call
```r
# Source: multivarious changelog
# In R/dkge-fit-core.R, replace:
multivarious::prep(multivarious::pass())
# With:
multivarious::fit(multivarious::pass())
```

### Add Missing Import
```r
# Source: NAMESPACE best practices
# In R/dkge-write.R, add to roxygen2:
#' @importFrom utils capture.output
# Or use explicit namespace:
utils::capture.output(...)
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `\dontrun{}` for slow examples | `\donttest{}` or `@examplesIf` | R 4.0.0 | `\donttest{}` now runs by default |
| `multivarious::prep()` | `multivarious::fit()` | multivarious 0.3.0 | Need to update call |
| Manual test shuffling | `testthat::test_dir(shuffle = TRUE)` | testthat 3.3.0 | Better isolation testing |
| S3 mocking workarounds | `local_mocked_s3_method()` | testthat 3.3.0 | Cleaner S3 testing |

**Deprecated/outdated:**
- `multivarious::prep()`: Deprecated, use `fit()` instead (affects dkge-fit-core.R)

## Open Questions

Things that couldn't be fully resolved:

1. **Test coverage 0% issue**
   - What we know: Tests pass (1177 PASS), covr reports 0%
   - What's unclear: Why covr doesn't trace execution
   - Recommendation: May be covr configuration issue; investigate or ignore for now if tests pass

2. **Vignette build failure**
   - What we know: `dkge-anchors.Rmd` fails with "non-conformable arrays"
   - What's unclear: Whether this is example code error or data dimension mismatch
   - Recommendation: Fix vignette code before R CMD check compliance

3. **GitHub-only dependencies**
   - What we know: Package depends on 4 GitHub packages (neuroim2, fmridesign, fmrireg, fmriAR)
   - What's unclear: How this affects CRAN submission
   - Recommendation: Document in Remotes: field, may need CRAN alternatives discussion

## R CMD Check Issues Summary

Current issues from `devtools::check(vignettes = FALSE)`:

### Warnings (6)
1. **Undocumented code object:** `dkge_jd_control` - add @export and documentation
2. **Codoc mismatch:** `dkge_fit` missing docs for `solver`, `jd_control`, `jd_mask`, `jd_init`
3. **Undocumented arguments:** `as.data.frame.dkge_inference` missing `row.names`, `optional`
4. **Rd syntax:** Lost braces in `dkge_fit.Rd`, `dkge_projector_K.Rd`, `dkge_regress.Rd`
5. **Unknown macro:** `\times` in `dkge_input_anchor.Rd`
6. **Documented arg not in usage:** `stringsAsFactors` in `as.data.frame.dkge_inference.Rd`

### Notes (3)
1. **Hidden files:** `.planning` directory - add to .Rbuildignore
2. **Non-standard file:** `repomix-output.txt` - add to .Rbuildignore
3. **Undefined global:** `capture.output` - add @importFrom

## Sources

### Primary (HIGH confidence)
- R Packages (2e) book - [Testing chapter](https://r-pkgs.org/testing-basics.html), [R CMD check appendix](https://r-pkgs.org/R-CMD-check.html), [Release chapter](https://r-pkgs.org/release.html)
- testthat documentation - [testthat.r-lib.org](https://testthat.r-lib.org/)
- CRAN submission checklist - [cran.r-project.org](https://cran.r-project.org/web/packages/submission_checklist.html)
- covr package - [covr.r-lib.org](https://covr.r-lib.org/)

### Secondary (MEDIUM confidence)
- [R-hub blog on examples](https://blog.r-hub.io/2020/01/27/examples/)
- [CRAN Cookbook](http://contributor.r-project.org/cran-cookbook/general_issues.html)
- [testthat 3.3.0 announcement](https://tidyverse.org/blog/2025/11/testthat-3-3-0/)

### Tertiary (LOW confidence)
- WebSearch results on coverage targets - general consensus ~80% is reasonable

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Standard R package tooling, well-documented
- Architecture patterns: HIGH - Based on existing codebase patterns and R Packages book
- Pitfalls: HIGH - Directly observed from running R CMD check
- R CMD check fixes: HIGH - Error messages are explicit about what to fix

**Research date:** 2026-01-20
**Valid until:** 60+ days (stable R package conventions)
