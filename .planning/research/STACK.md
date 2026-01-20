# Technology Stack: R Package Testing and Quality Assurance

**Project:** dkge - Design-Kernel Group Embedding
**Researched:** 2026-01-19
**Focus:** Making a scientific/numerical R package publication-ready

## Recommended Stack

### Core Testing Framework

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| testthat | 3.3.2 | Unit testing framework | Industry standard; Edition 3 with waldo-based comparisons; excellent numerical tolerance support; published 2026-01-11 |
| hedgehog | 0.2 | Property-based testing | Integrated shrinking finds minimal failing cases; essential for numerical code validation; published 2025-11-03 |
| quickcheck | 0.1.3 | Property-based testing (testthat integration) | Simpler API wrapping hedgehog; tibble generators useful for data-heavy tests |

**Recommendation:** Use testthat 3.3.2 as the primary framework (already in place). Add hedgehog for property-based testing of mathematical invariants (commutativity, associativity, idempotency). quickcheck is optional - hedgehog alone is sufficient for a scientific package.

**Confidence:** HIGH - verified via CRAN package pages (2026-01-11 for testthat, 2025-11-03 for hedgehog)

### Code Coverage

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| covr | 3.6.5 | Test coverage measurement | Standard for R; tracks R and C/C++ code; integrates with Codecov/Coveralls; published 2025-11-09 |

**Usage pattern:**
```r
# Local coverage report (requires DT package)
covr::report()

# Identify untested code
covr::zero_coverage(covr::package_coverage())

# For CI
covr::codecov()
```

**Target:** Aim for >80% coverage on core R/ files. Focus coverage on exported functions first.

**Confidence:** HIGH - verified via CRAN and official covr documentation

### Static Analysis

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| lintr | 3.3.0-1 | Static code analysis | Identifies syntax errors, style issues, potential bugs; integrates with CI; published 2025-11-27 |
| styler | 1.11.0 | Code formatting | Automatic tidyverse style enforcement; caching for speed; published 2025-10-13 |

**Recommendation:** Use lintr for linting (catches real bugs). styler is optional for formatting - useful but not critical for CRAN submission.

**Configuration (.lintr file):**
```
linters: linters_with_defaults(
  line_length_linter(120),
  object_name_linter(styles = c("snake_case", "symbols")),
  commented_code_linter = NULL
)
```

**Confidence:** HIGH - verified via CRAN package pages

### Documentation

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| roxygen2 | 7.3.3 | In-line documentation | Standard for R; generates .Rd files and NAMESPACE; markdown support; published 2025-09-03 |
| pkgdown | (latest) | Documentation website | Already configured in package; auto-generates from roxygen2 docs |

**Confidence:** HIGH - already in use per DESCRIPTION

### Visual Regression Testing

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| vdiffr | 1.0.8 | Plot snapshot testing | Tests ggplot2 output; integrates with testthat; published 2024-10-31 |

**Recommendation:** Add vdiffr ONLY if package has significant plotting functions. From DESCRIPTION, ggplot2 is an Import, so plot testing may be relevant. However, vdiffr snapshots can be fragile across platforms - use sparingly and skip on CRAN.

**Confidence:** MEDIUM - package has ggplot2 but unclear how central plotting is

### CI/CD Infrastructure

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| r-lib/actions | v2 | GitHub Actions workflows | Standard for R community; check-standard, test-coverage workflows |
| usethis | (latest) | Workflow setup | Generates proper workflow files |
| rcmdcheck | (latest) | R CMD check from R | Better error reporting than raw R CMD check |

**Essential Workflows:**

1. **check-standard** - R CMD check on Linux, macOS, Windows + R-devel
2. **test-coverage** - Coverage reporting to Codecov

**Setup commands:**
```r
usethis::use_github_action("check-standard")
usethis::use_github_action("test-coverage")
usethis::use_github_action("lint")  # optional but recommended
```

**Confidence:** HIGH - verified via r-lib/actions repository and usethis documentation

## Numerical Testing Strategy

### Tolerance Handling

testthat 3.3.2 provides two comparison levels:

| Function | Use Case | Default Tolerance |
|----------|----------|-------------------|
| `expect_equal()` | Numerical comparisons | `sqrt(.Machine$double.eps)` ~ 1.49e-8 |
| `expect_identical()` | Exact comparisons | None (bitwise identical) |

**For matrix/numerical code:**
```r
# Default tolerance - good for most numerical comparisons
expect_equal(computed_matrix, expected_matrix)

# Explicit tolerance for numerical algorithms
expect_equal(computed_matrix, expected_matrix, tolerance = 1e-10)

# Testing relative error (important for scaled quantities)
expect_equal(result, expected, tolerance = 1e-6, scale = expected)
```

**Best practices for dkge:**
- Use `expect_equal()` with explicit tolerance for all matrix comparisons
- Test mathematical properties (symmetry, positive definiteness) separately
- Use `tolerance = 1e-10` for decomposition comparisons where numerical stability matters

### Property-Based Testing for Mathematical Code

hedgehog enables testing mathematical invariants without manually specifying test cases:

```r
library(hedgehog)
library(testthat)

test_that("kernel is positive semidefinite", {
  forall(gen.element(2:10), function(n) {
    # Generate random positive definite matrix
    A <- matrix(rnorm(n * n), n, n)
    K <- A %*% t(A)  # guaranteed PSD

    # Test that eigenvalues are non-negative
    eigs <- eigen(K, symmetric = TRUE)$values
    expect_true(all(eigs >= -1e-10))
  })
})

test_that("Cholesky decomposition reconstructs original", {
  forall(gen.element(3:8), function(n) {
    A <- matrix(rnorm(n * n), n, n)
    K <- A %*% t(A) + diag(n) * 0.1  # ensure PD
    L <- chol(K)
    expect_equal(t(L) %*% L, K, tolerance = 1e-12)
  })
})
```

**Key properties to test for dkge:**
- Kernel matrices are symmetric and PSD
- Decompositions reconstruct originals
- Projections are idempotent
- Weights sum to expected values
- Rank constraints are respected

## Installation Commands

```bash
# Core testing (already have testthat)
Rscript -e "install.packages(c('covr', 'lintr'))"

# Property-based testing
Rscript -e "install.packages('hedgehog')"

# Optional - visual testing for plots
Rscript -e "install.packages('vdiffr')"

# Optional - code formatting
Rscript -e "install.packages('styler')"
```

## Alternatives Considered

| Category | Recommended | Alternative | Why Not Alternative |
|----------|-------------|-------------|---------------------|
| Testing | testthat | tinytest | testthat is industry standard; better tooling integration; snapshot testing |
| Coverage | covr | testCoverage | covr has better CI integration; actively maintained by r-lib |
| Linting | lintr | goodpractice | lintr is more focused; goodpractice is broader but less actionable |
| Property | hedgehog | quickcheck | hedgehog is the underlying engine; quickcheck adds little for numerical code |
| CI | GitHub Actions | Travis-CI | Travis-CI deprecated for open source; GHA is standard |

## What NOT to Use

| Tool | Why Avoid |
|------|-----------|
| Travis-CI | Deprecated for open source R packages; use GitHub Actions |
| AppVeyor | Superseded by GitHub Actions multi-platform support |
| tinytest | Lacks snapshot testing and tooling integration needed for complex packages |
| RUnit | Legacy; testthat is the modern standard |
| Manual R CMD check | Use rcmdcheck for better error capture and CI integration |

## CRAN Submission Checklist Integration

For R CMD check with 0 errors, 0 warnings, 0 notes:

```r
# Run comprehensive check
devtools::check(args = c("--as-cran"))

# Or via rcmdcheck for better output
rcmdcheck::rcmdcheck(args = c("--as-cran"), error_on = "warning")

# Check on multiple platforms via R-hub
rhub::check_for_cran()
```

**rcmdcheck configuration (tools/check.env):**
```
# Suppress inconsequential notes
Config/rcmdcheck/ignore-inconsequential-notes: true
```

## Current Package Status

Based on DESCRIPTION analysis:

| Item | Status | Action Needed |
|------|--------|---------------|
| testthat | Present in Suggests | Update to edition 3 explicitly (already have `Config/testthat/edition: 3`) |
| covr | Missing from Suggests | Add for coverage tracking |
| lintr | Missing | Add for static analysis |
| hedgehog | Missing | Add for property-based testing of numerical code |
| GitHub Actions | Not present | Set up check-standard and test-coverage workflows |

## Recommended DESCRIPTION Updates

Add to Suggests:
```
Suggests:
    testthat (>= 3.0.0),
    covr,
    lintr,
    hedgehog,
    ...existing suggests...
```

Add Config entries:
```
Config/testthat/edition: 3
Config/testthat/parallel: true
```

## Sources

- [testthat CRAN](https://cran.r-project.org/package=testthat) - Version 3.3.2, published 2026-01-11
- [testthat documentation](https://testthat.r-lib.org/)
- [covr CRAN](https://cran.r-project.org/package=covr) - Version 3.6.5, published 2025-11-09
- [covr documentation](https://covr.r-lib.org/)
- [lintr CRAN](https://cran.r-project.org/package=lintr) - Version 3.3.0-1, published 2025-11-27
- [lintr documentation](https://lintr.r-lib.org/)
- [hedgehog CRAN](https://cran.r-project.org/package=hedgehog) - Version 0.2, published 2025-11-03
- [quickcheck documentation](https://armcn.github.io/quickcheck/)
- [r-lib/actions GitHub](https://github.com/r-lib/actions) - v2 recommended
- [usethis documentation](https://usethis.r-lib.org/reference/github_actions.html)
- [rcmdcheck documentation](https://rcmdcheck.r-lib.org/)
- [R Packages (2e) - Testing](https://r-pkgs.org/testing-basics.html)
- [R Packages (2e) - R CMD check](https://r-pkgs.org/R-CMD-check.html)
- [CRAN Submission Checklist](https://cran.r-project.org/web/packages/submission_checklist.html)
