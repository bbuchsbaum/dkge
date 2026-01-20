# Testing Patterns

**Analysis Date:** 2026-01-19

## Test Framework

**Runner:**
- testthat version 3 (edition 3)
- Config: `Config/testthat/edition: 3` in `DESCRIPTION`

**Assertion Library:**
- testthat (built-in expectations)

**Run Commands:**
```bash
Rscript -e "devtools::test()"           # Run all tests
Rscript -e "devtools::test(filter='fit')" # Run tests matching 'fit'
R CMD check --as-cran dkge_*.tar.gz     # Full package check with tests
```

## Test File Organization

**Location:**
- Tests located in `tests/testthat/`
- Co-located with package (not separate directory)

**Naming:**
- Test files: `test-{module}.R` (matches source file naming)
- One exception: `test_dkge_analytic.R` (legacy naming with underscore)
- Helper files: `helper-{purpose}.R`

**Structure:**
```
tests/
  testthat/
    helper-fit-fixture.R      # Shared fixture generator
    helper-toy.R              # Toy data generators for weights
    helper-transport.R        # Transport test utilities
    helper-weights-null.R     # Null weight testing helpers
    test-fit.R                # Tests for dkge_fit()
    test-contrast.R           # Tests for dkge_contrast()
    test-weights-*.R          # Weight-related test suites
    test-transport-*.R        # Transport test suites
    test-toy-*.R              # Toy scenario integration tests
    ...
```

## Test Structure

**Suite Organization:**
```r
# test-fit.R
library(testthat)

# Local fixture generator
make_fit_fixture <- function(S = 3, q = 3, P = 4, T = 20, seed = 900) {
  set.seed(seed)
  betas <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(rnorm(T * q), T, q)
    qr.Q(qr(X))
  }, simplify = FALSE)
  list(betas = betas, designs = designs, K = diag(q))
}

test_that("pooled design Cholesky matches Gram matrix", {
  fixture <- make_fit_fixture()
  ruler <- dkge:::.dkge_compute_shared_ruler(fixture$designs)
  expect_equal(ruler$G_pool, Reduce(`+`, lapply(fixture$designs, crossprod)))
  expect_equal(t(ruler$R) %*% ruler$R, ruler$G_pool, tolerance = 1e-10)
})

test_that("dkge_fit accepts raw lists and dkge_data", {
  fixture <- make_fit_fixture()
  fit1 <- dkge_fit(fixture$betas, fixture$designs, K = fixture$K, rank = 2)
  data_bundle <- dkge_data(fixture$betas, designs = fixture$designs)
  fit2 <- dkge_fit(data_bundle, K = fixture$K, rank = 2)

  expect_equal(fit1$rank, 2)
  expect_equal(fit2$rank, 2)
})
```

**Patterns:**
- Each test file loads `library(testthat)` and optionally `library(dkge)`
- Fixture functions defined at file top with descriptive names
- `set.seed()` for reproducibility in stochastic tests
- Tests grouped by related functionality

## Mocking

**Framework:** No formal mocking framework; manual approaches used

**Patterns:**
```r
# Skip tests when optional package unavailable
skip_if_no_T4transport <- function() {
  if (!requireNamespace("T4transport", quietly = TRUE)) {
    skip("T4transport package not available")
  }
}

# Direct internal function testing via :::
ruler <- dkge:::.dkge_compute_shared_ruler(fixture$designs)
C_base <- dkge:::.dkge_cost_matrix(Aemb_s, Aemb_ref, ...)
```

**What to Mock:**
- Optional external packages (T4transport, neuroim2)
- Skip entire tests when dependencies unavailable

**What NOT to Mock:**
- Core linear algebra (eigen, crossprod, etc.)
- Internal helper functions (test directly via `:::`)

## Fixtures and Factories

**Test Data:**
```r
# Helper file: helper-toy.R
toy_kernel_info <- function() {
  levels <- list(A = 2L, B = 2L, time = 3L)
  Qcell  <- prod(unlist(levels))
  map    <- diag(Qcell)
  K      <- diag(Qcell)
  list(K = K, map = map, info = list(levels = levels))
}

toy_betas <- function(nsub = 3L, Q = 12L, V = 8L, seed = 123) {
  set.seed(seed)
  lapply(seq_len(nsub), function(s) {
    M <- matrix(rnorm(Q * V, mean = 0, sd = 1), nrow = Q, ncol = V)
    if (s == 1) M[1:2, ] <- M[1:2, ] + 0.7  # Signal injection
    M
  })
}

# Minimal dkge-like fit for fold tests
toy_fold_fit <- function(nsub = 3L, Q = 12L, V = 8L, seed = 123) {
  B_list <- toy_betas(nsub = nsub, Q = Q, V = V, seed = seed)
  fit <- list(
    Btil = B_list,
    weights = rep(1, nsub),
    K = diag(Q),
    U = diag(Q)[, seq_len(Q), drop = FALSE],
    ...
  )
  class(fit) <- c("dkge", "list")
  fit
}
```

**Location:**
- Helper files in `tests/testthat/helper-*.R`
- Auto-loaded by testthat before test execution
- Shared across test files

## Coverage

**Requirements:** None formally enforced

**View Coverage:**
```bash
# No explicit coverage commands in package
# Standard R approach:
Rscript -e "covr::package_coverage()"
```

## Test Types

**Unit Tests:**
- Internal helper function validation
- Input/output dimension checking
- Mathematical property verification (symmetry, orthonormality)
```r
test_that("kernel roots reconstruct original kernel", {
  K <- matrix(c(2, 0.5, 0.5, 1.5), 2, 2)
  roots <- dkge:::.dkge_kernel_roots(K)
  expect_equal(roots$Khalf %*% roots$Khalf, K, tolerance = 1e-10)
  expect_equal(roots$Khalf %*% roots$Kihalf, diag(2), tolerance = 1e-10)
})
```

**Integration Tests:**
- Full pipeline tests with toy data
- Cross-method consistency (LOSO vs analytic)
```r
test_that("Analytic approximation matches LOSO for small perturbations", {
  data <- create_toy_data(S = 8, q = 4, P = 20)
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 2)

  result_loso <- dkge_contrast(fit, c1, method = "loso")
  result_analytic <- dkge_contrast(fit, c1, method = "analytic")

  for (s in seq_len(data$S)) {
    cor_val <- cor(result_loso$values[[1]][[s]], result_analytic$values[[1]][[s]])
    expect_gt(cor_val, 0.95)
  }
})
```

**E2E Tests:**
- Not formally structured
- Toy scenario tests serve as workflow validation

## Common Patterns

**Async Testing:**
- Not applicable (R is single-threaded by default)
- Parallel tests via `future.apply` disabled in test suite

**Error Testing:**
```r
test_that("dkge_fit detects dimension mismatches", {
  fixture <- make_fit_fixture()
  bad_designs <- fixture$designs
  bad_designs[[1]] <- bad_designs[[1]][, -1]
  expect_error(dkge_fit(fixture$betas, bad_designs, K = fixture$K, rank = 2))
})

test_that("as.matrix requires transport when cluster counts differ", {
  data <- create_mismatched_data()
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 2)
  contrast <- suppressWarnings(dkge_contrast(fit, c(1, -1, 0), method = "loso"))

  expect_error(as.matrix(contrast, contrast = 1), "Subject cluster counts differ")
})
```

**Tolerance Testing:**
```r
# Numerical tolerance for floating-point comparisons
expect_equal(ruler$G_pool, expected_gram, tolerance = 1e-10)

# Mathematical property testing
gram <- t(U_loso) %*% fit$K %*% U_loso
expect_equal(gram, diag(r), tolerance = 1e-10)
```

**Warning Suppression:**
```r
# Suppress expected warnings in tests
contrast <- suppressWarnings(dkge_contrast(fit, c(1, -1, 0), method = "loso"))
res <- suppressWarnings(dkge_infer(fit, c(1, -1, 0), transport = transport_cfg))
```

**Custom Error Classes:**
```r
test_that("as.matrix hints about transport when cluster counts differ", {
  # ...
  expect_error(as.matrix(result),
               "dkge_transport_contrasts_to_medoid",
               class = "dkge_transport_needed")
})
```

## Test File Inventory

| File | Purpose | Tests |
|------|---------|-------|
| `test-fit.R` | Core fitting | 8 |
| `test-contrast.R` | Contrast engine | 12 |
| `test-inference.R` | Inference helpers | 3 |
| `test-weights-*.R` | Weight specifications | ~15 across files |
| `test-transport*.R` | Transport utilities | ~10 across files |
| `test-design-kernel.R` | Kernel construction | 6 |
| `test-toy-*.R` | Integration scenarios | ~12 across files |
| `test-classify.R` | Classification | ~6 |
| `test-*-*.R` | Various modules | ~50+ total |

## Test Helpers

**Shared Fixtures (`helper-*.R`):**
- `helper-fit-fixture.R`: `make_small_fit()` function
- `helper-toy.R`: `toy_kernel_info()`, `toy_betas()`, `toy_fold_fit()`, `toy_real_fit()`
- `helper-transport.R`: `create_mismatched_data()`, `make_simple_transforms()`, `skip_if_no_T4transport()`
- `helper-weights-null.R`: Permutation null testing utilities

**Context Setup:**
```r
# Edition specification in helper
testthat::local_edition(3)
```

---

*Testing analysis: 2026-01-19*
