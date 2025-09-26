# test-analytic.R
# Diagnostics for analytic LOSO approximation helpers

library(testthat)

make_analytic_toy <- function(S = 6, q = 4, P = 18, T = 60, seed = 2024) {
  set.seed(seed)
  betas <- replicate(S, matrix(rnorm(q * P, sd = 0.05), q, P), simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(rnorm(T * q), T, q)
    qr.Q(qr(X))
  }, simplify = FALSE)
  list(betas = betas, designs = designs, K = diag(q))
}

test_that("dkge_analytic_loso closely matches exact LOSO when perturbations are small", {
  toy <- make_analytic_toy()
  fit <- dkge_fit(toy$betas, toy$designs, K = toy$K,
                  rank = 2, w_method = "none")
  fit_mod <- fit
  fit_mod$weights <- fit$weights * 1e-3
  contrast <- c(1, -1, 0, 0)

  exact <- dkge_loso_contrast(fit_mod, s = 1, contrast)
  approx <- dkge_analytic_loso(fit_mod, s = 1, contrast, tol = 1e-6, fallback = TRUE)

  expect_equal(dim(approx$basis), dim(exact$basis))
  if (identical(approx$method, "analytic")) {
    expect_gt(cor(exact$v, approx$v), 0.95)
  } else {
    expect_equal(approx$v, exact$v)
  }
})

test_that("dkge_analytic_loso falls back to exact LOSO when tolerance is large", {
  toy <- make_analytic_toy()
  fit <- dkge_fit(toy$betas, toy$designs, K = toy$K, rank = 2)
  contrast <- c(1, -1, 0, 0)

  exact <- dkge_loso_contrast(fit, s = 2, contrast)
  fb <- dkge_analytic_loso(fit, s = 2, contrast, tol = 10, fallback = TRUE)

  expect_identical(fb$method, "fallback")
  expect_equal(fb$v, exact$v)
  expect_equal(fb$basis, exact$basis)
})

test_that("dkge_contrast analytic metadata records fallback usage", {
  toy <- make_analytic_toy(S = 5)
  fit <- dkge_fit(toy$betas, toy$designs, K = toy$K, rank = 2)
  contrast <- c(1, -1, 0, 0)

  res <- dkge_contrast(fit, contrast, method = "analytic", tol = 10, fallback = TRUE)

  expect_true(all(res$metadata$fallback_rates >= 0))
  expect_equal(res$metadata$fallback_rates[[1]], 1)
  expect_true(all(res$metadata$methods[[1]] == "fallback"))
})

test_that("analytic metadata includes diagnostic detail", {
  toy <- make_analytic_toy(S = 4)
  fit <- dkge_fit(toy$betas, toy$designs, K = toy$K, rank = 2)
  contrast <- c(1, -1, 0, 0)

  res <- dkge_contrast(fit, contrast, method = "analytic", tol = 1e-6, fallback = TRUE)
  detail <- res$metadata$fallback_detail
  expect_s3_class(detail, "data.frame")
  expect_true(all(c("contrast", "subject", "reason") %in% names(detail)))
})
