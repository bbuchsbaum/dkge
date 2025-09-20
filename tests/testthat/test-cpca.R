# test-cpca.R
# Minimal checks for CPCA tooling

library(testthat)

make_toy_fit <- function(S = 4, q = 4, P = 6, T = 15, rank = 3) {
  set.seed(21)
  betas <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(rnorm(T * q), T, q)
    qr.Q(qr(X))
  }, simplify = FALSE)
  dkge_fit(betas, designs, K = diag(q), rank = rank)
}

test_that("dkge_projector_K builds valid projectors", {
  fit <- make_toy_fit()
  Tmat <- diag(1, nrow(fit$U))[, 1:2, drop = FALSE]
  proj <- dkge_projector_K(Tmat, fit$K)

  expect_true(all(c("P_K", "P_hat") %in% names(proj)))
  expect_equal(dim(proj$P_K), c(nrow(fit$U), nrow(fit$U)))
  expect_equal(dim(proj$P_hat), c(nrow(fit$U), nrow(fit$U)))
})

test_that("dkge_fit_cpca returns design/residual bases", {
  fit <- make_toy_fit(rank = 2)
  res <- dkge_fit_cpca(fit, blocks = 1:2, part = "both", rank = 2)

  expect_true(all(c("U_design", "U_resid", "evals_design", "evals_resid") %in% names(res)))
  expect_equal(dim(res$U_design), dim(fit$U))
  expect_equal(dim(res$U_resid), dim(fit$U))
})

