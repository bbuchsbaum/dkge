# test-cpca-integration.R
# Tests for CPCA integration in dkge_fit/dkge

library(testthat)
library(dkge)

make_cpca_fixture <- function(S = 6, q = 4, P = 8, T = 40, seed = 2025) {
  set.seed(seed)
  effects <- paste0("eff", seq_len(q))
  betas <- replicate(S, {
    mat <- matrix(rnorm(q * P), q, P)
    rownames(mat) <- effects
    mat
  }, simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(rnorm(T * q), T, q)
    X <- qr.Q(qr(X))
    colnames(X) <- effects
    X
  }, simplify = FALSE)
  list(betas = betas, designs = designs, effects = effects, K = diag(q))
}

fixture <- make_cpca_fixture()

kernel_identity <- fixture$K

base_fit <- dkge(fixture$betas, designs = fixture$designs, kernel = kernel_identity,
                 keep_inputs = FALSE, rank = 3)

test_that("CPCA design basis replaces Chat and basis", {
  fit <- dkge(fixture$betas, designs = fixture$designs, kernel = kernel_identity,
              cpca_blocks = 1:2, cpca_part = "design", rank = 3)
  expect_equal(fit$cpca$part, "design")
  expect_true(!is.null(fit$cpca$Chat_design))
  expect_equal(fit$Chat, fit$cpca$Chat_design, tolerance = 1e-10)
  expect_equal(fit$cpca$U_design, fit$U)
  expect_equal(fit$cpca$evals_design, fit$evals)
})

test_that("CPCA residual basis is returned when requested", {
  fit_resid <- dkge(fixture$betas, designs = fixture$designs, kernel = kernel_identity,
                    cpca_blocks = 1:3, cpca_part = "resid", rank = 3)
  expect_equal(fit_resid$cpca$part, "resid")
  expect_true(!is.null(fit_resid$cpca$Chat_resid))
  expect_equal(fit_resid$cpca$U_resid, fit_resid$U)
})

test_that("CPCA both stores design and residual bases and they are K-orthogonal", {
  fit_both <- dkge(fixture$betas, designs = fixture$designs, kernel = kernel_identity,
                   cpca_blocks = 1:2, cpca_part = "both", rank = 3)
  expect_equal(fit_both$cpca$part, "both")
  Ud <- fit_both$cpca$U_design
  Ur <- fit_both$cpca$U_resid
  expect_true(!is.null(Ur))
  cross_term <- t(Ud) %*% fit_both$K %*% Ur
  expect_equal(cross_term, matrix(0, nrow = ncol(Ud), ncol = ncol(Ur)), tolerance = 1e-6)
})

test_that("CPCA ridge shifts filtered matrix", {
  lambda <- 1e-3
  fit_ridge <- dkge(fixture$betas, designs = fixture$designs, kernel = kernel_identity,
                    cpca_blocks = 1:2, cpca_part = "design", cpca_ridge = lambda, rank = 3)
  diff_mat <- fit_ridge$cpca$Chat_design - fit_ridge$cpca$Chat_design_raw
  expect_equal(diag(diff_mat), rep(lambda, length(diag(diff_mat))), tolerance = 1e-10)
})

test_that("dkge_cpca_fit wrapper matches dkge", {
  fit1 <- dkge(fixture$betas, designs = fixture$designs, kernel = kernel_identity,
               cpca_blocks = 1:2, cpca_part = "design", rank = 3)
  fit2 <- dkge_cpca_fit(betas = fixture$betas, designs = fixture$designs,
                        kernel = kernel_identity, cpca_blocks = 1:2,
                        cpca_part = "design", rank = 3)
  expect_equal(fit1$U, fit2$U)
  expect_equal(fit1$cpca$Chat_design, fit2$cpca$Chat_design, tolerance = 1e-5)
})
