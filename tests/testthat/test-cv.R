# test-cv.R
# Granular checks for dkge-cv helpers using toy data

library(testthat)

make_toy_inputs <- function(S = 4, q = 3, P = 5, T = 12) {
  set.seed(11)
  betas <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(rnorm(T * q), T, q)
    qr.Q(qr(X))
  }, simplify = FALSE)
  list(betas = betas, designs = designs, K = diag(q))
}

test_that("dkge_variance_explained returns normalized proportions", {
  toy <- make_toy_inputs()
  fit <- dkge_fit(toy$betas, toy$designs, K = toy$K, rank = 2)

  ve <- dkge_variance_explained(fit)
  expect_s3_class(ve, "data.frame")
  expect_equal(nrow(ve), 2)
  expect_true(all(c("component", "sdev", "variance", "prop_var", "cum_prop_var") %in% names(ve)))
  expect_equal(sum(ve$prop_var), 1, tolerance = 1e-6)
  expect_equal(ve$cum_prop_var[nrow(ve)], 1, tolerance = 1e-6)
})

test_that("dkge_cv_rank_loso evaluates candidate ranks", {
  toy <- make_toy_inputs(S = 5, q = 4, P = 6, T = 14)
  ranks <- 1:3
  cv <- dkge_cv_rank_loso(toy$betas, toy$designs, toy$K, ranks = ranks)

  expect_true(cv$pick %in% ranks)
  expect_true(cv$best %in% ranks)
  expect_s3_class(cv$table, "data.frame")
  expect_setequal(cv$table$param, ranks)
})

test_that("dkge_cv_kernel_grid scores candidate kernels", {
  toy <- make_toy_inputs()
  K_grid <- list(identity = diag(3), scaled = 1.2 * diag(3))
  cv <- dkge_cv_kernel_grid(toy$betas, toy$designs, K_grid, rank = 2)

  expect_setequal(cv$table$kernel, names(K_grid))
  expect_true(cv$pick %in% names(K_grid))
  expect_true(cv$best %in% names(K_grid))
})

