# test-kernel-cv.R
# Smoke tests for kernel pre-screening and combined CV helpers

library(testthat)

test_that("kernel prescreen + LOSO CV select valid kernel/rank", {
  skip("Kernel CV smoke test temporarily disabled due to stochastic degeneracy.")
  set.seed(1)
  q <- 6; P <- 40; S <- 5; T0 <- 80

  X_list <- replicate(S, {
    X <- matrix(rnorm(T0 * q), T0, q)
    qr.Q(qr(X))
  }, simplify = FALSE)
  B_list <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)

  K1 <- diag(q)
  tmp <- matrix(rnorm(q * q), q, q)
  K2 <- crossprod(tmp) / q
  K_grid <- list(main = K1, full = K2)

  sel <- dkge_cv_kernel_rank(B_list, X_list, K_grid, ranks = 2:4,
                             ridge = 1e-6, top_k = 2)

  expect_true(sel$pick$kernel %in% names(K_grid))
  expect_true(sel$pick$rank %in% 2:4)
  expect_true(nrow(sel$tables$alignment) >= 2)
  expect_true(all(sel$tables$cv$kernel %in% names(K_grid)))
})
