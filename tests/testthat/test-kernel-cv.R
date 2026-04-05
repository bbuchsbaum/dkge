# test-kernel-cv.R
# Smoke tests for kernel pre-screening and combined CV helpers

library(testthat)

make_kernel_cv_fixture <- function(q = 5L, P = 12L, S = 4L, seed = 42L) {
  set.seed(seed)
  design_base <- qr.Q(qr(matrix(rnorm(q * q), q, q)))
  X_list <- replicate(S, design_base, simplify = FALSE)
  B_list <- lapply(seq_len(S), function(s) {
    matrix(rnorm(q * P), q, P)
  })

  K_identity <- diag(q)
  K_ar1 <- stats::toeplitz(0.8 ^ (0:(q - 1)))
  list(B_list = B_list,
       X_list = X_list,
       K_grid = list(identity = K_identity, smooth = K_ar1))
}

test_that("kernel prescreen + LOSO CV select valid kernel/rank", {
  fixture <- make_kernel_cv_fixture()

  sel <- dkge_cv_kernel_rank(fixture$B_list, fixture$X_list, fixture$K_grid,
                             ranks = 2:4, ridge = 1e-6, top_k = 2)

  expect_true(sel$pick$kernel %in% names(fixture$K_grid))
  expect_true(sel$pick$rank %in% 2:4)
  expect_equal(nrow(sel$tables$alignment), length(fixture$K_grid))
  expect_true(all(sel$tables$cv$kernel %in% names(fixture$K_grid)))
  expect_true(all(sel$picks_per_kernel$kernel %in% names(fixture$K_grid)))
})
