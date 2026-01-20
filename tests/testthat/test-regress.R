skip_if_not_installed("testthat")

set.seed(99)

make_toy_fit <- function(Q = 4L, V = 6L, S = 3L, rank = 2L) {
  effects <- paste0("e", seq_len(Q))
  betas <- replicate(S, matrix(rnorm(Q * V), Q, V, dimnames = list(effects, NULL)), simplify = FALSE)
  designs <- replicate(S, diag(Q), simplify = FALSE)
  designs <- lapply(designs, function(X) { colnames(X) <- effects; X })
  dkge_fit(betas, designs = designs, K = diag(Q), rank = rank)
}

test_that("dkge_regress with lm engine reconstructs linear map", {
  fit <- make_toy_fit(rank = 2L)
  effects <- fit$effects
  expect_true(!is.null(fit$U))
  W <- matrix(rnorm(ncol(fit$U) * 2L), ncol = 2L)
  rownames(W) <- colnames(fit$U)
  Y <- fit$U %*% W
  rownames(Y) <- effects

  res <- dkge_regress(fit, Y, k = 2, engine = "lm", seed = 123)
  expect_s3_class(res, "dkge_regress")
  expect_equal(res$pred[effects, ], Y[effects, ], tolerance = 1e-6)
  expect_lt(res$metrics$rmse_micro, 1e-6)
  expect_gt(res$metrics$r2_micro, 0.99999)
})

test_that("dkge_regress works with glmnet engine", {
  skip_if_not_installed("glmnet")
  fit <- make_toy_fit(Q = 6L, rank = 2L)
  effects <- fit$effects
  W <- matrix(rnorm(ncol(fit$U) * 3L), ncol = 3L)
  rownames(W) <- colnames(fit$U)
  Y <- fit$U %*% W
  rownames(Y) <- effects

  res <- suppressWarnings(dkge_regress(fit, Y, k = 3, engine = "glmnet", seed = 42))
  expect_equal(dim(res$pred), dim(Y))
  expect_true(res$metrics$r2_micro > 0.9)
})
