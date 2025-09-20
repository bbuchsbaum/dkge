# test-components.R
# Tests for component-level helper

library(testthat)

make_component_fixture <- function(S = 4, q = 3, P = 5, T = 20, seed = 777) {
  set.seed(seed)
  betas <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(rnorm(T * q), T, q)
    qr.Q(qr(X))
  }, simplify = FALSE)
  centroids <- replicate(S, matrix(runif(P * 3), P, 3), simplify = FALSE)
  fit <- dkge_fit(betas, designs, K = diag(q), rank = 2)
  fit$centroids <- centroids
  list(fit = fit, centroids = centroids)
}

test_that("dkge_component_stats returns tidy summary", {
  fixture <- make_component_fixture()
  res <- dkge_component_stats(fixture$fit,
                              mapper = "ridge",
                              inference = list(type = "parametric"))
  expect_s3_class(res$summary, "data.frame")
  expect_true(all(c("component", "cluster", "stat", "p", "p_adj", "significant") %in% names(res$summary)))
  expect_equal(length(res$statistics), 2)
  expect_equal(length(res$transport), 2)
})

test_that("auto centroids and Sinkhorn mapper produce consensus summary", {
  fixture <- make_component_fixture(S = 5)
  res <- dkge_component_stats(fixture$fit,
                              mapper = list(strategy = "sinkhorn", epsilon = 0.05),
                              components = 1,
                              inference = list(type = "parametric"))
  expect_equal(unique(res$summary$component), 1)
  expect_equal(nrow(res$transport[[1]]), 5)
})
