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

test_that("dkge_component_stats analyzes all medoid clusters (not the first `rank`)", {
  # rank = 2, medoid parcellation has P = 5 clusters. The index-confusion bug
  # sliced each component matrix to its first `rank` columns, yielding rank*rank
  # rows and dropping clusters 3..5.
  fixture <- make_component_fixture(S = 4, q = 3, P = 5)
  res <- dkge_component_stats(fixture$fit,
                              mapper = "ridge",
                              inference = list(type = "parametric"))
  expect_equal(sort(unique(res$summary$component)), c(1, 2))
  expect_equal(sort(unique(res$summary$cluster)), 1:5)
  expect_equal(nrow(res$summary), 2L * 5L)
  expect_equal(dim(res$transport[[1]]), c(4L, 5L))
})

test_that("dkge_component_stats selects the requested component", {
  fixture <- make_component_fixture(S = 4, q = 3, P = 5)
  res <- dkge_component_stats(fixture$fit,
                              mapper = "ridge",
                              components = 2,
                              inference = list(type = "parametric"))
  expect_equal(unique(res$summary$component), 2)
  expect_equal(sort(unique(res$summary$cluster)), 1:5)
  expect_length(res$transport, 1L)
  expect_equal(dim(res$transport[[1]]), c(4L, 5L))
})

test_that("dkge_component_stats rejects out-of-range component indices", {
  fixture <- make_component_fixture(S = 4, q = 3, P = 5)  # rank = 2
  expect_error(
    dkge_component_stats(fixture$fit, mapper = "ridge",
                         components = 5, inference = list(type = "parametric")),
    "components"
  )
})
