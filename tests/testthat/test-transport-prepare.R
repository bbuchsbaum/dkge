# test-transport-prepare.R
# Coverage for dkge_prepare_transport helper

library(testthat)

set.seed(99)

test_that("dkge_prepare_transport builds mapper cache without prior operators", {
  fx <- make_small_fit(S = 3, q = 3, P = 5, T = 18, rank = 2, seed = 23)
  fit <- fx$fit
  betas <- fx$betas
  centroids <- replicate(length(betas), matrix(runif(ncol(betas[[1]]) * 3), ncol = 3), simplify = FALSE)
  loadings <- dkge_predict_loadings(fit, betas)

  mapper_spec <- dkge_mapper_spec("ridge", lambda = 1e-3)

  prep <- dkge_prepare_transport(fit,
                                 centroids = centroids,
                                 loadings = loadings,
                                 mapper = mapper_spec,
                                 medoid = 2L)

  expect_equal(length(prep$operators), length(loadings))
  expect_equal(prep$medoid, 2L)
  expect_true(is.matrix(prep$feature_ref))
  expect_identical(prep$feature_ref, prep$feature_list[[prep$medoid]])
  expect_true(all(vapply(prep$operators, is.matrix, logical(1))))

  ref_dim <- nrow(prep$feature_ref)
  expect_true(all(vapply(prep$operators, function(op) ncol(op) == ref_dim, logical(1))))
  expect_true(all(vapply(prep$operators, function(op) nrow(op) > 0, logical(1))))

  expect_equal(prep$mapper_spec$strategy, "ridge")
  expect_equal(prep$centroids, centroids)

  identity_ok <- prep$operators[[prep$medoid]]
  expect_true(isTRUE(all.equal(identity_ok, diag(1, ref_dim))))
})
