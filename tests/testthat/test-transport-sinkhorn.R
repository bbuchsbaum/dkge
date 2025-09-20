# test-transport-sinkhorn.R
# Cross-check Sinkhorn transport plans against T4transport reference

library(testthat)

skip_if_no_T4transport()

compare_sinkhorn <- function(C, mu, nu, epsilon, max_iter = 5000, tol = 1e-12) {
  ours <- dkge:::.dkge_sinkhorn_plan(C, mu, nu, epsilon = epsilon, max_iter = max_iter, tol = tol)
  ref <- T4transport::sinkhornD(C, p = 1, wx = mu, wy = nu, lambda = epsilon,
                                maxiter = max_iter, abstol = tol)
  list(ours = ours, ref = ref)
}

test_that("Sinkhorn plan matches T4transport for simple Gaussian blobs", {
  set.seed(123)
  X <- matrix(rnorm(3 * 2), 3, 2)
  Y <- matrix(rnorm(4 * 2) + 0.1, 4, 2)
  C <- as.matrix(dist(rbind(X, Y)))[seq_len(3), 3 + seq_len(4)]
  mu <- rep(1 / 3, 3)
  nu <- rep(1 / 4, 4)
  eps <- 0.1

  cmp <- compare_sinkhorn(C, mu, nu, epsilon = eps)

  expect_equal(cmp$ref$distance, sum(cmp$ours * C), tolerance = 1e-6)
  expect_equal(unname(cmp$ours), unname(cmp$ref$plan), tolerance = 1e-8)
  expect_equal(unname(rowSums(cmp$ours)), mu, tolerance = 1e-6)
  expect_equal(unname(colSums(cmp$ours)), nu, tolerance = 1e-6)
})

test_that("Sinkhorn plan matches T4transport with non-uniform weights", {
  set.seed(456)
  X <- matrix(runif(5 * 3), 5, 3)
  Y <- matrix(runif(6 * 3) + 0.3, 6, 3)
  D <- as.matrix(dist(rbind(X, Y)))[seq_len(5), 5 + seq_len(6)]
  mu <- runif(5); mu <- mu / sum(mu)
  nu <- runif(6); nu <- nu / sum(nu)
  eps <- 0.2

  cmp <- compare_sinkhorn(D, mu, nu, epsilon = eps)

  expect_equal(cmp$ref$distance, sum(cmp$ours * D), tolerance = 1e-6)
  expect_equal(unname(cmp$ours), unname(cmp$ref$plan), tolerance = 1e-8)
  expect_equal(unname(rowSums(cmp$ours)), mu, tolerance = 1e-6)
  expect_equal(unname(colSums(cmp$ours)), nu, tolerance = 1e-6)
})

test_that("spatial penalty biases Sinkhorn plan toward nearby targets", {
  source_feat <- matrix(0, 2, 2)
  target_feat <- matrix(0, 2, 2)
  source_xyz <- matrix(c(0, 0, 1, 0), ncol = 2, byrow = TRUE)
  target_xyz <- matrix(c(0, 0, 4, 0), ncol = 2, byrow = TRUE)

  spec <- dkge_mapper_spec("sinkhorn", lambda_emb = 0, lambda_spa = 1, sigma_mm = 1,
                           epsilon = 0.1, max_iter = 500)
  mapping <- fit_mapper(spec, source_feat = source_feat, target_feat = target_feat,
                        source_xyz = source_xyz, target_xyz = target_xyz)
  plan <- mapping$operator
  expect_gt(plan[1, 1], plan[1, 2])
  expect_gt(plan[2, 2], plan[2, 1])
})

test_that("CPP Sinkhorn wrapper honours return_plans flag", {
  v_list <- replicate(2, runif(3), simplify = FALSE)
  A_list <- replicate(2, matrix(runif(3 * 2), 3, 2), simplify = FALSE)
  centroids <- replicate(2, matrix(runif(3 * 3), 3, 3), simplify = FALSE)

  with_plans <- dkge_transport_to_medoid_sinkhorn_cpp(v_list, A_list, centroids,
                                                      medoid = 1, return_plans = TRUE)
  without_plans <- dkge_transport_to_medoid_sinkhorn_cpp(v_list, A_list, centroids,
                                                         medoid = 1, return_plans = FALSE)

  expect_false(is.null(with_plans$plans))
  expect_null(without_plans$plans)
})
