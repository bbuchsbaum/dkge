# test-transport.R
# Simple diagnostic tests for transport helpers

library(testthat)

make_transport_inputs <- function(S = 3, P = 4, r = 2) {
  set.seed(32)
  v_list <- replicate(S, runif(P), simplify = FALSE)
  A_list <- replicate(S, {
    mat <- matrix(rnorm(P * r), P, r)
    mat + 0.1  # avoid zero rows
  }, simplify = FALSE)
  centroids <- replicate(S, matrix(runif(P * 3), P, 3), simplify = FALSE)
  list(v_list = v_list, A_list = A_list, centroids = centroids)
}

test_that("sinkhorn transport preserves medoid values", {
  skip("Sinkhorn tests skipped until deterministic fixture implemented.")
  dat <- make_transport_inputs()
  res <- dkge_transport_to_medoid_sinkhorn(dat$v_list, dat$A_list, dat$centroids, medoid = 1)

  expect_equal(res$subj_values[1, ], dat$v_list[[1]])
  expect_equal(length(res$value), nrow(dat$A_list[[1]]))
})

test_that("dkge_transport_to_medoid dispatches correctly", {
  skip("Sinkhorn tests skipped until deterministic fixture implemented.")
  dat <- make_transport_inputs()
  res_r <- dkge_transport_to_medoid("sinkhorn", dat$v_list, dat$A_list, dat$centroids, medoid = 2)

  expect_equal(nrow(res_r$subj_values), length(dat$v_list))
  expect_equal(ncol(res_r$subj_values), nrow(dat$A_list[[2]]))
})
