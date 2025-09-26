# test-transport.R
# Simple diagnostic tests for transport helpers

library(testthat)

make_deterministic_sinkhorn_inputs <- function() {
  base_loadings <- matrix(c(1, 0,
                            0, 1,
                            1, 1), nrow = 3, byrow = TRUE)
  list(
    v_list = list(c(10, 20, 30), c(300, 600, 900), c(90, 30, 60)),
    A_list = list(
      base_loadings,
      base_loadings[c(3, 1, 2), , drop = FALSE],
      base_loadings
    ),
    centroids = list(
      matrix(c(0, 0, 0,
               1, 0, 0,
               0, 1, 0), nrow = 3, byrow = TRUE),
      matrix(c(0, 1, 0,
               0, 0, 0,
               1, 0, 0), nrow = 3, byrow = TRUE),
      matrix(c(0, 0, 0,
               1, 0, 0,
               0, 1, 0), nrow = 3, byrow = TRUE)
    )
  )
}

test_that("sinkhorn transport produces deterministic operators", {
  dat <- make_deterministic_sinkhorn_inputs()
  res <- dkge_transport_to_medoid_sinkhorn(
    dat$v_list,
    dat$A_list,
    dat$centroids,
    medoid = 1,
    epsilon = 1e-4,
    max_iter = 2000,
    tol = 1e-9
  )

  expect_equal(res$subj_values[1, ], dat$v_list[[1]])

  expected_plan_subject2 <- matrix(
    c(0, 0, 1 / 3,
      1 / 3, 0, 0,
      0, 1 / 3, 0),
    nrow = 3,
    byrow = TRUE
  )
  expected_plan_subject3 <- diag(1 / 3, 3)

  expect_equal(res$plans[[1]], diag(1, 3))
  expect_equal(res$plans[[2]], expected_plan_subject2, tolerance = 1e-6)
  expect_equal(res$plans[[3]], expected_plan_subject3, tolerance = 1e-6)

  expect_equal(res$subj_values[2, ], c(200, 300, 100), tolerance = 1e-6)
  expect_equal(res$subj_values[3, ], c(30, 10, 20), tolerance = 1e-6)
  expect_equal(res$value, c(30, 20, 30))
})

test_that("cpp sinkhorn wrapper matches deterministic R transport", {
  dat <- make_deterministic_sinkhorn_inputs()
  res_r <- dkge_transport_to_medoid_sinkhorn(
    dat$v_list,
    dat$A_list,
    dat$centroids,
    medoid = 2,
    epsilon = 1e-4,
    max_iter = 2000,
    tol = 1e-9
  )
  res_cpp <- dkge_transport_to_medoid_sinkhorn_cpp(
    dat$v_list,
    dat$A_list,
    dat$centroids,
    medoid = 2,
    epsilon = 1e-4,
    max_iter = 2000,
    tol = 1e-9,
    return_plans = TRUE
  )

  expect_equal(res_cpp$value, res_r$value, tolerance = 1e-6)
  expect_equal(res_cpp$subj_values, res_r$subj_values, tolerance = 1e-6)
  expect_equal(res_cpp$plans, res_r$plans, tolerance = 1e-6)
})

test_that("size penalty contributes to cost matrix for sub-unit masses", {
  Aemb_s <- diag(2)
  Aemb_ref <- diag(2)
  sizes_s <- c(0.2, 0.8)
  sizes_ref <- c(0.3, 0.6)

  C_base <- dkge:::.dkge_cost_matrix(Aemb_s, Aemb_ref,
                                     sizes_s = sizes_s,
                                     sizes_ref = sizes_ref,
                                     lambda_size = 0)
  C_pen <- dkge:::.dkge_cost_matrix(Aemb_s, Aemb_ref,
                                    sizes_s = sizes_s,
                                    sizes_ref = sizes_ref,
                                    lambda_size = 1)

  penalty <- C_pen - C_base
  expect_true(all(penalty > 0))
})
