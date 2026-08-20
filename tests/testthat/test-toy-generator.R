library(testthat)

set.seed(42)

test_that("dkge_sim_toy generates coherent shapes and PSD kernel", {
  factors <- list(A = list(L = 2, type = "nominal"))
  sim <- dkge_sim_toy(factors, active_terms = c("A"), S = 3, P = 10, snr = 5)

  expect_equal(length(sim$B_list), 3)
  q <- nrow(sim$K)
  expect_equal(nrow(sim$B_list[[1]]), q)
  expect_true(all(vapply(sim$B_list, ncol, integer(1)) == 10))

  Ksym <- (sim$K + t(sim$K)) / 2
  ev <- eigen(Ksym, symmetric = TRUE, only.values = TRUE)$values
  expect_gt(min(ev), -1e-8)

  UtKU <- t(sim$U_true) %*% sim$K %*% sim$U_true
  expect_equal(UtKU, diag(ncol(sim$U_true)), tolerance = 1e-8)
})

test_that("single-column term blocks are sampled by value, not as 1:n", {
  sim <- dkge_sim_toy(
    factors = list(condition = list(L = 2), phase = list(L = 2)),
    active_terms = c("condition", "phase"),
    S = 2, P = 3, seed = 9
  )

  expect_equal(sim$active_cols$condition, sim$info$blocks$condition)
  expect_equal(sim$active_cols$phase, sim$info$blocks$phase)
  expect_false(anyDuplicated(unlist(sim$active_cols)) > 0)
  expect_equal(crossprod(sim$U_true, sim$K %*% sim$U_true), diag(2),
               tolerance = 1e-10)
})
