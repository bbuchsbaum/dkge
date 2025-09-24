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
