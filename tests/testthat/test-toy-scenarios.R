library(testthat)

set.seed(101)

test_that("DKGE recovers a multi-component subspace", {
  skip_on_cran()
  factors <- list(A = list(L = 3), B = list(L = 2))
  sim <- dkge_sim_toy(factors,
                      active_terms = c("A", "B", "A:B"),
                      r_per_term = c("A" = 1, "B" = 1, "A:B" = 1),
                      S = 5,
                      P = 18,
                      snr = 12)

  fit <- dkge_fit(sim$B_list,
                  designs = sim$X_list,
                  K = sim$K,
                  w_method = "none",
                  rank = 3)

  cosines <- dkge_cosines_K(sim$U_true, fit$U[, 1:3, drop = FALSE], sim$K)
  expect_length(cosines, 3)
  expect_true(all(cosines > 0.95))
})

set.seed(202)

test_that("dkge_sim_toy allows heterogeneous cluster counts", {
  skip_on_cran()
  factors <- list(A = list(L = 2))
  P_vec <- c(8, 11, 9, 13)
  sim <- dkge_sim_toy(factors,
                      active_terms = "A",
                      S = length(P_vec),
                      P = P_vec,
                      snr = 7)

  ncols <- vapply(sim$B_list, ncol, integer(1))
  expect_equal(ncols, P_vec)
  expect_equal(length(sim$subject_ids), length(P_vec))
})

set.seed(303)

test_that("Simulated SNRs follow the requested subject profile", {
  skip_on_cran()
  factors <- list(A = list(L = 2))
  snr_target <- c(10, 4, 2)
  sim <- dkge_sim_toy(factors,
                      active_terms = "A",
                      S = length(snr_target),
                      P = 15,
                      snr = snr_target,
                      noise_scales = c(1, 1.5, 2))

  realised <- vapply(seq_along(sim$B_list), function(s) {
    signal <- sim$U_true %*% t(sim$M_list[[s]])
    resid <- sim$B_list[[s]] - signal
    fs <- sqrt(sum(signal^2))
    fn <- sqrt(sum(resid^2))
    fs / fn
  }, numeric(1))

  expect_equal(realised, snr_target, tolerance = 1.0)
})
