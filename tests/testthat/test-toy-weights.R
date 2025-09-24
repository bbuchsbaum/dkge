library(testthat)

set.seed(11)

test_that("MFA-type weighting downweights the noisiest subject", {
  skip_on_cran()
  factors <- list(A = list(L = 2))
  sim <- dkge_sim_toy(factors,
                      active_terms = "A",
                      S = 3,
                      P = 12,
                      snr = c(8, 8, 2),
                      noise_scales = c(1, 1, 4))

  fit <- dkge_fit(sim$B_list,
                  designs = sim$X_list,
                  K = sim$K,
                  w_method = "mfa_sigma1",
                  w_tau = 0,
                  rank = 1)

  expect_equal(which.min(fit$weights), 3L)
})
