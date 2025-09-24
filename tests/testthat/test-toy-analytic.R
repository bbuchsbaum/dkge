library(testthat)

set.seed(21)

test_that("Analytic LOSO agrees with exact LOSO in well-separated regime", {
  skip_on_cran()
  factors <- list(A = list(L = 2))
  sim <- dkge_sim_toy(factors, active_terms = "A", S = 4, P = 10, snr = 12)

  fit <- dkge_fit(sim$B_list,
                  designs = sim$X_list,
                  K = sim$K,
                  w_method = "none",
                  rank = 1)

  c_true <- as.numeric(sim$U_true[, 1])
  loso <- dkge_contrast(fit, c_true, method = "loso", align = FALSE)
  ana <- dkge_contrast(fit, c_true, method = "analytic", align = FALSE)

  cos_subj <- vapply(seq_along(sim$B_list), function(s) {
    v1 <- loso$values[[1]][[s]]
    v2 <- ana$values[[1]][[s]]
    sum(v1 * v2) / sqrt(sum(v1^2) * sum(v2^2))
  }, numeric(1))
  expect_gt(mean(cos_subj), 0.98)
})
