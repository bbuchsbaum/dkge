library(testthat)

set.seed(31)

test_that("K-fold with k = S matches LOSO", {
  skip_on_cran()
  factors <- list(A = list(L = 2))
  sim <- dkge_sim_toy(factors, active_terms = "A", S = 3, P = 10, snr = 9)

  fit <- dkge_fit(sim$B_list,
                  designs = sim$X_list,
                  K = sim$K,
                  w_method = "none",
                  rank = 1)

  c_true <- as.numeric(sim$U_true[, 1])
  loso <- dkge_contrast(fit, c_true, method = "loso", align = FALSE)
  kf <- dkge_contrast(fit, c_true, method = "kfold", folds = length(sim$B_list), align = FALSE)

  diffs <- vapply(seq_along(sim$B_list), function(s) {
    v1 <- loso$values[[1]][[s]]
    v2 <- kf$values[[1]][[s]]
    sqrt(mean((v1 - v2)^2))
  }, numeric(1))
  expect_lt(max(diffs), 1e-8)
})
