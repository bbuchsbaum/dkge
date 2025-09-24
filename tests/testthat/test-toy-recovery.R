library(testthat)

set.seed(1)

test_that("DKGE recovers a single main-effect component and contrasts line up", {
  skip_on_cran()
  factors <- list(A = list(L = 2, type = "nominal"))
  sim <- dkge_sim_toy(factors, active_terms = "A", S = 3, P = 12, snr = 10)

  fit <- dkge_fit(sim$B_list,
                  designs = sim$X_list,
                  K = sim$K,
                  w_method = "none",
                  rank = 1)

  cos1 <- dkge_cosines_K(sim$U_true[, 1, drop = FALSE],
                         fit$U[, 1, drop = FALSE],
                         sim$K)
  expect_gt(cos1[1], 0.98)

  c_true <- as.numeric(sim$U_true[, 1])
  cres <- dkge_contrast(fit, c_true, method = "loso", align = FALSE)
  cors <- vapply(seq_along(sim$B_list), function(s) {
    v_s <- cres$values[[1]][[s]]
    stats::cor(v_s, sim$M_list[[s]][, 1])
  }, numeric(1))
  expect_true(all(abs(cors) > 0.9))
})
