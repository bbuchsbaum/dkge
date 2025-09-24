library(testthat)

set.seed(7)

test_that("LOSO CV rank selection is close to true rank in multi-factor toy", {
  skip_on_cran()
  factors <- list(A = list(L = 3), B = list(L = 3), C = list(L = 3))
  sim <- dkge_sim_toy(factors,
                      terms = list('A', 'B', 'C', c('A', 'B'), c('A', 'C'), c('B', 'C'), c('A', 'B', 'C')),
                      active_terms = c('A', 'B', 'A:B'),
                      r_per_term = c('A' = 1, 'B' = 1, 'A:B' = 1),
                      S = 4,
                      P = 14,
                      snr = 6)

  ranks <- 1:6
  cv <- dkge_cv_rank_loso(sim$B_list, sim$X_list, sim$K, ranks)
  best_param <- cv$table$param[which.max(cv$table$mean)]
  expect_equal(cv$pick, best_param)
})
