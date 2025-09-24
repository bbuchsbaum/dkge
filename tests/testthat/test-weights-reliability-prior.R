testthat::local_edition(3)

test_that("Sign-invariant reliability prior maintains uniform null p-values", {
  set.seed(7)
  nsub <- 12
  V <- 96
  w_prior <- dkge_test_reliability_prior(Q = 2, V = V, nsub = nsub, seed = 77)

  pvals <- dkge_test_null_uniformity(nrep = 24,
                                     nsub = nsub,
                                     V = V,
                                     adapt = "kenergy_prec",
                                     w_prior = w_prior,
                                     mix = 0.5,
                                     n_perm = 149,
                                     seed = 407)

  expect_true(abs(mean(pvals) - 0.5) < 0.12)
  expect_true(abs(stats::median(pvals) - 0.5) < 0.12)
  ks <- suppressWarnings(stats::ks.test(pvals, "punif"))
  expect_gt(ks$p.value, 0.01)
})
