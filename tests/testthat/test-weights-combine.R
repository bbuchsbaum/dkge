test_that("combine product/sum/overrides are consistent and normalized", {
  norm_mean <- function(x) dkge:::.dkge_norm_vec(x, "mean")
  comb      <- dkge:::.dkge_combine_weights

  w_prior <- c(1, 2, 4, 0.5, 0.25)
  w_adapt <- c(2, 1, 0.5, 4, 0.25)

  shrink <- list(alpha = 1, winsor = 0.9999, normalize = "mean", roi_smooth = FALSE)

  w_prod <- comb(norm_mean(w_prior), norm_mean(w_adapt),
                 combine = "product", mix = 0.5, shrink = shrink)
  expect_equal(mean(w_prod), 1, tolerance = 1e-12)
  expected_raw <- sqrt(norm_mean(w_prior) * norm_mean(w_adapt))
  expected <- expected_raw / mean(expected_raw)
  expect_equal(w_prod, expected, tolerance = 1e-12)

  w_sum <- comb(norm_mean(w_prior), norm_mean(w_adapt),
                combine = "sum", mix = 0.3, shrink = shrink)
  expect_equal(mean(w_sum), 1, tolerance = 1e-12)

  w_ovrA <- comb(w_prior, w_adapt,
                 combine = "override_adapt", mix = 0.5, shrink = shrink)
  expect_equal(mean(w_ovrA), 1, tolerance = 1e-12)

  w_ovrP <- comb(w_prior, w_adapt,
                 combine = "override_prior", mix = 0.5, shrink = shrink)
  expect_equal(mean(w_ovrP), 1, tolerance = 1e-12)
})

test_that("winsor and shrink behave as expected", {
  comb  <- dkge:::.dkge_combine_weights
  w_prior <- c(rep(1, 9), 1000)
  w_adapt <- rep(1, 10)
  w <- comb(w_prior, w_adapt,
            combine = "product", mix = 0.5,
            shrink = list(alpha = 0.5, winsor = 0.9, normalize = "mean", roi_smooth = FALSE))
  expect_equal(mean(w), 1, tolerance = 1e-12)
  expect_lt(max(w), 5)
})

test_that("ROI smoothing produces piecewise-constant weights within ROIs", {
  comb  <- dkge:::.dkge_combine_weights
  labs  <- factor(rep(1:2, each = 5))
  w_prior <- c(1, 1, 2, 2, 10,  0.5, 0.5, 0.5, 2, 2)
  w_adapt <- rep(1, length(w_prior))
  w <- comb(w_prior, w_adapt,
            combine = "product", mix = 0.5,
            shrink = list(alpha = 1, winsor = 0.999, normalize = "mean", roi_smooth = TRUE),
            roi_labels = labs)
  expect_true(sd(w[labs == 1]) < 1e-12)
  expect_true(sd(w[labs == 2]) < 1e-12)
  expect_equal(mean(w), 1, tolerance = 1e-12)
})
