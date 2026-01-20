# test-plot-helpers.R
# Focused coverage for individual plotting utilities

library(testthat)

set.seed(314)

test_that("plot helpers return ggplot objects", {
  skip_if_not_installed("ggplot2")

  fx <- make_small_fit(S = 3, q = 4, P = 5, T = 20, rank = 3, seed = 7)
  fit <- fx$fit

  th <- theme_dkge()
  expect_s3_class(th, "theme")

  scree <- dkge_plot_scree(fit, one_se_pick = 1)
  expect_s3_class(scree, "ggplot")

  effects_plot <- dkge_plot_effect_loadings(fit, comps = 1:2)
  expect_s3_class(effects_plot, "ggplot")

  contrib <- dkge_plot_subject_contrib(fit, comps = 1:2)
  expect_true(is.list(contrib))
  expect_s3_class(contrib$weights, "ggplot")
  expect_s3_class(contrib$energy, "ggplot")

  bases <- list(
    fit$U,
    dkge_k_orthonormalize(fit$U + matrix(1e-4, nrow(fit$U), ncol(fit$U)), fit$K)
  )
  stability <- dkge_plot_subspace_stability(bases, fit$K)
  expect_s3_class(stability, "ggplot")

  info_panels <- dkge_plot_info_anchor(
    info_haufe = list(anchor = c(0.4, -0.2, 0.1)),
    info_loco = list(delta = c(0.05, 0, -0.03)),
    top = 0
  )
  expect_true(all(vapply(info_panels, inherits, logical(1), what = "ggplot")))
})
