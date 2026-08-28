library(testthat)

court_path <- dkge_validation_path("functional-alignment-court", "court.R")
if (!is.na(court_path)) source(court_path, local = TRUE)

skip_without_court <- function() {
  skip_if(is.na(court_path), "validation sources are not available in this layout")
}

# How much does the tested contrast move the transport operator?
#
# The frozen court answers this where it matters, in Type-I units. These tests
# probe the mechanism directly, by flipping the sign of an effect row and
# measuring how far the fitted operators travel. Two controls make the number
# interpretable:
#
#   * geometry-only features cannot see the contrast at all, so their operators
#     must be bit-identical -- a harness check, not a claim about DKGE;
#   * flipping a *nuisance* row instead of the tested row gives a same-scale
#     comparator, since residualisation removes the tested span but not the
#     nuisance content. The tested/nuisance ratio is therefore scale-free.
#
# Measured at S = 8 over seeds 5001-5006 (deterministic given those seeds):
#
#   legacy full-fit loadings   tested/nuisance ratio  mean 0.582
#   kernel-image residual      tested/nuisance ratio  mean 0.345
#
# Two of six seeds reverse individually, so only the mean-level claim is
# asserted. Note also what is NOT claimed: a companion measurement at S = 8 vs
# S = 24 (seeds 6001-6004) gave ratios 0.273 and 0.407, i.e. no decay with
# cohort size. Any "the residual dependence is O(1/S)" argument is unsupported
# by this evidence and must not be relied on without a properly powered study.

fa_leak_sensitivity <- function(cfg, seed, arm, rows, sign_seed = 7L) {
  dat <- fa_court_generate(cfg, seed)
  set.seed(sign_seed)
  sign <- sample(c(-1, 1), cfg$S, replace = TRUE)
  flipped <- dat
  flipped$B <- fa_court_apply_raw_sign(dat$B, rows, sign)
  map_one <- function(d) {
    obs <- fa_court_observed_objects(d, cfg)
    vals <- if (identical(arm, "legacy_fullfit_functional")) {
      obs$legacy_values
    } else {
      obs$values
    }
    fa_court_map(vals, obs$features[[arm]], d, cfg, arm,
                 geometry_only = identical(arm, "geometry_only"))
  }
  fa_court_operator_distance(map_one(flipped)$operators, map_one(dat)$operators)
}

fa_leak_ratio <- function(cfg, seed, arm) {
  fa_leak_sensitivity(cfg, seed, arm, rows = 1L) /
    fa_leak_sensitivity(cfg, seed, arm, rows = 2L)
}

test_that("geometry-only alignment is exactly invariant to the sign action", {
  skip_on_cran()
  skip_without_court()
  cfg <- fa_court_base_config()
  cfg$S <- 8L
  # Harness validity. Constant features cannot see the tested contrast, so the
  # operators must be bit-identical -- not merely close. If this ever drifts,
  # the sensitivities measured below are not attributable to features and the
  # whole comparison is uninterpretable.
  expect_equal(
    fa_leak_sensitivity(cfg, seed = 4101L, arm = "geometry_only", rows = 1L),
    0,
    tolerance = 0
  )
})

test_that("residualised features are less contrast-sensitive than full-fit loadings", {
  skip_on_cran()
  skip_without_court()
  cfg <- fa_court_base_config()
  cfg$S <- 8L
  seeds <- 5001:5006
  legacy <- vapply(seeds, function(s) {
    fa_leak_ratio(cfg, s, "legacy_fullfit_functional")
  }, numeric(1))
  residual <- vapply(seeds, function(s) {
    fa_leak_ratio(cfg, s, "kernel_image_residual_prototype")
  }, numeric(1))

  # Mean-level only: individual seeds reverse (2/6 in the calibration run).
  expect_lt(mean(residual), mean(legacy))
  expect_true(
    all(is.finite(c(legacy, residual))),
    info = paste("legacy:", paste(signif(legacy, 3), collapse = ", "),
                 "| residual:", paste(signif(residual, 3), collapse = ", "))
  )
})
