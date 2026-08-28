library(testthat)

benefit_path <- dkge_validation_path("functional-alignment-benefit", "benefit.R")
if (!is.na(benefit_path)) source(benefit_path, local = TRUE)

skip_without_benefit <- function() {
  skip_if(is.na(benefit_path), "validation sources are not available in this layout")
}

# Efficacy, not validity. The frozen Type-I court answers "does alignment
# leak?"; these tests answer "does alignment do anything?". They need a DGP the
# court deliberately does not have: one where anatomical position is genuinely
# ambiguous about functional correspondence.
#
# Gate constants were calibrated on seeds 201-215 (see fab_gates()) and are
# asserted here on disjoint seeds.

fab_test_seeds <- function() 301:304

# Each regime is exercised by several assertions; solving the transport once
# per regime keeps the file affordable without weakening any of them.
fab_cache <- new.env(parent = emptyenv())

fab_arms <- function(warp) {
  key <- paste0("warp_", warp)
  if (!is.null(fab_cache[[key]])) return(fab_cache[[key]])
  cfg <- fab_base_config()
  cfg$warp_amount <- warp
  runs <- lapply(fab_test_seeds(), function(seed) fab_compare(cfg, seed))
  fab_cache[[key]] <- runs
  runs
}

test_that("the corresponded regime really is solvable from geometry alone", {
  skip_on_cran()
  skip_without_benefit()
  # Guards the contrast the two-regime design rests on. If the unwarped DGP
  # were not near-monotone, "geometry should win here" would be untestable.
  dat <- fab_generate(
    utils::modifyList(fab_base_config(), list(warp_amount = 0)), 301
  )
  expect_equal(dat$monotone_fraction, 1)
  warped <- fab_generate(
    utils::modifyList(fab_base_config(), list(warp_amount = 0.25)), 301
  )
  # Anatomical order must actually break: a monotone warp is undone for free
  # by balanced OT on a line, which would make the study vacuous.
  expect_lt(warped$monotone_fraction, 0.85)
})

test_that("functional alignment recovers latent correspondence far better than chance", {
  skip_on_cran()
  skip_without_benefit()
  g <- fab_gates()
  runs <- fab_arms(g$warp_idiosyncratic)
  ratios <- vapply(runs, function(r) r$chance / r$functional$latent_error,
                   numeric(1))
  expect_true(
    all(ratios > g$min_chance_ratio),
    info = paste("chance/functional ratios:",
                 paste(signif(ratios, 3), collapse = ", "))
  )
})

test_that("functional alignment beats geometry when topography is idiosyncratic", {
  skip_on_cran()
  skip_without_benefit()
  g <- fab_gates()
  runs <- fab_arms(g$warp_idiosyncratic)

  error_ratio <- vapply(runs, function(r) {
    r$geometry$latent_error / r$functional$latent_error
  }, numeric(1))
  cor_gain <- vapply(runs, function(r) {
    r$functional$mean_subject_cor - r$geometry$mean_subject_cor
  }, numeric(1))

  expect_true(
    all(error_ratio > g$min_error_ratio),
    info = paste("latent-error ratios:",
                 paste(signif(error_ratio, 3), collapse = ", "))
  )
  expect_true(
    all(cor_gain > g$min_cor_gain),
    info = paste("reconstruction gains:",
                 paste(signif(cor_gain, 3), collapse = ", "))
  )
})

test_that("the advantage is correspondence, not extra smoothing", {
  skip_on_cran()
  skip_without_benefit()
  # A plan that simply averages more will reconstruct a smooth truth field
  # better without aligning anything. The functional arm must win while
  # diffusing no more than the geometry arm.
  g <- fab_gates()
  runs <- fab_arms(g$warp_idiosyncratic)
  excess <- vapply(runs, function(r) {
    r$functional$diffusion - r$geometry$diffusion
  }, numeric(1))
  expect_true(
    all(excess < g$max_diffusion_excess),
    info = paste("diffusion excess:", paste(signif(excess, 3), collapse = ", "))
  )
})

test_that("functional alignment does not manufacture an advantage when anatomy is right", {
  skip_on_cran()
  skip_without_benefit()
  # The decisive half of the design. A method that improves in both regimes is
  # measuring smoothing, not alignment.
  g <- fab_gates()
  runs <- fab_arms(g$warp_corresponded)
  cor_gain <- vapply(runs, function(r) {
    r$functional$mean_subject_cor - r$geometry$mean_subject_cor
  }, numeric(1))

  expect_true(
    all(cor_gain < g$max_cor_gain_corresponded),
    info = paste("gains where geometry suffices:",
                 paste(signif(cor_gain, 3), collapse = ", "))
  )
  # ...and the cost of aligning unnecessarily must stay bounded.
  expect_true(
    all(cor_gain > -g$max_cor_loss_corresponded),
    info = paste("losses where geometry suffices:",
                 paste(signif(cor_gain, 3), collapse = ", "))
  )
})
