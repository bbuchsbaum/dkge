library(testthat)

benefit_path <- dkge_validation_path("functional-alignment-benefit", "benefit.R")
if (!is.na(benefit_path)) source(benefit_path, local = TRUE)

skip_without_benefit <- function() {
  skip_if(is.na(benefit_path), "validation sources are not available in this layout")
}

fab_symmetry_fixture <- function(seed = 401L) {
  cfg <- fab_base_config()
  cfg$warp_amount <- 0.25
  dat <- fab_generate(cfg, seed)
  values <- lapply(seq_along(dat$ids), function(s) {
    as.numeric(dat$B[[s]][2L, ])
  })
  names(values) <- dat$ids
  transport <- fab_transport(dat, cfg, values)
  list(
    cfg = cfg,
    dat = dat,
    transport = transport,
    symmetry = fab_operator_symmetry(transport, dat, cfg)
  )
}

test_that("reference and non-reference subjects use one mapper policy", {
  skip_on_cran()
  skip_without_benefit()
  fixture <- fab_symmetry_fixture()
  tr <- fixture$transport
  ref <- fixture$cfg$reference_subject
  Q <- nrow(tr$feature_ref)

  # The reference is a genuine entropic self-map, not a hard-coded identity.
  expect_false(isTRUE(all.equal(tr$operators[[ref]], diag(Q), tolerance = 0)))
  expect_true(tr$diagnostics[[ref]]$self_map)
  expect_false(any(vapply(tr$diagnostics[-ref], `[[`, logical(1), "self_map")))

  epsilon <- vapply(tr$diagnostics, `[[`, numeric(1), "epsilon")
  expect_equal(epsilon, rep(fixture$cfg$epsilon, length(epsilon)))
  expect_true(all(vapply(tr$diagnostics, function(x) {
    is.list(x$point_spread) && is.list(x$amplitude_preservation)
  }, logical(1))))
})

test_that("effective point spread and amplitude diagnostics are auditable", {
  skip_on_cran()
  skip_without_benefit()
  fixture <- fab_symmetry_fixture()
  tr <- fixture$transport
  ref <- fixture$cfg$reference_subject
  Q <- nrow(tr$feature_ref)

  for (s in seq_along(tr$diagnostics)) {
    diagnostic <- tr$diagnostics[[s]]
    expect_length(diagnostic$point_spread$target_effective_points, Q)
    expect_true(all(diagnostic$point_spread$target_effective_points >= 1))
    expect_true(is.finite(
      diagnostic$amplitude_preservation$mapped_to_target_rms_ratio
    ))
    expect_true(isTRUE(diagnostic$converged))
    expect_lte(diagnostic$marginal_error, diagnostic$tolerance)
  }
  expect_gt(tr$diagnostics[[ref]]$point_spread$mean_effective_points, 1)
  expect_true(is.finite(tr$diagnostics[[ref]]$point_spread$self_mass))
  expect_true(all(is.finite(
    tr$diagnostics[[ref]]$point_spread$target_diagonal_mass
  )))
})

test_that("point-spread epsilon calibration is opt-in and fully recorded", {
  features <- list(
    rbind(c(1, 0), c(0, 1), c(1, 1)),
    rbind(c(1, 0), c(0, 1), c(1, 1))
  )
  centroids <- list(
    cbind(0:2, 0, 0),
    cbind(0:2, 0, 0)
  )
  sizes <- list(rep(1, 3), rep(1, 3))
  mapper <- dkge_mapper_spec(
    "sinkhorn",
    epsilon = 0.05,
    lambda_emb = 1,
    lambda_spa = 0.5,
    warm_start = FALSE,
    epsilon_calibration = list(
      target_effective_points = 1.5,
      epsilon_grid = c(0.01, 0.05, 0.2),
      calibration_data_hash = "independent-calibration-fixture",
      source = "test-heldout-calibration"
    )
  )

  tr <- dkge:::.dkge_transport_to_medoid(
    mapper, list(1:3, 1:3), features, centroids, sizes, 1L,
    subject_ids = c("s1", "s2"),
    preprocessing = list(source = "independent_test")
  )
  for (diagnostic in tr$diagnostics) {
    calibration <- diagnostic$epsilon_calibration
    expect_true(calibration$enabled)
    expect_equal(calibration$epsilon_grid, c(0.01, 0.05, 0.2))
    expect_true(calibration$selected_epsilon %in% calibration$epsilon_grid)
    expect_identical(
      calibration$calibration_data_hash,
      "independent-calibration-fixture"
    )
    expect_identical(calibration$source, "test-heldout-calibration")
  }

  invalid <- mapper
  invalid$params$epsilon_calibration$calibration_data_hash <- ""
  expect_error(
    dkge:::.dkge_transport_to_medoid(
      invalid, list(1:3, 1:3), features, centroids, sizes, 1L,
      subject_ids = c("s1", "s2")
    ),
    "calibration_data_hash",
    class = "dkge_alignment_calibration_error"
  )
})
