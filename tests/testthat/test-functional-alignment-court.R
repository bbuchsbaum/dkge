library(testthat)

court_path <- dkge_validation_path("functional-alignment-court", "court.R")
source(court_path, local = TRUE)

test_that("frozen court protocol and factor grid cover every declared axis", {
  protocol <- dkge_validation_path(
    "functional-alignment-court", "frozen-protocol.json"
  )
  expect_true(file.exists(protocol))
  expect_match(readLines(protocol, warn = FALSE),
               "full_permutation_reestimation", all = FALSE)
  grid <- fa_court_grid()
  expect_true(all(c(
    "S", "estimation_rank", "kernel_rank",
    "contrast_family_dimension", "eigengap", "nuisance_magnitude",
    "spatial_correlation", "spatial_heteroscedasticity",
    "parcel_count_heterogeneity", "epsilon", "functional_cost_weight",
    "dgp", "exact_oracle"
  ) %in% names(grid)))
  expect_gt(length(unique(grid$S)), 1L)
  expect_gt(length(unique(grid$estimation_rank)), 1L)
  expect_gt(length(unique(grid$kernel_rank)), 1L)
  expect_gt(length(unique(grid$contrast_family_dimension)), 1L)
  expect_true(any(grid$exact_oracle))
  expect_match(readLines(protocol, warn = FALSE),
               '"minimum": {"n_sim": 40', all = FALSE, fixed = TRUE)
})

test_that("court provenance hashes the live package source tree", {
  fixture <- tempfile("dkge-court-source-")
  dir.create(file.path(fixture, "R"), recursive = TRUE)
  writeLines("Package: dkge", file.path(fixture, "DESCRIPTION"))
  writeLines("exportPattern(\"^[[:alpha:]]+\")", file.path(fixture, "NAMESPACE"))
  writeLines("court_fixture <- 1", file.path(fixture, "R", "fixture.R"))
  on.exit(unlink(fixture, recursive = TRUE), add = TRUE)

  first <- fa_court_source_tree_hash(fixture)
  second <- fa_court_source_tree_hash(fixture)
  expect_identical(first, second)
  expect_match(first, "^[[:xdigit:]]{64}$")
  writeLines("court_fixture <- 2", file.path(fixture, "R", "fixture.R"))
  expect_false(identical(first, fa_court_source_tree_hash(fixture)))
})

test_that("court DGP and sign action are seeded, symmetric, and involutive", {
  cfg <- fa_court_base_config()
  first <- fa_court_generate(cfg, seed = 811)
  second <- fa_court_generate(cfg, seed = 811)
  expect_identical(first, second)
  sign <- rep(c(-1, 1), length.out = cfg$S)
  once <- fa_court_apply_raw_sign(first$B, first$tested_rows, sign)
  twice <- fa_court_apply_raw_sign(once, first$tested_rows, sign)
  expect_equal(twice, first$B, tolerance = 0)
  for (s in seq_len(cfg$S)) {
    nuisance <- setdiff(seq_len(cfg$q), first$tested_rows)
    expect_identical(once[[s]][nuisance, , drop = FALSE],
                     first$B[[s]][nuisance, , drop = FALSE])
  }
})

test_that("court exhaustive sign representatives exclude only identity", {
  signs <- fa_court_signs(5L, B = 15L, seed = 4L)
  expect_equal(dim(signs), c(5L, 15L))
  expect_true(all(signs[1L, ] == 1))
  expect_false(any(colSums(signs == 1) == nrow(signs)))
  expect_true(any(apply(signs, 2L, identical, c(1, rep(-1, 4L)))))
  complete <- cbind(identity = rep(1, 5L), signs)
  expect_equal(nrow(unique(t(complete))), 2^(5L - 1L))
})

test_that("kernel-image prototype reconstructs v and removes the family span", {
  cfg <- fa_court_base_config()
  cfg$contrast_family_dimension <- 2L
  dat <- fa_court_generate(cfg, seed = 812)
  fit <- fa_court_fit(dat, cfg)
  contrast <- dkge_contrast(fit, dat$contrasts, method = "loso", align = FALSE)
  features <- fa_court_residual_features(fit, contrast)
  diagnostics <- attr(features, "diagnostics")
  expect_true(all(vapply(diagnostics, function(x) {
    x$reconstruction_error < 1e-8
  }, logical(1))))
  expect_true(all(vapply(diagnostics, function(x) {
    x$orthogonality_error < 1e-8
  }, logical(1))))
  expect_true(all(vapply(diagnostics, function(x) {
    x$available_dimension == cfg$kernel_rank - cfg$contrast_family_dimension
  }, logical(1))))
})

test_that("tiny court run exposes every frozen arm and exact raw-refit arm", {
  cfg <- fa_court_base_config()
  cfg$S <- 6L
  cfg$base_parcels <- 8L
  cfg$cell <- "unit"
  cfg$cell_index <- 0L
  result <- fa_court_run_one(
    cfg, seed = 813, n_perm = 3L, include_full = TRUE, full_first = TRUE
  )
  expect_setequal(result$arm, c(
    "legacy_fullfit_functional", "l1_fixed_fullfit_functional",
    "l2_fold_loading", "geometry_only",
    "kernel_image_residual_prototype", "independent_alignment",
    "full_permutation_reestimation"
  ))
  expect_true(all(is.finite(result$p_fwer)))
  expect_true(all(result$convergence >= 0 & result$convergence <= 1))
  exact <- result[result$arm == "full_permutation_reestimation", ]
  expect_gt(exact$exact_elapsed_seconds, 0)
})
