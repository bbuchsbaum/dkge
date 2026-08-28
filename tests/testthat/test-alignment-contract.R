library(testthat)

make_contract_alignment <- function(
    feature_source = "independent",
    estimator_source = "same_data_rank_truncated",
    recompute_under_null = FALSE,
    verified = TRUE,
    S = 6L,
    P = 4L) {
  set.seed(260826)
  fx <- make_small_fit(S = S, q = 3L, P = P, T = 18L, rank = 2L,
                       seed = 260826)
  fit <- fx$fit
  ids <- fit$subject_ids
  contrast <- NULL
  features <- NULL
  centroids <- lapply(seq_len(S), function(s) {
    cbind(x = seq_len(P), y = s / 10, z = 0)
  })
  names(centroids) <- ids
  loadings <- dkge_predict_loadings(fit, fx$betas)
  if (identical(feature_source, "descriptive_adaptive") || !verified) {
    alignment <- dkge_prepare_transport(
      fit, centroids = centroids, loadings = loadings,
      mapper = dkge_mapper_spec("ridge", lambda = 1e-3), medoid = 2L
    )
  } else {
    stopifnot(identical(feature_source, "independent"),
              identical(estimator_source, "same_data_rank_truncated"),
              !isTRUE(recompute_under_null))
    contrast <- suppressWarnings(dkge_contrast(
      fit, c(1, -1, 0), method = "loso", align = FALSE
    ))
    independent <- lapply(seq_len(S), function(s) {
      matrix(
        rnorm(3L * P), 3L, P,
        dimnames = dimnames(fit$Braw[[s]])
      )
    })
    names(independent) <- ids
    features <- dkge_alignment_features(
      fit, contrast,
      independent_betas = independent,
      independent_data_hash = "independent-fixture-260826"
    )
    alignment <- dkge_prepare_alignment(
      fit, features, centroids,
      mapper = dkge_mapper_spec("ridge", lambda = 1e-3),
      reference_subject = ids[[2L]]
    )
  }
  list(
    fit = fit,
    contrast = contrast,
    features = features,
    loadings = loadings,
    centroids = centroids,
    alignment = alignment
  )
}

apply_contract_alignment <- function(fx, subject_weights = NULL) {
  dkge:::.dkge_apply_fitted_alignment(
    fx$contrast$values,
    fx$alignment,
    contrast_obj = fx$contrast,
    subject_weights = subject_weights,
    application_context = "alignment_contract_test"
  )
}

test_that("support and functional template are distinct immutable objects", {
  xyz <- cbind(x = 1:4, y = 0, z = 0)
  support <- dkge_reference_support(
    xyz, labels = paste0("a", 1:4),
    provenance = list(space = "MNI", content = "coordinates_only")
  )
  expect_s3_class(support, "dkge_reference_support")
  expect_false("features" %in% names(support))
  expect_match(capture.output(print(support)), "functional: none", all = FALSE)

  template <- dkge_functional_template(
    support,
    features = matrix(seq_len(12), 4, 3),
    masses = c(1, 2, 3, 4),
    provenance = list(independent_data_hash = "external-training-dataset")
  )
  expect_s3_class(template, "dkge_functional_template")
  expect_identical(template$support_hash, support$structural_hash)
  expect_true(template$source_verified)

  mutated_template <- template
  mutated_template$support_id <- "substituted-support"
  expect_error(
    print(mutated_template),
    "mutated",
    class = "dkge_functional_template_error"
  )

  mutated <- support
  mutated$coordinates[1, 1] <- -99
  expect_error(
    dkge_functional_template(mutated, matrix(0, 4, 1)),
    class = "dkge_reference_support_error"
  )
})

test_that("eligibility combines feature and estimator provenance", {
  independent_latent <- dkge_alignment_eligibility(
    "independent", "same_data_rank_truncated", source_verified = TRUE
  )
  expect_identical(independent_latent$status, "approximate")
  expect_false(independent_latent$exact)

  residual <- dkge_alignment_eligibility(
    "same_data_residualized", "same_data_rank_truncated",
    source_verified = TRUE
  )
  expect_identical(residual$status, "approximate")

  unverified <- dkge_alignment_eligibility(
    "independent", "independent", source_verified = FALSE
  )
  expect_identical(unverified$status, "ineligible")

  exact <- dkge_alignment_eligibility(
    "fully_recomputed", "fully_recomputed",
    recompute_under_null = TRUE, source_verified = TRUE
  )
  expect_true(exact$eligible)
  expect_true(exact$exact)
})

test_that("fitted alignment owns typed contract objects and solution hashes", {
  fx <- make_contract_alignment()
  alignment <- fx$alignment
  expect_s3_class(alignment, "dkge_fitted_alignment")
  expect_s3_class(alignment$reference_support, "dkge_reference_support")
  expect_s3_class(alignment$functional_template, "dkge_functional_template")
  expect_s3_class(alignment$eligibility, "dkge_alignment_eligibility")
  expect_identical(alignment$eligibility$status, "approximate")
  expect_true(all(c("operators", "plans", "diagnostics") %in%
                    names(alignment$solution_hashes)))

  bad <- alignment
  bad$operators[[1]][1, 1] <- bad$operators[[1]][1, 1] + 0.25
  values <- list(effect = matrix(rnorm(6 * 4), 6, 4))
  expect_error(
    dkge_aligned_maps(values, bad),
    "solution was mutated",
    class = "dkge_alignment_cache_mismatch"
  )
})

test_that("aligned maps preserve rows, support, estimand, and explicit weights", {
  fx <- make_contract_alignment()
  Y <- matrix(seq_len(24), 6, 4)
  aligned <- dkge_aligned_maps(
    list(c1 = Y), fx$alignment,
    estimand = list(target = "mean aligned activation")
  )
  expect_s3_class(aligned, "dkge_aligned_maps")
  expect_equal(dim(aligned$values$c1), c(6, 4))
  expect_identical(rownames(aligned$values$c1), fx$fit$subject_ids)
  expect_identical(aligned$subject_weighting, "equal_subject")
  expect_equal(unname(aligned$subject_weights), rep(1, 6))
  expect_identical(aligned$support_hash,
                   fx$alignment$reference_support$structural_hash)
  expect_identical(aligned$construction, "public_descriptive")
  expect_identical(aligned$eligibility$status, "ineligible")

  reversed <- Y[6:1, , drop = FALSE]
  rownames(reversed) <- rev(fx$fit$subject_ids)
  canonical <- dkge_aligned_maps(list(c1 = reversed), fx$alignment)
  expect_equal(unname(canonical$values$c1), unname(Y), tolerance = 0)
  expect_identical(rownames(canonical$values$c1), fx$fit$subject_ids)
  alien <- reversed
  rownames(alien) <- paste0("alien", seq_len(nrow(alien)))
  expect_error(
    dkge_aligned_maps(list(c1 = alien), fx$alignment),
    "different subjects",
    class = "dkge_aligned_maps_error"
  )

  weighted <- dkge_aligned_maps(
    list(c1 = Y), fx$alignment, subject_weights = 1:6
  )
  expect_identical(weighted$subject_weighting, "explicit_fixed")
  operator_bound_weighted <- apply_contract_alignment(
    fx, subject_weights = 1:6
  )
  expect_error(
    dkge_infer_aligned(
      operator_bound_weighted,
      n_perm = 100,
      allow_approximate_alignment = TRUE
    ),
    class = "dkge_alignment_weighting_error"
  )

  named_weights <- stats::setNames(1:6, fx$fit$subject_ids)
  named_canonical <- apply_contract_alignment(
    fx, subject_weights = named_weights
  )
  named_reversed <- apply_contract_alignment(
    fx, subject_weights = rev(named_weights)
  )
  expect_identical(named_reversed$subject_weights,
                   named_canonical$subject_weights)
  expect_error(
    apply_contract_alignment(
      fx,
      subject_weights = stats::setNames(1:6, paste0("alien", 1:6))
    ),
    "identify every fitted subject exactly once",
    class = "dkge_aligned_maps_error"
  )
  duplicate_names <- fx$fit$subject_ids
  duplicate_names[[6L]] <- duplicate_names[[1L]]
  expect_error(
    apply_contract_alignment(
      fx, subject_weights = stats::setNames(1:6, duplicate_names)
    ),
    "identify every fitted subject exactly once",
    class = "dkge_aligned_maps_error"
  )
})

test_that("subject IDs control cyclic row and typed-source ordering exactly", {
  fx <- make_contract_alignment()
  subject_ids <- fx$fit$subject_ids
  cycle <- c(2:6, 1L)
  canonical_values <- matrix(
    seq_len(6 * 4), 6, 4,
    dimnames = list(subject_ids, paste0("location", 1:4))
  )
  cycled_values <- canonical_values[cycle, , drop = FALSE]
  descriptive <- dkge_aligned_maps(
    list(c1 = cycled_values), fx$alignment,
    subject_ids = subject_ids
  )
  expect_equal(unname(descriptive$values$c1), unname(canonical_values),
               tolerance = 0)
  expect_identical(rownames(descriptive$values$c1), subject_ids)

  baseline <- apply_contract_alignment(fx)
  cycled_sources <- lapply(fx$contrast$values, function(values) {
    values[cycle]
  })
  observed <- dkge:::.dkge_apply_fitted_alignment(
    cycled_sources,
    fx$alignment,
    contrast_obj = fx$contrast,
    application_context = "cyclic_subject_order_test"
  )
  expect_identical(observed$values, baseline$values)

  invalid_name_sets <- list(
    partial = replace(subject_ids, 2L, ""),
    duplicate = replace(subject_ids, 6L, subject_ids[[1L]]),
    foreign = replace(subject_ids, 6L, "foreign-subject")
  )
  for (bad_names in invalid_name_sets) {
    bad_sources <- fx$contrast$values
    names(bad_sources[[1L]]) <- bad_names
    expect_error(
      dkge:::.dkge_apply_fitted_alignment(
        bad_sources,
        fx$alignment,
        contrast_obj = fx$contrast,
        application_context = "invalid_subject_names"
      ),
      "must match fitted subject IDs exactly",
      class = "dkge_alignment_feature_error"
    )
  }
})

test_that("aligned inference fails closed and labels any approximate override", {
  approx_fx <- make_contract_alignment(
    feature_source = "independent",
    estimator_source = "same_data_rank_truncated"
  )
  aligned <- apply_contract_alignment(approx_fx)
  expect_identical(aligned$construction, "operator_applied_internal")
  expect_identical(aligned$eligibility$status, "approximate")
  expect_error(
    dkge_infer_aligned(aligned, n_perm = 100),
    class = "dkge_alignment_ineligible_error"
  )
  result <- dkge_infer_aligned(
    aligned, n_perm = 100, allow_approximate_alignment = TRUE
  )
  expect_s3_class(result, "dkge_inference")
  expect_identical(result$metadata$alignment$status, "approximate")
  expect_true(result$metadata$alignment$approximate_override)
  mutated <- aligned
  mutated$values[[1L]][1L, 1L] <- mutated$values[[1L]][1L, 1L] + 1
  expect_error(
    dkge_infer_aligned(
      mutated, n_perm = 100, allow_approximate_alignment = TRUE
    ),
    "mutated",
    class = "dkge_aligned_maps_error"
  )
  substituted_support <- aligned
  substituted_support$support_id <- "substituted-support"
  expect_error(
    dkge_infer_aligned(
      substituted_support, n_perm = 100,
      allow_approximate_alignment = TRUE
    ),
    "mutated",
    class = "dkge_aligned_maps_error"
  )
  expect_error(
    dkge_infer_aligned(
      aligned,
      inference = "parametric",
      correction = "maxT",
      allow_approximate_alignment = TRUE
    ),
    "does not define a max-T null distribution",
    class = "dkge_inference_correction_error"
  )

  descriptive_fx <- make_contract_alignment(
    feature_source = "descriptive_adaptive",
    estimator_source = "descriptive",
    verified = FALSE
  )
  descriptive <- dkge_aligned_maps(
    list(c1 = matrix(rnorm(24), 6, 4)), descriptive_fx$alignment
  )
  expect_error(
    dkge_infer_aligned(
      descriptive, n_perm = 100, allow_approximate_alignment = TRUE
    ),
    class = "dkge_alignment_ineligible_error"
  )

  arbitrary <- dkge_aligned_maps(
    list(c1 = matrix(rnorm(24), 6, 4)), approx_fx$alignment
  )
  expect_identical(arbitrary$construction, "public_descriptive")
  expect_identical(arbitrary$eligibility$status, "ineligible")
  expect_null(arbitrary$application_receipt)
  expect_error(
    dkge_infer_aligned(
      arbitrary, n_perm = 100, allow_approximate_alignment = TRUE
    ),
    class = "dkge_alignment_ineligible_error"
  )
})

test_that("renderer consumes aligned values and never owns correspondence", {
  fx <- make_contract_alignment()
  aligned <- dkge_aligned_maps(
    list(c1 = matrix(seq_len(24), 6, 4)), fx$alignment
  )
  renderer <- dkge_renderer(fx$alignment$reference_support)
  expect_s3_class(renderer, "dkge_renderer")
  expect_false(any(c("operators", "plans", "feature_list", "mapper_fits") %in%
                     names(renderer)))
  rendered <- dkge_render_aligned(renderer, aligned, decode = FALSE)
  expect_equal(
    unname(rendered$support_values$c1),
    unname(colMeans(aligned$values$c1))
  )
  expect_identical(rendered$fitted_alignment_hash,
                   aligned$fitted_alignment_hash)

  mutated_renderer <- renderer
  mutated_renderer$support_id <- "substituted-support"
  expect_error(
    dkge_render_aligned(mutated_renderer, aligned),
    "identities do not match",
    class = "dkge_renderer_error"
  )

  other <- dkge_reference_support(
    fx$alignment$reference_support$coordinates + 1,
    labels = fx$alignment$reference_support$labels
  )
  expect_error(
    dkge_render_aligned(dkge_renderer(other), aligned),
    "different supports",
    class = "dkge_renderer_error"
  )
})
