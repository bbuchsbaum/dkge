library(testthat)

make_iterative_template_fixture <- local({
  cached <- NULL
  function(seed = 260850L) {
    if (!is.null(cached) && identical(cached$seed, seed)) return(cached)
    set.seed(seed)
    S <- 4L
    q <- 5L
    P <- 7L
    effects <- paste0("e", seq_len(q))
    parcels <- paste0("p", seq_len(P))
    K <- diag(seq(q, 1), q)
    dimnames(K) <- list(effects, effects)
    subjects <- lapply(seq_len(S), function(s) {
      B <- matrix(rnorm(q * P), q, P,
                  dimnames = list(effects, parcels))
      X <- matrix(rnorm((q + 6L) * q), q + 6L, q,
                  dimnames = list(NULL, effects))
      covariance <- diag(q)
      dimnames(covariance) <- list(effects, effects)
      dkge_subject(
        B, X, id = paste0("s", s), effect_noise_cov = covariance,
        residual_variance = stats::setNames(rep(1, P), parcels)
      )
    })
    fit <- suppressWarnings(dkge_fit(
      dkge_data(subjects), K = K, rank = 2,
      w_method = "none", effect_scaling = "pooled_design"
    ))
    contrast <- suppressWarnings(dkge_contrast(
      fit, c(1, -1, rep(0, q - 2L)), method = "loso", align = FALSE
    ))
    shared_channel <- matrix(
      rnorm(q * P), q, P,
      dimnames = dimnames(fit$Braw[[1L]])
    )
    independent <- lapply(seq_len(S), function(s) {
      shared_channel + matrix(
        rnorm(q * P, sd = 0.12), q, P,
        dimnames = dimnames(shared_channel)
      )
    })
    names(independent) <- fit$subject_ids
    features <- dkge_alignment_features(
      fit, contrast,
      independent_betas = independent,
      independent_data_hash = paste0("iterative-template-", seed)
    )
    centroids <- lapply(seq_len(S), function(s) {
      cbind(
        x = 10 * seq_len(P),
        y = 5 * sin(seq_len(P) / 2) + s / 2,
        z = 0
      )
    })
    names(centroids) <- fit$subject_ids
    sizes <- stats::setNames(lapply(seq_len(S), function(s) {
      seq(1, 2, length.out = P)
    }), fit$subject_ids)
    support <- dkge_reference_support(
      centroids[[2]],
      provenance = list(kind = "fixed_mni_grid", source = "test_fixture")
    )
    mapper <- dkge_mapper_spec(
      "sinkhorn", epsilon = 0.1, lambda_emb = 0.8, lambda_spa = 0.2,
      sigma_mm = 15, warm_start = FALSE
    )
    template <- dkge_fit_functional_template(
      support, features, centroids, sizes = sizes, mapper = mapper,
      max_iter = 30L, tolerance = 1e-3, objective_tolerance = 1e-3
    )
    cached <<- list(
      seed = seed, fit = fit, contrast = contrast, features = features,
      centroids = centroids, sizes = sizes, support = support,
      mapper = mapper, template = template
    )
    cached
  }
})

test_that("iterative template exposes convergence, scale, rank, and plan diagnostics", {
  fx <- make_iterative_template_fixture()
  template <- fx$template
  expect_s3_class(template, "dkge_functional_template")
  expect_true(template$fitting$converged)
  expect_null(template$fitting$failure_reason)
  expect_true(template$fitting$all_subjects_refit_after_final_update)
  expect_true(template$fitting$no_identity_branch)
  expect_identical(template$fitting$subject_weighting, "equal_subject")
  expect_equal(template$fitting$subject_weights,
               stats::setNames(rep(1, 4), fx$fit$subject_ids))
  expect_identical(template$alignment_features_hash,
                   fx$features$structural_hash)
  expect_true(template$provenance$functional_features_learned_on_support)
  expect_true(template$provenance$bare_coordinates_are_not_correspondence)

  trajectory <- template$fitting$trajectory
  expect_true(all(is.finite(trajectory$objective)))
  expect_true(all(is.finite(trajectory$template_delta)))
  expect_equal(sqrt(rowSums(template$features^2)),
               rep(1, nrow(template$features)), tolerance = 1e-12)
  expect_true(all(trajectory$normalized_feature_rms > 0))
  expect_true(all(trajectory$effective_rank > 1))
  expect_true(template$fitting$final_feature_stats$effective_rank > 1)
  expect_true(template$fitting$rank_gate$passed)
  expect_gte(
    template$fitting$rank_gate$rank_tolerance,
    fx$features$control$rank_tolerance
  )

  diagnostics <- template$fitting$diagnostics
  expect_length(diagnostics, length(fx$fit$subject_ids))
  expect_true(all(vapply(diagnostics, function(x) {
    isTRUE(x$converged) && is.finite(x$marginal_error) &&
      is.finite(x$plan_entropy) &&
      is.finite(x$point_spread$mean_effective_points) &&
      is.finite(x$amplitude_preservation$mapped_to_target_rms_ratio)
  }, logical(1))))
  # Subject 2 supplies the display coordinates but receives no identity branch.
  expect_false(isTRUE(all.equal(
    template$fitting$operators[[2]], diag(fx$support$n_locations),
    tolerance = 0
  )))
})

test_that("template iteration and neighbourhood counts require whole numbers", {
  fx <- make_iterative_template_fixture()
  expect_error(
    dkge_fit_functional_template(
      fx$support, fx$features, fx$centroids, sizes = fx$sizes,
      mapper = fx$mapper, max_iter = 1.5
    ),
    "max_iter.*whole number"
  )
  expect_error(
    dkge_fit_functional_template(
      fx$support, fx$features, fx$centroids, sizes = fx$sizes,
      mapper = fx$mapper, knn_k = 2.5, max_iter = 1L
    ),
    "positive `k`"
  )
})

test_that("selected support provenance survives template fitting", {
  fx <- make_iterative_template_fixture(seed = 260851L)
  selection <- dkge_select_reference_subject(
    fx$centroids, method = "explicit", reference_subject = "s2",
    sizes = fx$sizes,
    provenance = list(fixed_before_template_fit = TRUE)
  )
  support <- dkge_reference_support(
    fx$centroids[[2L]],
    provenance = list(kind = "selected_subject_support")
  )
  template <- dkge_fit_functional_template(
    support, fx$features, fx$centroids, sizes = fx$sizes,
    mapper = fx$mapper, reference_selection = selection,
    max_iter = 2L, tolerance = 1, objective_tolerance = 1
  )
  expect_identical(template$reference_selection$structural_hash,
                   selection$structural_hash)
  aligned <- dkge_align_to_template(
    template, support, fx$contrast
  )
  fitted <- attr(aligned, "fitted_alignment")
  expect_identical(fitted$reference_selection$structural_hash,
                   selection$structural_hash)

  expect_error(
    dkge_fit_functional_template(
      support, fx$features, fx$centroids, sizes = fx$sizes,
      mapper = fx$mapper, max_iter = 2L
    ),
    "reference_selection|selected.*support",
    class = "dkge_template_reference_error"
  )
  moved <- support
  moved$coordinates[1L, 1L] <- moved$coordinates[1L, 1L] + 1
  moved$structural_hash <- dkge:::.dkge_object_hash(
    dkge:::.dkge_reference_support_payload(moved)
  )
  expect_error(
    dkge_fit_functional_template(
      moved, fx$features, fx$centroids, sizes = fx$sizes,
      mapper = fx$mapper, reference_selection = selection,
      max_iter = 2L
    ),
    "selected subject|coordinates|support",
    class = "dkge_template_reference_error"
  )

  initializer_without_selection <- dkge_functional_template(
    support, fx$template$features, feature_source = "independent",
    provenance = list(independent_data_hash = "separate-initializer")
  )
  expect_error(
    dkge_fit_functional_template(
      support, fx$features, fx$centroids, sizes = fx$sizes,
      mapper = fx$mapper, reference_selection = selection,
      initialization = "supplied_independent",
      initial_template = initializer_without_selection, max_iter = 2L
    ),
    "different reference-selection provenance",
    class = "dkge_template_reference_error"
  )
})

test_that("descriptive reference selection cannot be laundered by a template", {
  fx <- make_iterative_template_fixture(seed = 260852L)
  selection <- dkge_select_reference_subject(
    fx$centroids, alignment_features = fx$features,
    sizes = fx$sizes, mapper = fx$mapper,
    method = "descriptive_training"
  )
  support <- dkge_reference_support(
    fx$centroids[[selection$reference_subject]],
    provenance = list(kind = "selected_subject_support")
  )
  template <- dkge_fit_functional_template(
    support, fx$features, fx$centroids, sizes = fx$sizes,
    mapper = fx$mapper, reference_selection = selection,
    max_iter = 2L, tolerance = 1, objective_tolerance = 1
  )
  expect_identical(template$eligibility$status, "ineligible")
  values <- lapply(fx$features$features, function(x) rep(1, nrow(x)))
  aligned <- dkge_align_to_template(template, support, values)
  expect_error(
    dkge_infer_aligned(
      aligned, n_perm = 100L, allow_approximate_alignment = TRUE
    ),
    class = "dkge_alignment_ineligible_error"
  )
})

test_that("template fitting fails closed on Sinkhorn and rank collapse", {
  fx <- make_iterative_template_fixture(seed = 260853L)
  failed_mapper <- fx$mapper
  failed_mapper$params$epsilon <- 1e-4
  failed_mapper$params$max_iter <- 1L
  failed_mapper$params$tol <- 1e-14
  expect_error(
    suppressWarnings(dkge_fit_functional_template(
      fx$support, fx$features, fx$centroids, sizes = fx$sizes,
      mapper = failed_mapper, max_iter = 2L
    )),
    "did not converge|numerical",
    class = "dkge_template_numerical_error"
  )

  rank_one <- outer(
    seq_len(fx$support$n_locations),
    rep(1, ncol(fx$template$features))
  )
  rank_one <- rank_one / sqrt(rowSums(rank_one^2))
  initializer <- dkge_functional_template(
    fx$support, rank_one, feature_source = "independent",
    provenance = list(independent_data_hash = "rank-one-initializer")
  )
  expect_error(
    dkge_fit_functional_template(
      fx$support, fx$features, fx$centroids, sizes = fx$sizes,
      mapper = fx$mapper, initialization = "supplied_independent",
      initial_template = initializer, max_iter = 2L
    ),
    "rank|effective rank",
    class = "dkge_template_rank_error"
  )
})

test_that("mass-aware template update matches a direct plan oracle", {
  fx <- make_iterative_template_fixture()
  fitting <- fx$template$fitting
  cycle <- dkge:::.dkge_template_cycle(
    fx$template$features,
    fitting$source_features,
    fitting$source_centroids,
    fitting$source_sizes,
    fx$support,
    fitting$support_masses,
    fitting$subject_weights,
    fitting$mapper_spec
  )
  numerator <- matrix(0, nrow(fx$template$features),
                      ncol(fx$template$features))
  denominator <- numeric(nrow(fx$template$features))
  for (s in seq_along(cycle$plans)) {
    plan <- cycle$plans[[s]]
    weight <- fitting$subject_weights[[s]]
    numerator <- numerator + weight *
      crossprod(plan, fitting$source_features[[s]])
    denominator <- denominator + weight * colSums(plan)
  }
  expected <- sweep(numerator, 1L, denominator, "/")
  expect_equal(cycle$update, expected, tolerance = 1e-12)
})

test_that("fitted template creates aligned subject rows on the fixed support", {
  fx <- make_iterative_template_fixture()
  values <- lapply(seq_along(fx$features$features), function(s) {
    rep(7, nrow(fx$features$features[[s]]))
  })
  names(values) <- fx$fit$subject_ids
  aligned <- dkge_align_to_template(fx$template, fx$support, values)
  expect_s3_class(aligned, "dkge_aligned_maps")
  expect_identical(aligned$support_hash, fx$support$structural_hash)
  expect_identical(aligned$subject_weighting, "equal_subject")
  expect_equal(unname(aligned$values[[1]]),
               matrix(7, 4, fx$support$n_locations),
               tolerance = 1e-3)
  expect_identical(aligned$construction, "public_descriptive")
  expect_identical(aligned$eligibility$status, "ineligible")
  expect_null(aligned$application_receipt)
  expect_null(attr(aligned, "fitted_alignment"))
  expect_error(
    dkge_infer_aligned(
      aligned, n_perm = 100L, allow_approximate_alignment = TRUE
    ),
    class = "dkge_alignment_ineligible_error"
  )

  typed <- dkge_align_to_template(
    fx$template, fx$support, fx$contrast
  )
  expect_identical(typed$construction, "operator_applied_internal")
  expect_identical(typed$eligibility$status, "approximate")
  expect_identical(
    typed$application_receipt$contrast_binding,
    dkge:::.dkge_alignment_contrast_binding(fx$contrast)
  )
  expect_s3_class(
    dkge_infer_aligned(
      typed,
      inference = "parametric",
      correction = "none",
      allow_approximate_alignment = TRUE
    ),
    "dkge_inference"
  )
  tampered_family <- typed
  tampered_family$application_receipt$contrast_family_binding <-
    paste0(tampered_family$application_receipt$contrast_family_binding, "-tampered")
  expect_error(
    dkge_infer_aligned(
      tampered_family,
      n_perm = 100L,
      allow_approximate_alignment = TRUE
    ),
    "mutated|provenance",
    class = "dkge_aligned_maps_error"
  )

  fitted <- attr(typed, "fitted_alignment")
  expect_s3_class(fitted, "dkge_fitted_alignment")
  expect_null(fitted$reference_subject)
  expect_null(fitted$medoid)
  expect_identical(fitted$reference_method, "iterative_group_template")
  expect_false(fitted$reference_is_medoid)
  expect_output(print(fitted), "group functional template", fixed = TRUE)
})

test_that("independent template application retains its exact contrast binding", {
  fx <- make_iterative_template_fixture()
  other_contrast <- suppressWarnings(dkge_contrast(
    fx$fit, c(0, 1, -1, rep(0, nrow(fx$fit$K) - 3L)),
    method = "loso", align = FALSE
  ))
  expect_error(
    dkge_align_to_template(fx$template, fx$support, other_contrast),
    "different typed contrast result or family",
    class = "dkge_alignment_cache_mismatch"
  )

  fitted <- dkge:::.dkge_new_template_alignment(fx$template, fx$support)
  expect_error(
    dkge:::.dkge_apply_fitted_alignment(
      other_contrast$values,
      fitted,
      contrast_obj = other_contrast,
      application_context = "binding_adversary"
    ),
    "different typed contrast result or family",
    class = "dkge_alignment_cache_mismatch"
  )
  happy <- dkge:::.dkge_apply_fitted_alignment(
    fx$contrast$values,
    fitted,
    contrast_obj = fx$contrast,
    application_context = "binding_happy_path"
  )
  expect_s3_class(happy, "dkge_aligned_maps")

  laundered <- happy
  receipt <- laundered$application_receipt
  receipt$contrast_binding <-
    dkge:::.dkge_alignment_contrast_binding(other_contrast)
  receipt$contrast_family_binding <-
    dkge:::.dkge_alignment_contrast_family_binding(other_contrast)
  receipt$receipt_hash <- NULL
  receipt$receipt_hash <- dkge:::.dkge_object_hash(receipt)
  laundered$application_receipt <- receipt
  laundered$structural_hash <- dkge:::.dkge_object_hash(
    dkge:::.dkge_aligned_maps_payload(laundered)
  )
  expect_error(
    dkge_infer_aligned(
      laundered,
      inference = "parametric",
      correction = "none",
      allow_approximate_alignment = TRUE
    ),
    "bindings do not match",
    class = "dkge_aligned_maps_error"
  )
})

test_that("template subject weights are matched by exact subject identity", {
  fx <- make_iterative_template_fixture(seed = 260857L)
  weights <- stats::setNames(seq_along(fx$fit$subject_ids),
                             fx$fit$subject_ids)
  fit_with_weights <- function(subject_weights) {
    dkge_fit_functional_template(
      fx$support, fx$features, fx$centroids, sizes = fx$sizes,
      mapper = fx$mapper, subject_weights = subject_weights,
      max_iter = 2L, tolerance = 1, objective_tolerance = 1
    )
  }
  canonical <- fit_with_weights(weights)
  reversed <- fit_with_weights(rev(weights))
  expect_identical(reversed$fitting$subject_weights,
                   canonical$fitting$subject_weights)
  expect_equal(reversed$features, canonical$features, tolerance = 0)
  expect_identical(reversed$fitting$subject_weighting, "explicit_fixed")

  expect_error(
    fit_with_weights(stats::setNames(
      weights, paste0("alien", seq_along(weights))
    )),
    "identify every fitted subject exactly once",
    class = "dkge_template_error"
  )
  duplicate_names <- fx$fit$subject_ids
  duplicate_names[[length(duplicate_names)]] <- duplicate_names[[1L]]
  expect_error(
    fit_with_weights(stats::setNames(weights, duplicate_names)),
    "identify every fitted subject exactly once",
    class = "dkge_template_error"
  )
})

test_that("supplied initialization must be verified and independently identified", {
  fx <- make_iterative_template_fixture()
  supplied <- dkge_functional_template(
    fx$support,
    fx$template$features,
    feature_source = "independent",
    provenance = list(independent_data_hash = "separate-initial-template")
  )
  one_step <- dkge_fit_functional_template(
    fx$support, fx$features, fx$centroids, sizes = fx$sizes,
    mapper = fx$mapper,
    initialization = "supplied_independent",
    initial_template = supplied,
    max_iter = 1L
  )
  expect_identical(
    one_step$fitting$initialization$method,
    "supplied_verified_independent_template"
  )
  expect_identical(one_step$fitting$initialization$template_hash,
                   supplied$structural_hash)

  unverified <- dkge_functional_template(
    fx$support,
    fx$template$features,
    feature_source = "independent"
  )
  expect_error(
    dkge_fit_functional_template(
      fx$support, fx$features, fx$centroids, sizes = fx$sizes,
      mapper = fx$mapper,
      initialization = "supplied_independent",
      initial_template = unverified,
      max_iter = 1L
    ),
    "independent data|provenance receipt",
    class = "dkge_template_initialization_error"
  )
})

test_that("non-convergence is durable and fails closed for aligned maps", {
  fx <- make_iterative_template_fixture()
  unfinished <- dkge_fit_functional_template(
    fx$support, fx$features, fx$centroids, sizes = fx$sizes,
    mapper = fx$mapper, max_iter = 1L
  )
  expect_false(unfinished$fitting$converged)
  expect_identical(unfinished$fitting$failure_reason,
                   "maximum_iterations_reached")
  expect_identical(unfinished$eligibility$status, "ineligible")
  values <- lapply(fx$features$features, function(x) rep(1, nrow(x)))
  expect_error(
    dkge_align_to_template(unfinished, fx$support, values),
    "did not converge",
    class = "dkge_alignment_ineligible_error"
  )
  inspected <- dkge_align_to_template(
    unfinished, fx$support, values, allow_nonconverged = TRUE
  )
  expect_identical(inspected$eligibility$status, "ineligible")
})

test_that("same-data residualized templates retain their approximate label", {
  fx <- make_iterative_template_fixture()
  same_data <- dkge_alignment_features(
    fx$fit, fx$contrast, feature_source = "same_data_residualized"
  )
  approximate <- dkge_fit_functional_template(
    fx$support, same_data, fx$centroids, sizes = fx$sizes,
    mapper = fx$mapper,
    max_iter = 2L, tolerance = 1, objective_tolerance = 1
  )
  expect_true(approximate$fitting$converged)
  expect_identical(approximate$feature_source, "same_data_residualized")
  expect_identical(approximate$eligibility$status, "approximate")
  expect_false(approximate$eligibility$eligible)

  aligned <- dkge_align_to_template(
    approximate, fx$support, fx$contrast
  )
  expect_s3_class(aligned, "dkge_aligned_maps")
  other_contrast <- suppressWarnings(dkge_contrast(
    fx$fit, c(0, 1, -1, rep(0, nrow(fx$fit$K) - 3L)),
    method = "loso", align = FALSE
  ))
  expect_error(
    dkge_align_to_template(approximate, fx$support, other_contrast),
    "different typed contrast result or family",
    class = "dkge_alignment_cache_mismatch"
  )
  expect_error(
    dkge_align_to_template(
      approximate, fx$support, other_contrast$values
    ),
    "requires its bound",
    class = "dkge_alignment_cache_mismatch"
  )
})

test_that("public template constructor cannot mint fitted eligibility", {
  fx <- make_iterative_template_fixture(seed = 260856L)
  expect_error(
    dkge_functional_template(
      fx$support,
      fx$template$features,
      masses = fx$template$masses,
      feature_source = "independent",
      provenance = fx$template$provenance,
      fitting = fx$template$fitting,
      eligibility = fx$template$eligibility,
      alignment_features_hash = fx$template$alignment_features_hash,
      subject_ids = fx$template$subject_ids
    ),
    "internal receipts",
    class = "dkge_functional_template_error"
  )
})

test_that("a bare MNI grid is support, not functional correspondence", {
  fx <- make_iterative_template_fixture()
  bare <- dkge_functional_template(
    fx$support,
    matrix(0, fx$support$n_locations, 1L),
    feature_source = "geometry_only",
    provenance = list(kind = "coordinates_only")
  )
  values <- lapply(fx$features$features, function(x) rep(1, nrow(x)))
  expect_error(
    dkge_align_to_template(bare, fx$support, values),
    "operators are unavailable",
    class = "dkge_template_error"
  )
  expect_identical(fx$template$provenance$support_kind, "fixed_mni_grid")
  expect_true(fx$template$provenance$functional_features_learned_on_support)
})

test_that("template state is immutable", {
  fx <- make_iterative_template_fixture()
  mutated <- fx$template
  mutated$fitting$trajectory$objective[[1]] <-
    mutated$fitting$trajectory$objective[[1]] + 1
  expect_error(
    print(mutated), "mutated", class = "dkge_functional_template_error"
  )
})
