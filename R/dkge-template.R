# dkge-template.R
# Iterative group functional templates on an identified display support.

.dkge_template_abort <- function(message, class = "dkge_template_error") {
  .dkge_abort(message, class)
}

.dkge_template_feature_stats <- function(features, masses,
                                         rank_tolerance = 1e-9) {
  features <- as.matrix(features)
  weights <- as.numeric(masses) / sum(masses)
  center <- colSums(features * weights)
  centered <- sweep(features, 2L, center, "-")
  scale <- sqrt(sum(weights * rowSums(centered^2)) / ncol(centered))
  covariance <- crossprod(centered * sqrt(weights), centered)
  eigenvalues <- pmax(
    eigen(covariance, symmetric = TRUE, only.values = TRUE)$values,
    0
  )
  threshold <- rank_tolerance * max(eigenvalues, 1)
  positive <- eigenvalues[eigenvalues > threshold]
  effective_rank <- if (!length(positive)) {
    0
  } else {
    sum(positive)^2 / sum(positive^2)
  }
  list(
    center = center,
    rms = scale,
    rank = length(positive),
    effective_rank = effective_rank,
    eigenvalues = eigenvalues
  )
}

.dkge_normalize_template_features <- function(features, masses,
                                               target_row_norm = 1,
                                               rank_tolerance = 1e-9) {
  raw_stats <- .dkge_template_feature_stats(
    features, masses, rank_tolerance = rank_tolerance
  )
  row_norms <- sqrt(rowSums(as.matrix(features)^2))
  if (any(!is.finite(row_norms)) || any(row_norms <= 1e-12)) {
    .dkge_template_abort(
      "A template feature row collapsed during normalization.",
      "dkge_template_rank_error"
    )
  }
  normalized <- as.matrix(features) * (target_row_norm / row_norms)
  normalized_stats <- .dkge_template_feature_stats(
    normalized, masses, rank_tolerance = rank_tolerance
  )
  list(
    features = normalized,
    raw = raw_stats,
    normalized = normalized_stats,
    scale_multiplier = target_row_norm / row_norms,
    rule = paste0(
      "target-row L2 normalization; every functional signature has norm ",
      format(target_row_norm, digits = 8)
    )
  )
}

.dkge_template_rank_requirements <- function(alignment_features) {
  control <- alignment_features$control %||%
    .dkge_resolve_alignment_feature_control(NULL)
  minimum_rank <- if (identical(
    alignment_features$feature_source, "same_data_residualized"
  )) {
    control$min_residual_rank
  } else {
    control$min_feature_rank
  }
  list(
    minimum_rank = as.integer(minimum_rank),
    minimum_effective_rank = as.numeric(control$min_effective_rank),
    rank_tolerance = as.numeric(control$rank_tolerance)
  )
}

.dkge_assert_template_feature_rank <- function(stats, requirements, stage) {
  passed <- is.finite(stats$rank) && is.finite(stats$effective_rank) &&
    stats$rank >= requirements$minimum_rank &&
    stats$effective_rank >= requirements$minimum_effective_rank
  if (!passed) {
    .dkge_template_abort(
      sprintf(
        paste0(
          "Template %s collapsed below its declared functional-rank gate: ",
          "rank %d (minimum %d), effective rank %.3f (minimum %.3f)."
        ),
        stage, stats$rank, requirements$minimum_rank,
        stats$effective_rank, requirements$minimum_effective_rank
      ),
      "dkge_template_rank_error"
    )
  }
  invisible(list(
    passed = TRUE,
    stage = stage,
    numerical_rank = stats$rank,
    effective_rank = stats$effective_rank,
    minimum_rank = requirements$minimum_rank,
    minimum_effective_rank = requirements$minimum_effective_rank,
    rank_tolerance = requirements$rank_tolerance
  ))
}

.dkge_validate_template_reference <- function(
    support, reference_selection, alignment_features,
    centroids, sizes, mapper_spec) {
  kind <- support$provenance$kind %||% ""
  selected_kinds <- c("selected_subject_support", "reference_subject_support")
  if (is.null(reference_selection)) {
    if (kind %in% selected_kinds) {
      .dkge_template_abort(
        paste0(
          "A selected/reference-subject support requires its typed ",
          "`reference_selection` receipt."
        ),
        "dkge_template_reference_error"
      )
    }
    fixed <- isTRUE(support$provenance$fixed_by_caller) ||
      startsWith(kind, "fixed_") ||
      kind %in% c("external_fixed_support", "fixed_mni_grid",
                  "fixed_group_support")
    if (!fixed) {
      .dkge_template_abort(
        paste0(
          "Selection-free template fitting requires explicit fixed-support ",
          "provenance (`kind = 'fixed_*'` or `fixed_by_caller = TRUE`)."
        ),
        "dkge_template_reference_error"
      )
    }
    return(invisible(NULL))
  }
  .dkge_validate_reference_selection(
    reference_selection,
    subject_ids = alignment_features$subject_ids,
    centroids = centroids,
    sizes = sizes,
    alignment_features = alignment_features,
    mapper_spec = mapper_spec
  )
  index <- reference_selection$reference_subject
  selected_coordinates <- as.matrix(centroids[[index]])
  same_coordinates <- identical(dim(support$coordinates),
                                dim(selected_coordinates)) &&
    isTRUE(all.equal(unname(support$coordinates),
                     unname(selected_coordinates), tolerance = 0))
  if (!same_coordinates) {
    .dkge_template_abort(
      "Template support coordinates do not match the selected reference subject.",
      "dkge_template_reference_error"
    )
  }
  declared_id <- support$provenance$subject_id %||%
    support$provenance$reference_subject_id %||% NULL
  if (!is.null(declared_id) &&
      !identical(as.character(declared_id),
                 reference_selection$reference_subject_id)) {
    .dkge_template_abort(
      "Template support subject provenance conflicts with reference selection.",
      "dkge_template_reference_error"
    )
  }
  invisible(reference_selection)
}

.dkge_standardize_template_sources <- function(features, sizes,
                                               subject_weights,
                                               rank_tolerance = 1e-9) {
  row_normalized <- lapply(features, function(x) {
    x <- as.matrix(x)
    norms <- sqrt(rowSums(x^2))
    if (any(!is.finite(norms)) || any(norms <= 1e-12)) {
      .dkge_template_abort(
        "Alignment features contain a collapsed source row.",
        "dkge_template_rank_error"
      )
    }
    x / norms
  })
  names(row_normalized) <- names(features)
  list(
    features = row_normalized,
    center = rep(0, ncol(row_normalized[[1]])),
    rms = 1 / sqrt(ncol(row_normalized[[1]])),
    target_row_norm = 1,
    rank_tolerance = rank_tolerance,
    rule = "source-row L2 normalization in the shared kernel-image gauge",
    structural_hash = .dkge_object_hash(list(
      subject_weights = subject_weights, sizes = sizes,
      rule = "source_row_l2"
    ))
  )
}

.dkge_initialize_template_knn <- function(features, centroids, sizes,
                                          support, subject_weights,
                                          k = 8L, sigma = 5) {
  Q <- support$n_locations
  valid_k <- is.numeric(k) && is.null(dim(k)) && length(k) == 1L &&
    !is.na(k) && is.finite(k) && k == floor(k) && k >= 1 &&
    k <= .Machine$integer.max
  if (!valid_k ||
      length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
    .dkge_template_abort("Spatial initialization requires positive `k` and `sigma`.")
  }
  k <- min(as.integer(k), Q)
  numerator <- matrix(0, Q, ncol(features[[1]]))
  denominator <- numeric(Q)
  mapper <- dkge_mapper("knn", k = k, sigx = sigma)
  for (s in seq_along(features)) {
    fitted <- fit_mapper(
      mapper,
      subj_points = centroids[[s]],
      anchor_points = support$coordinates
    )
    weighted_features <- features[[s]] * sizes[[s]]
    numerator <- numerator + subject_weights[[s]] * apply_mapper(
      fitted, weighted_features, normalize_by_reliab = FALSE
    )
    denominator <- denominator + subject_weights[[s]] * apply_mapper(
      fitted, sizes[[s]], normalize_by_reliab = FALSE
    )
  }
  covered <- denominator > 1e-12
  if (!any(covered)) {
    .dkge_template_abort(
      "Spatial initialization did not cover any target support row.",
      "dkge_template_initialization_error"
    )
  }
  initialized <- matrix(0, Q, ncol(numerator))
  initialized[covered, ] <- numerator[covered, , drop = FALSE] /
    denominator[covered]
  uncovered <- which(!covered)
  if (length(uncovered)) {
    nearest <- FNN::get.knnx(
      support$coordinates[covered, , drop = FALSE],
      support$coordinates[uncovered, , drop = FALSE],
      k = 1L
    )$nn.index[, 1L]
    initialized[uncovered, ] <- initialized[which(covered)[nearest], , drop = FALSE]
  }
  list(
    features = initialized,
    provenance = list(
      method = "one_pass_spatial_knn_pooling",
      k = k,
      sigma = sigma,
      mass_aware = TRUE,
      covered_rows = sum(covered),
      uncovered_rows_filled_from_nearest_covered = uncovered,
      mapper = mapper
    )
  )
}

.dkge_template_cycle <- function(template_features, source_features,
                                 centroids, sizes, support, support_masses,
                                 subject_weights, mapper_spec) {
  S <- length(source_features)
  Q <- nrow(template_features)
  d <- ncol(template_features)
  numerator <- matrix(0, Q, d)
  denominator <- numeric(Q)
  mappings <- vector("list", S)
  plans <- vector("list", S)
  operators <- vector("list", S)
  diagnostics <- vector("list", S)
  regularized_objective <- feature_distortion <- numeric(S)
  names(mappings) <- names(plans) <- names(operators) <- names(diagnostics) <-
    names(source_features)
  epsilon <- mapper_spec$params$epsilon %||% 0.05

  for (s in seq_len(S)) {
    mapping <- tryCatch(
      .dkge_fit_mapper_policy(
        mapper_spec,
        source_feat = source_features[[s]],
        target_feat = template_features,
        source_weights = sizes[[s]],
        target_weights = support_masses,
        source_xyz = centroids[[s]],
        target_xyz = support$coordinates,
        self_map = FALSE
      ),
      dkge_alignment_numerical_error = function(e) {
        .dkge_abort(
          sprintf("Template mapping for subject '%s' is numerically invalid: %s",
                  names(source_features)[[s]] %||% s, conditionMessage(e)),
          "dkge_template_numerical_error"
        )
      }
    )
    .dkge_assert_mapping_numerically_valid(
      mapping, mapper_spec,
      context = sprintf("Template mapping for subject '%s'",
                        names(source_features)[[s]] %||% s),
      class = "dkge_template_numerical_error"
    )
    plan <- as.matrix(mapping$plan)
    operator <- as.matrix(mapping$operator)
    mapped <- as.matrix(.dkge_apply_operator(operator, source_features[[s]]))
    target_mass <- colSums(plan)
    numerator <- numerator + subject_weights[[s]] *
      crossprod(plan, source_features[[s]])
    denominator <- denominator + subject_weights[[s]] * target_mass

    cost <- mapping$fit_info$cost
    positive <- plan > 0
    entropy_term <- sum(plan[positive] * (log(plan[positive]) - 1))
    regularized_objective[[s]] <- sum(plan * cost) + epsilon * entropy_term
    feature_distortion[[s]] <- sum(
      plan * .dkge_pairwise_sqdist(source_features[[s]], template_features)
    )
    numerical <- mapping$fit_info$diagnostics %||% list()
    operator_diagnostics <- .dkge_operator_diagnostics(
      operator, plan, source_features[[s]], template_features,
      self_map = FALSE
    )
    diagnostics[[s]] <- utils::modifyList(
      numerical,
      c(operator_diagnostics, list(
        plan_entropy = -sum(plan[positive] * log(plan[positive])),
        regularized_objective = regularized_objective[[s]],
        feature_distortion = feature_distortion[[s]],
        mapped_feature_rms = sqrt(mean(mapped^2))
      )),
      keep.null = TRUE
    )
    mappings[[s]] <- mapping
    plans[[s]] <- plan
    operators[[s]] <- operator
  }
  if (any(!is.finite(denominator)) || any(denominator <= 0)) {
    .dkge_template_abort(
      "A template update received an empty target marginal.",
      "dkge_template_numerical_error"
    )
  }
  update <- sweep(numerator, 1L, denominator, "/")
  weights <- subject_weights / sum(subject_weights)
  list(
    update = update,
    mappings = mappings,
    plans = plans,
    operators = operators,
    diagnostics = diagnostics,
    objective = sum(weights * regularized_objective),
    feature_distortion = sum(weights * feature_distortion)
  )
}

.dkge_template_eligibility <- function(alignment_features, converged,
                                       failure_reason = NULL,
                                       reference_selection = NULL,
                                       rank_gate = NULL) {
  eligibility <- dkge_alignment_eligibility(
    feature_source = alignment_features$feature_source,
    estimator_source = alignment_features$estimator_source,
    recompute_under_null = FALSE,
    source_verified = alignment_features$source_verified,
    reference_selection_status = reference_selection$eligibility$status %||%
      "eligible"
  )
  eligibility$solver_status <- "converged"
  eligibility$solver_converged <- TRUE
  eligibility$solver_diagnostics <- list(
    required = TRUE, passed = TRUE,
    context = "all iterative and final template mappings"
  )
  if (!isTRUE(converged)) {
    eligibility <- .dkge_invalidate_alignment_eligibility(
      eligibility,
      paste0("template optimization failed: ", failure_reason)
    )
  }
  if (!is.null(rank_gate) && !isTRUE(rank_gate$passed)) {
    eligibility <- .dkge_invalidate_alignment_eligibility(
      eligibility, "template functional-rank gate failed"
    )
  }
  eligibility
}

#' Fit an iterative group functional template on fixed support
#'
#' Fits all subjects to one identified display support while learning the
#' target's functional features from a typed alignment-feature channel. The
#' default initialization is the existing one-pass spatial kNN pooling. Each
#' iteration alternates fixed-policy Sinkhorn fits with the mass-aware
#' barycentric update that minimizes feature distortion conditional on the
#' plans. Every updated target signature is then returned to unit L2 norm,
#' preventing entropic averaging from silently turning off the functional cost
#' while preserving the shared kernel-image gauge. Every subject is refitted
#' once more against the converged template.
#'
#' A fixed anatomical or MNI support supplies coordinates only. This function
#' makes it a functional target by learning and recording template features on
#' that support; passing a bare grid directly is not functional alignment.
#'
#' @param support A [dkge_reference_support()] identifying the display target.
#' @param alignment_features Typed features from [dkge_alignment_features()].
#' @param reference_selection Optional typed selection receipt when `support`
#'   was derived from a cohort subject. Selection-free supports must declare
#'   fixed ancillary provenance.
#' @param centroids Named list of subject parcel coordinates.
#' @param sizes Optional positive subject parcel masses.
#' @param support_masses Optional positive target-support masses.
#' @param mapper Fixed Sinkhorn transport specification. Functional cost weight
#'   must be positive and epsilon calibration is not permitted inside template
#'   fitting.
#' @param initialization Either one-pass spatial kNN pooling or a supplied
#'   template on the same support with separately identified independent-data
#'   provenance.
#' @param initial_template Required for `"supplied_independent"`.
#' @param subject_weights Optional fixed positive subject weights. Named weights
#'   are reordered by exact subject ID; unnamed weights are positional in fitted
#'   subject order. Equal subject weighting is the default; fit-level MFA
#'   weights are never imported.
#' @param knn_k,knn_sigma Spatial initialization controls.
#' @param max_iter Maximum number of alternating updates.
#' @param tolerance Relative target-feature change required for convergence.
#' @param objective_tolerance Relative objective change required jointly with
#'   `tolerance` after at least two iterations.
#' @param update_rate Under-relaxation in `(0, 1]` toward each exact
#'   plan-conditional barycentric update. The default damps small deterministic
#'   Sinkhorn/template oscillations without changing the fixed point.
#' @param rank_tolerance Relative eigenvalue threshold for rank diagnostics.
#'   The alignment-feature object's declared tolerance is a non-relaxable floor;
#'   this argument may make the template gate stricter but not weaker.
#' @return An immutable `dkge_functional_template` containing final operators,
#'   plans, objective/scale/rank trajectories, numerical diagnostics,
#'   convergence status, and inferential eligibility.
#' @export
dkge_fit_functional_template <- function(
    support,
    alignment_features,
    centroids,
    reference_selection = NULL,
    sizes = NULL,
    support_masses = NULL,
    mapper = dkge_mapper_spec(
      "sinkhorn", epsilon = 0.05, lambda_emb = 1, lambda_spa = 0.5,
      sigma_mm = 15, warm_start = FALSE
    ),
    initialization = c("spatial_knn", "supplied_independent"),
    initial_template = NULL,
    subject_weights = NULL,
    knn_k = 8L,
    knn_sigma = 5,
    max_iter = 25L,
    tolerance = 1e-4,
    objective_tolerance = 1e-4,
    update_rate = 0.5,
    rank_tolerance = 1e-9) {
  .dkge_validate_reference_support(support)
  .dkge_validate_alignment_features(alignment_features)
  initialization <- match.arg(initialization)
  if (!inherits(mapper, "dkge_mapper_spec") ||
      !identical(mapper$strategy, "sinkhorn")) {
    .dkge_template_abort(
      "Iterative functional templates currently require a Sinkhorn mapper.",
      "dkge_template_mapper_error"
    )
  }
  if ((mapper$params$lambda_emb %||% 1) <= 0) {
    .dkge_template_abort(
      "A functional template requires a positive functional cost weight.",
      "dkge_template_mapper_error"
    )
  }
  if (!is.null(.dkge_epsilon_calibration(mapper))) {
    .dkge_template_abort(
      "Calibrate epsilon separately; template fitting requires one fixed policy.",
      "dkge_template_mapper_error"
    )
  }
  valid_max_iter <- is.numeric(max_iter) && is.null(dim(max_iter)) &&
    length(max_iter) == 1L && !is.na(max_iter) && is.finite(max_iter) &&
    max_iter == floor(max_iter) && max_iter >= 1 &&
    max_iter <= .Machine$integer.max
  if (!valid_max_iter) {
    .dkge_template_abort("`max_iter` must be one finite positive whole number.")
  }
  scalar_positive <- list(
    tolerance = tolerance,
    objective_tolerance = objective_tolerance,
    rank_tolerance = rank_tolerance
  )
  if (any(!vapply(scalar_positive, function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0
  }, logical(1)))) {
    .dkge_template_abort("Iteration controls must be finite and positive.")
  }
  if (!is.numeric(update_rate) || length(update_rate) != 1L ||
      !is.finite(update_rate) || update_rate <= 0 || update_rate > 1) {
    .dkge_template_abort("`update_rate` must lie in (0, 1].")
  }
  max_iter <- as.integer(max_iter)

  subject_ids <- alignment_features$subject_ids
  source_features <- .dkge_order_subject_list(
    alignment_features$features, subject_ids, "alignment features"
  )
  canonical <- .dkge_reference_centroids(centroids, subject_ids)
  centroids <- canonical$centroids
  sizes <- .dkge_reference_sizes(sizes, centroids, subject_ids)
  .dkge_validate_template_reference(
    support, reference_selection, alignment_features,
    centroids, sizes, mapper
  )
  dimensions <- vapply(source_features, ncol, integer(1))
  if (length(unique(dimensions)) != 1L ||
      any(vapply(seq_along(source_features), function(s) {
        nrow(source_features[[s]]) != nrow(centroids[[s]])
      }, logical(1)))) {
    .dkge_template_abort(
      "Feature dimensions or parcel rows do not match the supplied cohort geometry."
    )
  }
  if (ncol(support$coordinates) != ncol(centroids[[1]])) {
    .dkge_template_abort(
      "Source centroids and target support use different coordinate dimensions."
    )
  }
  if (is.null(support_masses)) support_masses <- rep(1, support$n_locations)
  support_masses <- as.numeric(support_masses)
  if (length(support_masses) != support$n_locations ||
      any(!is.finite(support_masses)) || any(support_masses <= 0)) {
    .dkge_template_abort("`support_masses` must be finite and positive.")
  }
  resolved_weights <- .dkge_resolve_subject_weights(
    subject_weights, subject_ids, "dkge_template_error"
  )
  subject_weights <- resolved_weights$weights

  rank_requirements <- .dkge_template_rank_requirements(alignment_features)
  rank_requirements$rank_tolerance <- max(
    rank_tolerance, rank_requirements$rank_tolerance
  )
  gate_rank_tolerance <- rank_requirements$rank_tolerance

  standardized <- .dkge_standardize_template_sources(
    source_features, sizes, subject_weights,
    rank_tolerance = gate_rank_tolerance
  )
  source_features <- standardized$features

  if (identical(initialization, "supplied_independent")) {
    if (is.null(initial_template)) {
      .dkge_template_abort(
        "Supplied initialization requires `initial_template`.",
        "dkge_template_initialization_error"
      )
    }
    .dkge_validate_functional_template(initial_template, support)
    initializer_selection_hash <-
      initial_template$reference_selection$structural_hash %||% NULL
    fitted_selection_hash <- reference_selection$structural_hash %||% NULL
    if (!identical(initializer_selection_hash, fitted_selection_hash)) {
      .dkge_template_abort(
        paste0(
          "The supplied initializer and template fit carry different ",
          "reference-selection provenance."
        ),
        "dkge_template_reference_error"
      )
    }
    if (!identical(initial_template$feature_source, "independent") ||
        !isTRUE(initial_template$source_verified)) {
      .dkge_template_abort(
        paste0(
          "A supplied initializer must declare independent data through a ",
          "valid provenance receipt."
        ),
        "dkge_template_initialization_error"
      )
    }
    if (ncol(initial_template$features) != dimensions[[1]]) {
      .dkge_template_abort(
        "The supplied template and source features have different dimensions.",
        "dkge_template_initialization_error"
      )
    }
    initialized <- list(
      features = initial_template$features,
      provenance = list(
        method = "supplied_verified_independent_template",
        template_hash = initial_template$structural_hash
      )
    )
  } else {
    if (!is.null(initial_template)) {
      .dkge_template_abort(
        "`initial_template` is used only with supplied-independent initialization."
      )
    }
    initialized <- .dkge_initialize_template_knn(
      source_features, centroids, sizes, support, subject_weights,
      k = knn_k, sigma = knn_sigma
    )
  }
  normalized <- .dkge_normalize_template_features(
    initialized$features, support_masses,
    target_row_norm = standardized$target_row_norm,
    rank_tolerance = gate_rank_tolerance
  )
  .dkge_assert_template_feature_rank(
    normalized$normalized, rank_requirements, "initialization"
  )
  current <- normalized$features
  initial_stats <- normalized$normalized
  trajectory <- vector("list", max_iter)
  previous_objective <- NA_real_
  converged <- FALSE
  iterations <- max_iter

  for (iteration in seq_len(max_iter)) {
    cycle <- .dkge_template_cycle(
      current, source_features, centroids, sizes, support, support_masses,
      subject_weights, mapper
    )
    candidate <- (1 - update_rate) * current + update_rate * cycle$update
    updated <- .dkge_normalize_template_features(
      candidate, support_masses,
      target_row_norm = standardized$target_row_norm,
      rank_tolerance = gate_rank_tolerance
    )
    .dkge_assert_template_feature_rank(
      updated$normalized, rank_requirements,
      sprintf("iteration %d", iteration)
    )
    target_mass <- support_masses / sum(support_masses)
    delta <- sqrt(sum(target_mass * rowSums(
      (updated$features - current)^2
    ))) / max(
      sqrt(sum(target_mass * rowSums(current^2))),
      .Machine$double.xmin
    )
    objective_change <- if (is.finite(previous_objective)) {
      abs(cycle$objective - previous_objective) /
        max(1, abs(previous_objective))
    } else {
      NA_real_
    }
    trajectory[[iteration]] <- data.frame(
      iteration = iteration,
      objective = cycle$objective,
      feature_distortion = cycle$feature_distortion,
      relative_objective_change = objective_change,
      template_delta = delta,
      raw_feature_rms = updated$raw$rms,
      normalized_feature_rms = updated$normalized$rms,
      normalized_row_norm_min = min(sqrt(rowSums(updated$features^2))),
      normalized_row_norm_mean = mean(sqrt(rowSums(updated$features^2))),
      normalized_row_norm_max = max(sqrt(rowSums(updated$features^2))),
      numerical_rank = updated$normalized$rank,
      effective_rank = updated$normalized$effective_rank
    )
    current <- updated$features
    if (iteration >= 2L && delta <= tolerance &&
        objective_change <= objective_tolerance) {
      converged <- TRUE
      iterations <- iteration
      trajectory <- trajectory[seq_len(iteration)]
      break
    }
    previous_objective <- cycle$objective
  }
  trajectory <- do.call(rbind, trajectory)
  final <- .dkge_template_cycle(
    current, source_features, centroids, sizes, support, support_masses,
    subject_weights, mapper
  )
  final_stats <- .dkge_template_feature_stats(
    current, support_masses, rank_tolerance = gate_rank_tolerance
  )
  rank_gate <- .dkge_assert_template_feature_rank(
    final_stats, rank_requirements, "final state"
  )
  failure_reason <- if (converged) NULL else "maximum_iterations_reached"
  eligibility <- .dkge_template_eligibility(
    alignment_features, converged, failure_reason,
    reference_selection = reference_selection,
    rank_gate = rank_gate
  )
  fitting <- list(
    schema_version = "1.0.0",
    algorithm = "mass_aware_entropic_barycentric_template_v1",
    objective_definition = paste0(
      "subject-weighted entropic Sinkhorn objective under a fixed mapper; ",
      "the conditional template step minimizes plan-weighted feature ",
      "distortion, is under-relaxed toward that minimizer, and is projected ",
      "to the fixed unit-row-norm constraint"
    ),
    convergence_rule = paste0(
      "after at least two iterations: relative template change <= ",
      format(tolerance), " and relative objective change <= ",
      format(objective_tolerance), "; update rate = ", format(update_rate)
    ),
    converged = converged,
    iterations = iterations,
    failure_reason = failure_reason,
    trajectory = trajectory,
    initial_feature_stats = initial_stats,
    final_feature_stats = final_stats,
    rank_gate = rank_gate,
    final_objective = final$objective,
    final_feature_distortion = final$feature_distortion,
    operators = final$operators,
    plans = final$plans,
    diagnostics = final$diagnostics,
    mapper_spec = mapper,
    update_rate = update_rate,
    initialization = initialized$provenance,
    normalization = list(
      rule = normalized$rule,
      target_row_norm = standardized$target_row_norm,
      source_standardization = standardized
    ),
    subject_weights = subject_weights,
    subject_weighting = if (all(subject_weights == 1)) {
      "equal_subject"
    } else {
      "explicit_fixed"
    },
    source_features = source_features,
    source_sizes = sizes,
    source_centroids = centroids,
    support_masses = support_masses,
    reference_selection_hash = reference_selection$structural_hash %||% NULL,
    all_subjects_refit_after_final_update = TRUE,
    no_identity_branch = TRUE
  )
  provenance <- list(
    source = "iterative_group_functional_template",
    alignment_features_hash = alignment_features$structural_hash,
    fit_binding = alignment_features$fit_binding,
    contrast_binding = alignment_features$contrast_binding,
    contrast_family_binding = alignment_features$contrast_family_binding,
    alignment_feature_source = alignment_features$feature_source,
    estimator_source = alignment_features$estimator_source,
    independent_data_hash = alignment_features$provenance$independent_data_hash %||%
      NULL,
    support_kind = support$provenance$kind %||% "identified_coordinates",
    functional_features_learned_on_support = TRUE,
    bare_coordinates_are_not_correspondence = TRUE,
    initialization = initialized$provenance
  )
  provenance$reference_selection_hash <-
    reference_selection$structural_hash %||% NULL
  provenance$reference_subject_id <-
    reference_selection$reference_subject_id %||% NULL
  .dkge_construct_functional_template(
    support = support,
    features = current,
    masses = support_masses,
    feature_source = alignment_features$feature_source,
    provenance = provenance,
    fitting = fitting,
    eligibility = eligibility,
    alignment_features_hash = alignment_features$structural_hash,
    subject_ids = subject_ids,
    reference_selection = reference_selection,
    construction = "fitted_internal"
  )
}

.dkge_new_template_alignment <- function(template, support) {
  .dkge_validate_functional_template(template, support)
  fitting <- template$fitting
  mapper_spec <- fitting$mapper_spec
  preprocessing <- list(
    source = "iterative_group_functional_template",
    feature_source = template$feature_source,
    estimator_source = template$provenance$estimator_source %||%
      "same_data_rank_truncated",
    recompute_under_null = FALSE,
    fit_binding = template$provenance$fit_binding %||% NULL,
    contrast_binding = template$provenance$contrast_binding %||% NULL,
    contrast_family_binding = template$provenance$contrast_family_binding %||%
      NULL,
    feature_provenance = template$provenance,
    template_hash = template$structural_hash
  )
  receipt <- .dkge_alignment_structural_receipt(
    mapper_spec,
    fitting$source_features,
    fitting$source_sizes,
    fitting$source_centroids,
    reference_subject = NULL,
    subject_ids = template$subject_ids,
    preprocessing = preprocessing,
    reference_support = support,
    functional_template = template,
    eligibility = template$eligibility,
    reference_selection = template$reference_selection
  )
  solution_hashes <- vapply(
    list(
      operators = fitting$operators,
      plans = fitting$plans,
      diagnostics = fitting$diagnostics
    ),
    .dkge_object_hash,
    character(1)
  )
  out <- list(
    schema_version = "1.0.0",
    operators = fitting$operators,
    plans = fitting$plans,
    diagnostics = fitting$diagnostics,
    mapper_spec = mapper_spec,
    feature_list = fitting$source_features,
    size_list = fitting$source_sizes,
    feature_ref = template$features,
    size_ref = template$masses,
    centroids = fitting$source_centroids,
    target_centroids = support$coordinates,
    medoid = NULL,
    reference_subject = NULL,
    subject_ids = template$subject_ids,
    preprocessing = preprocessing,
    reference_support = support,
    functional_template = template,
    feature_source = template$feature_source,
    estimator_source = preprocessing$estimator_source,
    eligibility = template$eligibility,
    reference_selection = template$reference_selection,
    reference_method = "iterative_group_template",
    reference_is_medoid = !is.null(template$reference_selection$medoid),
    structural_receipt = receipt,
    solution_hashes = solution_hashes
  )
  out$fitted_hash <- .dkge_object_hash(c(
    receipt$structural_hash,
    solution_hashes,
    eligibility = .dkge_object_hash(template$eligibility)
  ))
  structure(out, class = c("dkge_fitted_alignment", "list"))
}

#' Apply a fitted functional template to subject-level values
#'
#' @param template A template returned by [dkge_fit_functional_template()].
#' @param support The same identified [dkge_reference_support()] used to fit the
#'   template.
#' @param values Either a [dkge_contrast()] result or one list of subject vectors
#'   (or a named list of such lists) per contrast. A typed contrast result binds
#'   estimator, family, operator, and output provenance. Raw values can be
#'   transported for description, but their result is always inferentially
#'   ineligible. Typed input must match the exact contrast result and family
#'   bound when the template features were built. Same-data residualized
#'   templates additionally reject raw values.
#' @param contrast_ids Optional contrast identifiers.
#' @param subject_weights Optional fixed subject aggregation weights. Named
#'   weights are reordered by exact subject ID; unnamed weights are positional
#'   in fitted subject order. Equal subject weighting is the default and is
#'   distinct from DKGE fit-level MFA weights.
#' @param allow_nonconverged Logical; permit construction of explicitly
#'   ineligible aligned maps from a non-converged template.
#' @return A `dkge_aligned_maps` object. Typed contrast input produces an
#'   operator-bound object with its fitted alignment attached; raw input produces
#'   a descriptive/ineligible object without an inferential application receipt.
#' @export
dkge_align_to_template <- function(template, support, values,
                                   contrast_ids = NULL,
                                   subject_weights = NULL,
                                   allow_nonconverged = FALSE) {
  .dkge_validate_functional_template(template, support)
  fitting <- template$fitting %||% NULL
  if (is.null(fitting) || is.null(fitting$operators)) {
    .dkge_template_abort("Template correspondence operators are unavailable.")
  }
  typed_contrasts <- inherits(values, "dkge_contrasts")
  contrast_obj <- if (typed_contrasts) values else NULL
  if (identical(template$feature_source, "same_data_residualized") &&
      !typed_contrasts) {
    .dkge_template_abort(
      paste0(
        "A residualized template requires its bound `dkge_contrasts` object; ",
        "raw values cannot prove which contrast family was projected out."
      ),
      "dkge_alignment_cache_mismatch"
    )
  }
  if (typed_contrasts) {
    typed_ids <- names(values$values) %||% names(values$contrasts)
    if (!is.null(contrast_ids) &&
        !identical(as.character(contrast_ids), as.character(typed_ids))) {
      .dkge_template_abort(
        "`contrast_ids` conflict with the typed contrast result.",
        "dkge_alignment_cache_mismatch"
      )
    }
    contrast_ids <- typed_ids
    values <- values$values
  }
  if (!isTRUE(fitting$converged) && !isTRUE(allow_nonconverged)) {
    .dkge_template_abort(
      paste0(
        "Template fitting did not converge: ", fitting$failure_reason,
        ". Set `allow_nonconverged = TRUE` only to inspect ineligible maps."
      ),
      "dkge_alignment_ineligible_error"
    )
  }
  S <- length(template$subject_ids)
  is_subject_list <- is.list(values) && length(values) == S &&
    all(vapply(values, is.numeric, logical(1)))
  if (is_subject_list) values <- list(contrast1 = values)
  if (!is.list(values) || !length(values) ||
      !all(vapply(values, is.list, logical(1)))) {
    .dkge_template_abort(
      "`values` must contain one numeric source vector per subject and contrast."
    )
  }
  contrast_ids <- as.character(
    contrast_ids %||% names(values) %||% paste0("contrast", seq_along(values))
  )
  if (length(contrast_ids) != length(values) || anyNA(contrast_ids) ||
      any(!nzchar(contrast_ids)) || anyDuplicated(contrast_ids)) {
    .dkge_template_abort("`contrast_ids` must be unique and non-empty.")
  }
  names(values) <- contrast_ids
  source_values <- lapply(values, function(contrast) {
    contrast <- .dkge_order_subject_list(
      contrast, template$subject_ids, "subject values"
    )
    if (any(vapply(seq_len(S), function(s) {
      !is.numeric(contrast[[s]]) ||
        length(contrast[[s]]) != nrow(fitting$operators[[s]]) ||
        any(!is.finite(contrast[[s]]))
    }, logical(1)))) {
      .dkge_template_abort(
        "Subject values do not match the fitted source supports."
      )
    }
    names(contrast) <- template$subject_ids
    contrast
  })
  alignment <- .dkge_new_template_alignment(template, support)
  aligned <- .dkge_apply_fitted_alignment(
    source_values,
    fitted_alignment = alignment,
    contrast_obj = contrast_obj,
    contrast_ids = contrast_ids,
    subject_weights = subject_weights,
    estimand = list(
      population = "analysis subjects represented on one learned group functional template",
      correspondence = "fixed fitted template operators",
      alignment_status = if (typed_contrasts) {
        template$eligibility$status
      } else {
        "descriptive"
      },
      estimator_provenance = if (typed_contrasts) {
        "typed_dkge_contrasts"
      } else {
        "unverified_raw_subject_values"
      },
      aggregation = if (is.null(subject_weights)) {
        "equal_subject"
      } else {
        "explicit_fixed"
      },
      value_semantics = fitting$mapper_spec$params$value_type %||% "intensive"
    ),
    application_context = "dkge_align_to_template"
  )
  aligned
}
