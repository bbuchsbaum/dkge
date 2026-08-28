# dkge-alignment-contract.R
# Typed public objects for support, functional correspondence, aligned rows,
# and inferential eligibility.

.dkge_alignment_feature_sources <- c(
  "independent",
  "geometry_only",
  "same_data_residualized",
  "fully_recomputed",
  "descriptive_adaptive"
)

.dkge_alignment_estimator_sources <- c(
  "independent",
  "fixed_ancillary",
  "same_data_rank_truncated",
  "fully_recomputed",
  "descriptive"
)

.dkge_reference_support_payload <- function(x) {
  list(
    schema_version = x$schema_version,
    support_id = x$support_id,
    coordinates = x$coordinates,
    labels = x$labels,
    topology = x$topology,
    decoder = x$decoder,
    provenance = x$provenance,
    n_locations = x$n_locations
  )
}

.dkge_validate_reference_support <- function(x) {
  if (!inherits(x, "dkge_reference_support")) {
    .dkge_abort("Expected a `dkge_reference_support` object.",
                "dkge_reference_support_error")
  }
  expected <- .dkge_object_hash(.dkge_reference_support_payload(x))
  if (!identical(x$structural_hash, expected)) {
    .dkge_abort("Reference support was mutated after construction.",
                "dkge_reference_support_error")
  }
  invisible(x)
}

#' Define the identified support on which aligned maps are represented
#'
#' A reference support owns only anatomical/display information: coordinates,
#' labels or topology, and an optional decoder. Functional features belong to
#' [dkge_functional_template()] and are deliberately excluded here. A bare MNI
#' lattice is therefore a valid support but is not functional correspondence.
#'
#' @param coordinates Finite numeric matrix with one row per target location.
#' @param labels Optional target labels in the same row order.
#' @param topology Optional topology or adjacency metadata.
#' @param decoder Optional pre-fitted support-to-display decoder.
#' @param support_id Optional stable identifier. A content-derived identifier is
#'   used by default.
#' @param provenance Optional list describing how the support was obtained.
#' @return An immutable `dkge_reference_support` object.
#' @export
dkge_reference_support <- function(coordinates,
                                   labels = NULL,
                                   topology = NULL,
                                   decoder = NULL,
                                   support_id = NULL,
                                   provenance = NULL) {
  coordinates <- as.matrix(coordinates)
  if (!is.numeric(coordinates) || nrow(coordinates) < 1L ||
      ncol(coordinates) < 1L || any(!is.finite(coordinates))) {
    .dkge_abort("`coordinates` must be a non-empty finite numeric matrix.",
                "dkge_reference_support_error")
  }
  if (is.null(labels)) {
    labels <- rownames(coordinates) %||% seq_len(nrow(coordinates))
  }
  if (length(labels) != nrow(coordinates) || anyNA(labels)) {
    .dkge_abort("`labels` must identify every support row without missing values.",
                "dkge_reference_support_error")
  }
  if (anyDuplicated(as.character(labels))) {
    .dkge_abort("Reference-support labels must be unique.",
                "dkge_reference_support_error")
  }
  if (!is.null(provenance) && !is.list(provenance)) {
    .dkge_abort("`provenance` must be a list or NULL.",
                "dkge_reference_support_error")
  }
  identity_payload <- list(
    coordinates = coordinates,
    labels = labels,
    topology = topology,
    decoder = decoder
  )
  support_id <- support_id %||% paste0(
    "support-", substr(.dkge_object_hash(identity_payload), 1L, 12L)
  )
  if (length(support_id) != 1L || is.na(support_id) || !nzchar(support_id)) {
    .dkge_abort("`support_id` must be one non-empty string.",
                "dkge_reference_support_error")
  }
  out <- list(
    schema_version = "1.0.0",
    support_id = as.character(support_id),
    coordinates = coordinates,
    labels = labels,
    topology = topology,
    decoder = decoder,
    provenance = provenance,
    n_locations = nrow(coordinates)
  )
  out$structural_hash <- .dkge_object_hash(.dkge_reference_support_payload(out))
  structure(out, class = c("dkge_reference_support", "list"))
}

.dkge_functional_template_payload <- function(x) {
  list(
    schema_version = x$schema_version,
    support_id = x$support_id,
    support_hash = x$support_hash,
    features = x$features,
    masses = x$masses,
    normalized_masses = x$normalized_masses,
    feature_source = x$feature_source,
    provenance = x$provenance,
    source_verified = x$source_verified,
    construction = x$construction %||% "legacy_unknown",
    fitting = x$fitting %||% NULL,
    eligibility = x$eligibility %||% NULL,
    alignment_features_hash = x$alignment_features_hash %||% NULL,
    subject_ids = x$subject_ids %||% NULL,
    reference_selection = x$reference_selection %||% NULL
  )
}

.dkge_validate_functional_template <- function(x, support = NULL) {
  if (!inherits(x, "dkge_functional_template")) {
    .dkge_abort("Expected a `dkge_functional_template` object.",
                "dkge_functional_template_error")
  }
  expected <- .dkge_object_hash(.dkge_functional_template_payload(x))
  if (!identical(x$structural_hash, expected)) {
    .dkge_abort("Functional template was mutated after construction.",
                "dkge_functional_template_error")
  }
  if (!is.null(support)) {
    .dkge_validate_reference_support(support)
    if (!identical(x$support_id, support$support_id) ||
        !identical(x$support_hash, support$structural_hash)) {
      .dkge_abort("Functional template and reference support do not match.",
                  "dkge_functional_template_error")
    }
  }
  if (!is.null(x$reference_selection)) {
    .dkge_validate_reference_selection(x$reference_selection)
  }
  construction <- x$construction %||% "legacy_unknown"
  if (!construction %in% c("public_descriptive", "fitted_internal")) {
    .dkge_abort("Functional template has no trusted construction receipt.",
                "dkge_functional_template_error")
  }
  if (identical(construction, "public_descriptive")) {
    if (!is.null(x$fitting) ||
        !inherits(x$eligibility, "dkge_alignment_eligibility") ||
        !identical(x$eligibility$status, "ineligible")) {
      .dkge_abort(
        "A public descriptive template cannot carry fitted or eligible correspondence.",
        "dkge_functional_template_error"
      )
    }
  } else if (is.null(x$fitting) ||
             !inherits(x$eligibility, "dkge_alignment_eligibility") ||
             is.null(x$alignment_features_hash) ||
             is.null(x$subject_ids)) {
    .dkge_abort("Fitted template construction receipt is incomplete.",
                "dkge_functional_template_error")
  }
  invisible(x)
}

#' Attach functional template features to a reference support
#'
#' The template owns target functional features and masses. It does not own
#' coordinates or rendering topology. `feature_source = "independent"` is the
#' recommended functional source. DKGE verifies that its provenance receipt
#' includes a non-empty `independent_data_hash`; the caller remains responsible
#' for the scientific claim that the identified acquisition/training channel is
#' statistically independent. Missing receipt evidence fails closed.
#'
#' @param support A [dkge_reference_support()] object.
#' @param features Finite target feature matrix with one row per support point.
#' @param masses Positive target masses. Defaults to equal masses.
#' @param feature_source Provenance category for correspondence features.
#' @param provenance Optional evidence, including `independent_data_hash` for an
#'   independently trained template.
#' @param fitting,eligibility Reserved internal receipts. Public callers must
#'   leave them `NULL`; fitted templates are created only by
#'   [dkge_fit_functional_template()].
#' @param alignment_features_hash Optional structural hash of the typed source
#'   feature object from which the template was learned.
#' @param subject_ids Optional identities of subjects used to learn the
#'   template.
#' @param reference_selection Optional immutable reference-selection receipt
#'   identifying how a subject-derived display support was chosen.
#' @return An immutable descriptive template/initializer. It cannot mint fitted
#'   correspondence or inferential eligibility.
#' @export
dkge_functional_template <- function(
    support,
    features,
    masses = NULL,
    feature_source = c(
      "independent", "geometry_only", "same_data_residualized",
      "fully_recomputed", "descriptive_adaptive"
    ),
    provenance = NULL,
    fitting = NULL,
    eligibility = NULL,
    alignment_features_hash = NULL,
    subject_ids = NULL,
    reference_selection = NULL) {
  if (!is.null(fitting) || !is.null(eligibility)) {
    .dkge_abort(
      paste0(
        "`fitting` and `eligibility` are internal receipts and cannot be ",
        "supplied to `dkge_functional_template()`. Use ",
        "`dkge_fit_functional_template()` to fit correspondence."
      ),
      "dkge_functional_template_error"
    )
  }
  descriptive_eligibility <- dkge_alignment_eligibility(
    feature_source = "descriptive_adaptive",
    estimator_source = "descriptive",
    source_verified = FALSE
  )
  .dkge_construct_functional_template(
    support = support,
    features = features,
    masses = masses,
    feature_source = feature_source,
    provenance = provenance,
    fitting = NULL,
    eligibility = descriptive_eligibility,
    alignment_features_hash = alignment_features_hash,
    subject_ids = subject_ids,
    reference_selection = reference_selection,
    construction = "public_descriptive"
  )
}

.dkge_construct_functional_template <- function(
    support,
    features,
    masses = NULL,
    feature_source = c(
      "independent", "geometry_only", "same_data_residualized",
      "fully_recomputed", "descriptive_adaptive"
    ),
    provenance = NULL,
    fitting = NULL,
    eligibility = NULL,
    alignment_features_hash = NULL,
    subject_ids = NULL,
    reference_selection = NULL,
    construction = c("public_descriptive", "fitted_internal")) {
  .dkge_validate_reference_support(support)
  feature_source <- match.arg(feature_source)
  construction <- match.arg(construction)
  features <- as.matrix(features)
  if (!is.numeric(features) || nrow(features) != support$n_locations ||
      ncol(features) < 1L || any(!is.finite(features))) {
    .dkge_abort(
      "`features` must be a finite numeric matrix with one row per support point.",
      "dkge_functional_template_error"
    )
  }
  if (is.null(masses)) masses <- rep(1, nrow(features))
  masses <- as.numeric(masses)
  if (length(masses) != nrow(features) || any(!is.finite(masses)) ||
      any(masses <= 0)) {
    .dkge_abort("`masses` must contain one finite positive value per support row.",
                "dkge_functional_template_error")
  }
  if (!is.null(provenance) && !is.list(provenance)) {
    .dkge_abort("`provenance` must be a list or NULL.",
                "dkge_functional_template_error")
  }
  if (!is.null(fitting) && !is.list(fitting)) {
    .dkge_abort("`fitting` must be a list or NULL.",
                "dkge_functional_template_error")
  }
  if (!is.null(eligibility) &&
      !inherits(eligibility, "dkge_alignment_eligibility")) {
    .dkge_abort("`eligibility` must be a typed alignment-eligibility record.",
                "dkge_functional_template_error")
  }
  if (!is.null(alignment_features_hash) &&
      (!is.character(alignment_features_hash) ||
       length(alignment_features_hash) != 1L ||
       is.na(alignment_features_hash) || !nzchar(alignment_features_hash))) {
    .dkge_abort("`alignment_features_hash` must be one non-empty string.",
                "dkge_functional_template_error")
  }
  if (!is.null(subject_ids)) {
    subject_ids <- as.character(subject_ids)
    if (!length(subject_ids) || anyNA(subject_ids) ||
        any(!nzchar(subject_ids)) || anyDuplicated(subject_ids)) {
      .dkge_abort("`subject_ids` must be unique and non-empty.",
                  "dkge_functional_template_error")
    }
  }
  if (!is.null(reference_selection)) {
    .dkge_validate_reference_selection(reference_selection)
  }
  out <- list(
    schema_version = "1.0.0",
    support_id = support$support_id,
    support_hash = support$structural_hash,
    features = features,
    masses = masses,
    normalized_masses = masses / sum(masses),
    feature_source = feature_source,
    provenance = provenance,
    source_verified = if (identical(feature_source, "independent")) {
      is.character(provenance$independent_data_hash %||% NULL) &&
        length(provenance$independent_data_hash) == 1L &&
        nzchar(provenance$independent_data_hash)
    } else {
      !identical(feature_source, "descriptive_adaptive")
    },
    construction = construction,
    fitting = fitting,
    eligibility = eligibility,
    alignment_features_hash = alignment_features_hash,
    subject_ids = subject_ids,
    reference_selection = reference_selection
  )
  out$structural_hash <- .dkge_object_hash(.dkge_functional_template_payload(out))
  structure(out, class = c("dkge_functional_template", "list"))
}

#' Classify inferential eligibility of a fitted alignment
#'
#' Eligibility is composite. Independent correspondence features are not an
#' exactness guarantee when contrast values still depend on a same-data,
#' rank-truncated latent span. The frozen v1 court therefore makes
#' `same_data_rank_truncated` approximate even with independent features.
#'
#' @param feature_source Correspondence-feature provenance category.
#' @param estimator_source Provenance of the contrast-generating latent span.
#' @param recompute_under_null Logical; whether every beta-dependent quantity and
#'   correspondence operator is recomputed under each valid null action.
#' @param source_verified Logical; whether the claimed feature source has the
#'   receipt evidence required by its contract. This validates provenance, not
#'   the scientific independence of an external data-generating process.
#' @param reference_selection_status Inferential status of the reference-support
#'   selection channel.
#' @return A typed `dkge_alignment_eligibility` record.
#' @export
dkge_alignment_eligibility <- function(
    feature_source = .dkge_alignment_feature_sources,
    estimator_source = .dkge_alignment_estimator_sources,
    recompute_under_null = FALSE,
    source_verified = FALSE,
    reference_selection_status = c("eligible", "approximate", "ineligible")) {
  feature_source <- match.arg(feature_source, .dkge_alignment_feature_sources)
  estimator_source <- match.arg(estimator_source,
                                .dkge_alignment_estimator_sources)
  source_verified <- isTRUE(source_verified)
  recompute_under_null <- isTRUE(recompute_under_null)
  reference_selection_status <- match.arg(reference_selection_status)

  feature_status <- switch(
    feature_source,
    independent = if (source_verified) "eligible" else "ineligible",
    geometry_only = "eligible",
    same_data_residualized = "approximate",
    fully_recomputed = if (recompute_under_null) "eligible" else "ineligible",
    descriptive_adaptive = "ineligible"
  )
  estimator_status <- switch(
    estimator_source,
    independent = "eligible",
    fixed_ancillary = "eligible",
    same_data_rank_truncated = "approximate",
    fully_recomputed = if (recompute_under_null) "eligible" else "ineligible",
    descriptive = "ineligible"
  )
  statuses <- c(feature_status, estimator_status, reference_selection_status)
  status <- if ("ineligible" %in% statuses) {
    "ineligible"
  } else if ("approximate" %in% statuses) {
    "approximate"
  } else {
    "eligible"
  }
  reasons <- character()
  if (identical(feature_source, "independent") && !source_verified) {
    reasons <- c(reasons, "independent feature provenance is unverified")
  }
  if (identical(feature_source, "descriptive_adaptive")) {
    reasons <- c(reasons, "correspondence was fitted adaptively for description")
  }
  if (identical(feature_source, "same_data_residualized")) {
    reasons <- c(reasons, "residualized same-data correspondence is approximate")
  }
  if (identical(estimator_source, "same_data_rank_truncated")) {
    reasons <- c(
      reasons,
      "contrast values depend on an estimated same-data rank-truncated latent span"
    )
  }
  if (identical(feature_source, "fully_recomputed") &&
      !recompute_under_null) {
    reasons <- c(reasons, "the fitted operator is frozen rather than recomputed under the null")
  }
  if (identical(estimator_source, "fully_recomputed") &&
      !recompute_under_null) {
    reasons <- c(reasons, "the estimator is not actually recomputed under the null")
  }
  if (identical(reference_selection_status, "approximate")) {
    reasons <- c(reasons, "reference selection is approximate")
  } else if (identical(reference_selection_status, "ineligible")) {
    reasons <- c(reasons, "reference selection reused an ineligible data channel")
  }
  if (!length(reasons) && identical(status, "eligible")) {
    reasons <- "correspondence and estimator provenance satisfy the declared null action"
  }
  structure(
    list(
      schema_version = "1.0.0",
      status = status,
      eligible = identical(status, "eligible"),
      exact = identical(status, "eligible"),
      feature_source = feature_source,
      feature_status = feature_status,
      estimator_source = estimator_source,
      estimator_status = estimator_status,
      reference_selection_status = reference_selection_status,
      recompute_under_null = recompute_under_null,
      source_verified = source_verified,
      solver_status = "not_evaluated",
      solver_converged = NA,
      solver_diagnostics = NULL,
      reason = paste(reasons, collapse = "; "),
      court_protocol = "dkge-functional-alignment-type1-v1"
    ),
    class = c("dkge_alignment_eligibility", "list")
  )
}

.dkge_invalidate_alignment_eligibility <- function(eligibility, reason,
                                                   solver_status = NULL,
                                                   solver_converged = NULL,
                                                   solver_diagnostics = NULL) {
  if (!inherits(eligibility, "dkge_alignment_eligibility")) {
    .dkge_abort("Cannot invalidate an untyped alignment eligibility record.",
                "dkge_alignment_ineligible_error")
  }
  eligibility$status <- "ineligible"
  eligibility$eligible <- FALSE
  eligibility$exact <- FALSE
  existing <- eligibility$reason %||% ""
  eligibility$reason <- paste(unique(c(
    existing[nzchar(existing)], as.character(reason)
  )), collapse = "; ")
  if (!is.null(solver_status)) eligibility$solver_status <- solver_status
  if (!is.null(solver_converged)) {
    eligibility$solver_converged <- isTRUE(solver_converged)
  }
  if (!is.null(solver_diagnostics)) {
    eligibility$solver_diagnostics <- solver_diagnostics
  }
  eligibility
}

.dkge_alignment_solver_eligibility <- function(
    eligibility,
    mapper_spec,
    diagnostics,
    operators,
    plans,
    subject_ids = NULL) {
  if (!inherits(eligibility, "dkge_alignment_eligibility")) {
    .dkge_abort("Alignment solver gating requires typed eligibility.",
                "dkge_alignment_ineligible_error")
  }
  if (!identical(mapper_spec$strategy, "sinkhorn")) {
    eligibility$solver_status <- "not_applicable"
    eligibility$solver_converged <- NA
    eligibility$solver_diagnostics <- list(required = FALSE, passed = TRUE)
    return(eligibility)
  }
  S <- length(operators)
  subject_ids <- as.character(subject_ids %||% paste0("subject", seq_len(S)))
  statuses <- lapply(seq_len(S), function(s) {
    mapping <- list(
      strategy = "sinkhorn",
      operator = operators[[s]] %||% NULL,
      plan = plans[[s]] %||% NULL,
      fit_info = list(
        tol = mapper_spec$params$tol %||% 1e-4,
        diagnostics = diagnostics[[s]] %||% list()
      )
    )
    .dkge_mapping_numerical_status(mapping, mapper_spec)
  })
  names(statuses) <- subject_ids
  passed <- vapply(statuses, `[[`, logical(1), "passed")
  summary <- list(
    required = TRUE,
    passed = all(passed),
    subject_status = statuses,
    failed_subjects = subject_ids[!passed]
  )
  if (all(passed)) {
    eligibility$solver_status <- "converged"
    eligibility$solver_converged <- TRUE
    eligibility$solver_diagnostics <- summary
    return(eligibility)
  }
  .dkge_invalidate_alignment_eligibility(
    eligibility,
    paste0(
      "Sinkhorn numerical contract failed for subject(s): ",
      paste(subject_ids[!passed], collapse = ", ")
    ),
    solver_status = "failed",
    solver_converged = FALSE,
    solver_diagnostics = summary
  )
}

.dkge_resolve_alignment_feature_source <- function(preprocessing = NULL) {
  explicit <- preprocessing$feature_source %||% NULL
  if (is.character(explicit) && length(explicit) == 1L &&
      explicit %in% .dkge_alignment_feature_sources) {
    return(explicit)
  }
  source <- preprocessing$source %||% ""
  if (grepl("independent", source, fixed = TRUE)) return("independent")
  if (grepl("geometry", source, fixed = TRUE)) return("geometry_only")
  if (grepl("residual", source, fixed = TRUE)) return("same_data_residualized")
  if (grepl("fully_recomputed", source, fixed = TRUE)) return("fully_recomputed")
  "descriptive_adaptive"
}

.dkge_alignment_source_verified <- function(feature_source, preprocessing) {
  if (identical(feature_source, "independent")) {
    hash <- preprocessing$feature_provenance$independent_data_hash %||%
      preprocessing$independent_data_hash %||% NULL
    return(is.character(hash) && length(hash) == 1L && nzchar(hash))
  }
  !identical(feature_source, "descriptive_adaptive")
}

.dkge_alignment_objects_from_inputs <- function(feature_list, size_list,
                                                centroids,
                                                reference_subject,
                                                preprocessing = NULL,
                                                reference_selection = NULL) {
  if (!is.null(reference_selection)) {
    .dkge_validate_reference_selection(reference_selection)
  }
  feature_source <- .dkge_resolve_alignment_feature_source(preprocessing)
  support <- dkge_reference_support(
    centroids[[reference_subject]],
    labels = rownames(feature_list[[reference_subject]]) %||%
      seq_len(nrow(feature_list[[reference_subject]])),
    provenance = list(
      kind = "reference_subject_support",
      reference_subject = as.integer(reference_subject),
      reference_subject_id = reference_selection$reference_subject_id %||%
        NULL,
      reference_selection_hash = reference_selection$structural_hash %||% NULL
    )
  )
  feature_provenance <- preprocessing$feature_provenance %||% list(
    source = preprocessing$source %||% feature_source
  )
  template <- dkge_functional_template(
    support,
    feature_list[[reference_subject]],
    masses = size_list[[reference_subject]],
    feature_source = feature_source,
    provenance = feature_provenance,
    reference_selection = reference_selection
  )
  estimator_source <- preprocessing$estimator_source %||%
    "same_data_rank_truncated"
  recompute <- isTRUE(preprocessing$recompute_under_null)
  eligibility <- dkge_alignment_eligibility(
    feature_source,
    estimator_source,
    recompute_under_null = recompute,
    source_verified = .dkge_alignment_source_verified(
      feature_source, preprocessing %||% list()
    ),
    reference_selection_status = reference_selection$eligibility$status %||%
      "eligible"
  )
  list(
    support = support,
    template = template,
    feature_source = feature_source,
    estimator_source = estimator_source,
    eligibility = eligibility,
    reference_selection = reference_selection
  )
}

.dkge_validate_fitted_alignment_object <- function(x) {
  if (!inherits(x, "dkge_fitted_alignment")) {
    .dkge_abort("Expected a typed `dkge_fitted_alignment` object.",
                "dkge_alignment_cache_mismatch")
  }
  .dkge_validate_reference_support(x$reference_support)
  .dkge_validate_functional_template(x$functional_template,
                                     x$reference_support)
  if (!is.null(x$reference_selection)) {
    .dkge_validate_reference_selection(x$reference_selection)
  }
  if (!is.null(x$fit_binding) || !is.null(x$contrast_binding) ||
      !is.null(x$contrast_family_binding)) {
    .dkge_abort(
      paste0(
        "Fitted-alignment bindings must live only in the immutable ",
        "preprocessing receipt."
      ),
      "dkge_alignment_cache_mismatch"
    )
  }
  if (x$feature_source %in% c("independent", "same_data_residualized") &&
      (is.null(x$preprocessing$fit_binding) ||
       is.null(x$preprocessing$contrast_binding) ||
       is.null(x$preprocessing$contrast_family_binding))) {
    .dkge_abort("Typed fitted alignment has incomplete fit/contrast bindings.",
                "dkge_alignment_cache_mismatch")
  }
  structural_now <- .dkge_alignment_structural_receipt(
    x$mapper_spec,
    x$feature_list,
    x$size_list,
    x$centroids,
    x$reference_subject %||% x$medoid,
    subject_ids = x$subject_ids,
    preprocessing = x$preprocessing,
    reference_support = x$reference_support,
    functional_template = x$functional_template,
    eligibility = x$eligibility,
    reference_selection = x$reference_selection
  )
  if (!identical(x$structural_receipt$structural_hash,
                 structural_now$structural_hash)) {
    bad <- names(structural_now$component_hashes)[
      x$structural_receipt$component_hashes != structural_now$component_hashes
    ]
    .dkge_abort(
      sprintf("Fitted alignment state was mutated after fitting: %s.",
              paste(bad, collapse = ", ")),
      "dkge_alignment_cache_mismatch"
    )
  }
  solution_hashes <- vapply(
    list(operators = x$operators, plans = x$plans,
         diagnostics = x$diagnostics),
    .dkge_object_hash,
    character(1)
  )
  if (!identical(x$solution_hashes, solution_hashes)) {
    bad <- names(solution_hashes)[x$solution_hashes != solution_hashes]
    .dkge_abort(
      sprintf("Fitted alignment solution was mutated after fitting: %s.",
              paste(bad, collapse = ", ")),
      "dkge_alignment_cache_mismatch"
    )
  }
  expected <- .dkge_object_hash(c(
    x$structural_receipt$structural_hash,
    solution_hashes,
    eligibility = .dkge_object_hash(x$eligibility)
  ))
  if (!identical(x$fitted_hash, expected)) {
    .dkge_abort("Fitted-alignment provenance hash is invalid.",
                "dkge_alignment_cache_mismatch")
  }
  invisible(x)
}

.dkge_aligned_maps_payload <- function(x) {
  list(
    schema_version = x$schema_version,
    construction = x$construction,
    values = x$values,
    subject_ids = x$subject_ids,
    contrast_ids = x$contrast_ids,
    support_id = x$support_id,
    support_hash = x$support_hash,
    fitted_alignment_hash = x$fitted_alignment_hash,
    operators_hash = x$operators_hash,
    application_receipt = x$application_receipt,
    subject_weights = x$subject_weights,
    subject_weighting = x$subject_weighting,
    estimand = x$estimand,
    eligibility = x$eligibility,
    feature_source = x$feature_source
  )
}

.dkge_resolve_subject_weights <- function(subject_weights, subject_ids,
                                          error_class) {
  subject_ids <- as.character(subject_ids)
  S <- length(subject_ids)
  if (is.null(subject_weights)) {
    weights <- stats::setNames(rep(1, S), subject_ids)
    return(list(weights = weights, weighting = "equal_subject"))
  }

  weight_names <- names(subject_weights)
  if (!is.null(weight_names)) {
    weight_names <- as.character(weight_names)
    if (length(weight_names) != S || anyNA(weight_names) ||
        any(!nzchar(weight_names)) || anyDuplicated(weight_names) ||
        !setequal(weight_names, subject_ids)) {
      .dkge_abort(
        paste0(
          "Named `subject_weights` must identify every fitted subject exactly ",
          "once; unnamed weights are interpreted in fitted subject order."
        ),
        error_class
      )
    }
  }

  weights <- as.numeric(subject_weights)
  if (length(weights) != S || any(!is.finite(weights)) ||
      any(weights <= 0)) {
    .dkge_abort(
      "`subject_weights` must be finite, positive, and aligned to subjects.",
      error_class
    )
  }
  if (!is.null(weight_names)) {
    weights <- weights[match(subject_ids, weight_names)]
  }
  weights <- weights / mean(weights)
  names(weights) <- subject_ids
  list(weights = weights, weighting = "explicit_fixed")
}

.dkge_assert_fitted_alignment_contrast_binding <- function(
    fitted_alignment, contrast_obj) {
  if (!inherits(contrast_obj, "dkge_contrasts")) {
    .dkge_abort(
      "A typed fitted correspondence requires a `dkge_contrasts` result.",
      "dkge_alignment_cache_mismatch"
    )
  }
  preprocessing <- fitted_alignment$preprocessing %||% list()
  expected_contrast <- preprocessing$contrast_binding %||% NULL
  expected_family <- preprocessing$contrast_family_binding %||% NULL
  observed_contrast <- .dkge_alignment_contrast_binding(contrast_obj)
  observed_family <- .dkge_alignment_contrast_family_binding(contrast_obj)
  if ((!is.null(expected_contrast) &&
       !identical(expected_contrast, observed_contrast)) ||
      (!is.null(expected_family) &&
       !identical(expected_family, observed_family))) {
    .dkge_abort(
      paste0(
        "Fitted correspondence was constructed for a different typed ",
        "contrast result or family."
      ),
      "dkge_alignment_cache_mismatch"
    )
  }
  invisible(fitted_alignment)
}

.dkge_validate_aligned_maps <- function(x) {
  if (!inherits(x, "dkge_aligned_maps")) {
    .dkge_abort("Expected a `dkge_aligned_maps` object.",
                "dkge_aligned_maps_error")
  }
  expected <- .dkge_object_hash(.dkge_aligned_maps_payload(x))
  if (!identical(x$structural_hash, expected)) {
    .dkge_abort("Aligned subject maps were mutated after construction.",
                "dkge_aligned_maps_error")
  }
  if (!x$construction %in% c("public_descriptive", "operator_applied_internal")) {
    .dkge_abort("Aligned maps have an unknown construction route.",
                "dkge_aligned_maps_error")
  }
  if (identical(x$construction, "public_descriptive")) {
    if (!is.null(x$application_receipt) ||
        !identical(x$eligibility$status, "ineligible") ||
        isTRUE(x$eligibility$eligible)) {
      .dkge_abort(
        "Public raw aligned maps cannot carry operator-application eligibility.",
        "dkge_aligned_maps_error"
      )
    }
    return(invisible(x))
  }

  receipt <- x$application_receipt
  required <- c(
    "schema_version", "application_context", "source_value_kind",
    "source_values_hash", "contrast_family_binding", "contrast_binding",
    "fitted_alignment_hash", "operators_hash", "output_values_hash",
    "receipt_hash"
  )
  if (!is.list(receipt) || !all(required %in% names(receipt))) {
    .dkge_abort("Operator-applied maps have an incomplete application receipt.",
                "dkge_aligned_maps_error")
  }
  receipt_payload <- receipt
  receipt_payload$receipt_hash <- NULL
  if (!identical(receipt$receipt_hash, .dkge_object_hash(receipt_payload)) ||
      !identical(receipt$output_values_hash, .dkge_object_hash(x$values)) ||
      !identical(receipt$fitted_alignment_hash, x$fitted_alignment_hash) ||
      !identical(receipt$operators_hash, x$operators_hash) ||
      !is.character(receipt$source_values_hash) ||
      length(receipt$source_values_hash) != 1L ||
      is.na(receipt$source_values_hash) || !nzchar(receipt$source_values_hash) ||
      !is.character(receipt$contrast_family_binding) ||
      length(receipt$contrast_family_binding) != 1L ||
      is.na(receipt$contrast_family_binding) ||
      !nzchar(receipt$contrast_family_binding)) {
    .dkge_abort("Aligned-map operator-application provenance is invalid.",
                "dkge_aligned_maps_error")
  }
  fitted_alignment <- attr(x, "fitted_alignment", exact = TRUE)
  if (!inherits(fitted_alignment, "dkge_fitted_alignment")) {
    .dkge_abort("Operator-applied maps lost their fitted alignment receipt.",
                "dkge_aligned_maps_error")
  }
  .dkge_validate_fitted_alignment_object(fitted_alignment)
  if (!identical(fitted_alignment$fitted_hash, x$fitted_alignment_hash) ||
      !identical(fitted_alignment$reference_support$support_id,
                 x$support_id) ||
      !identical(fitted_alignment$reference_support$structural_hash,
                 x$support_hash) ||
      !identical(fitted_alignment$solution_hashes[["operators"]],
                 x$operators_hash) ||
      !identical(fitted_alignment$eligibility, x$eligibility)) {
    .dkge_abort("Aligned maps do not match their fitted alignment receipt.",
                "dkge_aligned_maps_error")
  }
  expected_contrast <-
    fitted_alignment$preprocessing$contrast_binding %||% NULL
  expected_family <-
    fitted_alignment$preprocessing$contrast_family_binding %||% NULL
  if ((!is.null(expected_contrast) &&
       !identical(receipt$contrast_binding, expected_contrast)) ||
      (!is.null(expected_family) &&
       !identical(receipt$contrast_family_binding, expected_family))) {
    .dkge_abort(
      paste0(
        "Aligned-map contrast bindings do not match their fitted ",
        "correspondence receipt."
      ),
      "dkge_aligned_maps_error"
    )
  }
  invisible(x)
}

.dkge_construct_aligned_maps <- function(
    values,
    fitted_alignment,
    subject_ids = NULL,
    contrast_ids = NULL,
    subject_weights = NULL,
    estimand = NULL,
    construction = c("public_descriptive", "operator_applied_internal"),
    source_values = NULL,
    contrast_binding = NULL,
    contrast_family_binding = NULL,
    source_value_kind = NULL,
    application_context = NULL) {
  construction <- match.arg(construction)
  .dkge_validate_fitted_alignment_object(fitted_alignment)
  if (is.matrix(values)) values <- list(contrast1 = values)
  if (!is.list(values) || !length(values)) {
    .dkge_abort("`values` must be a non-empty list of subject-by-support matrices.",
                "dkge_aligned_maps_error")
  }
  values <- lapply(values, as.matrix)
  S <- nrow(values[[1]])
  Q <- fitted_alignment$reference_support$n_locations
  if (S < 1L || any(vapply(values, nrow, integer(1)) != S) ||
      any(vapply(values, ncol, integer(1)) != Q) ||
      any(!vapply(values, function(x) is.numeric(x) && all(is.finite(x)),
                  logical(1)))) {
    .dkge_abort(
      "Every aligned map must be a finite S-by-Q numeric matrix on the identified support.",
      "dkge_aligned_maps_error"
    )
  }
  row_id_list <- lapply(values, rownames)
  has_row_ids <- !vapply(row_id_list, is.null, logical(1))
  if (any(has_row_ids) && !all(has_row_ids)) {
    .dkge_abort("Aligned matrices mix named and unnamed subject rows.",
                "dkge_aligned_maps_error")
  }
  value_ids <- NULL
  if (all(has_row_ids)) {
    value_ids <- as.character(row_id_list[[1L]])
    valid_rows <- function(ids) {
      length(ids) == S && !anyNA(ids) && all(nzchar(ids)) &&
        !anyDuplicated(ids) && setequal(ids, value_ids)
    }
    if (!all(vapply(row_id_list, valid_rows, logical(1)))) {
      .dkge_abort("Aligned matrices do not identify the same subject cohort.",
                  "dkge_aligned_maps_error")
    }
    values <- Map(function(Y, ids) {
      Y[match(value_ids, ids), , drop = FALSE]
    }, values, row_id_list)
  }
  subject_ids <- as.character(subject_ids %||% value_ids %||%
                                fitted_alignment$subject_ids %||%
                                paste0("subject", seq_len(S)))
  if (length(subject_ids) != S || anyNA(subject_ids) ||
      any(!nzchar(subject_ids)) || anyDuplicated(subject_ids)) {
    .dkge_abort("`subject_ids` must uniquely identify every aligned row.",
                "dkge_aligned_maps_error")
  }
  if (!is.null(value_ids) && !setequal(subject_ids, value_ids)) {
    .dkge_abort("`subject_ids` conflict with aligned matrix row names.",
                "dkge_aligned_maps_error")
  }
  fitted_subject_ids <- fitted_alignment$subject_ids %||% NULL
  if (!is.null(fitted_subject_ids)) {
    fitted_subject_ids <- as.character(fitted_subject_ids)
    if (!setequal(subject_ids, fitted_subject_ids)) {
      .dkge_abort(
        "Aligned maps and fitted correspondence identify different subjects.",
        "dkge_aligned_maps_error"
      )
    }
  }
  source_subject_ids <- value_ids %||% subject_ids
  final_subject_ids <- fitted_subject_ids %||% subject_ids
  order_index <- match(final_subject_ids, source_subject_ids)
  values <- lapply(values, function(Y) Y[order_index, , drop = FALSE])
  subject_ids <- final_subject_ids
  contrast_ids <- as.character(contrast_ids %||% names(values) %||%
                                 paste0("contrast", seq_along(values)))
  if (length(contrast_ids) != length(values) || anyNA(contrast_ids) ||
      anyDuplicated(contrast_ids)) {
    .dkge_abort("`contrast_ids` must uniquely identify every aligned matrix.",
                "dkge_aligned_maps_error")
  }
  names(values) <- contrast_ids
  values <- lapply(values, function(x) {
    dimnames(x) <- list(subject_ids,
                        as.character(fitted_alignment$reference_support$labels))
    x
  })
  resolved_weights <- .dkge_resolve_subject_weights(
    subject_weights, subject_ids, "dkge_aligned_maps_error"
  )
  subject_weights <- resolved_weights$weights
  subject_weighting <- resolved_weights$weighting
  estimand <- estimand %||% list(
    population = "analysis subjects represented on one identified reference support",
    location = "subject-level aligned contrast value",
    aggregation = subject_weighting,
    value_semantics = fitted_alignment$mapper_spec$params$value_type %||%
      "intensive"
  )
  application_receipt <- NULL
  eligibility <- fitted_alignment$eligibility
  feature_source <- fitted_alignment$feature_source
  operators_hash <- fitted_alignment$solution_hashes[["operators"]]
  if (identical(construction, "public_descriptive")) {
    if (!is.null(source_values) || !is.null(contrast_binding) ||
        !is.null(contrast_family_binding) || !is.null(source_value_kind) ||
        !is.null(application_context)) {
      .dkge_abort("Public descriptive maps cannot receive internal application receipts.",
                  "dkge_aligned_maps_error")
    }
    eligibility <- dkge_alignment_eligibility(
      feature_source = "descriptive_adaptive",
      estimator_source = "descriptive",
      source_verified = FALSE,
      reference_selection_status = fitted_alignment$reference_selection$eligibility$status %||%
        "ineligible"
    )
    feature_source <- "descriptive_adaptive"
    operators_hash <- NULL
  } else {
    if (is.null(source_values) || !is.character(source_value_kind) ||
        length(source_value_kind) != 1L || is.na(source_value_kind) ||
        !nzchar(source_value_kind) || !is.character(application_context) ||
        length(application_context) != 1L || is.na(application_context) ||
        !nzchar(application_context) ||
        !is.character(contrast_family_binding) ||
        length(contrast_family_binding) != 1L ||
        is.na(contrast_family_binding) || !nzchar(contrast_family_binding)) {
      .dkge_abort("Internal aligned maps require complete source/application bindings.",
                  "dkge_aligned_maps_error")
    }
    application_receipt <- list(
      schema_version = "1.0.0",
      application_context = application_context,
      source_value_kind = source_value_kind,
      source_values_hash = .dkge_object_hash(source_values),
      contrast_family_binding = contrast_family_binding,
      contrast_binding = contrast_binding,
      fitted_alignment_hash = fitted_alignment$fitted_hash,
      operators_hash = operators_hash,
      output_values_hash = .dkge_object_hash(values)
    )
    application_receipt$receipt_hash <-
      .dkge_object_hash(application_receipt)
  }
  out <- list(
    schema_version = "1.0.0",
    construction = construction,
    values = values,
    subject_ids = subject_ids,
    contrast_ids = contrast_ids,
    support_id = fitted_alignment$reference_support$support_id,
    support_hash = fitted_alignment$reference_support$structural_hash,
    fitted_alignment_hash = fitted_alignment$fitted_hash,
    operators_hash = operators_hash,
    application_receipt = application_receipt,
    subject_weights = subject_weights,
    subject_weighting = subject_weighting,
    estimand = estimand,
    eligibility = eligibility,
    feature_source = feature_source
  )
  out$structural_hash <- .dkge_object_hash(.dkge_aligned_maps_payload(out))
  out <- structure(out, class = c("dkge_aligned_maps", "list"))
  if (identical(construction, "operator_applied_internal")) {
    attr(out, "fitted_alignment") <- fitted_alignment
  }
  out
}

.dkge_apply_fitted_alignment <- function(
    source_values,
    fitted_alignment,
    contrast_obj = NULL,
    contrast_ids = NULL,
    subject_weights = NULL,
    estimand = NULL,
    application_context = "fixed_fitted_alignment") {
  .dkge_validate_fitted_alignment_object(fitted_alignment)
  subject_ids <- as.character(fitted_alignment$subject_ids)
  S <- length(subject_ids)
  typed_contrasts <- inherits(contrast_obj, "dkge_contrasts")
  if (!is.null(contrast_obj) && !typed_contrasts) {
    .dkge_abort("`contrast_obj` must be a typed `dkge_contrasts` object.",
                "dkge_aligned_maps_error")
  }
  if (typed_contrasts) {
    .dkge_assert_fitted_alignment_contrast_binding(
      fitted_alignment, contrast_obj
    )
    typed_ids <- names(contrast_obj$values) %||% names(contrast_obj$contrasts)
    if (!is.null(contrast_ids) &&
        !identical(as.character(contrast_ids), as.character(typed_ids))) {
      .dkge_abort("`contrast_ids` conflict with the typed contrast result.",
                  "dkge_alignment_cache_mismatch")
    }
    contrast_ids <- typed_ids
  }
  is_subject_list <- is.list(source_values) && length(source_values) == S &&
    all(vapply(source_values, is.numeric, logical(1)))
  if (is_subject_list) source_values <- list(contrast1 = source_values)
  if (!is.list(source_values) || !length(source_values) ||
      !all(vapply(source_values, is.list, logical(1)))) {
    .dkge_abort(
      "Source values must contain one numeric vector per subject and contrast.",
      "dkge_aligned_maps_error"
    )
  }
  contrast_ids <- as.character(
    contrast_ids %||% names(source_values) %||%
      paste0("contrast", seq_along(source_values))
  )
  if (length(contrast_ids) != length(source_values) || anyNA(contrast_ids) ||
      any(!nzchar(contrast_ids)) || anyDuplicated(contrast_ids)) {
    .dkge_abort("`contrast_ids` must be unique and non-empty.",
                "dkge_aligned_maps_error")
  }
  names(source_values) <- contrast_ids
  source_values <- lapply(source_values, function(contrast) {
    contrast <- .dkge_order_subject_list(
      contrast, subject_ids, "subject source values"
    )
    bad <- vapply(seq_len(S), function(s) {
      !is.numeric(contrast[[s]]) ||
        length(contrast[[s]]) != nrow(fitted_alignment$operators[[s]]) ||
        any(!is.finite(contrast[[s]]))
    }, logical(1))
    if (any(bad)) {
      .dkge_abort("Subject values do not match the fitted source supports.",
                  "dkge_aligned_maps_error")
    }
    names(contrast) <- subject_ids
    contrast
  })
  if (typed_contrasts) {
    bound_values <- contrast_obj$values
    bound_is_subject_list <- is.list(bound_values) &&
      length(bound_values) == S &&
      all(vapply(bound_values, is.numeric, logical(1)))
    if (bound_is_subject_list) bound_values <- list(contrast1 = bound_values)
    names(bound_values) <- contrast_ids
    bound_values <- lapply(bound_values, function(contrast) {
      contrast <- .dkge_order_subject_list(
        contrast, subject_ids, "bound contrast values"
      )
      names(contrast) <- subject_ids
      contrast
    })
    if (!identical(.dkge_object_hash(source_values),
                   .dkge_object_hash(bound_values))) {
      .dkge_abort("Source values do not match the bound contrast result.",
                  "dkge_alignment_cache_mismatch")
    }
  }
  if (identical(fitted_alignment$feature_source, "same_data_residualized") &&
      !typed_contrasts) {
    .dkge_abort(
      "Residualized correspondence requires its exact bound contrast result and family.",
      "dkge_alignment_cache_mismatch"
    )
  }
  mapped <- lapply(source_values, function(contrast) {
    do.call(rbind, lapply(seq_len(S), function(s) {
      as.numeric(.dkge_apply_operator(
        fitted_alignment$operators[[s]], contrast[[s]]
      ))
    }))
  })
  if (!typed_contrasts) {
    return(.dkge_construct_aligned_maps(
      mapped,
      fitted_alignment = fitted_alignment,
      subject_ids = subject_ids,
      contrast_ids = contrast_ids,
      subject_weights = subject_weights,
      estimand = estimand,
      construction = "public_descriptive"
    ))
  }
  family_binding <- .dkge_alignment_contrast_family_binding(contrast_obj)
  .dkge_construct_aligned_maps(
    mapped,
    fitted_alignment = fitted_alignment,
    subject_ids = subject_ids,
    contrast_ids = contrast_ids,
    subject_weights = subject_weights,
    estimand = estimand,
    construction = "operator_applied_internal",
    source_values = source_values,
    contrast_binding = .dkge_alignment_contrast_binding(contrast_obj),
    contrast_family_binding = family_binding,
    source_value_kind = "typed_dkge_contrasts",
    application_context = application_context
  )
}

#' Store descriptive subject-by-support maps
#'
#' This constructor accepts values that are already on a reference support, so
#' it cannot prove that a fitted correspondence produced them. Its result is
#' always descriptive and inferentially ineligible, even when
#' `fitted_alignment` itself is eligible. Use
#' [dkge_transport_contrasts_to_reference()] or [dkge_align_to_template()] to
#' apply fitted operators and obtain operator-bound aligned maps.
#'
#' @param values Named list of matrices, one per contrast. Every matrix must be
#'   finite with one row per subject and one column per reference location.
#' @param fitted_alignment A typed `dkge_fitted_alignment`; used to identify the
#'   support, not to certify that it produced `values`.
#' @param subject_ids Optional subject identifiers.
#' @param contrast_ids Optional contrast identifiers.
#' @param subject_weights Optional fixed group weights. Named weights are
#'   reordered by exact subject ID; unnamed weights are positional in fitted
#'   subject order. Equal subject weighting is the default; DKGE fit-level MFA
#'   weights are never inherited silently.
#' @param estimand Optional explicit estimand description.
#' @return An immutable, descriptive `dkge_aligned_maps` object. It can be
#'   rendered but is refused by group-inference functions.
#' @export
dkge_aligned_maps <- function(values,
                              fitted_alignment,
                              subject_ids = NULL,
                              contrast_ids = NULL,
                              subject_weights = NULL,
                              estimand = NULL) {
  .dkge_construct_aligned_maps(
    values,
    fitted_alignment = fitted_alignment,
    subject_ids = subject_ids,
    contrast_ids = contrast_ids,
    subject_weights = subject_weights,
    estimand = estimand,
    construction = "public_descriptive"
  )
}

.dkge_assert_alignment_eligible <- function(x, allow_approximate = FALSE) {
  eligibility <- if (inherits(x, "dkge_aligned_maps")) {
    x$eligibility
  } else if (inherits(x, "dkge_fitted_alignment")) {
    x$eligibility
  } else {
    x
  }
  if (!inherits(eligibility, "dkge_alignment_eligibility")) {
    .dkge_abort("Alignment has no typed inferential-eligibility record.",
                "dkge_alignment_ineligible_error")
  }
  if (isTRUE(eligibility$eligible)) return(invisible(eligibility))
  if (isTRUE(allow_approximate) && identical(eligibility$status, "approximate")) {
    return(invisible(eligibility))
  }
  .dkge_abort(
    paste0(
      "Group inference refuses alignment status '", eligibility$status,
      "' by default: ", eligibility$reason,
      if (identical(eligibility$status, "approximate")) {
        ". Set `allow_approximate_alignment = TRUE` only for an explicitly labelled approximate analysis"
      } else {
        ""
      }, "."
    ),
    "dkge_alignment_ineligible_error"
  )
}

#' @export
print.dkge_reference_support <- function(x, ...) {
  cat("<dkge_reference_support>\n")
  cat("  id        :", x$support_id, "\n")
  cat("  locations :", x$n_locations, "\n")
  cat("  functional: none (support only)\n")
  invisible(x)
}

#' @export
print.dkge_functional_template <- function(x, ...) {
  .dkge_validate_functional_template(x)
  cat("<dkge_functional_template>\n")
  cat("  support :", x$support_id, "\n")
  cat("  features:", ncol(x$features), "\n")
  cat("  source  :", x$feature_source, "\n")
  cat("  provenance receipt:", x$source_verified, "\n")
  if (!is.null(x$fitting)) {
    cat("  fit     :", if (isTRUE(x$fitting$converged)) "converged" else
      paste0("failed (", x$fitting$failure_reason, ")"), "\n")
    cat("  iterations:", x$fitting$iterations, "\n")
    cat("  eligibility:", x$eligibility$status, "\n")
  }
  invisible(x)
}

#' @export
print.dkge_fitted_alignment <- function(x, ...) {
  .dkge_validate_fitted_alignment_object(x)
  cat("<dkge_fitted_alignment>\n")
  cat("  support       :", x$reference_support$support_id, "\n")
  reference_label <- x$reference_selection$reference_subject_id %||%
    if (!is.null(x$reference_subject %||% x$medoid)) {
      x$subject_ids[[x$reference_subject %||% x$medoid]]
    } else {
      "group functional template"
    }
  cat("  reference     :", reference_label, "\n")
  cat("  selection     :", x$reference_method %||% "explicit_legacy", "\n")
  cat("  subjects      :", length(x$subject_ids), "\n")
  cat("  feature source:", x$feature_source, "\n")
  cat("  eligibility   :", x$eligibility$status, "\n")
  invisible(x)
}

#' @export
print.dkge_aligned_maps <- function(x, ...) {
  cat("<dkge_aligned_maps>\n")
  cat("  support   :", x$support_id, "\n")
  cat("  subjects  :", length(x$subject_ids), "\n")
  cat("  contrasts :", length(x$contrast_ids), "\n")
  cat("  weighting :", x$subject_weighting, "\n")
  cat("  eligibility:", x$eligibility$status, "\n")
  invisible(x)
}
