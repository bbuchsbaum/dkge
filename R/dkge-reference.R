# dkge-reference.R
# Auditable reference-subject selection for functional alignment.

.dkge_reference_selection_payload <- function(x) {
  list(
    schema_version = x$schema_version,
    method = x$method,
    criterion = x$criterion,
    reference_subject = x$reference_subject,
    reference_subject_id = x$reference_subject_id,
    medoid = x$medoid,
    subject_ids = x$subject_ids,
    scores = x$scores,
    pairwise_loss = x$pairwise_loss,
    directed_loss = x$directed_loss,
    diagnostics = x$diagnostics,
    mapper_spec = x$mapper_spec,
    fixed_scaling = x$fixed_scaling,
    alignment_features_hash = x$alignment_features_hash,
    validation_features_hash = x$validation_features_hash,
    centroids_hash = x$centroids_hash,
    sizes_hash = x$sizes_hash,
    eligibility = x$eligibility,
    provenance = x$provenance
  )
}

.dkge_validate_reference_selection <- function(x,
                                               subject_ids = NULL,
                                               centroids = NULL,
                                               sizes = NULL,
                                               alignment_features = NULL,
                                               mapper_spec = NULL) {
  if (!inherits(x, "dkge_reference_selection")) {
    .dkge_abort("Expected a typed `dkge_reference_selection` object.",
                "dkge_reference_selection_error")
  }
  expected <- .dkge_object_hash(.dkge_reference_selection_payload(x))
  if (!identical(expected, x$structural_hash)) {
    .dkge_abort("Reference selection was mutated after fitting.",
                "dkge_reference_selection_error")
  }
  if (!is.null(subject_ids) &&
      !identical(as.character(subject_ids), x$subject_ids)) {
    .dkge_abort("Reference selection and subject order have mixed provenance.",
                "dkge_reference_selection_error")
  }
  if (!is.null(centroids)) {
    canonical <- .dkge_reference_centroids(
      centroids, subject_ids %||% x$subject_ids
    )
    canonical_sizes <- .dkge_reference_sizes(
      sizes, canonical$centroids, canonical$subject_ids
    )
    if (!identical(.dkge_object_hash(canonical$centroids), x$centroids_hash)) {
      .dkge_abort("Reference selection and centroids have mixed provenance.",
                  "dkge_reference_selection_error")
    }
    if (!identical(.dkge_object_hash(canonical_sizes), x$sizes_hash)) {
      .dkge_abort("Reference selection and masses have mixed provenance.",
                  "dkge_reference_selection_error")
    }
  }
  if (!is.null(alignment_features)) {
    .dkge_validate_alignment_features(alignment_features)
    if (!is.null(x$alignment_features_hash) &&
        !identical(alignment_features$structural_hash,
                   x$alignment_features_hash)) {
      .dkge_abort(
        "Reference selection was fitted from a different alignment-feature object.",
        "dkge_reference_selection_error"
      )
    }
  }
  if (!is.null(mapper_spec) && !is.null(x$mapper_spec)) {
    expected_mapper <- if (identical(x$method, "geometry_only")) {
      .dkge_reference_mapper(mapper_spec, "geometry_only")
    } else {
      mapper_spec
    }
    mapper_matches <- identical(
      .dkge_object_hash(expected_mapper), .dkge_object_hash(x$mapper_spec)
    )
  } else {
    mapper_matches <- TRUE
  }
  if (!mapper_matches) {
    .dkge_abort(
      "Reference selection and transport use different mapper scaling.",
      "dkge_reference_selection_error"
    )
  }
  invisible(x)
}

.dkge_reference_subject_index <- function(reference_subject, subject_ids) {
  if (is.character(reference_subject)) {
    if (length(reference_subject) != 1L || is.na(reference_subject) ||
        !nzchar(reference_subject)) {
      .dkge_abort("Character `reference_subject` must be one non-empty ID.",
                  "dkge_reference_selection_error")
    }
    index <- match(reference_subject, subject_ids)
    if (is.na(index)) {
      .dkge_abort("The requested reference subject is not in the cohort.",
                  "dkge_reference_selection_error")
    }
    return(as.integer(index))
  }
  valid_index <- is.numeric(reference_subject) && is.null(dim(reference_subject)) &&
    length(reference_subject) == 1L && !is.na(reference_subject) &&
    is.finite(reference_subject) &&
    reference_subject == floor(reference_subject) &&
    reference_subject >= 1 && reference_subject <= length(subject_ids) &&
    reference_subject <= .Machine$integer.max
  if (!valid_index) {
    .dkge_abort("`reference_subject` is not a valid subject index or ID.",
                "dkge_reference_selection_error")
  }
  as.integer(reference_subject)
}

.dkge_reference_sizes <- function(sizes, centroids, subject_ids) {
  if (is.null(sizes)) {
    sizes <- lapply(centroids, function(x) rep(1, nrow(x)))
  } else {
    sizes <- .dkge_order_subject_list(sizes, subject_ids, "sizes")
  }
  sizes <- Map(function(mass, xyz, id) {
    mass <- as.numeric(mass)
    if (length(mass) != nrow(xyz) || any(!is.finite(mass)) ||
        any(mass <= 0)) {
      .dkge_abort(
        sprintf("Reference-selection masses are invalid for subject '%s'.", id),
        "dkge_reference_selection_error"
      )
    }
    mass
  }, sizes, centroids, subject_ids)
  names(sizes) <- subject_ids
  sizes
}

.dkge_reference_centroids <- function(centroids, subject_ids = NULL) {
  if (!is.list(centroids) || length(centroids) < 2L) {
    .dkge_abort("Reference selection requires centroids for at least two subjects.",
                "dkge_reference_selection_error")
  }
  subject_ids <- as.character(subject_ids %||% names(centroids) %||%
                                paste0("subject", seq_along(centroids)))
  if (length(subject_ids) != length(centroids) || anyNA(subject_ids) ||
      any(!nzchar(subject_ids)) || anyDuplicated(subject_ids)) {
    .dkge_abort("Reference-selection subject IDs must be unique and non-empty.",
                "dkge_reference_selection_error")
  }
  centroids <- .dkge_order_subject_list(centroids, subject_ids, "centroids")
  centroids <- lapply(centroids, function(x) {
    x <- as.matrix(x)
    if (!is.numeric(x) || nrow(x) < 1L || ncol(x) < 1L ||
        any(!is.finite(x))) {
      .dkge_abort("Every centroid block must be a finite numeric matrix.",
                  "dkge_reference_selection_error")
    }
    x
  })
  names(centroids) <- subject_ids
  list(centroids = centroids, subject_ids = subject_ids)
}

.dkge_reference_global_scale <- function(values) {
  dimensions <- vapply(values, ncol, integer(1))
  if (length(unique(dimensions)) != 1L) {
    .dkge_abort("Held-out reference features must share one column dimension.",
                "dkge_reference_selection_error")
  }
  pooled <- do.call(rbind, lapply(values, as.matrix))
  center <- colMeans(pooled)
  centered <- sweep(pooled, 2L, center, "-")
  scale <- sqrt(colMeans(centered^2))
  scale[!is.finite(scale) | scale <= 1e-12] <- 1
  transformed <- lapply(values, function(x) {
    sweep(sweep(as.matrix(x), 2L, center, "-"), 2L, scale, "/")
  })
  list(values = transformed, center = center, scale = scale,
       rule = "global pooled-column RMS; fixed across every candidate/direction")
}

.dkge_reference_feature_channel <- function(alignment_features,
                                            validation_features,
                                            method) {
  if (!inherits(alignment_features, "dkge_alignment_features")) {
    .dkge_abort(
      sprintf("`%s` reference selection requires typed alignment features.",
              method),
      "dkge_reference_selection_error"
    )
  }
  .dkge_validate_alignment_features(alignment_features)
  if (identical(method, "functional_heldout")) {
    if (!inherits(validation_features, "dkge_alignment_features")) {
      .dkge_abort(
        paste0("Functional reference selection requires a second, held-out ",
               "`validation_features` object."),
        "dkge_reference_selection_error"
      )
    }
    .dkge_validate_alignment_features(validation_features)
    eligible_feature <- function(x) {
      identical(x$feature_source, "independent") &&
        isTRUE(x$source_verified) &&
        identical(x$eligibility$feature_status, "eligible")
    }
    if (!eligible_feature(alignment_features) ||
        !eligible_feature(validation_features)) {
      .dkge_abort(
        paste0(
          "Held-out functional selection requires two separately identified ",
          "independent channels with valid provenance receipts."
        ),
        "dkge_reference_selection_error"
      )
    }
    train_hash <- alignment_features$provenance$independent_data_hash
    validation_hash <- validation_features$provenance$independent_data_hash
    if (identical(train_hash, validation_hash) ||
        identical(alignment_features$provenance$independent_content_hash,
                  validation_features$provenance$independent_content_hash)) {
      .dkge_abort(
        "Reference fitting and validation channels are not independently identified.",
        "dkge_reference_selection_error"
      )
    }
    if (!identical(alignment_features$fit_binding,
                   validation_features$fit_binding) ||
        !identical(alignment_features$contrast_binding,
                   validation_features$contrast_binding) ||
        !identical(alignment_features$subject_ids,
                   validation_features$subject_ids) ||
        !identical(alignment_features$cluster_ids,
                   validation_features$cluster_ids)) {
      .dkge_abort(
        "Held-out functional channels do not describe the same fitted cohort/support.",
        "dkge_reference_selection_error"
      )
    }
    return(list(
      training = alignment_features$features,
      validation = validation_features$features,
      training_hash = alignment_features$structural_hash,
      validation_hash = validation_features$structural_hash,
      channel = "verified_independent_holdout",
      eligible = TRUE
    ))
  }
  if (!is.null(validation_features)) {
    .dkge_abort(
      "Descriptive training-feature selection does not accept a held-out channel.",
      "dkge_reference_selection_error"
    )
  }
  list(
    training = alignment_features$features,
    validation = alignment_features$features,
    training_hash = alignment_features$structural_hash,
    validation_hash = alignment_features$structural_hash,
    channel = "training_features_reused_descriptively",
    eligible = FALSE
  )
}

.dkge_reference_mapper <- function(mapper, method) {
  mapper_spec <- .dkge_resolve_mapper_spec(mapper, method = NULL, dots = list())
  if (!is.null(.dkge_epsilon_calibration(mapper_spec))) {
    .dkge_abort(
      paste0("Reference selection requires one fixed mapper scaling; point-",
             "spread calibration must be performed separately after selection."),
      "dkge_reference_selection_error"
    )
  }
  if (identical(method, "geometry_only")) {
    if (!identical(mapper_spec$strategy, "sinkhorn")) {
      .dkge_abort("Geometry-only reference selection requires the Sinkhorn mapper.",
                  "dkge_reference_selection_error")
    }
    mapper_spec$params$lambda_emb <- 0
    if ((mapper_spec$params$lambda_spa %||% 0.5) <= 0) {
      .dkge_abort("Geometry-only selection requires positive spatial cost weight.",
                  "dkge_reference_selection_error")
    }
  }
  mapper_spec
}

.dkge_reference_pairwise_losses <- function(training, validation,
                                            centroids, sizes, mapper_spec,
                                            subject_ids) {
  prepared <- .dkge_prepare_mapping_inputs(training, centroids, sizes)
  training <- prepared$features
  validation_scale <- .dkge_reference_global_scale(validation)
  validation <- validation_scale$values
  S <- length(subject_ids)
  directed <- matrix(NA_real_, S, S,
                     dimnames = list(subject_ids, subject_ids))
  spread <- matrix(NA_real_, S, S,
                   dimnames = list(subject_ids, subject_ids))
  convergence <- matrix(NA, S, S,
                        dimnames = list(subject_ids, subject_ids))
  marginal_error <- matrix(NA_real_, S, S,
                           dimnames = list(subject_ids, subject_ids))
  for (source in seq_len(S)) {
    for (target in setdiff(seq_len(S), source)) {
      mapping <- tryCatch(
        .dkge_fit_mapper_policy(
          mapper_spec,
          source_feat = training[[source]],
          target_feat = training[[target]],
          source_weights = sizes[[source]],
          target_weights = sizes[[target]],
          source_xyz = centroids[[source]],
          target_xyz = centroids[[target]],
          self_map = FALSE
        ),
        dkge_alignment_numerical_error = function(e) {
          .dkge_abort(
            sprintf("Reference-selection mapping %s -> %s is numerically invalid: %s",
                    subject_ids[[source]], subject_ids[[target]],
                    conditionMessage(e)),
            "dkge_reference_numerical_error"
          )
        }
      )
      .dkge_assert_mapping_numerically_valid(
        mapping, mapper_spec,
        context = sprintf("Reference-selection mapping %s -> %s",
                          subject_ids[[source]], subject_ids[[target]]),
        class = "dkge_reference_numerical_error"
      )
      predicted <- as.matrix(.dkge_apply_operator(
        mapping$operator, validation[[source]]
      ))
      target_values <- as.matrix(validation[[target]])
      target_energy <- mean(target_values^2)
      directed[source, target] <- mean((predicted - target_values)^2) /
        max(target_energy, .Machine$double.xmin)
      operator_diagnostics <- .dkge_operator_diagnostics(
        mapping$operator, mapping$plan %||% NULL,
        training[[source]], training[[target]], self_map = FALSE
      )
      spread[source, target] <-
        operator_diagnostics$point_spread$mean_effective_points
      diagnostics <- mapping$fit_info$diagnostics %||% list()
      convergence[source, target] <- isTRUE(diagnostics$converged %||% TRUE)
      marginal_error[source, target] <- diagnostics$marginal_error %||% NA_real_
    }
  }
  symmetric <- (directed + t(directed)) / 2
  diag(symmetric) <- NA_real_
  scores <- vapply(seq_len(S), function(target) {
    mean(symmetric[-target, target], na.rm = TRUE)
  }, numeric(1))
  names(scores) <- subject_ids
  list(
    scores = scores,
    pairwise = symmetric,
    directed = directed,
    diagnostics = list(
      mean_effective_points = spread,
      converged = convergence,
      marginal_error = marginal_error
    ),
    validation_scaling = validation_scale[c("center", "scale", "rule")]
  )
}

#' Select an auditable reference subject for functional alignment
#'
#' A reference subject supplies the target support; it is a *medoid* only when
#' selected by a stated cohort criterion. The recommended functional criterion
#' fits each directed mapper on one provenance-declared independent feature
#' channel and scores it on a second, separately identified channel. DKGE
#' validates the receipts and nonidentity of those channels; the study owner is
#' responsible for the scientific independence claim. Pair losses are symmetrized and
#' the subject with the smallest mean held-out normalized reconstruction error
#' is selected. Mapper cost weights, spatial scale, epsilon, and feature
#' normalization are fixed across every candidate.
#'
#' Geometry-only selection is ancillary and uses the same symmetric
#' reconstruction criterion on coordinates with the functional cost disabled.
#' Explicit selection is recorded as fixed by the caller. Reusing the training
#' functional features for both fit and score is available only through the
#' clearly ineligible `"descriptive_training"` method.
#'
#' @param centroids Named list of subject coordinate matrices.
#' @param alignment_features Typed training features from
#'   [dkge_alignment_features()].
#' @param validation_features A second typed independent feature object for
#'   held-out functional selection.
#' @param sizes Optional positive parcel masses.
#' @param mapper Fixed mapper specification used for every candidate/direction.
#' @param method Selection channel.
#' @param reference_subject Required subject index or ID for explicit mode.
#' @param subject_ids Optional subject IDs, primarily for geometry/explicit
#'   selection when no feature object supplies them.
#' @param provenance Optional explicit-selection provenance.
#' @param tie_tolerance Relative tolerance for deterministic score ties; ties
#'   are broken by lexicographic subject ID, never input position.
#' @return An immutable `dkge_reference_selection` object.
#' @export
dkge_select_reference_subject <- function(
    centroids,
    alignment_features = NULL,
    validation_features = NULL,
    sizes = NULL,
    mapper = dkge_mapper_spec("sinkhorn"),
    method = c(
      "functional_heldout", "geometry_only", "explicit",
      "descriptive_training"
    ),
    reference_subject = NULL,
    subject_ids = NULL,
    provenance = NULL,
    tie_tolerance = 1e-10) {
  method <- match.arg(method)
  if (!is.numeric(tie_tolerance) || length(tie_tolerance) != 1L ||
      !is.finite(tie_tolerance) || tie_tolerance < 0) {
    .dkge_abort("`tie_tolerance` must be one finite non-negative number.",
                "dkge_reference_selection_error")
  }
  if (!is.null(alignment_features)) {
    .dkge_validate_alignment_features(alignment_features)
    subject_ids <- alignment_features$subject_ids
  }
  centroid_info <- .dkge_reference_centroids(centroids, subject_ids)
  centroids <- centroid_info$centroids
  subject_ids <- centroid_info$subject_ids
  sizes <- .dkge_reference_sizes(sizes, centroids, subject_ids)
  if (!is.null(provenance) && !is.list(provenance)) {
    .dkge_abort("`provenance` must be a list or NULL.",
                "dkge_reference_selection_error")
  }

  scores <- stats::setNames(rep(NA_real_, length(subject_ids)), subject_ids)
  pairwise <- directed <- matrix(
    NA_real_, length(subject_ids), length(subject_ids),
    dimnames = list(subject_ids, subject_ids)
  )
  diagnostics <- list()
  feature_hash <- validation_hash <- NULL
  mapper_spec <- NULL
  fixed_scaling <- list()

  if (identical(method, "explicit")) {
    if (is.null(reference_subject)) {
      .dkge_abort("Explicit selection requires `reference_subject`.",
                  "dkge_reference_selection_error")
    }
    reference_index <- .dkge_reference_subject_index(
      reference_subject, subject_ids
    )
    eligible <- TRUE
    channel <- "caller_fixed_explicit_reference"
    criterion <- "explicit_fixed_reference"
    provenance <- utils::modifyList(
      list(fixed_by_caller = TRUE,
           inferential_contract = paste0(
             "Reference is conditioned on as a fixed caller input; DKGE cannot ",
             "verify choices made outside this object."
           )),
      provenance %||% list(), keep.null = TRUE
    )
  } else {
    mapper_spec <- .dkge_reference_mapper(mapper, method)
    if (identical(method, "geometry_only")) {
      if (!is.null(validation_features)) {
        .dkge_abort("Geometry-only selection does not use functional validation data.",
                    "dkge_reference_selection_error")
      }
      training <- lapply(centroids, function(x) matrix(0, nrow(x), 1L))
      sigma <- mapper_spec$params$sigma_mm %||% 15
      validation <- lapply(centroids, function(x) as.matrix(x) / sigma)
      channel <- "ancillary_geometry_only"
      eligible <- TRUE
      feature_hash <- NULL
      validation_hash <- NULL
    } else {
      channel_info <- .dkge_reference_feature_channel(
        alignment_features, validation_features, method
      )
      training <- channel_info$training
      validation <- channel_info$validation
      feature_hash <- channel_info$training_hash
      validation_hash <- channel_info$validation_hash
      channel <- channel_info$channel
      eligible <- channel_info$eligible
    }
    loss <- .dkge_reference_pairwise_losses(
      training, validation, centroids, sizes, mapper_spec, subject_ids
    )
    scores <- loss$scores
    pairwise <- loss$pairwise
    directed <- loss$directed
    diagnostics <- loss$diagnostics
    fixed_scaling <- list(
      training_features = "subject-row L2 normalization",
      validation_features = loss$validation_scaling,
      mapper = mapper_spec,
      candidate_specific_rescaling = FALSE,
      epsilon_calibration = FALSE
    )
    minimum <- min(scores)
    tolerance <- tie_tolerance * max(1, abs(minimum))
    tied <- which(abs(scores - minimum) <= tolerance)
    winner_id <- sort(subject_ids[tied], method = "radix")[[1L]]
    reference_index <- match(winner_id, subject_ids)
    criterion <- "symmetric_normalized_reconstruction_loss"
    provenance <- utils::modifyList(
      list(
        channel = channel,
        heldout = identical(method, "functional_heldout"),
        inferential_contract = if (eligible) {
          paste0(
            "Reference selection depends only on an ancillary geometry channel ",
            "or two independently identified functional channels."
          )
        } else {
          paste0(
            "Training features were reused to score reference candidates; this ",
            "selection is descriptive and not inferentially eligible."
          )
        }
      ),
      provenance %||% list(), keep.null = TRUE
    )
  }

  selection_status <- if (eligible) "eligible" else "ineligible"
  out <- list(
    schema_version = "1.0.0",
    method = method,
    criterion = criterion,
    reference_subject = as.integer(reference_index),
    reference_subject_id = subject_ids[[reference_index]],
    medoid = if (method %in% c(
      "functional_heldout", "geometry_only", "descriptive_training"
    )) as.integer(reference_index) else NULL,
    subject_ids = subject_ids,
    scores = scores,
    pairwise_loss = pairwise,
    directed_loss = directed,
    diagnostics = diagnostics,
    mapper_spec = mapper_spec,
    fixed_scaling = fixed_scaling,
    alignment_features_hash = feature_hash,
    validation_features_hash = validation_hash,
    centroids_hash = .dkge_object_hash(centroids),
    sizes_hash = .dkge_object_hash(sizes),
    eligibility = list(
      status = selection_status,
      eligible = eligible,
      channel = channel,
      reason = provenance$inferential_contract
    ),
    provenance = provenance
  )
  out$structural_hash <- .dkge_object_hash(
    .dkge_reference_selection_payload(out)
  )
  structure(out, class = c("dkge_reference_selection", "list"))
}

#' @export
print.dkge_reference_selection <- function(x, ...) {
  .dkge_validate_reference_selection(x)
  cat("<dkge_reference_selection>\n")
  cat("  subject    :", x$reference_subject_id,
      "(index", x$reference_subject, ")\n")
  cat("  method     :", x$method, "\n")
  cat("  criterion  :", x$criterion, "\n")
  cat("  eligibility:", x$eligibility$status, "\n")
  invisible(x)
}
