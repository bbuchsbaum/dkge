# dkge-alignment-api.R
# Reference-oriented public workflow and medoid lifecycle shims.

.dkge_resolve_public_reference <- function(
    fit,
    alignment_features,
    centroids,
    sizes,
    mapper_spec,
    reference_selection = NULL,
    validation_features = NULL,
    reference_subject = NULL,
    selection_method = c(
      "auto", "functional_heldout", "geometry_only", "explicit",
      "descriptive_training"
    )) {
  selection_method <- match.arg(selection_method)
  .dkge_validate_alignment_features(alignment_features, fit = fit)
  subject_ids <- alignment_features$subject_ids
  canonical <- .dkge_reference_centroids(centroids, subject_ids)
  centroids <- canonical$centroids
  sizes <- .dkge_reference_sizes(sizes, centroids, subject_ids)

  if (!is.null(reference_selection)) {
    if (!identical(selection_method, "auto") ||
        !is.null(reference_subject) || !is.null(validation_features)) {
      .dkge_abort(
        paste0(
          "A fitted `reference_selection` cannot be combined with a new ",
          "selection method, subject, or validation channel."
        ),
        "dkge_reference_selection_error"
      )
    }
    .dkge_validate_reference_selection(
      reference_selection,
      subject_ids = subject_ids,
      centroids = centroids,
      sizes = sizes,
      alignment_features = alignment_features,
      mapper_spec = mapper_spec
    )
    return(list(
      selection = reference_selection,
      centroids = centroids,
      sizes = sizes,
      selection_method = reference_selection$method
    ))
  }

  if (identical(selection_method, "auto")) {
    selection_method <- if (!is.null(reference_subject)) {
      "explicit"
    } else if (!is.null(validation_features)) {
      "functional_heldout"
    } else {
      "geometry_only"
    }
  }
  if (!identical(selection_method, "explicit") &&
      !is.null(reference_subject)) {
    .dkge_abort(
      "`reference_subject` is used only by explicit reference selection.",
      "dkge_reference_selection_error"
    )
  }
  selection <- dkge_select_reference_subject(
    centroids = centroids,
    alignment_features = if (selection_method %in% c(
      "functional_heldout", "descriptive_training"
    )) alignment_features else NULL,
    validation_features = if (identical(
      selection_method, "functional_heldout"
    )) validation_features else NULL,
    sizes = sizes,
    mapper = mapper_spec,
    method = selection_method,
    reference_subject = reference_subject,
    subject_ids = subject_ids,
    provenance = if (identical(selection_method, "explicit")) {
      list(source = "dkge_reference_oriented_api")
    } else {
      NULL
    }
  )
  list(
    selection = selection,
    centroids = centroids,
    sizes = sizes,
    selection_method = selection_method
  )
}

#' Prepare a reference-oriented fitted alignment
#'
#' This is the strict high-level replacement for constructing transport state
#' around an unnamed or default subject index. It requires a typed functional
#' feature object. With a second independent feature channel it selects a
#' functional medoid by held-out symmetric reconstruction; with no held-out
#' channel it selects a geometry-only medoid. A caller may instead supply a
#' fixed `reference_subject` or a previously fitted selection object.
#'
#' @param fit A fitted `dkge` object.
#' @param alignment_features Typed functional features from
#'   [dkge_alignment_features()]. Independent features are recommended.
#' @param centroids Named subject centroid matrices, or `NULL` to use geometry
#'   stored on `fit`.
#' @param sizes Optional positive subject parcel masses.
#' @param mapper Fixed mapper specification.
#' @param reference_selection Optional immutable result from
#'   [dkge_select_reference_subject()].
#' @param validation_features Optional second, separately identified independent
#'   channel with a validated provenance receipt, used for held-out functional
#'   reference selection.
#' @param reference_subject Optional fixed subject index or ID. This creates an
#'   explicit reference, not a medoid.
#' @param selection_method Reference-selection channel. `"auto"` chooses
#'   held-out functional selection when `validation_features` is present,
#'   explicit selection when `reference_subject` is present, and geometry-only
#'   selection otherwise. Training-feature selection is opt-in and ineligible.
#' @param ... Mapper parameters used only when `mapper` is a shorthand.
#' @return A typed `dkge_fitted_alignment` with immutable reference-selection,
#'   feature, support, template, mapper, and numerical receipts.
#' @export
dkge_prepare_alignment <- function(
    fit,
    alignment_features,
    centroids = NULL,
    sizes = NULL,
    mapper = dkge_mapper_spec("sinkhorn"),
    reference_selection = NULL,
    validation_features = NULL,
    reference_subject = NULL,
    selection_method = c(
      "auto", "functional_heldout", "geometry_only", "explicit",
      "descriptive_training"
    ),
    ...) {
  if (!inherits(fit, "dkge")) {
    .dkge_abort("`fit` must be a fitted DKGE object.", "dkge_alignment_error")
  }
  .dkge_validate_alignment_features(alignment_features, fit = fit)
  centroids <- centroids %||% fit$centroids %||% fit$input$centroids %||%
    .dkge_abort("Centroids are required for functional alignment.",
                "dkge_alignment_error")
  mapper_spec <- .dkge_resolve_mapper_spec(
    mapper, method = NULL, dots = list(...)
  )
  resolved <- .dkge_resolve_public_reference(
    fit, alignment_features, centroids, sizes, mapper_spec,
    reference_selection = reference_selection,
    validation_features = validation_features,
    reference_subject = reference_subject,
    selection_method = selection_method
  )
  dkge_prepare_transport(
    fit,
    centroids = resolved$centroids,
    sizes = resolved$sizes,
    mapper = mapper_spec,
    medoid = resolved$selection$reference_subject,
    reference_selection = resolved$selection,
    alignment_features = alignment_features
  )
}

#' Align cross-fitted contrasts to an identified reference support
#'
#' Strict reference-oriented counterpart to the deprecated
#' [dkge_transport_contrasts_to_medoid()]. Typed alignment features are
#' mandatory, preventing an inferential call from silently falling back to
#' full-fit loadings. Reference selection follows the same held-out,
#' geometry-only, or explicit rules as [dkge_prepare_alignment()].
#'
#' @param fit A fitted `dkge` object.
#' @param contrast_obj Cross-fitted contrasts from [dkge_contrast()].
#' @param alignment_features Typed functional features.
#' @param centroids Named subject centroid matrices, or `NULL` to use `fit`.
#' @param sizes Optional subject parcel masses.
#' @param mapper Fixed mapper specification.
#' @param reference_selection,validation_features,reference_subject Reference
#'   selection object, held-out validation channel, or fixed subject; see
#'   [dkge_prepare_alignment()].
#' @param selection_method Reference-selection channel; see
#'   [dkge_prepare_alignment()].
#' @param transport_cache Optional exactly matching fitted alignment.
#' @param ... Mapper parameters used only when `mapper` is shorthand.
#' @return A named transport result with attached `dkge_fitted_alignment`,
#'   `dkge_aligned_maps`, and `dkge_reference_selection` objects.
#' @export
dkge_transport_contrasts_to_reference <- function(
    fit,
    contrast_obj,
    alignment_features,
    centroids = NULL,
    sizes = NULL,
    mapper = dkge_mapper_spec("sinkhorn"),
    reference_selection = NULL,
    validation_features = NULL,
    reference_subject = NULL,
    selection_method = c(
      "auto", "functional_heldout", "geometry_only", "explicit",
      "descriptive_training"
    ),
    transport_cache = NULL,
    ...) {
  if (!inherits(fit, "dkge") || !inherits(contrast_obj, "dkge_contrasts")) {
    .dkge_abort("`fit` and `contrast_obj` must be typed DKGE objects.",
                "dkge_alignment_error")
  }
  .dkge_validate_alignment_features(
    alignment_features, fit = fit, contrast_obj = contrast_obj
  )
  centroids <- centroids %||% fit$centroids %||% fit$input$centroids %||%
    .dkge_abort("Centroids are required for functional alignment.",
                "dkge_alignment_error")
  mapper_spec <- if (!is.null(transport_cache) && missing(mapper) &&
                     !length(list(...))) {
    transport_cache$mapper_spec
  } else {
    .dkge_resolve_mapper_spec(mapper, method = NULL, dots = list(...))
  }
  if (is.null(reference_selection) && !is.null(transport_cache)) {
    reference_selection <- transport_cache$reference_selection %||% NULL
  }
  resolved <- .dkge_resolve_public_reference(
    fit, alignment_features, centroids, sizes, mapper_spec,
    reference_selection = reference_selection,
    validation_features = validation_features,
    reference_subject = reference_subject,
    selection_method = selection_method
  )
  alignment_mode <- switch(
    alignment_features$feature_source,
    independent = "independent",
    same_data_residualized = "contrast_orthogonal",
    .dkge_abort(
      "The typed feature source is unsupported by the reference workflow.",
      "dkge_alignment_feature_error"
    )
  )
  result <- .dkge_transport_contrasts_to_reference_core(
    fit = fit,
    contrast_obj = contrast_obj,
    medoid = resolved$selection$reference_subject,
    centroids = resolved$centroids,
    sizes = resolved$sizes,
    mapper = mapper_spec,
    transport_cache = transport_cache,
    reference_selection = resolved$selection,
    alignment_features = alignment_features,
    alignment_mode = alignment_mode,
    .method_missing = FALSE,
    .alignment_mode_missing = FALSE
  )
  attr(result, "selection_method") <- resolved$selection_method
  result
}

#' Paint values defined on an identified reference support
#'
#' @param values Vector or matrix of values defined on the reference labels.
#' @param labels Labels describing the identified target support.
#' @param out_file Optional output path.
#' @return A `BrainVolume` or output path, depending on `out_file`.
#' @export
dkge_paint_reference_map <- function(values, labels, out_file = NULL) {
  dkge_write_group_map(values, labels, out_file = out_file)
}

#' Transport subject contrasts to a medoid parcellation (deprecated)
#'
#' `dkge_transport_contrasts_to_medoid()` is retained for one migration cycle.
#' It preserves the legacy descriptive/fold-safe behavior where safe, while
#' current cache-provenance checks still fail closed. New analyses should use
#' [dkge_transport_contrasts_to_reference()], which requires typed features and
#' an auditable reference-selection channel.
#'
#' @inheritParams dkge_transport_contrasts_to_reference
#' @param medoid Legacy integer reference-subject index.
#' @param loadings Optional legacy loose loading matrices.
#' @param betas Optional legacy beta matrices used to derive loadings.
#' @param method Legacy mapper strategy.
#' @param alignment_mode Legacy feature-provenance mode.
#' @param reference_selection Optional typed selection object supported by the
#'   migration shim.
#' @export
dkge_transport_contrasts_to_medoid <- function(
    fit, contrast_obj, medoid, centroids = NULL,
    loadings = NULL, betas = NULL, sizes = NULL, mapper = NULL,
    method = c("sinkhorn", "ridge", "ols", "sinkhorn_cpp"),
    transport_cache = NULL, reference_selection = NULL,
    alignment_features = NULL,
    alignment_mode = c(
      "fold_safe", "independent", "contrast_orthogonal", "descriptive"
    ),
    ...) {
  method_missing <- missing(method)
  mode_missing <- missing(alignment_mode)
  .Deprecated(
    "dkge_transport_contrasts_to_reference",
    package = "dkge",
    old = "dkge_transport_contrasts_to_medoid"
  )
  .dkge_transport_contrasts_to_reference_core(
    fit = fit,
    contrast_obj = contrast_obj,
    medoid = medoid,
    centroids = centroids,
    loadings = loadings,
    betas = betas,
    sizes = sizes,
    mapper = mapper,
    method = method,
    transport_cache = transport_cache,
    reference_selection = reference_selection,
    alignment_features = alignment_features,
    alignment_mode = alignment_mode,
    .method_missing = method_missing,
    .alignment_mode_missing = mode_missing,
    ...
  )
}

#' Paint medoid values back to a label volume (deprecated)
#'
#' @inheritParams dkge_paint_reference_map
#' @export
dkge_paint_medoid_map <- function(values, labels, out_file = NULL) {
  .Deprecated(
    "dkge_paint_reference_map",
    package = "dkge",
    old = "dkge_paint_medoid_map"
  )
  dkge_paint_reference_map(values, labels, out_file = out_file)
}
