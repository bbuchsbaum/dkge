# dkge-components.R
# Convenience helpers for component-level inference and transport

#' Deprecated descriptive component consensus
#'
#' This legacy helper transports full-fit component loadings with correspondence
#' learned from those same loadings. It is retained only for descriptive
#' summaries and cannot compute p-values, confidence claims, or significance.
#' Use [dkge_prepare_alignment()], apply correspondence with
#' [dkge_transport_contrasts_to_reference()] or [dkge_align_to_template()], and
#' then call [dkge_infer_aligned()] for the typed inferential workflow.
#'
#' @param fit A fitted `dkge` object.
#' @param mapper Mapper strategy (string or [dkge_mapper_spec()]). Defaults to
#'   "sinkhorn".
#' @param centroids Optional list of subject centroid matrices; defaults to
#'   centroids stored in `fit` if available.
#' @param sizes Optional list of cluster masses (one vector per subject).
#' @param inference Deprecated. Must be `NULL`; legacy full-fit correspondence
#'   is descriptive/ineligible and cannot enter an inference boundary.
#' @param medoid Reference subject index for the descriptive display.
#' @param components Optional vector of component indices or names; default is
#'   all components.
#' @param adjust Deprecated and ignored because inferential p-values are no
#'   longer produced by this helper.
#' @param ... Additional mapper-specific parameters (e.g. `epsilon`).
#'
#' @return A list with fields:
#'   - `summary`: tidy data frame of descriptive means and standard deviations.
#'   - `statistics`: per-component mean vectors.
#'   - `transport`: per-component transported subject matrices.
#'   - `eligibility`: an ineligible/descriptive alignment receipt.
#' @examples
#' \donttest{
#' toy <- dkge_sim_toy(
#'   factors = list(A = list(L = 2), B = list(L = 3)),
#'   active_terms = c("A", "B"), S = 3, P = 15, snr = 5
#' )
#' fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 2)
#' centroids <- lapply(toy$B_list, function(B) matrix(rnorm(ncol(B) * 3), ncol(B), 3))
#' res <- suppressWarnings(dkge_component_stats(fit,
#'                             centroids = centroids,
#'                             mapper = "ridge",
#'                             inference = NULL,
#'                             components = 1))
#' head(res$summary)
#' }
#' @export
dkge_component_stats <- function(fit,
                                 mapper = "sinkhorn",
                                 centroids = NULL,
                                 sizes = NULL,
                                 inference = NULL,
                                 medoid = 1L,
                                 components = NULL,
                                 adjust = "fdr",
                                 ...) {
  stopifnot(inherits(fit, "dkge"))
  .Deprecated(
    "dkge_prepare_alignment",
    package = "dkge",
    msg = paste0(
      "`dkge_component_stats()` is deprecated and descriptive only; use the ",
      "typed alignment and `dkge_infer_aligned()` workflow for inference."
    )
  )
  if (!is.null(inference)) {
    .dkge_abort(
      paste0(
        "Legacy component correspondence is fitted from full-fit loadings and ",
        "is descriptive/ineligible. `dkge_component_stats()` cannot compute ",
        "inferential p-values or significance."
      ),
      "dkge_alignment_ineligible_error"
    )
  }

  centroids <- centroids %||% fit$centroids %||% fit$input$centroids %||%
    stop("Centroids must be supplied or stored in the fit object.")

  # Build mapper specification
  mapper_spec <- .dkge_resolve_mapper_spec(mapper, method = NULL, dots = list(...))

  loadings <- .dkge_fit_subject_loadings(fit)
  rank <- ncol(loadings[[1]])

  if (is.null(components)) {
    comp_idx <- seq_len(rank)
  } else if (is.numeric(components)) {
    comp_idx <- components
  } else {
    comp_idx <- match(components, colnames(fit$U))
  }

  transport <- suppressWarnings(dkge_transport_loadings_to_medoid(fit,
                                                 medoid = medoid,
                                                 centroids = centroids,
                                                 loadings = loadings,
                                                 sizes = sizes,
                                                 mapper = mapper_spec))

  # transport$subjects is a length-`rank` list (one S x Q matrix per component);
  # select the requested components, not columns (clusters) of every component.
  if (any(is.na(comp_idx)) || any(comp_idx < 1L) || any(comp_idx > length(transport$subjects))) {
    stop("`components` must index existing components (1..rank).", call. = FALSE)
  }
  subj_mats <- transport$subjects[comp_idx]
  alignment <- transport$cache
  .dkge_validate_fitted_alignment_object(alignment)
  if (!isTRUE(alignment$eligibility$solver_converged) &&
      !is.na(alignment$eligibility$solver_converged)) {
    .dkge_abort(
      "Legacy component transport did not satisfy its numerical contract.",
      "dkge_alignment_numerical_error"
    )
  }
  descriptive <- .dkge_component_descriptive(subj_mats, comp_idx)

  list(summary = descriptive$summary,
       statistics = descriptive$means,
       transport = subj_mats,
       eligibility = alignment$eligibility,
       metadata = list(
         status = "descriptive",
         inferential = FALSE,
         subject_weighting = "equal_subject",
         reference_subject = alignment$reference_subject %||% alignment$medoid
       ))
}

#' @rdname dkge_component_stats
#' @param file Path to the CSV file where component statistics will be written.
#' @export
dkge_write_component_stats <- function(fit, file, ...) {
  res <- dkge_component_stats(fit, ...)
  utils::write.csv(res$summary, file = file, row.names = FALSE)
  invisible(res)
}

.dkge_component_descriptive <- function(subj_mats, comp_idx) {
  means <- lapply(subj_mats, colMeans)
  sds <- lapply(subj_mats, function(Y) apply(Y, 2, stats::sd))
  out <- Map(function(mean_value, sd_value, Y, comp) {
    data.frame(component = comp,
               cluster = seq_along(mean_value),
               mean = mean_value,
               sd = sd_value,
               n_subjects = nrow(Y),
               stringsAsFactors = FALSE)
  }, means, sds, subj_mats, comp_idx)
  list(summary = do.call(rbind, out), means = means)
}
