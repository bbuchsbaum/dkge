# dkge-voxel.R
# Utilities for voxel-level consensus mapping.

#' Descriptively transport DKGE quantities directly to voxel space
#'
#' This legacy convenience path learns voxel correspondence from full-fit
#' loadings and returns untyped display values. It is descriptive only and its
#' output is not accepted by [dkge_infer_aligned()]. For group functional
#' alignment, fit typed features and use
#' [dkge_transport_contrasts_to_reference()] or [dkge_align_to_template()].
#'
#' @param fit A `dkge` object containing subject loadings and centroids.
#' @param values List of subject value vectors (one per subject, length P_s).
#' @param voxels List of voxel feature matrices.
#' @param coords Optional list of voxel coordinate matrices.
#' @param mapper Mapper specification or shorthand (defaults to `"ridge"`).
#' @param sizes Optional list of cluster masses (one vector per subject).
#' @param ... Additional mapper parameters.
#'
#' @return Descriptive list with `subj_values` (S x V matrix), `value` (mean
#'   across subjects), and explicit non-inferential metadata.
#' @examples
#' \donttest{
#' toy <- dkge_sim_toy(
#'   factors = list(A = list(L = 2), B = list(L = 3)),
#'   active_terms = c("A", "B"), S = 3, P = 15, snr = 5
#' )
#' fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 2)
#' fit$centroids <- lapply(toy$B_list, function(B) matrix(rnorm(ncol(B) * 3), ncol(B), 3))
#' values <- lapply(toy$B_list, function(B) rnorm(ncol(B)))
#' voxels <- lapply(toy$B_list, function(B) matrix(rnorm(10 * fit$rank), 10, fit$rank))
#' out <- dkge_transport_to_voxels(fit, values = values, voxels = voxels, mapper = "ridge")
#' dim(out$subj_values)
#' }
#' @export
dkge_transport_to_voxels <- function(fit,
                                     values,
                                     voxels,
                                     coords = NULL,
                                     mapper = "ridge",
                                     sizes = NULL,
                                     ...) {
  stopifnot(inherits(fit, "dkge"))
  mapper_spec <- .dkge_resolve_mapper_spec(mapper, method = NULL, dots = list(...))

  loadings <- .dkge_fit_subject_loadings(fit)
  if (is.null(coords)) {
    coords <- vector("list", length(voxels))
  }

  mapped <- vector("list", length(values))
  for (s in seq_along(values)) {
    sw <- if (is.null(sizes)) NULL else sizes[[s]]
    map_fit <- fit_mapper(mapper_spec,
                          source_feat = loadings[[s]],
                          target_feat = voxels[[s]],
                          source_weights = sw,
                          source_xyz = fit$centroids[[s]],
                          target_xyz = coords[[s]])
    mapped[[s]] <- predict_mapper(map_fit, values[[s]])
  }
  subj_mat <- do.call(rbind, mapped)
  list(
    subj_values = subj_mat,
    value = colMeans(subj_mat),
    metadata = list(
      status = "descriptive",
      inferential = FALSE,
      correspondence_source = "full_fit_loadings"
    )
  )
}
