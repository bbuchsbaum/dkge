# dkge-voxel.R
# Utilities for voxel-level consensus mapping.

#' Transport DKGE quantities directly to voxel space
#'
#' @param fit A `dkge` object containing subject loadings and centroids.
#' @param values List of subject value vectors (one per subject, length P_s).
#' @param voxels List of voxel feature matrices.
#' @param coords Optional list of voxel coordinate matrices.
#' @param mapper Mapper specification or shorthand (defaults to `"ridge"`).
#' @param sizes Optional list of cluster masses (one vector per subject).
#' @param ... Additional mapper parameters.
#'
#' @return List with `subj_values` (S × V matrix) and `value` (mean across subjects).
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

  loadings <- lapply(fit$Btil, function(Bts) t(Bts) %*% fit$K %*% fit$U)
  if (is.null(coords)) {
    coords <- vector("list", length(voxels))
  }

  mapped <- vector("list", length(values))
  for (s in seq_along(values)) {
    map_fit <- fit_mapper(mapper_spec,
                          source_feat = loadings[[s]],
                          target_feat = voxels[[s]],
                          source_weights = sizes[[s]],
                          source_xyz = fit$centroids[[s]],
                          target_xyz = coords[[s]])
    mapped[[s]] <- predict_mapper(map_fit, values[[s]])
  }
  subj_mat <- do.call(rbind, mapped)
  list(subj_values = subj_mat,
       value = colMeans(subj_mat))
}
