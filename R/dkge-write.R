#' Write a group map as NIfTI using a medoid label image
#'
#' Paints medoid-level values onto the reference parcellation and either
#' returns a `neuroim2::BrainVolume` or writes it to disk.
#'
#' @param group_values Numeric vector of medoid-cluster values (length Q).
#' @param medoid_labels A `neuroim2::BrainVolume` containing integer medoid labels.
#' @param label_table Optional data frame with cluster metadata (currently unused).
#' @param out_file Optional output path (`.nii` or `.nii.gz`). When `NULL`, the
#'   painted volume is returned without writing to disk.
#' @return Either the output path (when `out_file` is supplied) or a
#'   `neuroim2::BrainVolume`.
#' @export
dkge_write_group_map <- function(group_values, medoid_labels, label_table = NULL, out_file = NULL) {
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    stop("Install 'neuroim2' to write NIfTI group maps.")
  }
  lab_img <- medoid_labels
  neuro_ns <- asNamespace("neuroim2")
  values_fun <- utils::getFromNamespace("values", neuro_ns)
  space_fun <- utils::getFromNamespace("space", neuro_ns)
  brain_volume <- utils::getFromNamespace("BrainVolume", neuro_ns)
  write_nifti <- utils::getFromNamespace("write_nifti", neuro_ns)

  lab_mat <- values_fun(lab_img)
  ids <- sort(unique(as.vector(lab_mat)))
  ids <- ids[ids > 0]
  Q <- length(ids)
  stopifnot(length(group_values) == Q)

  vox <- array(0, dim = dim(lab_mat))
  idx <- match(lab_mat, ids, nomatch = 0L)
  nz <- which(idx > 0L)
  vox[nz] <- group_values[idx[nz]]

  out_img <- brain_volume(vox, space = space_fun(medoid_labels))
  if (is.null(out_file)) {
    return(out_img)
  }
  write_nifti(out_img, out_file)
  out_file
}
