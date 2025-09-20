
# dkge-write.R (v0.4)
# Writers for neuroim2: paint group maps on medoid parcellation and save NIfTI.

#' Write a group map to NIfTI using a medoid label image (neuroim2)
#'
#' @param group_values Q-vector of medoid-cluster values
#' @param medoid_labels a neuroim2::BrainVolume (integer labels) for the medoid parcellation
#' @param label_table data.frame with columns id (cluster id in labels), name (optional)
#' @param out_file path to write (.nii or .nii.gz); if NULL, returns a BrainVolume
#' @return path or BrainVolume
#' @export
dkge_write_group_map <- function(group_values, medoid_labels, label_table = NULL, out_file = NULL) {
  if (!requireNamespace("neuroim2", quietly = TRUE)) {
    stop("Install 'neuroim2' to write NIfTI group maps.")
  }
  lab_img <- medoid_labels
  lab_mat <- neuroim2::values(lab_img)
  ids <- sort(unique(as.vector(lab_mat)))
  ids <- ids[ids > 0]
  Q <- length(ids)
  stopifnot(length(group_values) == Q)

  # build voxel image by mapping ids -> values
  val_map <- setNames(as.numeric(group_values), ids)
  vox <- array(0, dim = dim(lab_mat))
  idx <- match(lab_mat, ids, nomatch = 0L)
  nz <- which(idx > 0L)
  vox[nz] <- group_values[idx[nz]]

  out_img <- neuroim2::BrainVolume(vox, space = neuroim2::space(medoid_labels))
  if (is.null(out_file)) return(out_img)
  neuroim2::write_nifti(out_img, out_file)
  out_file
}
