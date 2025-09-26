#' Write a group map as NIfTI using a medoid label image
#'
#' Paints medoid-level values onto the reference parcellation and either
#' returns a `neuroim2::BrainVolume` or writes it to disk.
#'
#' @param group_values Numeric vector of medoid-cluster values (length Q). When
#'   named, entries are matched to label IDs before fallback to positional order.
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
  if (exists("BrainVolume", envir = neuro_ns, inherits = FALSE)) {
    volume_ctor <- utils::getFromNamespace("BrainVolume", neuro_ns)
  } else if (exists("NeuroVol", envir = neuro_ns, inherits = FALSE)) {
    volume_ctor <- utils::getFromNamespace("NeuroVol", neuro_ns)
  } else {
    stop("neuroim2: neither 'BrainVolume' nor 'NeuroVol' constructors are available.")
  }
  write_nifti <- NULL
  if (!is.null(out_file)) {
    if (exists("write_nifti", envir = neuro_ns, inherits = FALSE)) {
      write_nifti <- utils::getFromNamespace("write_nifti", neuro_ns)
    } else if (exists("write_vol", envir = neuro_ns, inherits = FALSE)) {
      write_nifti <- utils::getFromNamespace("write_vol", neuro_ns)
    } else {
      stop("neuroim2: volume writing helper ('write_nifti' or 'write_vol') not available.")
    }
  }

  lab_mat <- values_fun(lab_img)
  ids <- sort(unique(as.vector(lab_mat)))
  ids <- ids[ids > 0]
  Q <- length(ids)
  stopifnot(length(group_values) == Q)

  gv_names <- names(group_values)
  use_names <- !is.null(gv_names) && any(nzchar(gv_names))
  mapping_idx <- NULL
  if (use_names) {
    if (any(!nzchar(gv_names))) {
      stop("group_values names must be non-empty when provided.")
    }
    if (anyDuplicated(gv_names)) {
      dup_table <- as.data.frame(table(gv_names), stringsAsFactors = FALSE)
      stop(sprintf("group_values names must be unique. Duplicate counts:\n%s",
                   paste(capture.output(dup_table), collapse = "\n")))
    }
    mapping_idx <- match(as.character(ids), gv_names)
    if (anyNA(mapping_idx)) {
      warning("Some label IDs not found among group_values names; using positional mapping.")
      mapping_idx <- NULL
      use_names <- FALSE
    }
  }

  if (!use_names) {
    group_values_mapped <- group_values
  } else {
    group_values_mapped <- group_values[mapping_idx]
  }

  vox <- array(0, dim = dim(lab_mat))
  idx <- match(lab_mat, ids, nomatch = 0L)
  nz <- which(idx > 0L)
  vox[nz] <- group_values_mapped[idx[nz]]

  out_img <- volume_ctor(vox, space = space_fun(medoid_labels))
  if (is.null(out_file)) {
    return(out_img)
  }
  stopifnot(!is.null(write_nifti))
  write_nifti(out_img, out_file)
  out_file
}
