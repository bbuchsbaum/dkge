# dkge-anchor-targets.R
# Utilities for building anchor-space classification targets.

#' Assemble anchor targets from prototype sets
#'
#' Convenience helper that stacks calls to
#' [dkge_anchor_contrast_from_prototypes()] so the resulting matrix can be fed
#' directly to [dkge_classify()] as a weight specification.
#'
#' @param anchors Matrix of anchor coordinates (`L x d`).
#' @param prototypes Named list whose elements are matrices (rows = prototypes
#'   in the same feature space as `anchors`). List names become class labels.
#' @param negatives Optional named list of matrices providing negative
#'   prototypes per class (matched by name). When `NULL`, classes are contrasted
#'   against the origin.
#' @param sigma Optional bandwidth passed to
#'   [dkge_anchor_contrast_from_prototypes()]. Defaults to the per-class median
#'   heuristic.
#' @param normalize Logical indicator forwarded to
#'   [dkge_anchor_contrast_from_prototypes()].
#'
#' @return Matrix with one row per class and `nrow(anchors)` columns.
#' @export
dkge_anchor_targets_from_prototypes <- function(anchors,
                                                prototypes,
                                                negatives = NULL,
                                                sigma = NULL,
                                                normalize = TRUE) {
  stopifnot(is.matrix(anchors))
  stopifnot(is.list(prototypes), length(prototypes) > 0L)
  class_names <- names(prototypes)
  if (is.null(class_names) || any(!nzchar(class_names))) {
    stop("`prototypes` must be a named list.")
  }

  neg_list <- list()
  if (!is.null(negatives)) {
    stopifnot(is.list(negatives))
    for (nm in class_names) {
      neg_list[[nm]] <- negatives[[nm]] %||% NULL
    }
  }

  weights <- lapply(class_names, function(nm) {
    pos <- prototypes[[nm]]
    if (is.vector(pos)) pos <- matrix(pos, nrow = 1)
    stopifnot(is.matrix(pos), ncol(pos) == ncol(anchors))
    neg <- neg_list[[nm]]
    if (!is.null(neg) && is.vector(neg)) neg <- matrix(neg, nrow = 1)
    if (!is.null(neg)) {
      stopifnot(is.matrix(neg), ncol(neg) == ncol(anchors))
    }
    dkge_anchor_contrast_from_prototypes(anchors,
                                         positives = pos,
                                         negatives = neg,
                                         sigma = sigma,
                                         normalize = normalize)
  })
  res <- do.call(rbind, weights)
  rownames(res) <- class_names
  res
}

#' Assemble anchor targets from feature-space directions
#'
#' Converts named direction vectors into anchor weight rows using
#' [dkge_anchor_contrast_from_direction()].
#'
#' @param anchors Matrix of anchor coordinates (`L x d`).
#' @param directions Either a named list of numeric vectors (length `d`) or a
#'   matrix whose rows are named directions in the same feature space as
#'   `anchors`.
#' @param sigma Optional bandwidth forwarded to
#'   [dkge_anchor_contrast_from_direction()].
#' @param normalize Logical; when `TRUE` (default) the resulting weights are
#'   L2-normalised.
#'
#' @return Matrix with one row per supplied direction and `nrow(anchors)` columns.
#' @export
dkge_anchor_targets_from_directions <- function(anchors,
                                                directions,
                                                sigma = NULL,
                                                normalize = TRUE) {
  stopifnot(is.matrix(anchors))
  if (is.matrix(directions)) {
    if (is.null(rownames(directions))) {
      stop("Direction matrix must have row names identifying classes.")
    }
    dir_list <- lapply(seq_len(nrow(directions)), function(i) directions[i, ])
    names(dir_list) <- rownames(directions)
  } else if (is.list(directions)) {
    if (is.null(names(directions)) || any(!nzchar(names(directions)))) {
      stop("Named directions are required when supplying a list.")
    }
    dir_list <- directions
  } else {
    stop("`directions` must be a named list or a matrix with row names.")
  }

  weights <- lapply(names(dir_list), function(nm) {
    vec <- dir_list[[nm]]
    dkge_anchor_contrast_from_direction(anchors,
                                        direction = vec,
                                        sigma = sigma,
                                        normalize = normalize)
  })
  res <- do.call(rbind, weights)
  rownames(res) <- names(dir_list)
  res
}
