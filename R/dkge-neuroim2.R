# dkge-neuroim2.R
# Helpers integrating DKGE with neuroim2 objects and streaming loaders.

#' Aggregate voxel time series into cluster means
#'
#' @param bv A `neuroim2::NeuroVec` object (4D time-series data, TxXxYxZ).
#' @param labels A `neuroim2::NeuroVol` with integer cluster identifiers
#'   (0 indicates background).
#' @param ids Optional subset of cluster IDs to retain.
#' @param chunker Optional function returning list(mat = TxV_block, vox_idx).
#' @return Matrix of dimension TxP where P is the number of clusters retained.
#' @export
dkge_cluster_ts <- function(bv, labels, ids = NULL, chunker = NULL) {
  lab <- as.vector(neuroim2::values(labels))
  if (is.null(ids)) {
    ids <- sort(unique(lab))
    ids <- ids[ids > 0]
  }
  p <- length(ids)
  id_to_col <- setNames(seq_len(p), ids)

  if (is.null(chunker)) {
    # Convert NeuroVec to matrix (TxV format)
    y_mat <- neuroim2::as.matrix(bv)
    # y_mat is now voxelsxtime, need to transpose for TxV
    y_mat <- t(y_mat)
    keep <- which(lab > 0)
    grp <- id_to_col[as.character(lab[keep])]
    y_k <- y_mat[, keep, drop = FALSE]
    z_t <- rowsum(t(y_k), group = grp, reorder = FALSE)
    return(t(z_t))
  }

  t_n <- dim(bv)[1]
  z_mat <- matrix(0, t_n, p)
  counts <- numeric(p)
  i <- 1L
  repeat {
    chunk <- chunker(i)
    if (is.null(chunk)) break
    y_b <- chunk$mat
    vox <- chunk$vox_idx
    valid <- which(lab[vox] > 0)
    if (!length(valid)) {
      i <- i + 1L
      next
    }
    grp <- id_to_col[as.character(lab[vox[valid]])]
    z_t <- rowsum(t(y_b[, valid, drop = FALSE]), group = grp,
                  reorder = FALSE)
    z_mat <- z_mat + t(z_t)
    counts_tab <- as.numeric(
      rowsum(rep(1, length(valid)), group = grp, reorder = FALSE)
    )
    unique_grp <- sort(unique(grp))
    counts[unique_grp] <- counts[unique_grp] + counts_tab
    i <- i + 1L
  }

  nz <- counts > 0
  z_mat[, nz] <- sweep(z_mat[, nz, drop = FALSE], 2, counts[nz], "/")
  z_mat
}

#' Compute cluster-level betas from neuroim2 objects
#'
#' @param bv A `neuroim2::NeuroVec` object containing 4D time-series data.
#' @param x_mat Subject design matrix (Txq).
#' @param labels A `neuroim2::NeuroVol` object with cluster label assignments.
#' @return Matrix of GLM betas (qxP) where P is the number of clusters.
#' @export
dkge_cluster_betas <- function(bv, x_mat, labels) {
  z_mat <- dkge_cluster_ts(bv, labels)
  fit <- fmrireg::fmri_ols_fit(x_mat, z_mat)
  fit$betas
}

#' Build a streaming loader backed by neuroim2 objects
#'
#' The returned loader exposes `n()`, `X(s)`, `B(s)`, and `Omega(s)` methods
#' compatible with streaming DKGE fits.
#'
#' @param design_objs List of design matrices or `fmridesign` objects.
#' @param bv_list List of `neuroim2::NeuroVec` objects, `ClusteredNeuroVec`
#'   objects, or file paths readable by `neuroim2::read_vec()`.
#' @param labels_list List of `neuroim2::NeuroVol` label volumes (same length
#'   as `bv_list`). Can be NULL if bv_list contains ClusteredNeuroVec objects.
#' @param omega_fun Optional function mapping a label volume or
#'   ClusteredNeuroVec to cluster weights. Default uses cluster sizes as
#'   weights.
#' @return Loader list suitable for streaming DKGE fits.
#' @export
dkge_neuro_loader <- function(design_objs, bv_list,
                               labels_list = NULL,
                               omega_fun = NULL) {
  # Check if all objects are ClusteredNeuroVec
  all_clustered <- all(vapply(
    bv_list,
    function(x) inherits(x, "ClusteredNeuroVec"),
    logical(1)
  ))

  if (!all_clustered) {
    stopifnot(!is.null(labels_list),
              length(bv_list) == length(labels_list))
  }

  stopifnot(length(design_objs) == length(bv_list))
  s_count <- length(bv_list)

  get_bv <- function(i) {
    bv <- bv_list[[i]]
    if (inherits(bv, c("NeuroVec", "ClusteredNeuroVec"))) return(bv)
    neuroim2::read_vec(bv)
  }

  get_x <- function(i) {
    x_i <- design_objs[[i]]
    if (is.matrix(x_i)) return(x_i)
    if (inherits(x_i, "fmridesign")) {
      return(fmridesign::design_matrix(x_i))
    }
    stop("Provide design matrices or fmridesign objects.")
  }

  if (is.null(omega_fun)) {
    omega_fun <- function(labels) {
      lab <- as.vector(neuroim2::values(labels))
      tab <- table(lab[lab > 0])
      # Return cluster sizes in order of cluster IDs
      cluster_ids <- as.integer(names(tab))
      as.numeric(tab[order(cluster_ids)])
    }
  }

  list(
    n = function() s_count,
    X = function(s) get_x(s),
    B = function(s) {
      bv <- get_bv(s)
      x_mat <- get_x(s)

      if (inherits(bv, "ClusteredNeuroVec")) {
        # For ClusteredNeuroVec, compute betas directly
        cluster_ts <- neuroim2::as.matrix(bv)  # TxK matrix
        fit <- fmrireg::fmri_ols_fit(x_mat, cluster_ts)
        fit$betas
      } else {
        # For regular NeuroVec, use labels for clustering
        labs <- labels_list[[s]]
        dkge_cluster_betas(bv, x_mat, labs)
      }
    },
    Omega = function(s) {
      if (all_clustered) {
        bv <- get_bv(s)
        if (!is.null(omega_fun)) {
          omega_fun(bv)
        } else {
          # Use cluster sizes from ClusteredNeuroVec
          cluster_sizes <- table(bv@cl_map[bv@cl_map > 0])
          as.numeric(cluster_sizes)
        }
      } else {
        labs <- labels_list[[s]]
        omega_fun(labs)
      }
    }
  )
}