
# dkge-neuroim2-stream.R (v0.5 add-on)
# Subject/run-level streaming readers for neuroim2 on-disk series.

#' Compute cluster-level time series from a NeuroVec and a label image
#'
#' By default uses neuroim2::as.matrix() and rowsum aggregation.
#' If your NeuroVec is memory-mapped (the usual case), this is reasonably memory-light.
#' For very large data, provide a custom 'chunker' that returns a T×V_block matrix per call.
#'
#' @param bv a neuroim2::NeuroVec (T×V)
#' @param labels a neuroim2::NeuroVol of integer cluster IDs (V voxels; 0 = background)
#' @param ids optional vector of cluster IDs to keep (defaults to positive unique labels)
#' @param chunker optional function(i) -> list(mat=T×V_block, vox_idx=integer V_block)
#' @return T×P matrix (P = number of clusters)
#' @export
dkge_cluster_ts <- function(bv, labels, ids = NULL, chunker = NULL) {
  lab <- as.vector(neuroim2::values(labels))
  if (is.null(ids)) {
    ids <- sort(unique(lab)); ids <- ids[ids > 0]
  }
  P <- length(ids)
  id_to_col <- setNames(seq_len(P), ids)

  if (is.null(chunker)) {
    # default: aggregate with rowsum over full matrix
    Y <- t(neuroim2::as.matrix(bv))         # T×V (may be memory-mapped)
    keep <- which(lab > 0)
    grp  <- id_to_col[as.character(lab[keep])]
    # rowsum on transposed for speed (aggregate voxels → clusters)
    Yk <- Y[, keep, drop = FALSE]
    Zt <- rowsum(t(Yk), group = grp, reorder = FALSE) # P×T
    Z  <- t(Zt)                                       # T×P
    return(Z)
  }

  # chunked path: initialize accumulators
  Tn <- dim(bv)[1]  # get time dimension
  Z <- matrix(0, Tn, P)
  counts <- numeric(P)

  i <- 1L
  repeat {
    ch <- chunker(i); if (is.null(ch)) break
    Yb <- ch$mat; vox <- ch$vox_idx
    idx <- which(lab[vox] > 0)
    if (length(idx) == 0) { i <- i + 1L; next }
    grp <- id_to_col[as.character(lab[vox[idx]])]
    Zt  <- rowsum(t(Yb[, idx, drop = FALSE]), group = grp, reorder = FALSE) # P×T
    Z   <- Z + t(Zt)
    counts_tab <- as.numeric(rowsum(rep(1, length(idx)), group = grp, reorder = FALSE))
    # align counts_tab to P
    ugrp <- sort(unique(grp))
    counts[ugrp] <- counts[ugrp] + counts_tab
    i <- i + 1L
  }
  # average by counts where counts > 0
  nz <- counts > 0
  Z[, nz] <- sweep(Z[, nz, drop = FALSE], 2, counts[nz], "/")
  Z
}

#' Compute cluster-level betas for a subject from neuroim2 data + fmrireg
#' @param bv NeuroVec (T×V)
#' @param X T×q design matrix
#' @param labels NeuroVol of clusters
#' @return q×P matrix of betas
#' @export
dkge_cluster_betas <- function(bv, X, labels) {
  Z <- dkge_cluster_ts(bv, labels)     # T×P
  fit <- fmrireg::fmri_ols_fit(X, Z)   # returns list with $betas (q×P)
  fit$betas
}

#' Build a subject/run-level streaming loader backed by neuroim2 objects
#'
#' @param design_objs list of fmridesign objects (or T×q matrices) per subject
#' @param bv_list list of NeuroVec objects or file paths (neuroim2 can read)
#' @param labels_list list of NeuroVol label images (cluster IDs) per subject
#' @param omega_fun optional function(ids, labels) -> vector of cluster weights; default uses cluster sizes
#' @return loader with n(), X(s), B(s), Omega(s)
#' @export
dkge_neuro_loader <- function(design_objs, bv_list, labels_list, omega_fun = NULL) {
  stopifnot(length(design_objs) == length(bv_list), length(bv_list) == length(labels_list))
  S <- length(bv_list)

  # normalize inputs
  get_bv <- function(i) {
    bv <- bv_list[[i]]
    if (inherits(bv, "NeuroVec")) return(bv)
    # try reading from path
    neuroim2::read_vec(bv)
  }
  get_X <- function(i) {
    Xi <- design_objs[[i]]
    if (is.matrix(Xi)) return(Xi)
    if (inherits(Xi, "fmridesign")) {
      return(fmridesign::design_matrix(Xi))
    }
    stop("Provide design matrices or fmridesign objects.")
  }
  if (is.null(omega_fun)) {
    omega_fun <- function(labels) {
      lab <- as.vector(neuroim2::values(labels))
      ids <- sort(unique(lab)); ids <- ids[ids > 0]
      tab <- table(lab[lab > 0])
      as.numeric(tab[as.character(ids)])
    }
  }

  list(
    n = function() S,
    X = function(s) get_X(s),
    B = function(s) {
      bv <- get_bv(s)
      X  <- get_X(s)
      labs <- labels_list[[s]]
      dkge_cluster_betas(bv, X, labs)
    },
    Omega = function(s) {
      labs <- labels_list[[s]]
      omega_fun(labs)
    }
  )
}
