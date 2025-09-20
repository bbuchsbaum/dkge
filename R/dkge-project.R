# dkge-project.R
# Projection helpers wrapping the multivarious interface for DKGE objects.

#' Preprocess a subject block into DKGE training space
#'
#' Applies the same transformations used during fitting (shared ruler,
#' kernel whitening, optional spatial weights, and subject weighting) to a
#' subject's beta matrix.
#'
#' @param fit A `dkge` object.
#' @param B_s q×P matrix of subject betas.
#' @param Omega_s Optional weights (vector length P or P×P matrix) matching the
#'   columns of `B_s`.
#' @param w_s Optional subject-level weight (defaults to 1 when omitted).
#' @return q×P matrix in the DKGE training space.
#' @export
dkge_transform_block <- function(fit, B_s, Omega_s = NULL, w_s = NULL) {
  stopifnot(inherits(fit, "dkge"))
  effects <- fit$effects
  B_s <- as.matrix(B_s)
  if (!is.null(effects)) {
    stopifnot(nrow(B_s) == length(effects))
    if (!is.null(rownames(B_s))) {
      idx <- match(effects, rownames(B_s))
      if (anyNA(idx)) {
        stop("Row names of new betas must cover all training effects.")
      }
      B_s <- B_s[idx, , drop = FALSE]
    }
  }
  R <- fit$R
  Khalf <- fit$Khalf
  stopifnot(nrow(B_s) == nrow(R))
  Btil <- t(R) %*% B_s
  if (!is.null(Omega_s)) {
    if (is.vector(Omega_s)) {
      stopifnot(length(Omega_s) == ncol(Btil))
      Btil <- Btil * rep(sqrt(as.numeric(Omega_s)), each = nrow(Btil))
    } else {
      Omega_s <- as.matrix(Omega_s)
      stopifnot(nrow(Omega_s) == ncol(Btil), ncol(Omega_s) == ncol(Btil))
      Btil <- Btil %*% sqrtm_sym(Omega_s)
    }
  }
  weight <- if (is.null(w_s)) 1 else sqrt(max(as.numeric(w_s), 0))
  weight * (Khalf %*% Btil)
}

#' Preprocess multiple blocks into the DKGE training space
#'
#' @param fit A `dkge` object returned by [dkge_fit()].
#' @param B_list List of subject beta matrices.
#' @param Omega_list Optional list of spatial weights aligned with `B_list`.
#' @param w Optional vector of subject weights for the new data.
#' @return q×(sum P_s) matrix matching the training block layout.
#' @export
dkge_preprocess_blocks <- function(fit, B_list, Omega_list = NULL, w = NULL) {
  stopifnot(inherits(fit, "dkge"))
  S <- length(B_list)
  if (is.null(Omega_list)) Omega_list <- vector("list", S)
  if (is.null(w)) w <- rep(1, S)
  stopifnot(length(Omega_list) == S, length(w) == S)
  total_cols <- nrow(fit$v)
  q <- nrow(fit$U)
  X <- matrix(0, q, total_cols)
  for (s in seq_len(S)) {
    idx <- fit$block_indices[[s]]
    Xs <- dkge_transform_block(fit, B_list[[s]], Omega_list[[s]], w[s])
    if (ncol(Xs) != length(idx)) {
      stop(sprintf("Block %d produced %d columns but %d were expected.",
                   s, ncol(Xs), length(idx)))
    }
    X[, idx] <- Xs
  }
  X
}

#' Project new blocks into DKGE score space
#'
#' @param fit A `dkge` object.
#' @inheritParams dkge_preprocess_blocks
#' @return Matrix of projected scores (q×rank).
#' @export
dkge_project_blocks <- function(fit, B_list, Omega_list = NULL, w = NULL) {
  Xnew <- dkge_preprocess_blocks(fit, B_list, Omega_list, w)
  multivarious::project(fit, Xnew)
}

#' @title Project DKGE data into component space
#'
#' @description Functions for projecting subject-standardised betas or new blocks into DKGE component coordinates.
#' @param fit A `dkge` object.
#' @param s Block index (subject) to project against.
#' @param B_s Beta matrix for the new data block.
#' @inheritParams dkge_transform_block
#' @param least_squares Logical; pass to [multivarious::project_block()].

#' @description Convenience wrapper for projecting row-standardised betas onto DKGE components.
#' @param fit A `dkge` object.
#' @param Btil Either a q×P matrix or a list of such matrices (e.g. `fit$Btil`).
#' @return List of P×rank matrices; returns a single matrix when `Btil` is a matrix.
#' @describeIn dkge_project_block Project subject-standardised betas into component space
#' @export
dkge_project_btil <- function(fit, Btil) {
  stopifnot(inherits(fit, "dkge"))
  KsU <- fit$K %*% fit$U
  project_one <- function(mat) {
    mat <- as.matrix(mat)
    stopifnot(nrow(mat) == nrow(fit$U))
    t(mat) %*% KsU
  }
  if (is.list(Btil)) {
    lapply(Btil, project_one)
  } else {
    project_one(Btil)
  }
}

#' Project a single block against its training block index
#'
#' @param fit A `dkge` object.
#' @param s Block index (subject) to project against.
#' @param B_s Beta matrix for the new data block.
#' @inheritParams dkge_transform_block
#' @param least_squares Logical; pass to [multivarious::project_block()].
#' @return Projection scores (q×rank) restricted to block `s`.
#' @rdname dkge_project_block
#' @export
dkge_project_block <- function(fit, s, B_s, Omega_s = NULL, w_s = NULL,
                               least_squares = TRUE) {
  stopifnot(inherits(fit, "dkge"))
  stopifnot(s >= 1L, s <= length(fit$block_indices))
  Xs <- dkge_transform_block(fit, B_s, Omega_s, w_s)
  if (ncol(Xs) != length(fit$block_indices[[s]])) {
    stop("New block has different width than training block.")
  }
  multivarious::project_block(fit, new_data = Xs, block = s,
                              least_squares = least_squares)
}

#' Project a new cluster/voxel vector onto DKGE components
#'
#' @param fit A `dkge` object.
#' @param b Numeric vector of length q (effects).
#' @param omega Optional scalar or matrix weight.
#' @param w Optional scalar subject weight.
#' @return Numeric vector of length `rank` (component scores).
#' @export
dkge_project_cluster <- function(fit, b, omega = 1, w = 1) {
  stopifnot(inherits(fit, "dkge"))
  b <- as.numeric(b)
  stopifnot(length(b) == nrow(fit$U))
  ctil <- t(fit$R) %*% matrix(b, ncol = 1)
  if (is.matrix(omega)) {
    ctil <- ctil %*% sqrtm_sym(omega)
  } else {
    ctil <- ctil * sqrt(as.numeric(omega))
  }
  x <- sqrt(max(as.numeric(w), 0)) * (fit$Khalf %*% ctil)
  as.numeric(multivarious::project_vars(fit, x))
}

#' Project multiple cluster/voxel vectors
#'
#' @param fit A `dkge` object.
#' @param B q×P matrix of cluster betas.
#' @param omega_vec Optional vector of per-cluster weights.
#' @param w Optional subject weight.
#' @return P×rank matrix of projected coordinates.
#' @export
dkge_project_clusters <- function(fit, B, omega_vec = NULL, w = 1) {
  stopifnot(inherits(fit, "dkge"))
  B <- as.matrix(B)
  stopifnot(nrow(B) == nrow(fit$U))
  ctil <- t(fit$R) %*% B
  if (!is.null(omega_vec)) {
    stopifnot(length(omega_vec) == ncol(B))
    ctil <- ctil * rep(sqrt(as.numeric(omega_vec)), each = nrow(ctil))
  }
  x <- sqrt(max(as.numeric(w), 0)) * (fit$Khalf %*% ctil)
  t(multivarious::project_vars(fit, x))
}
