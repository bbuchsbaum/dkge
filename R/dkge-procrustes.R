# dkge-procrustes.R (robust K-Procrustes utilities)
# Provides numerically stable alignment/consensus helpers for DKGE bases.

`%||%` <- function(a, b) if (is.null(a)) b else a

# Internal helper: kernel roots with jitter (exported via design-kernel)
.dkge_kernel_roots <- function(K, jitter = 1e-10) {
  Ksym <- (K + t(K))/2
  ee <- eigen(Ksym, symmetric = TRUE)
  vals <- pmax(ee$values, jitter)
  V <- ee$vectors
  Khalf  <- V %*% diag(sqrt(vals),  length(vals)) %*% t(V)
  Kihalf <- V %*% diag(1/sqrt(vals), length(vals)) %*% t(V)
  list(Khalf = Khalf, Kihalf = Kihalf, evals = vals)
}

#' Robust K-orthonormalization
#'
#' Ensures the columns of `W` are orthonormal with respect to the
#' design kernel metric: U^T K U = I.
#'
#' @param W q×r matrix (columns = basis vectors)
#' @param K q×q design kernel (PSD)
#' @param Kroots optional precomputed kernel roots from `.dkge_kernel_roots`
#' @return q×r matrix with K-orthonormal columns
#' @export
dkge_k_orthonormalize <- function(W, K, Kroots = NULL) {
  stopifnot(is.matrix(W), is.matrix(K), nrow(K) == nrow(W))
  Kroots <- Kroots %||% .dkge_kernel_roots(K)
  # Use Euclidean QR on K^{1/2} W for stability
  Q <- qr.Q(qr(Kroots$Khalf %*% W))
  Kroots$Kihalf %*% Q
}

#' K-orthogonal Procrustes alignment
#'
#' Aligns basis `U` to reference `Uref` by solving
#' max_R tr((U_ref^T K U) R) subject to R^T R = I.
#'
#' @param Uref reference basis (q×r)
#' @param U basis to align (q×r)
#' @param K q×q design kernel
#' @param allow_reflection logical; if FALSE, forces det(R)=+1
#' @return list(U_aligned, R, d=sum(singular values), cosines, det)
#' @export
dkge_procrustes_K <- function(Uref, U, K, allow_reflection = TRUE) {
  stopifnot(is.matrix(Uref), is.matrix(U), is.matrix(K))
  stopifnot(nrow(Uref) == nrow(U), ncol(Uref) == ncol(U), nrow(K) == nrow(U))
  C <- t(Uref) %*% K %*% U
  sv <- svd(C)
  R <- sv$u %*% t(sv$v)
  if (!allow_reflection && det(R) < 0) {
    Uu <- sv$u; Vv <- sv$v
    Uu[, ncol(Uu)] <- -Uu[, ncol(Uu)]
    R <- Uu %*% t(Vv)
  }
  list(U_aligned = U %*% R,
       R = R,
       d = sum(sv$d),
       cosines = pmin(sv$d, 1),
       determinant = det(R))
}

#' Align multiple bases to a reference
#'
#' @param U_list list of q×r bases
#' @param K q×q design kernel
#' @param ref reference basis or index (default 1)
#' @param allow_reflection logical passed to `dkge_procrustes_K`
#' @param weights optional numeric weights (stored for convenience)
#' @return list(U_aligned, R, Uref, score, weights)
#' @export
dkge_align_bases_K <- function(U_list, K, ref = 1L,
                               allow_reflection = TRUE,
                               weights = NULL) {
  stopifnot(length(U_list) >= 1)
  Uref <- if (is.matrix(ref)) ref else U_list[[ref]]
  aligned <- vector("list", length(U_list))
  Rlist   <- vector("list", length(U_list))
  scores  <- numeric(length(U_list))
  for (i in seq_along(U_list)) {
    pr <- dkge_procrustes_K(Uref, U_list[[i]], K, allow_reflection)
    aligned[[i]] <- pr$U_aligned
    Rlist[[i]] <- pr$R
    scores[i] <- pr$d
  }
  list(U_aligned = aligned,
       R = Rlist,
       Uref = Uref,
       score = scores,
       weights = weights)
}

#' Consensus K-orthonormal basis (K-Procrustes mean)
#'
#' Iteratively aligns bases, takes a weighted average, and retraction via
#' K-orthonormalization until convergence.
#'
#' @param U_list list of q×r K-orthonormal bases
#' @param K q×q design kernel
#' @param weights optional numeric weights (default equal)
#' @param Kroots optional precomputed kernel roots
#' @param max_iter maximum iterations
#' @param tol convergence tolerance on 1 - mean principal cosines
#' @param allow_reflection passed to alignment step
#' @return list(U, iters, converged, gaps, scores)
#' @export
dkge_consensus_basis_K <- function(U_list, K,
                                   weights = NULL,
                                   Kroots = NULL,
                                   max_iter = 50,
                                   tol = 1e-6,
                                   allow_reflection = TRUE) {
  stopifnot(length(U_list) >= 1)
  Kroots <- Kroots %||% .dkge_kernel_roots(K)
  Uc <- U_list[[1]]
  r <- ncol(Uc)
  S <- length(U_list)
  if (is.null(weights)) weights <- rep(1, S)
  weights <- as.numeric(weights)
  weights <- weights / sum(weights)

  gaps <- numeric(0)
  scores_hist <- numeric(0)
  for (it in seq_len(max_iter)) {
    aligned <- dkge_align_bases_K(U_list, K, ref = Uc,
                                  allow_reflection = allow_reflection,
                                  weights = weights)
    Abar <- Reduce(`+`, Map(function(Ui, w) w * Ui, aligned$U_aligned, weights))
    Unew <- dkge_k_orthonormalize(Abar, K, Kroots)
    svals <- svd(t(Uc) %*% K %*% Unew)$d
    gap <- 1 - mean(pmin(svals, 1))
    gaps <- c(gaps, gap)
    scores_hist <- c(scores_hist, mean(aligned$score))
    Uc <- Unew
    if (gap < tol) {
      return(list(U = Uc,
                  iters = it,
                  converged = TRUE,
                  gaps = gaps,
                  scores = scores_hist))
    }
  }
  list(U = Uc,
       iters = max_iter,
       converged = FALSE,
       gaps = gaps,
       scores = scores_hist)
}
