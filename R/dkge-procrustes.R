# dkge-procrustes.R (robust K-Procrustes utilities)
# Provides numerically stable alignment/consensus helpers for DKGE bases.

# Internal helper: kernel roots with jitter (exported via design-kernel)
.dkge_kernel_roots <- function(K, jitter = 1e-10) {
  Ksym <- .dkge_validate_kernel(K)
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
#' @param W qxr matrix (columns = basis vectors)
#' @param K qxq design kernel (PSD)
#' @param Kroots Optional precomputed kernel roots retained for API
#'   compatibility. Orthonormalization is computed from the exact Gram matrix
#'   `t(W) %*% K %*% W`, so null directions of a PSD kernel are not jittered
#'   into artificial metric dimensions.
#' @return qxr matrix with K-orthonormal columns
#' @export
#' @examples
#' K <- diag(5)
#' W <- matrix(rnorm(10), 5, 2)
#' U <- dkge_k_orthonormalize(W, K)
#' # Verify K-orthonormality
#' round(t(U) %*% K %*% U, 10)
dkge_k_orthonormalize <- function(W, K, Kroots = NULL) {
  if (!is.matrix(W) || !is.numeric(W) || any(!is.finite(W))) {
    stop("`W` must be a finite numeric matrix.", call. = FALSE)
  }
  K <- .dkge_validate_kernel(K)
  if (nrow(K) != nrow(W) || ncol(W) < 1L) {
    stop("`W` must have as many rows as `K` and at least one column.", call. = FALSE)
  }

  # Whitening the small r-by-r Gram matrix preserves span(W) and works for
  # semidefinite K exactly when W has full column rank in the K metric.
  gram <- crossprod(W, K %*% W)
  gram <- (gram + t(gram)) / 2
  ee <- eigen(gram, symmetric = TRUE)
  scale <- max(1, max(abs(ee$values)))
  rank_tol <- 1e-10 * scale
  if (min(ee$values) <= rank_tol) {
    stop("`W` must have full column rank in the K metric; its Gram matrix is singular.",
         call. = FALSE)
  }
  inv_sqrt <- ee$vectors %*%
    diag(1 / sqrt(ee$values), length(ee$values)) %*%
    t(ee$vectors)
  U <- W %*% inv_sqrt
  post_error <- max(abs(crossprod(U, K %*% U) - diag(ncol(U))))
  if (!is.finite(post_error) || post_error > 1e-7) {
    stop(sprintf("K-orthonormalization failed its postcondition (Gram error %.3e).",
                 post_error), call. = FALSE)
  }
  U
}

.dkge_validate_basis_K <- function(U, K, label, tol = 1e-6) {
  if (!is.matrix(U) || !is.numeric(U) || any(!is.finite(U))) {
    stop(sprintf("`%s` must be a finite numeric matrix.", label), call. = FALSE)
  }
  if (nrow(U) != nrow(K) || ncol(U) < 1L) {
    stop(sprintf("`%s` must have as many rows as `K` and at least one column.", label),
         call. = FALSE)
  }
  KU <- K %*% U
  err <- max(abs(crossprod(U, KU) - diag(ncol(U))))
  if (!is.finite(err) || err > tol) {
    stop(sprintf("`%s` must be K-orthonormal (maximum Gram error %.3e).",
                 label, err), call. = FALSE)
  }
  KU
}

.dkge_procrustes_K_validated <- function(Uref, U, K,
                                         allow_reflection = TRUE,
                                         KUref = NULL,
                                         KU = NULL) {
  if (!is.matrix(Uref) || !is.matrix(U) ||
      nrow(Uref) != nrow(U) || ncol(Uref) != ncol(U) ||
      nrow(K) != nrow(U)) {
    stop("`Uref` and `U` must be conformable q-by-r matrices for q-by-q `K`.",
         call. = FALSE)
  }
  if (is.null(KUref)) KUref <- .dkge_validate_basis_K(Uref, K, "Uref")
  if (is.null(KU)) KU <- .dkge_validate_basis_K(U, K, "U")
  C <- crossprod(Uref, KU)
  sv <- svd(C)
  R <- sv$v %*% t(sv$u)
  if (!allow_reflection && det(R) < 0) {
    Vv <- sv$v
    Vv[, ncol(Vv)] <- -Vv[, ncol(Vv)]
    R <- Vv %*% t(sv$u)
  }
  objective <- sum(diag(C %*% R))
  list(U_aligned = U %*% R,
       R = R,
       d = objective,
       unconstrained_d = sum(sv$d),
       cosines = pmax(0, pmin(sv$d, 1)),
       determinant = det(R))
}

#' K-orthogonal Procrustes alignment
#'
#' Aligns basis `U` to reference `Uref` by solving
#' max_R tr((U_ref^T K U) R) subject to R^T R = I.
#'
#' @param Uref K-orthonormal reference basis (qxr)
#' @param U K-orthonormal basis to align (qxr)
#' @param K qxq finite, symmetric, positive-semidefinite design kernel
#' @param allow_reflection logical; if FALSE, forces det(R)=+1
#' @return A list with the aligned basis `U_aligned`, orthogonal rotation `R`,
#'   achieved objective `d`, unconstrained objective `unconstrained_d`,
#'   principal `cosines`, and rotation `determinant`. When reflection is
#'   forbidden, `d` can be smaller than `unconstrained_d`.
#' @export
dkge_procrustes_K <- function(Uref, U, K, allow_reflection = TRUE) {
  K <- .dkge_validate_kernel(K)
  .dkge_procrustes_K_validated(Uref, U, K, allow_reflection)
}

#' Align multiple bases to a reference
#'
#' @param U_list list of qxr bases
#' @param K qxq design kernel
#' @param ref reference basis or index (default 1)
#' @param allow_reflection logical passed to `dkge_procrustes_K`
#' @param weights optional numeric weights (stored for convenience)
#' @return list(U_aligned, R, Uref, score, weights)
#' @export
dkge_align_bases_K <- function(U_list, K, ref = 1L,
                               allow_reflection = TRUE,
                               weights = NULL) {
  K <- .dkge_validate_kernel(K)
  .dkge_align_bases_K_validated(U_list, K, ref, allow_reflection, weights)
}

.dkge_align_bases_K_validated <- function(U_list, K, ref = 1L,
                                          allow_reflection = TRUE,
                                          weights = NULL) {
  if (!is.list(U_list) || length(U_list) < 1L) {
    stop("`U_list` must be a non-empty list of basis matrices.", call. = FALSE)
  }
  if (is.matrix(ref)) {
    Uref <- ref
  } else {
    if (length(ref) != 1L || !is.finite(ref) || ref != as.integer(ref) ||
        ref < 1L || ref > length(U_list)) {
      stop("`ref` must be a basis matrix or a valid index into `U_list`.",
           call. = FALSE)
    }
    Uref <- U_list[[as.integer(ref)]]
  }
  if (!is.null(weights)) {
    if (!is.numeric(weights) || length(weights) != length(U_list) ||
        any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) {
      stop("`weights` must be finite, non-negative, match `U_list`, and have positive sum.",
           call. = FALSE)
    }
  }
  aligned <- vector("list", length(U_list))
  Rlist   <- vector("list", length(U_list))
  scores  <- numeric(length(U_list))
  KUref <- .dkge_validate_basis_K(Uref, K, "Uref")
  for (i in seq_along(U_list)) {
    pr <- .dkge_procrustes_K_validated(Uref, U_list[[i]], K,
                                       allow_reflection, KUref = KUref)
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
#' @param U_list list of qxr K-orthonormal bases
#' @param K qxq design kernel
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
  if (!is.list(U_list) || length(U_list) < 1L) {
    stop("`U_list` must be a non-empty list of basis matrices.", call. = FALSE)
  }
  K <- .dkge_validate_kernel(K)
  if (length(max_iter) != 1L || !is.finite(max_iter) ||
      max_iter < 1 || max_iter != as.integer(max_iter)) {
    stop("`max_iter` must be a positive integer.", call. = FALSE)
  }
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`tol` must be finite and positive.", call. = FALSE)
  }
  max_iter <- as.integer(max_iter)
  Uc <- U_list[[1]]
  S <- length(U_list)
  if (is.null(weights)) weights <- rep(1, S)
  weights <- as.numeric(weights)
  if (length(weights) != S || any(!is.finite(weights)) ||
      any(weights < 0) || sum(weights) <= 0) {
    stop("`weights` must be finite, non-negative, match `U_list`, and have positive sum.",
         call. = FALSE)
  }
  weights <- weights / sum(weights)

  gaps <- numeric(0)
  scores_hist <- numeric(0)
  for (it in seq_len(max_iter)) {
    aligned <- .dkge_align_bases_K_validated(
      U_list, K, ref = Uc,
      allow_reflection = allow_reflection,
      weights = weights
    )
    Abar <- Reduce(`+`, Map(function(Ui, w) w * Ui, aligned$U_aligned, weights))
    Unew <- dkge_k_orthonormalize(Abar, K, Kroots)
    svals <- svd(t(Uc) %*% K %*% Unew)$d
    gap <- 1 - mean(pmin(svals, 1))
    gaps <- c(gaps, gap)
    scores_hist <- c(scores_hist, sum(weights * aligned$score))
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
