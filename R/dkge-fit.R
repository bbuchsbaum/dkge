
# dkge-fit.R
# Core DKGE batch fit, LOSO contrasts, and user-facing helpers for multiple subjects.


# -------------------------------------------------------------------------
# Internal helpers used by the fitter ------------------------------------
# -------------------------------------------------------------------------

#' Compute pooled Gram matrix and Cholesky ruler
#'
#' @param X_list List of subject design matrices (`T_sxq`).
#' @param jitter Numerical jitter added to the diagonal for stability.
#' @return List with the pooled Gram matrix and its upper-triangular Cholesky factor.
#' @keywords internal
#' @noRd
.dkge_compute_shared_ruler <- function(X_list, jitter = 1e-10) {
  q <- ncol(X_list[[1]])
  G_pool <- matrix(0, q, q)
  for (Xs in X_list) {
    G_pool <- G_pool + crossprod(Xs)
  }
  diag(G_pool) <- diag(G_pool) + jitter
  list(R = chol(G_pool), G_pool = G_pool)
}

#' Row-standardise subject betas using pooled design Cholesky factor
#'
#' @param B_list List of qxP subject beta matrices.
#' @param R Upper-triangular Cholesky factor from the pooled Gram matrix.
#' @return List of row-standardised betas (`R^T B_s`).
#' @keywords internal
#' @noRd
.dkge_row_standardize <- function(B_list, R) {
  lapply(B_list, function(B) {
    stopifnot(nrow(B) == nrow(R))
    t(R) %*% B
  })
}

#' Compute symmetric square roots of the design kernel
#'
#' @param K qxq positive semi-definite design kernel.
#' @return List containing `Khalf`, `Kihalf`, and the raw eigen decomposition.
#' @keywords internal
#' @noRd
.dkge_kernel_roots <- function(K) {
  Ksym <- (K + t(K)) / 2
  eigK <- eigen(Ksym, symmetric = TRUE)
  vals <- pmax(eigK$values, 1e-10)
  V <- eigK$vectors
  list(
    Khalf = V %*% (diag(sqrt(vals), length(vals))) %*% t(V),
    Kihalf = V %*% (diag(1 / sqrt(vals), length(vals))) %*% t(V),
    eigen = eigK
  )
}

#' Power iteration helper for leading singular value
#'
#' @param X Numeric matrix.
#' @param tol Convergence tolerance.
#' @param max_iter Maximum iterations.
#' @return Approximation to the squared leading singular value of `X`.
#' @keywords internal
#' @noRd
.dkge_leading_sv_squared <- function(X, tol = 1e-6, max_iter = 50) {
  X <- as.matrix(X)
  if (!all(is.finite(X))) {
    X[!is.finite(X)] <- 0
  }
  n <- nrow(X)
  if (n == 0L || ncol(X) == 0L) {
    return(0)
  }
  v <- stats::rnorm(n)
  v_norm <- sqrt(sum(v * v))
  if (!is.finite(v_norm) || v_norm == 0) {
    return(0)
  }
  v <- v / v_norm
  sigma_sq <- 0
  for (iter in seq_len(max_iter)) {
    w <- X %*% (t(X) %*% v)
    w_norm <- sqrt(sum(w * w))
    if (!is.finite(w_norm) || w_norm == 0) {
      break
    }
    v <- w / w_norm
    s_sq_new <- sum((t(X) %*% v)^2)
    if (abs(s_sq_new - sigma_sq) <= tol * max(1, sigma_sq)) {
      sigma_sq <- s_sq_new
      break
    }
    sigma_sq <- s_sq_new
  }
  sigma_sq
}

#' Derive optional subject-level weights
#'
#' @param Btil List of row-standardised betas.
#' @param Omega_list Optional per-subject spatial weights.
#' @param Khalf Kernel square root used for energy computations.
#' @param w_method Weighting scheme (`"mfa_sigma1"`, `"energy"`, or `"none"`).
#' @param w_tau Shrinkage parameter toward equal weights.
#' @return Numeric vector of subject weights.
#' @keywords internal
#' @noRd
.dkge_subject_weights <- function(Btil, Omega_list, Khalf, w_method, w_tau) {
  S <- length(Btil)
  if (w_method == "none") {
    return(rep(1, S))
  }
  # Compute weights on X_s = K^{1/2} Btil_s Omega_s^{1/2} so that the block
  # normalisation matches the metric used in the eigensolve. A numeric vector
  # Omega_s is treated as diagonal reliabilities; a full matrix allows custom
  # covariance or smoothing.
  q <- nrow(Btil[[1]])
  weights <- numeric(S)
  for (s in seq_len(S)) {
    Bts <- Btil[[s]]
    Omega <- Omega_list[[s]]
    KsBt <- if (is.null(Omega)) {
      Khalf %*% Bts
    } else if (is.vector(Omega)) {
      stopifnot(length(Omega) == ncol(Bts))
      Khalf %*% (Bts * rep(sqrt(Omega), each = q))
    } else {
      Omega <- as.matrix(Omega)
      stopifnot(nrow(Omega) == ncol(Bts), ncol(Omega) == ncol(Bts))
      Khalf %*% Bts %*% sqrtm_sym(Omega)
    }
    if (w_method == "mfa_sigma1") {
      sigma_sq <- .dkge_leading_sv_squared(KsBt)
      weights[s] <- 1 / (sigma_sq + 1e-12)
    } else {
      weights[s] <- 1 / (sum(KsBt * KsBt) + 1e-12)
    }
  }
  w_norm <- weights / mean(weights)
  weights <- (1 - w_tau) * w_norm + w_tau * 1
  weights
}

#' Accumulate compressed covariance in the K-metric
#'
#' @param Btil List of row-standardised betas.
#' @param Omega_list Optional spatial weights aligned with `Btil`.
#' @param Khalf Kernel square root used to project into the K-metric.
#' @param weights Subject weights applied during accumulation.
#' @return List with the symmetrised compressed covariance (`Chat`) and per-subject contributions.
#' @keywords internal
#' @noRd
.dkge_accumulate_chat <- function(Btil, Omega_list, Khalf, weights, voxel_weights = NULL) {
  S <- length(Btil)
  q <- nrow(Btil[[1]])
  Chat <- matrix(0, q, q)
  contribs <- vector("list", S)

  scale_columns <- function(B, w) {
    if (is.null(w) || length(w) == 0L) return(B)
    if (length(w) != ncol(B)) {
      w <- rep(w, length.out = ncol(B))
    }
    sweep(B, 2L, sqrt(pmax(w, 0)), "*")
  }

  for (s in seq_len(S)) {
    Bts <- Btil[[s]]
    w_s <- if (is.list(voxel_weights)) voxel_weights[[s]] else voxel_weights
    Bw <- scale_columns(Bts, w_s)
    Omega <- Omega_list[[s]]
    right <- if (is.null(Omega)) {
      Bw %*% t(Bw)
    } else if (is.vector(Omega)) {
      stopifnot(length(Omega) == ncol(Bw))
      W <- sweep(Bw, 2L, sqrt(pmax(Omega, 0)), "*")
      W %*% t(W)
    } else {
      Omega <- as.matrix(Omega)
      stopifnot(nrow(Omega) == ncol(Bw), ncol(Omega) == ncol(Bw))
      Bw %*% Omega %*% t(Bw)
    }
    S_s <- Khalf %*% right %*% Khalf
    contribs[[s]] <- S_s
    Chat <- Chat + weights[s] * S_s
  }
  list(Chat = (Chat + t(Chat)) / 2, contribs = contribs)
}

# -------------------------------------------------------------------------
# Core fit ----------------------------------------------------------------
# -------------------------------------------------------------------------

#' Fit a Design-Kernel Group Embedding (DKGE) model
#'
#' @param data `dkge_data` bundle (preferred) or a list of subject beta matrices.
#' @param designs Per-subject design matrices (ignored when `data` already carries
#'   them). Columns must align with the effects encoded by `K`.
#' @param K qxq design kernel in effect space. Supply either a raw matrix or the
#'   list returned by [design_kernel()].
#' @param Omega_list Optional per-subject spatial reliabilities/weights. Each
#'   element may be `NULL` (no extra weighting), a numeric vector of length
#'   `P_s` (diagonal weights for clusters/voxels), or a full `P_s x P_s` matrix
#'   specifying custom covariance. These weights are applied both when
#'   accumulating the compressed covariance and when computing MFA/energy block
#'   normalisation.
#' @param w_method Subject-level weighting scheme.
#'   * `"mfa_sigma1"` (default): inverse squared leading singular value of
#'     `K^{1/2} Btil_s Omega_s^{1/2}` (Multiple Factor Analysis scaling).
#'   * `"energy"`: inverse Frobenius norm squared of the same block.
#'   * `"none"`: disable block scaling (all weights = 1).
#' @param w_tau Shrinkage parameter (0..1) toward equal weights. 0 keeps the raw
#'   weights, 1 forces equal weighting.
#' @param ridge Ridge regularization parameter (default 0).
#' @param rank Desired rank for the decomposition. If NULL, uses full rank.
#' @param keep_X Logical; when `TRUE`, store the concatenated training matrix
#'   used to build the multiblock projection (can be large).
#' @param cpca_blocks Optional integer vector specifying which effect rows span
#'   the CPCA design subspace. Ignored when `cpca_part = "none"` or when
#'   `cpca_T` is provided.
#' @param cpca_T Optional qxq0 matrix giving the CPCA design basis explicitly.
#'   When supplied it overrides `cpca_blocks`.
#' @param cpca_part Which CPCA-filtered component to estimate: `"none"` (default)
#'   performs the standard DKGE fit; `"design"` uses only the CPCA design part;
#'   `"resid"` uses the residual part; `"both"` fits the design part but also
#'   stores the residual basis for inspection.
#' @param cpca_ridge Optional ridge added to the CPCA-filtered matrices before
#'   eigen-decomposition.
#' @param weights Either `NULL` (default) to apply the weighting scheme specified
#'   by `w_method`, or a [`dkge_weights()`] specification controlling additional
#'   voxel/anchor-level weighting. When supplied, it must inherit from
#'   `dkge_weights` and is resolved via [dkge_weights()].
#' @return A fitted `dkge` object. When the \pkg{multivarious} package is installed
#'   the return value additionally inherits from `multiblock_biprojector`, making
#'   it compatible with the multivarious multiblock interface.
#' @export
dkge_fit <- function(data, designs = NULL, K = NULL, Omega_list = NULL,
                     w_method = c("mfa_sigma1", "energy", "none"),
                     w_tau = 0.3,
                     ridge = 0,
                     rank = NULL,
                     keep_X = FALSE,
                     cpca_blocks = NULL,
                     cpca_T = NULL,
                     cpca_part = c("none", "design", "resid", "both"),
                     cpca_ridge = 0,
                     weights = NULL) {
  w_method <- match.arg(w_method)
  cpca_part <- match.arg(cpca_part)

  prepped <- .dkge_fit_prepare(data,
                               designs = designs,
                               K = K,
                               Omega_list = Omega_list,
                               weights = weights,
                               rank = rank)

  accum <- .dkge_fit_accumulate(prepped,
                                w_method = w_method,
                                w_tau = w_tau)

  solved <- .dkge_fit_solve(prepped,
                            accum,
                            rank = prepped$rank,
                            cpca_part = cpca_part,
                            cpca_blocks = cpca_blocks,
                            cpca_T = cpca_T,
                            cpca_ridge = cpca_ridge,
                            ridge = ridge)

  .dkge_fit_assemble(prepped,
                     accum,
                     solved,
                     keep_X = keep_X,
                     w_method = w_method,
                     w_tau = w_tau,
                     ridge = ridge)
}

#' Symmetric matrix square root helper
#'
#' @param M Symmetric positive semi-definite matrix.
#' @return Matrix square root obtained from the eigen decomposition.
#' @keywords internal
#' @noRd
sqrtm_sym <- function(M) {
  M <- as.matrix(M)
  eig <- eigen((M + t(M)) / 2, symmetric = TRUE)
  vals <- pmax(eig$values, 0)
  eig$vectors %*% (diag(sqrt(vals), length(vals))) %*% t(eig$vectors)
}
