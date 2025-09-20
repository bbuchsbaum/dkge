
# dkge-fit.R
# Core DKGE batch fit, LOSO contrasts, and user-facing helpers for multiple subjects.


# -------------------------------------------------------------------------
# Internal helpers used by the fitter ------------------------------------
# -------------------------------------------------------------------------

#' Compute pooled Gram matrix and Cholesky ruler
#'
#' @param X_list List of subject design matrices (`T_s×q`).
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

#' Row-standardise subject betas using shared ruler
#'
#' @param B_list List of q×P subject beta matrices.
#' @param R Upper-triangular shared ruler from the pooled Gram matrix.
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
#' @param K q×q positive semi-definite design kernel.
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
    d <- svd(KsBt, nu = 0, nv = 0)$d
    if (w_method == "mfa_sigma1") {
      weights[s] <- 1 / (d[1]^2 + 1e-12)
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
.dkge_accumulate_chat <- function(Btil, Omega_list, Khalf, weights) {
  S <- length(Btil)
  q <- nrow(Btil[[1]])
  # Each subject contributes S_s = K^{1/2} (Btil_s Omega_s Btil_s^T) K^{1/2}.
  # Omega_s may be NULL, a length-P vector (diagonal weights), or a P×P matrix
  # (full covariance). Contributions are cached for LOSO refits.
  Chat <- matrix(0, q, q)
  contribs <- vector("list", S)
  for (s in seq_len(S)) {
    Bts <- Btil[[s]]
    Omega <- Omega_list[[s]]
    right <- if (is.null(Omega)) {
      Bts %*% t(Bts)
    } else if (is.vector(Omega)) {
      stopifnot(length(Omega) == ncol(Bts))
      W <- Bts * rep(sqrt(Omega), each = q)
      W %*% t(W)
    } else {
      Omega <- as.matrix(Omega)
      stopifnot(nrow(Omega) == ncol(Bts), ncol(Omega) == ncol(Bts))
      Bts %*% Omega %*% t(Bts)
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
#' @param K q×q design kernel in effect space. Supply either a raw matrix or the
#'   list returned by [design_kernel()].
#' @param Omega_list Optional per-subject spatial reliabilities/weights. Each
#'   element may be `NULL` (no extra weighting), a numeric vector of length
#'   `P_s` (diagonal weights for clusters/voxels), or a full `P_s × P_s` matrix
#'   specifying custom covariance. These weights are applied both when
#'   accumulating the compressed covariance and when computing MFA/energy block
#'   normalisation.
#' @inheritParams dkge
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
#' @param cpca_T Optional q×q0 matrix giving the CPCA design basis explicitly.
#'   When supplied it overrides `cpca_blocks`.
#' @param cpca_part Which CPCA-filtered component to estimate: `"none"` (default)
#'   performs the standard DKGE fit; `"design"` uses only the CPCA design part;
#'   `"resid"` uses the residual part; `"both"` fits the design part but also
#'   stores the residual basis for inspection.
#' @param cpca_ridge Optional ridge added to the CPCA-filtered matrices before
#'   eigen-decomposition.
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
                     cpca_ridge = 0) {
  w_method <- match.arg(w_method)
  cpca_part <- match.arg(cpca_part)

  if (inherits(data, "dkge_data")) {
    dataset <- data
    if (!is.null(Omega_list)) {
      stopifnot(length(Omega_list) == dataset$n_subjects)
      dataset$omega <- Omega_list
    }
  } else {
    dataset <- dkge_data(data, designs = designs, omega = Omega_list)
  }

  betas <- dataset$betas
  designs <- dataset$designs
  Omega_list <- dataset$omega
  subject_ids <- dataset$subject_ids
  effects <- dataset$effects
  q <- dataset$q
  S <- dataset$n_subjects

  kernel_info <- NULL
  if (is.list(K) && !is.null(K$K)) {
    kernel_info <- if (!is.null(K$info)) K$info else NULL
    K <- K$K
  }

  stopifnot(is.matrix(K), nrow(K) == q, ncol(K) == q)
  if (is.null(rank)) rank <- min(q, 10L)
  rank <- max(1L, min(rank, q))

  if (is.null(Omega_list)) {
    Omega_list <- vector("list", S)
  }

  ruler <- .dkge_compute_shared_ruler(designs)
  Btil <- .dkge_row_standardize(betas, ruler$R)
  kernels <- .dkge_kernel_roots(K)
  weights <- .dkge_subject_weights(Btil, Omega_list, kernels$Khalf, w_method, w_tau)
  accum <- .dkge_accumulate_chat(Btil, Omega_list, kernels$Khalf, weights)
  Chat_sym <- accum$Chat
  Chat <- Chat_sym

  cpca_info <- NULL
  if (cpca_part != "none") {
    if (!is.null(cpca_T)) {
      T_mat <- as.matrix(cpca_T)
      stopifnot(nrow(T_mat) == q)
    } else {
      stopifnot(!is.null(cpca_blocks), length(cpca_blocks) >= 1)
      T_mat <- diag(1, q)[, unique(cpca_blocks), drop = FALSE]
    }
    split <- dkge_cpca_split_chat(Chat, T_mat, K)
    cpca_info <- list(
      part = cpca_part,
      blocks = cpca_blocks,
      T = T_mat,
      ridge = cpca_ridge,
      P_hat = split$P_hat,
      Chat_design_raw = split$Chat_design,
      Chat_resid_raw = split$Chat_resid
    )
    if (cpca_part %in% c("design", "both")) {
      Chat_design <- split$Chat_design
      if (cpca_ridge > 0) Chat_design <- Chat_design + cpca_ridge * diag(q)
      Chat_design <- (Chat_design + t(Chat_design)) / 2
      cpca_info$Chat_design <- Chat_design
    }
    if (cpca_part %in% c("resid", "both")) {
      Chat_resid <- split$Chat_resid
      if (cpca_ridge > 0) Chat_resid <- Chat_resid + cpca_ridge * diag(q)
      Chat_resid <- (Chat_resid + t(Chat_resid)) / 2
      cpca_info$Chat_resid <- Chat_resid
    }
    if (cpca_part == "design") {
      Chat <- cpca_info$Chat_design
    } else if (cpca_part == "resid") {
      Chat <- cpca_info$Chat_resid
    } else if (cpca_part == "both") {
      Chat <- cpca_info$Chat_design
    }
  }

  if (ridge > 0) Chat <- Chat + ridge * diag(q)

  eigChat <- eigen((Chat + t(Chat))/2, symmetric = TRUE)
  eig_vectors_full <- eigChat$vectors
  eig_values_full <- eigChat$values

  eig_vectors <- eig_vectors_full[, seq_len(rank), drop = FALSE]
  eig_values <- eig_values_full[seq_len(rank)]
  pos_idx <- eig_values > 1e-12
  if (!all(pos_idx)) {
    eig_vectors <- eig_vectors[, pos_idx, drop = FALSE]
    eig_values <- eig_values[pos_idx]
    rank <- length(eig_values)
  }
  sdev <- sqrt(pmax(eig_values, 0))
  U_hat <- eig_vectors
  U <- kernels$Kihalf %*% U_hat

  if (!is.null(cpca_info)) {
    if (cpca_part %in% c("design", "both")) {
      cpca_info$U_design <- U
      cpca_info$evals_design <- eigChat$values
    }
    if (cpca_part == "resid") {
      cpca_info$U_resid <- U
      cpca_info$evals_resid <- eigChat$values
    } else if (cpca_part == "both") {
      eg_resid <- eigen((cpca_info$Chat_resid + t(cpca_info$Chat_resid)) / 2, symmetric = TRUE)
      cpca_info$evals_resid <- eg_resid$values
      cpca_info$U_resid <- kernels$Kihalf %*% eg_resid$vectors[, seq_len(rank), drop = FALSE]
    }
  }

  # Assemble multiblock representation (variables = subject blocks)
  total_clusters <- 0L
  block_indices <- vector("list", S)
  X_blocks <- vector("list", S)
  for (s in seq_len(S)) {
    Bts <- Btil[[s]]
    q <- nrow(Bts)
    P_s <- ncol(Bts)
    idx <- (total_clusters + 1L):(total_clusters + P_s)
    block_indices[[s]] <- idx
    total_clusters <- total_clusters + P_s

    Omega <- Omega_list[[s]]
    if (is.null(Omega)) {
      block <- Bts
    } else if (is.vector(Omega)) {
      block <- Bts * rep(sqrt(Omega), each = q)
    } else {
      block <- Bts %*% sqrtm_sym(Omega)
    }
    block <- kernels$Khalf %*% block
    block <- sqrt(max(weights[s], 0)) * block
    X_blocks[[s]] <- block
  }
  X_concat <- do.call(cbind, X_blocks)
  if (rank > 0) {
    safe_sdev <- ifelse(sdev > 0, sdev, 1)
    V <- t(X_concat) %*% U_hat %*% diag(1 / safe_sdev, nrow = rank)
    zero_cols <- which(sdev <= 0)
    if (length(zero_cols) > 0) {
      V[, zero_cols] <- 0
    }
    scores <- U_hat %*% diag(sdev, nrow = rank)
  } else {
    V <- matrix(0, nrow = total_clusters, ncol = 0)
    scores <- matrix(0, nrow = nrow(Btil[[1]]), ncol = 0)
  }

  preproc_obj <- multivarious::prep(multivarious::pass())
  multivarious_obj <- multivarious::multiblock_biprojector(
    v = V,
    s = scores,
    sdev = sdev,
    preproc = preproc_obj,
    block_indices = block_indices,
    classes = "dkge_core"
  )

  X_store <- if (keep_X) X_concat else NULL

  fit <- list(
    v = V,
    s = scores,
    sdev = sdev,
    U = U,
    evals = eigChat$values,
    R = ruler$R,
    K = K,
    Khalf = kernels$Khalf,
    Kihalf = kernels$Kihalf,
    Chat = Chat,
    contribs = accum$contribs,
    weights = weights,
    Btil = Btil,
    Omega = Omega_list,
    subject_ids = subject_ids,
    effects = effects,
    kernel_info = kernel_info,
    block_indices = block_indices,
    X_concat = X_store,
    V_full = eig_vectors_full,
    eig_vectors_full = eig_vectors_full,
    eig_values_full = eig_values_full,
    rank = rank,
    cpca = cpca_info
  )

  fit$Chat_sym <- Chat_sym
  fit$KU <- fit$K %*% fit$U

  fit$variables <- fit$v
  fit$scores_matrix <- fit$s

  for (nm in names(fit)) {
    multivarious_obj[[nm]] <- fit[[nm]]
  }
  fit <- multivarious_obj
  class(fit) <- unique(c("dkge", class(fit)))

  fit
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
