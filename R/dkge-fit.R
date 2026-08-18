
# dkge-fit.R
# Core DKGE batch fit, LOSO contrasts, and user-facing helpers for multiple subjects.


# -------------------------------------------------------------------------
# Internal helpers used by the fitter ------------------------------------
# -------------------------------------------------------------------------

#' Compute pooled Gram matrix and Cholesky ruler
#'
#' @param X_list List of subject design matrices (`T_sxq`).
#' @param information_list Optional list of q-by-q effect-information matrices.
#'   A non-NULL entry replaces `crossprod(X_s)` for that subject, allowing GLS
#'   trial models to contribute `X_s' W_s X_s` to the shared ruler.
#' @param jitter Numerical jitter added to the diagonal for stability.
#' @return List with the pooled Gram matrix and its upper-triangular Cholesky factor.
#' @keywords internal
#' @noRd
.dkge_compute_shared_ruler <- function(X_list, information_list = NULL,
                                       jitter = 1e-10) {
  q <- ncol(X_list[[1]])
  G_pool <- matrix(0, q, q)
  if (is.null(information_list)) {
    information_list <- vector("list", length(X_list))
  }
  if (length(information_list) != length(X_list)) {
    stop("Effect-information list must match the number of subject designs.",
         call. = FALSE)
  }
  source <- character(length(X_list))
  for (s in seq_along(X_list)) {
    Xs <- X_list[[s]]
    info <- information_list[[s]]
    if (is.null(info)) {
      info <- crossprod(Xs)
      source[[s]] <- "design_crossproduct"
    } else {
      info <- as.matrix(info)
      if (!identical(dim(info), c(q, q)) || any(!is.finite(info))) {
        stop("Each effect-information matrix must be finite and q x q.",
             call. = FALSE)
      }
      asym <- max(abs(info - t(info)))
      if (asym > 1e-8 * max(1, max(abs(info)))) {
        stop("Effect-information matrices must be symmetric.", call. = FALSE)
      }
      info <- (info + t(info)) / 2
      source[[s]] <- "effect_information"
    }
    G_pool <- G_pool + info
  }
  diag(G_pool) <- diag(G_pool) + jitter

  # Check condition number before Cholesky
  .dkge_check_condition(G_pool, threshold = 1e8, name = "pooled Gram matrix")

  list(
    R = chol(G_pool),
    G_pool = G_pool,
    source = if (length(unique(source)) == 1L) unique(source) else "mixed"
  )
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
  Khalf <- V %*% (diag(sqrt(vals), length(vals))) %*% t(V)
  Kihalf <- V %*% (diag(1 / sqrt(vals), length(vals))) %*% t(V)
  list(
    Khalf = Khalf,
    Kihalf = Kihalf,
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
  # Deterministic data-driven start (row sums), so subject weights are
  # reproducible and do not perturb the caller's RNG stream. Power iteration
  # converges to the same leading value for any start not orthogonal to the top
  # singular vector; fall back to a constant vector if the row sums vanish.
  v <- rowSums(X)
  v_norm <- sqrt(sum(v * v))
  if (!is.finite(v_norm) || v_norm == 0) {
    v <- rep(1, n)
    v_norm <- sqrt(n)
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
#' @param Btil List of beta matrices in the raw effect coordinate. Historical
#'   callers may pass already-standardised matrices by leaving `R = NULL`.
#' @param Omega_list Optional per-subject spatial weights.
#' @param Khalf Kernel square root used for energy computations.
#' @param w_method Weighting scheme (`"mfa_sigma1"`, `"energy"`, or `"none"`).
#' @param w_tau Shrinkage parameter toward equal weights.
#' @return Numeric vector of subject weights.
#' @keywords internal
#' @noRd
.dkge_subject_weights <- function(Btil, Omega_list, Khalf, w_method, w_tau,
                                  obs_masks = NULL, R = NULL,
                                  voxel_weights = NULL) {
  S <- length(Btil)
  if (w_method == "none") {
    return(rep(1, S))
  }
  # Compute weights on the exact physical block used by the eigensolve. Masks
  # act on raw effects before R or K can mix coordinates.
  q <- nrow(Btil[[1]])
  if (is.null(R)) R <- diag(q)
  weights <- numeric(S)
  for (s in seq_len(S)) {
    Bts <- Btil[[s]]
    mask_s <- if (is.null(obs_masks) || length(obs_masks) < s) NULL else obs_masks[[s]]
    obs <- .dkge_observed_rows(mask_s, q)
    if (!length(obs)) {
      weights[s] <- 1 / 1e-12
      next
    }
    Omega <- Omega_list[[s]]
    vw <- if (is.list(voxel_weights)) voxel_weights[[s]] else voxel_weights
    Bwork <- .dkge_right_weighted_beta(Bts, Omega, vw, mask_s)
    KsBt <- Khalf %*% t(R) %*% Bwork
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

#' Resolve observed effect rows for a subject
#'
#' @param mask Optional logical observation mask.
#' @param q Number of global effect rows.
#' @return Integer indices of observed rows.
#' @keywords internal
#' @noRd
.dkge_observed_rows <- function(mask, q) {
  if (is.null(mask)) {
    return(seq_len(q))
  }
  mask <- as.logical(mask)
  if (length(mask) != q) {
    return(seq_len(q))
  }
  which(mask)
}

#' Reorder provenance masks to the fit subject order
#'
#' @param provenance Fit/data provenance.
#' @param subject_ids Subject IDs in fit order.
#' @param q Number of effect rows.
#' @return List of logical masks or `NULL` when unavailable.
#' @keywords internal
#' @noRd
.dkge_obs_masks_from_provenance <- function(provenance, subject_ids, q) {
  masks <- provenance$obs_mask %||% NULL
  if (is.null(masks) || !length(masks)) {
    return(NULL)
  }
  if (!is.null(names(masks)) && !is.null(subject_ids)) {
    order_idx <- match(subject_ids, names(masks))
    if (!anyNA(order_idx)) {
      masks <- masks[order_idx]
    }
  }
  lapply(masks, function(mask) {
    mask <- as.logical(mask)
    if (length(mask) != q) {
      rep(TRUE, q)
    } else {
      mask
    }
  })
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
.dkge_accumulate_chat <- function(Btil, Omega_list, Khalf, weights,
                                  voxel_weights = NULL,
                                  obs_masks = NULL) {
  S <- length(Btil)
  q <- nrow(Btil[[1]])
  Chat <- matrix(0, q, q)
  contribs <- vector("list", S)

  for (s in seq_len(S)) {
    Bts <- Btil[[s]]
    mask_s <- if (is.null(obs_masks) || length(obs_masks) < s) NULL else obs_masks[[s]]
    obs <- .dkge_observed_rows(mask_s, q)
    if (!length(obs)) {
      contribs[[s]] <- matrix(0, q, q)
      next
    }
    w_s <- if (is.list(voxel_weights)) voxel_weights[[s]] else voxel_weights
    Omega <- Omega_list[[s]]
    Bw <- .dkge_right_weighted_beta(Bts, Omega, w_s, mask_s)
    block <- Khalf %*% Bw
    S_s <- tcrossprod(block)
    S_s <- (S_s + t(S_s)) / 2
    contribs[[s]] <- S_s
    Chat <- Chat + weights[s] * S_s
  }
  list(Chat = (Chat + t(Chat)) / 2, contribs = contribs)
}

#' Reconstruct pairwise effect coverage from subject observation masks
#'
#' @param provenance Data provenance carrying `obs_mask`.
#' @param subject_ids Subject order used by the fit.
#' @param indices Integer subject indices to include.
#' @param q Number of effect rows.
#' @param dimnames Optional dimnames for the returned matrix.
#' @return Integer qxq pair-count matrix, or `NULL` when masks are unavailable.
#' @keywords internal
#' @noRd
.dkge_pair_counts_from_provenance <- function(provenance,
                                              subject_ids,
                                              indices,
                                              q,
                                              dimnames = NULL) {
  masks <- .dkge_obs_masks_from_provenance(provenance, subject_ids, q)
  if (is.null(masks) || !length(masks) || !length(indices)) {
    return(NULL)
  }

  if (length(masks) < max(indices)) {
    return(NULL)
  }

  effect_ids <- provenance$effect_ids %||% NULL
  if (!is.null(effect_ids) && length(effect_ids) == q) {
    dimnames <- list(effect_ids, effect_ids)
  }

  pair_counts <- matrix(0L, nrow = q, ncol = q)
  if (!is.null(dimnames)) {
    dimnames(pair_counts) <- dimnames
  }

  for (idx in indices) {
    mask <- masks[[idx]]
    if (is.null(mask)) next
    mask <- as.logical(mask)
    if (length(mask) != q) next
    sel <- which(mask)
    if (!length(sel)) next
    pair_counts[sel, sel] <- pair_counts[sel, sel] + 1L
  }

  pair_counts
}

#' Apply a partial-effect missingness policy to compressed covariance
#'
#' @param Chat qxq compressed covariance matrix.
#' @param pair_counts qxq integer matrix of observed effect-pair counts.
#' @param missingness One of `"none"`, `"rescale"`, `"mask"`, or `"shrink"`.
#' @param miss_args Policy-specific arguments.
#' @return Symmetric qxq covariance matrix after policy application.
#' @keywords internal
#' @noRd
.dkge_apply_missingness <- function(Chat,
                                    pair_counts = NULL,
                                    missingness = c("none", "rescale", "mask", "shrink"),
                                    miss_args = list()) {
  missingness <- match.arg(missingness)
  Chat <- as.matrix(Chat)
  if (is.null(pair_counts) || identical(missingness, "none")) {
    return((Chat + t(Chat)) / 2)
  }
  if (is.null(dimnames(Chat)) && !is.null(dimnames(pair_counts))) {
    dimnames(Chat) <- dimnames(pair_counts)
  }

  pair_counts <- as.matrix(pair_counts)
  if (!identical(dim(Chat), dim(pair_counts))) {
    stop("`pair_counts` must have the same dimensions as `Chat`.", call. = FALSE)
  }

  if (identical(missingness, "rescale")) {
    pc_safe <- pmax(pair_counts, 1L)
    Chat <- Chat / pc_safe
    Chat[pair_counts == 0L] <- 0
  } else if (identical(missingness, "mask")) {
    threshold <- miss_args$min_pairs %||% 1L
    if (!is.numeric(threshold) || length(threshold) != 1L || threshold < 0) {
      stop("`miss_args$min_pairs` must be a non-negative numeric scalar.", call. = FALSE)
    }
    Chat[pair_counts < threshold] <- 0
  } else if (identical(missingness, "shrink")) {
    pc_safe <- pmax(pair_counts, 1L)
    rescaled <- Chat / pc_safe
    rescaled[pair_counts == 0L] <- 0

    max_pc <- max(pair_counts)
    gamma <- miss_args$gamma %||% 1
    if (!is.numeric(gamma) || length(gamma) != 1L || gamma < 0) {
      stop("`miss_args$gamma` must be a non-negative numeric scalar.", call. = FALSE)
    }
    if (max_pc <= 0) {
      weights_mat <- matrix(0, nrow = nrow(Chat), ncol = ncol(Chat))
    } else {
      weights_mat <- (pair_counts / max_pc)^gamma
    }
    if (!is.null(dimnames(Chat))) {
      dimnames(weights_mat) <- dimnames(Chat)
    }

    diag_part <- diag(diag(Chat), nrow = nrow(Chat), ncol = ncol(Chat))
    if (!is.null(dimnames(Chat))) {
      dimnames(diag_part) <- dimnames(Chat)
    }
    Chat <- weights_mat * rescaled + (1 - weights_mat) * diag_part
  }

  (Chat + t(Chat)) / 2
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
#'     \eqn{K^{1/2} Btil_s \Omega_s^{1/2}} (Multiple Factor Analysis scaling).
#'   * `"energy"`: inverse Frobenius norm squared of the same block.
#'   * `"none"`: disable block scaling (all weights = 1).
#' @param w_tau Shrinkage parameter (0..1) toward equal weights. 0 keeps the raw
#'   weights, 1 forces equal weighting.
#' @param ridge Ridge regularization parameter (default 0).
#' @param rank Desired rank for the decomposition. If NULL, uses full rank.
#' @param effect_scaling Effect-space scaling applied before accumulating the
#'   compressed covariance. `"pooled_design"` (default) multiplies each subject
#'   beta/effect matrix by the pooled design Cholesky factor, preserving the
#'   original GLM-beta behaviour. `"none"` leaves rows on their input scale and
#'   uses an identity contrast ruler, which is appropriate when rows already have
#'   a common absolute interpretation such as AUC minus chance.
#' @param keep_X Logical; when `TRUE`, store the concatenated training matrix
#'   used to build the exact multiblock projection (can be large). This is only
#'   available when the fitted moment has a physical block factor; pairwise
#'   reliability normalization, debiasing, ridge, CPCA, and JD fits fail closed.
#' @param force_qspace Logical; force the fitted object to use the q-space-only
#'   representation even when an exact physical block factor exists. This is
#'   used by resampling refits, which need only the effect-space eigensystem.
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
#' @param effect_weights Optional [dkge_effect_weights()] specification. Count
#'   explicit precision, or random-effects weighting is applied to raw q-space
#'   effect moments before the pooled ruler and design kernel mix effect rows.
#' @param debias Finite-trial noise treatment. `"none"` preserves the observed
#'   second moment. `"analytic"` subtracts
#'   `noise_trace * effect_noise_cov` for each subject before pooling; these
#'   sufficient statistics are produced by [dkge_trial_subject()].
#'   `"split_half"` uses the symmetrized cross-product of two independently
#'   estimated trial halves and requires `dkge_trial_subject(split = ...)`.
#'   The fitted `$moment_estimator` distinguishes IID analytic,
#'   covariance-aware analytic, independent split, and exploratory split
#'   estimators; `debias = "analytic"` by itself does not imply IID errors.
#' @param solver Solver used for the q-space problem. `"pooled"` keeps the
#'   original eigen-decomposition; `"jd"` performs joint diagonalisation using
#'   [dkge_jd_control()] settings.
#' @param missingness Partial-effect coverage policy used when subjects do not
#'   all observe the same effect rows. `"none"` keeps the raw zero-filled union
#'   accumulation, `"rescale"` divides entries by observed pair counts,
#'   `"mask"` zeros under-covered entries, and `"shrink"` blends rescaled entries
#'   toward the diagonal according to coverage. These policies operate on raw
#'   q-space effect moments before the pooled ruler and design kernel mix rows.
#' @param miss_args List of arguments for `missingness`; currently
#'   `min_pairs` for `"mask"` and `gamma` for `"shrink"`.
#' @param jd_control Control parameters for the JD solver.
#' @param jd_mask Optional mask (matrix or list of matrices) applied to the
#'   off-diagonal penalty during JD. When supplied as a single matrix it is
#'   recycled across subjects.
#' @param jd_init Optional orthogonal initialiser for the JD solver expressed in
#'   the whitened \eqn{K^{1/2}} metric (qxq matrix).
#' @return A fitted `dkge` object. Exact unregularized pooled block moments
#'   additionally inherit from `multiblock_biprojector` and advertise physical
#'   block loadings in `$v`. Pair-normalized, debiased, ridge, CPCA, and JD fits
#'   instead inherit from `dkge_qspace`, set `$v` and `$X_concat` to `NULL`, and
#'   expose only the coherent q-space eigensystem. `$representation` records the
#'   active contract.
#'   The object also records `$effect_moment`, `$pair_weight`, `$pair_ess`,
#'   `$effect_precision_diagnostics`, and `$moment_diagnostics` for auditing
#'   reliability weighting, subject dominance, and any negative spectral mass
#'   introduced by debiasing. [summary.dkge()] presents these diagnostics in a
#'   compact report.
#' @examples
#' # Simulate toy data with known structure
#' toy <- dkge_sim_toy(
#'   factors = list(A = list(L = 2), B = list(L = 3)),
#'   active_terms = c("A", "B"), S = 4, P = 30, snr = 5
#' )
#'
#' # Bundle into dkge_data and fit
#' data <- dkge_data(toy$B_list, toy$X_list)
#' fit <- dkge_fit(data, K = toy$K, rank = 2)
#' fit$rank
#' @export
dkge_fit <- function(data, designs = NULL, K = NULL, Omega_list = NULL,
                     w_method = c("mfa_sigma1", "energy", "none"),
                     w_tau = 0.3,
                     ridge = 0,
                     rank = NULL,
                     effect_scaling = c("pooled_design", "none"),
                     keep_X = FALSE,
                     force_qspace = FALSE,
                     cpca_blocks = NULL,
                     cpca_T = NULL,
                     cpca_part = c("none", "design", "resid", "both"),
                     cpca_ridge = 0,
                     weights = NULL,
                     effect_weights = NULL,
                     debias = c("none", "analytic", "split_half"),
                     solver = c("pooled", "jd"),
                     missingness = c("none", "rescale", "mask", "shrink"),
                     miss_args = list(),
                     jd_control = dkge_jd_control(),
                     jd_mask = NULL,
                     jd_init = NULL) {
  w_method <- match.arg(w_method)
  effect_scaling <- match.arg(effect_scaling)
  cpca_part <- match.arg(cpca_part)
  solver <- match.arg(solver)
  missingness <- match.arg(missingness)
  debias <- match.arg(debias)
  if (!is.logical(force_qspace) || length(force_qspace) != 1L ||
      is.na(force_qspace)) {
    stop("`force_qspace` must be TRUE or FALSE.", call. = FALSE)
  }

  effect_weight_spec <- effect_weights %||% dkge_effect_weights("none")
  if (identical(solver, "jd") &&
      (!identical(effect_weight_spec$method, "none") ||
       !identical(missingness, "none") ||
       !identical(debias, "none"))) {
    stop("The JD solver does not yet support effect weighting, missingness normalization, or debiasing; use solver='pooled'.",
         call. = FALSE)
  }

  prepped <- .dkge_fit_prepare(data,
                               designs = designs,
                               K = K,
                               Omega_list = Omega_list,
                               weights = weights,
                               effect_weights = effect_weight_spec,
                               rank = rank,
                               effect_scaling = effect_scaling)

  accum <- .dkge_fit_accumulate(prepped,
                                w_method = w_method,
                                w_tau = w_tau,
                                missingness = missingness,
                                miss_args = miss_args,
                                debias = debias)

  solved <- .dkge_fit_solve(prepped,
                            accum,
                            rank = prepped$rank,
                            cpca_part = cpca_part,
                            cpca_blocks = cpca_blocks,
                            cpca_T = cpca_T,
                            cpca_ridge = cpca_ridge,
                            ridge = ridge,
                            solver = solver,
                            jd_control = jd_control,
                            jd_mask = jd_mask,
                            jd_init = jd_init)

  fit <- .dkge_fit_assemble(prepped,
                            accum,
                            solved,
                            keep_X = keep_X,
                            w_method = w_method,
                            w_tau = w_tau,
                            ridge = ridge,
                            force_qspace = force_qspace,
                            missingness = missingness,
                            miss_args = miss_args)
  fit$refit_spec <- list(
    w_method = w_method,
    w_tau = w_tau,
    ridge = ridge,
    rank = rank,
    effect_scaling = effect_scaling,
    force_qspace = force_qspace,
    cpca_blocks = cpca_blocks,
    cpca_T = cpca_T,
    cpca_part = cpca_part,
    cpca_ridge = cpca_ridge,
    weights = prepped$weight_spec,
    effect_weights = effect_weight_spec,
    debias = debias,
    solver = solver,
    missingness = missingness,
    miss_args = miss_args,
    jd_control = jd_control,
    jd_mask = jd_mask,
    jd_init = jd_init
  )
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
