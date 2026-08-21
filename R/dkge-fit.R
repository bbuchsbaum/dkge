
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

  # Check condition number before Cholesky
  .dkge_check_condition(G_pool, threshold = 1e8, name = "pooled Gram matrix")

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

#' Squared leading singular value of a q x P matrix
#'
#' Computed exactly from the q x q Gram matrix `X X'` (O(q^2 P + q^3)), which
#' is the same cost as a single power-iteration step at DKGE sizes. A
#' randomly-seeded power iteration was used previously; it made the default
#' `w_method = "mfa_sigma1"` subject weights -- and therefore `fit$U` and the
#' sign of every downstream component -- depend on ambient RNG state.
#'
#' @param X Numeric matrix (q x P).
#' @return The squared leading singular value of `X` (0 for empty input).
#' @keywords internal
#' @noRd
.dkge_leading_sv_squared <- function(X) {
  X <- as.matrix(X)
  if (!all(is.finite(X))) {
    X[!is.finite(X)] <- 0
  }
  if (nrow(X) == 0L || ncol(X) == 0L) {
    return(0)
  }
  G <- tcrossprod(X)
  G <- (G + t(G)) / 2
  max(0, eigen(G, symmetric = TRUE, only.values = TRUE)$values[1L])
}

#' Frozen MFA power-iteration approximation
#'
#' The original public default used this 50-step approximation. Its random
#' initial vector was formerly drawn from the caller's RNG stream. The enclosing
#' subject-weight calculation now runs it in one frozen, private RNG scope, which
#' preserves the historical numerical result for the canonical compatibility
#' seed without perturbing or depending on caller state.
#'
#' @param X Numeric matrix.
#' @param tol Relative convergence tolerance.
#' @param max_iter Maximum number of iterations.
#' @return Approximate squared leading singular value.
#' @keywords internal
#' @noRd
.dkge_mfa_leading_sv_squared <- function(X, tol = 1e-6, max_iter = 50L) {
  X <- as.matrix(X)
  if (!all(is.finite(X))) X[!is.finite(X)] <- 0
  n <- nrow(X)
  if (n == 0L || ncol(X) == 0L) return(0)

  v <- stats::rnorm(n)
  v_norm <- sqrt(sum(v * v))
  if (!is.finite(v_norm) || v_norm == 0) return(0)
  v <- v / v_norm
  sigma_sq <- 0
  for (iter in seq_len(max_iter)) {
    w <- X %*% (t(X) %*% v)
    w_norm <- sqrt(sum(w * w))
    if (!is.finite(w_norm) || w_norm == 0) break
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

.dkge_mfa_compatibility_seed <- 8172026L

#' Derive optional subject-level weights
#'
#' Weights are computed on the same quantity the eigensolve sees: the subject's
#' contribution to `Chat` is \eqn{w_s X_s X_s^\top} with
#' \eqn{X_s = K^{1/2} \tilde B_s \Omega_s^{1/2}}. Rows a subject does not
#' observe are zero-filled in *raw* effect space by [dkge_data()] before the
#' pooled ruler `R` is applied, so no row selection happens here: the full
#' `Khalf` congruence is applied, exactly as in `.dkge_transform_effect_moment()`.
#'
#' @param Btil List of row-standardised betas.
#' @param Omega_list Optional per-subject spatial weights.
#' @param Khalf Kernel square root used for energy computations.
#' @param w_method Weighting scheme (`"mfa_sigma1"`, `"energy"`, or `"none"`).
#' @param w_tau Shrinkage parameter toward equal weights.
#' @param obs_masks Optional per-subject observation masks, used only to detect
#'   subjects that observe no effect row at all (they receive weight 0).
#'   Voxel/parcel weights are ignored here: this uses the non-debiased `Btil`
#'   matrices (row-standardised betas) and optional `Omega_list` only.
#' @return Numeric vector of subject weights.
#' @keywords internal
#' @noRd
.dkge_subject_weights <- function(Btil, Omega_list, Khalf, w_method, w_tau, obs_masks = NULL) {
  S <- length(Btil)
  if (w_method == "none") {
    return(rep(1, S))
  }
  if (identical(w_method, "mfa_sigma1")) {
    rng_state <- .dkge_seed_enter(.dkge_mfa_compatibility_seed)
    on.exit(.dkge_seed_exit(rng_state), add = TRUE)
  }
  q <- nrow(Btil[[1]])
  weights <- numeric(S)
  usable <- rep(TRUE, S)
  for (s in seq_len(S)) {
    Bts <- Btil[[s]]
    mask_s <- if (is.null(obs_masks) || length(obs_masks) < s) NULL else obs_masks[[s]]
    if (!length(.dkge_observed_rows(mask_s, q))) {
      usable[s] <- FALSE
      weights[s] <- 0
      next
    }
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
      sigma_sq <- .dkge_mfa_leading_sv_squared(KsBt)
      weights[s] <- 1 / (sigma_sq + 1e-12)
    } else {
      weights[s] <- 1 / (sum(KsBt * KsBt) + 1e-12)
    }
  }
  if (!any(usable)) {
    stop("No subject observes any effect row; cannot derive subject weights.",
         call. = FALSE)
  }
  if (!all(usable)) {
    warning(sprintf(
      "Subject(s) %s observe no effect rows and receive weight 0.",
      paste(which(!usable), collapse = ", ")
    ), call. = FALSE)
  }
  w_norm <- weights
  w_norm[usable] <- weights[usable] / mean(weights[usable])
  (1 - w_tau) * w_norm + w_tau * as.numeric(usable)
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
    stop(sprintf(
      "Observation mask must have one entry per effect row (got %d, expected %d).",
      length(mask), q
    ), call. = FALSE)
  }
  which(mask)
}

#' Reorder provenance masks to the fit subject order
#'
#' Reordering by name is only safe when the stored names are unique; duplicate
#' subject IDs would otherwise make several subjects share one mask, so the
#' stored (positional) order is kept in that case.
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
  mask_names <- names(masks)
  if (!is.null(mask_names) && !anyDuplicated(mask_names) &&
      !is.null(subject_ids) && !anyDuplicated(subject_ids)) {
    order_idx <- match(subject_ids, mask_names)
    if (!anyNA(order_idx)) {
      masks <- masks[order_idx]
    }
  }
  lapply(masks, function(mask) {
    mask <- as.logical(mask)
    if (length(mask) != q) {
      stop(sprintf(
        "Provenance observation mask has length %d but the fit has %d effect rows.",
        length(mask), q
      ), call. = FALSE)
    }
    mask
  })
}

#' Apply a partial-effect missingness policy to compressed covariance
#'
#' `pair_counts` is a *sample-weighted mass*, not an integer tally: it is
#' \eqn{\sum_s a_s 1[\text{pair observed by } s]} for sample weights `a_s`.
#' Those weights are 1 in a batch fit (so the mass is the plain number of
#' contributing subjects) but are fractional under the multiplier bootstrap
#' schemes (`"exp"`, `"bayes"`), which reach this code through
#' `.dkge_repool_fit(fit, sample_weights = ...)`. The policies therefore gate on
#' positive mass rather than on a count of at least one, and `min_pairs` is
#' compared against the mass expressed on the unweighted (subject-count) scale,
#' `pair_counts / weight_scale`.
#'
#' @param Chat qxq raw effect-space pooled moment (not the compressed
#'   covariance).
#' @param pair_counts qxq matrix of sample-weighted observed effect-pair mass.
#' @param missingness One of `"none"`, `"rescale"`, `"mask"`, or `"shrink"`.
#' @param miss_args Policy-specific arguments.
#' @param weight_scale Mean sample weight behind `pair_counts`; `min_pairs` is
#'   applied to `pair_counts / weight_scale` so the threshold keeps its
#'   subject-count meaning under fractional multiplier weights.
#' @param pair_weight Optional qxq summed subject-weight mass per pair. When
#'   supplied, `"rescale"` and `"shrink"` divide by this instead of
#'   `pair_counts` so unequal subject weights do not leak into the per-pair
#'   mean. `pair_counts` still drives the `"mask"` threshold and the shrink
#'   reliability weights.
#' @return Symmetric qxq moment after policy application.
#' @keywords internal
#' @noRd
.dkge_apply_missingness <- function(Chat,
                                    pair_counts = NULL,
                                    missingness = c("none", "rescale", "mask", "shrink"),
                                    miss_args = list(),
                                    weight_scale = 1,
                                    pair_weight = NULL) {
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
  storage.mode(pair_counts) <- "double"
  weight_scale <- suppressWarnings(as.numeric(weight_scale)[1L])
  if (!length(weight_scale) || is.na(weight_scale) || !is.finite(weight_scale) ||
      weight_scale <= 0) {
    weight_scale <- 1
  }
  if (is.null(pair_weight)) {
    pair_weight <- pair_counts
  } else {
    pair_weight <- as.matrix(pair_weight)
    if (!identical(dim(Chat), dim(pair_weight))) {
      stop("`pair_weight` must have the same dimensions as `Chat`.",
           call. = FALSE)
    }
    storage.mode(pair_weight) <- "double"
  }
  covered <- pair_weight > 0

  if (identical(missingness, "rescale")) {
    rescaled <- matrix(0, nrow(Chat), ncol(Chat), dimnames = dimnames(Chat))
    rescaled[covered] <- Chat[covered] / pair_weight[covered]
    Chat <- rescaled
  } else if (identical(missingness, "mask")) {
    threshold <- .dkge_miss_scalar(miss_args$min_pairs, "min_pairs", 1L)
    Chat[(pair_counts / weight_scale) < threshold] <- 0
  } else if (identical(missingness, "shrink")) {
    rescaled <- matrix(0, nrow(Chat), ncol(Chat), dimnames = dimnames(Chat))
    rescaled[covered] <- Chat[covered] / pair_weight[covered]

    max_pc <- max(pair_counts)
    gamma <- .dkge_miss_scalar(miss_args$gamma, "gamma", 1)
    if (max_pc <= 0) {
      weights_mat <- matrix(0, nrow = nrow(Chat), ncol = ncol(Chat))
    } else {
      weights_mat <- (pair_counts / max_pc)^gamma
    }
    if (!is.null(dimnames(Chat))) {
      dimnames(weights_mat) <- dimnames(Chat)
    }

    # Blend the *rescaled* matrix toward its own diagonal: mixing rescaled
    # off-diagonals with un-rescaled diagonal entries would compare quantities
    # on different scales.
    diag_part <- diag(diag(rescaled), nrow = nrow(Chat), ncol = ncol(Chat))
    if (!is.null(dimnames(Chat))) {
      dimnames(diag_part) <- dimnames(Chat)
    }
    Chat <- weights_mat * rescaled + (1 - weights_mat) * diag_part
  }

  (Chat + t(Chat)) / 2
}

#' Accumulate compressed covariance in the K-metric
#'
#' Used by fold re-pooling ([.dkge_fold_weight_context()]) to rebuild `Chat`
#' from a training subset of `Btil` without going through the full fit
#' pipeline. Subject contributions are `K^{1/2} B_s B_s' K^{1/2}` (or the
#' Omega-weighted analogue), then pooled with the supplied subject weights.
#'
#' @param Btil List of row-standardised betas.
#' @param Omega_list Optional spatial weights aligned with `Btil`.
#' @param Khalf Kernel square root used to project into the K-metric.
#' @param weights Subject weights applied during accumulation.
#' @param voxel_weights Optional per-subject or shared voxel/parcel weights.
#' @return List with the symmetrised compressed covariance (`Chat`) and
#'   per-subject contributions.
#' @keywords internal
#' @noRd
.dkge_accumulate_chat <- function(Btil, Omega_list, Khalf, weights,
                                  voxel_weights = NULL) {
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
    S_s <- (S_s + t(S_s)) / 2
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
#'   normalization.
#' @param w_method Subject-level weighting scheme.
#'   * `"mfa_sigma1"` (default): inverse squared leading singular value of
#'     \eqn{K^{1/2} Btil_s \Omega_s^{1/2}} (Multiple Factor Analysis scaling).
#'   * `"energy"`: inverse Frobenius norm squared of the same block.
#'   * `"none"`: disable block scaling (all weights = 1).
#' @param w_tau Shrinkage parameter (0..1) toward equal weights. 0 keeps the raw
#'   weights, 1 forces equal weighting.
#' @param ridge Ridge regularization parameter (default 0).
#' @param rank Desired rank for the decomposition. If NULL, uses full rank.
#' @param keep_X Logical; when `TRUE`, store the concatenated training matrix
#'   used to build the multiblock projection (can be large). Only available
#'   when the fit is an exact `block_biprojector`; q-space moments have no
#'   matching subject-by-voxel factor.
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
#' @param solver Solver used for the q-space problem. `"pooled"` keeps the
#'   original eigen-decomposition; `"jd"` performs joint diagonalisation using
#'   [dkge_jd_control()] settings.
#' @param jd_control Control parameters for the JD solver.
#' @param jd_mask Optional mask (matrix or list of matrices) applied to the
#'   off-diagonal penalty during JD. When supplied as a single matrix it is
#'   recycled across subjects.
#' @param jd_init Optional orthogonal initialiser for the JD solver expressed in
#'   the whitened \eqn{K^{1/2}} metric (qxq matrix).
#' @param effect_scaling Effect-space scaling applied before accumulating the
#'   compressed covariance. `"pooled_design"` (default) multiplies each subject
#'   beta/effect matrix by the pooled design Cholesky factor, preserving the
#'   original GLM-beta behavior. `"none"` leaves rows on their input scale and
#'   uses an identity contrast ruler, which is appropriate when rows already have
#'   a common absolute interpretation such as AUC minus chance.
#' @param effect_weights Optional [dkge_effect_weights()] specification. Count,
#'   explicit precision, or DerSimonian-Laird `"random_effects"` weighting is
#'   applied to raw q-space effect moments before the pooled ruler and design
#'   kernel mix effect rows. Random-effects fits are q-space-only.
#' @param debias Finite-trial noise treatment. `"none"` preserves the observed
#'   second moment. `"analytic"` subtracts
#'   `noise_trace * effect_noise_cov` for each subject before pooling; these
#'   sufficient statistics are produced by [dkge_trial_subject()].
#'   `"split_half"` uses the symmetrized cross-product of two independently
#'   estimated trial halves and requires `dkge_trial_subject(split = ...)`.
#' @param missingness Partial-effect coverage policy used when subjects do not
#'   all observe the same effect rows. `"none"` keeps the raw zero-filled union
#'   accumulation, `"rescale"` divides entries by observed pair counts,
#'   `"mask"` zeros under-covered entries, and `"shrink"` blends rescaled entries
#'   toward the diagonal according to coverage. These policies operate on raw
#'   q-space effect moments before the pooled ruler and design kernel mix rows.
#'   When `effect_weights` is active, `"none"`/`"mask"` scale the
#'   precision-weighted per-pair mean by observed subject mass
#'   \eqn{\sum_s a_s w_s 1[\mathrm{obs}_s]} (equal to the full cohort mass
#'   only under complete coverage), `"rescale"` returns that per-pair mean,
#'   and `"mask"`/`"shrink"` threshold on the effective pair count `$pair_ess`
#'   instead of the raw `$pair_counts`. Under unit precision, equal subject
#'   weights, **and full coverage** these coincide with the unweighted
#'   policies; under partial coverage they keep the effect-weighted mean on
#'   the observed-mass scale rather than mean-imputing absent subjects.
#' @param miss_args List of arguments for `missingness`. Known fields are
#'   `min_pairs` (used by `"mask"` and `"shrink"`) and `gamma` (used by
#'   `"shrink"`). Unknown names are an error.
#' @return A fitted `dkge` object. Exact unregularized pooled moments
#'   additionally inherit from `multiblock_biprojector` and advertise physical
#'   block loadings `$v`. Pair-normalized, missingness-transformed, debiased,
#'   ridged, CPCA, and JD fits inherit from `dkge_qspace` instead, set `$v`
#'   and `$X_concat` to `NULL`, and expose only the coherent q-space
#'   eigensystem. `$representation` records which contract applies.
#'   The object also records `$effect_moment`, `$pair_weight`, `$pair_ess`, and
#'   `$moment_diagnostics` for auditing reliability weighting and any negative
#'   spectral mass introduced by debiasing.
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
                     keep_X = FALSE,
                     cpca_blocks = NULL,
                     cpca_T = NULL,
                     cpca_part = c("none", "design", "resid", "both"),
                     cpca_ridge = 0,
                     weights = NULL,
                     solver = c("pooled", "jd"),
                     jd_control = dkge_jd_control(),
                     jd_mask = NULL,
                     jd_init = NULL,
                     effect_scaling = c("pooled_design", "none"),
                     effect_weights = NULL,
                     debias = c("none", "analytic", "split_half"),
                     missingness = c("none", "rescale", "mask", "shrink"),
                     miss_args = list()) {
  w_method <- match.arg(w_method)
  effect_scaling <- match.arg(effect_scaling)
  cpca_part <- match.arg(cpca_part)
  solver <- match.arg(solver)
  missingness <- match.arg(missingness)
  debias <- match.arg(debias)
  miss_args <- .dkge_validate_miss_args(missingness, miss_args)

  effect_weight_spec <- if (is.null(effect_weights)) {
    dkge_effect_weights("none")
  } else {
    effect_weights
  }
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

  .dkge_fit_assemble(prepped,
                     accum,
                     solved,
                     keep_X = keep_X,
                     w_method = w_method,
                     w_tau = w_tau,
                     ridge = ridge,
                     missingness = missingness,
                     miss_args = miss_args)
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
