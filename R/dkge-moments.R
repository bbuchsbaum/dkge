# dkge-moments.R
# Central q-space moment construction, debiasing, and reliability pooling.

#' Apply voxel/parcel weights to the columns of an effect matrix
#'
#' Strict: the caller is responsible for handing over weights that already match
#' this subject's column count (see `.dkge_subject_voxel_weights()`). Recycling
#' here would apply one subject's spatial profile to another's parcels.
#'
#' @noRd
.dkge_scale_effect_columns <- function(B, weights) {
  if (is.null(weights) || !length(weights)) return(B)
  w <- as.numeric(weights)
  if (length(w) != ncol(B)) {
    stop(sprintf(
      "Voxel weights must have one entry per beta column (got %d, expected %d).",
      length(w), ncol(B)
    ), call. = FALSE)
  }
  sweep(B, 2L, sqrt(pmax(w, 0)), "*")
}

#' Select the voxel weights that apply to one subject
#'
#' Weight specifications resolve either to one vector per subject or to a single
#' cohort-level vector sized to whichever subject was measured first. Subjects
#' may have different numbers of clusters (they are reconciled later by
#' transport), so a cohort vector only applies to a subject when the lengths
#' agree. When they do not:
#'
#' * a constant vector carries no per-column information, so it is re-expanded
#'   to the subject's width (identical to no weighting for the usual all-ones
#'   case, and it does not silently drop a uniform scale);
#' * a vector that varies across columns is rejected, because stretching one
#'   subject's spatial profile onto another's parcels is silently wrong.
#'
#' @noRd
.dkge_active_voxel_weights <- function(voxel_weights) {
  if (is.null(voxel_weights) || !length(voxel_weights)) return(FALSE)
  w <- if (is.list(voxel_weights)) {
    unlist(voxel_weights, use.names = FALSE)
  } else {
    as.numeric(voxel_weights)
  }
  if (!length(w)) return(FALSE)
  any(abs(w - 1) > 1e-12)
}

.dkge_resolve_subject <- function(fit, subject) {
  ids <- fit$subject_ids
  S <- length(fit$Btil %||% ids)
  if (is.character(subject) && length(subject) == 1L) {
    idx <- match(subject, ids)
    if (is.na(idx)) {
      stop(sprintf("Unknown subject id '%s'.", subject), call. = FALSE)
    }
    return(as.integer(idx))
  }
  if (is.numeric(subject) && length(subject) == 1L &&
      is.finite(subject) && subject == as.integer(subject)) {
    idx <- as.integer(subject)
    if (idx < 1L || idx > S) {
      stop(sprintf("`subject` index %d is out of range [1, %d].", idx, S),
           call. = FALSE)
    }
    return(idx)
  }
  stop("`subject` must be a single training-subject index or id.",
       call. = FALSE)
}

.dkge_subject_voxel_weights <- function(voxel_weights, s, n_cols,
                                        subject_id = NULL) {
  w <- if (is.list(voxel_weights)) voxel_weights[[s]] else voxel_weights
  if (is.null(w) || !length(w)) return(NULL)
  w <- as.numeric(w)
  if (length(w) == n_cols) return(w)
  if (isTRUE(all(w == w[[1]]))) return(rep(w[[1]], n_cols))
  stop(sprintf(
    paste0("Voxel weights of length %d do not apply to subject %s, which has ",
           "%d columns. Supply per-subject weights when cluster counts differ."),
    length(w), subject_id %||% paste0("#", s), n_cols
  ), call. = FALSE)
}

#' Zero q x q effect moment carrying the effect labels of `B`
#'
#' @noRd
.dkge_empty_effect_moment <- function(B) {
  q <- nrow(B)
  nms <- rownames(B)
  matrix(0, q, q,
         dimnames = if (is.null(nms)) NULL else list(nms, nms))
}

#' Validate a scalar entry of `miss_args`
#'
#' Shared by the pre-transform (`.dkge_pool_effect_moments`) and
#' raw-effect-space (`.dkge_apply_missingness`) policy implementations so
#' that both reject the same inputs.
#'
#' @noRd
.dkge_validate_miss_args <- function(policy, miss_args) {
  policy <- match.arg(policy, c("none", "rescale", "mask", "shrink"))
  if (is.null(miss_args)) miss_args <- list()
  if (!is.list(miss_args)) {
    stop("`miss_args` must be a list.", call. = FALSE)
  }
  # Names used by any policy. Extra known fields are ignored (a shared
  # `list(min_pairs=..., gamma=...)` is valid for every policy); typos
  # such as `gama` are not.
  known <- c("min_pairs", "gamma")
  nms <- names(miss_args)
  if (is.null(nms)) nms <- character(0)
  unknown <- setdiff(nms[nzchar(nms)], known)
  if (length(unknown)) {
    stop(sprintf(
      "`miss_args` has unknown field(s): %s. Known: %s.",
      paste(unknown, collapse = ", "),
      paste(known, collapse = ", ")
    ), call. = FALSE)
  }
  miss_args
}

.dkge_miss_scalar <- function(value, name, default) {
  value <- value %||% default
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 0) {
    stop(sprintf("`miss_args$%s` must be a non-negative numeric scalar.", name),
         call. = FALSE)
  }
  as.numeric(value)
}

#' @noRd
.dkge_effect_moment <- function(B, Omega = NULL, voxel_weights = NULL,
                                obs_mask = NULL) {
  B <- as.matrix(B)
  q <- nrow(B)
  obs <- .dkge_observed_rows(obs_mask, q)
  if (!length(obs)) return(.dkge_empty_effect_moment(B))
  Bwork <- .dkge_scale_effect_columns(B, voxel_weights)
  if (length(obs) < q) Bwork[-obs, ] <- 0

  M <- if (is.null(Omega)) {
    tcrossprod(Bwork)
  } else if (is.vector(Omega)) {
    if (length(Omega) != ncol(Bwork)) {
      stop("Diagonal Omega must have one entry per beta column.", call. = FALSE)
    }
    tcrossprod(sweep(Bwork, 2L, sqrt(pmax(as.numeric(Omega), 0)), "*"))
  } else {
    Omega <- as.matrix(Omega)
    if (!identical(dim(Omega), c(ncol(Bwork), ncol(Bwork)))) {
      stop("Full Omega must be square with one row per beta column.",
           call. = FALSE)
    }
    Bwork %*% Omega %*% t(Bwork)
  }
  M <- (M + t(M)) / 2
  if (!is.null(rownames(B))) dimnames(M) <- list(rownames(B), rownames(B))
  M
}

#' @noRd
.dkge_split_effect_moment <- function(split_betas, Omega = NULL,
                                      voxel_weights = NULL,
                                      obs_mask = NULL,
                                      subject_id = NULL) {
  if (is.null(split_betas) || !is.list(split_betas) || length(split_betas) != 2L) {
    who <- if (is.null(subject_id)) "every subject" else
      sprintf("subject %s", subject_id)
    stop(sprintf(
      "Split-half debiasing requires two stored split beta matrices for %s.",
      who
    ), call. = FALSE)
  }
  B1 <- .dkge_scale_effect_columns(as.matrix(split_betas[[1]]), voxel_weights)
  B2 <- .dkge_scale_effect_columns(as.matrix(split_betas[[2]]), voxel_weights)
  if (!identical(dim(B1), dim(B2))) {
    stop("Stored split beta matrices must have identical dimensions.", call. = FALSE)
  }
  q <- nrow(B1)
  obs <- .dkge_observed_rows(obs_mask, q)
  if (!length(obs)) return(.dkge_empty_effect_moment(B1))
  if (length(obs) < q) {
    B1[-obs, ] <- 0
    B2[-obs, ] <- 0
  }
  cross <- if (is.null(Omega)) {
    B1 %*% t(B2)
  } else if (is.vector(Omega)) {
    if (length(Omega) != ncol(B1)) {
      stop("Diagonal Omega must have one entry per beta column.", call. = FALSE)
    }
    w <- sqrt(pmax(as.numeric(Omega), 0))
    sweep(B1, 2L, w, "*") %*% t(sweep(B2, 2L, w, "*"))
  } else {
    Omega <- as.matrix(Omega)
    if (!identical(dim(Omega), c(ncol(B1), ncol(B1)))) {
      stop("Full Omega must be square with one row per beta column.",
           call. = FALSE)
    }
    B1 %*% Omega %*% t(B2)
  }
  out <- (cross + t(cross)) / 2
  if (!is.null(rownames(B1))) dimnames(out) <- list(rownames(B1), rownames(B1))
  out
}

#' @noRd
.dkge_noise_trace <- function(subject, Omega = NULL, voxel_weights = NULL) {
  # A precomputed trace is only safe when it already includes the same
  # spatial/voxel weights the caller is asking for, or when the caller
  # supplied an explicit override and no extra spatial weights are applied.
  # Auto-computed unweighted caches are ignored so voxel/Omega weights can
  # recompute `sum(vw * omega * sigma2)`.
  cached <- subject$noise_trace
  if (!is.null(cached) &&
      is.null(Omega) &&
      !.dkge_active_voxel_weights(voxel_weights) &&
      (identical(subject$noise_trace_scope, "weighted") ||
       identical(subject$noise_trace_scope, "supplied"))) {
    return(as.numeric(cached))
  }
  sigma2 <- subject$residual_variance
  if (is.null(sigma2)) {
    stop(sprintf(
      "Subject '%s' lacks residual variance information for analytic debiasing.",
      subject$id %||% "(unnamed)"
    ), call. = FALSE)
  }
  sigma2 <- as.numeric(sigma2)
  if (anyNA(sigma2) || any(!is.finite(sigma2)) || any(sigma2 < 0)) {
    stop(sprintf(
      "Subject '%s' has unavailable residual variances; analytic debiasing requires positive residual degrees of freedom.",
      subject$id %||% "(unnamed)"
    ), call. = FALSE)
  }
  P <- length(sigma2)
  vw <- if (is.null(voxel_weights) || !length(voxel_weights)) rep(1, P) else {
    if (length(voxel_weights) != P) {
      stop("Voxel weights must match residual variances.", call. = FALSE)
    }
    as.numeric(voxel_weights)
  }
  omega_diag <- if (is.null(Omega)) {
    rep(1, P)
  } else if (is.vector(Omega)) {
    if (length(Omega) != P) stop("Omega must match residual variances.", call. = FALSE)
    as.numeric(Omega)
  } else {
    Omega <- as.matrix(Omega)
    if (!identical(dim(Omega), c(P, P))) {
      stop("Omega must match residual variances.", call. = FALSE)
    }
    diag(Omega)
  }
  sum(vw * omega_diag * sigma2)
}

#' @noRd
.dkge_effect_noise_moment <- function(subject, Omega = NULL,
                                      voxel_weights = NULL,
                                      obs_mask = NULL) {
  Lambda <- subject$effect_noise_cov
  if (is.null(Lambda)) {
    stop(sprintf(
      "Subject '%s' lacks `effect_noise_cov` for analytic debiasing.",
      subject$id %||% "(unnamed)"
    ), call. = FALSE)
  }
  Lambda <- as.matrix(Lambda)
  q <- nrow(Lambda)
  obs <- .dkge_observed_rows(obs_mask, q)
  if (!length(obs)) return(.dkge_empty_effect_moment(Lambda))
  if (length(obs) < q) {
    Lambda[-obs, ] <- 0
    Lambda[, -obs] <- 0
  }
  out <- .dkge_noise_trace(subject, Omega, voxel_weights) * Lambda
  (out + t(out)) / 2
}

#' Map a raw q-space effect moment into the K-metric
#'
#' Computes \eqn{K^{1/2} R^\top M R K^{1/2}}. `A = R K^{1/2}` may be supplied by
#' the caller so that a pool of moments shares one product; the two forms agree
#' to machine precision because `Khalf` is symmetric.
#'
#' @noRd
.dkge_transform_effect_moment <- function(M, R, Khalf, A = NULL) {
  if (is.null(A)) A <- R %*% Khalf
  out <- crossprod(A, M %*% A)
  out <- (out + t(out)) / 2
  # Effect labels come from the moment (which inherits the beta row names);
  # the kernel roots themselves are unnamed.
  dn <- dimnames(M) %||% dimnames(Khalf)
  if (!is.null(dn)) dimnames(out) <- dn
  out
}

#' Cache the sample-independent pair matrices used when pooling
#'
#' `structural` is the observed-pair indicator and `reliability` the pairwise
#' precision weight; both depend only on the subject's mask and precision
#' vector, so resampling loops can build them once and reuse them.
#'
#' @noRd
.dkge_pool_cache <- function(obs_masks, effect_precision, q) {
  S <- length(effect_precision)
  structural <- vector("list", S)
  reliability <- vector("list", S)
  for (s in seq_len(S)) {
    mask <- if (is.null(obs_masks)) rep(TRUE, q) else as.logical(obs_masks[[s]])
    if (length(mask) != q) {
      stop(sprintf("Observation mask for subject %d must have length %d.", s, q),
           call. = FALSE)
    }
    p <- effect_precision[[s]]
    if (length(p) != q) stop("Effect precision must have length q.", call. = FALSE)
    m <- as.numeric(mask)
    structural[[s]] <- tcrossprod(m)
    reliability[[s]] <- tcrossprod(sqrt(pmax(as.numeric(p), 0)) * m)
  }
  list(structural = structural, reliability = reliability, q = q)
}

#' Pool subject effect moments
#'
#' This is the single pooling engine used by batch and resampled fits. With no
#' effect weighting it preserves the historical missingness policies, now in
#' unmixed effect space. With count or explicit precision weighting it computes
#' a pair-normalized precision-weighted mean and rescales it by the observed
#' subject mass \eqn{\sum_s a_s w_s 1[\mathrm{obs}_s]} so that unit precision
#' recovers the unweighted sum (and equals the legacy cohort-mass scale under
#' full coverage).
#'
#' Two coverage tallies are returned and the missingness policies use whichever
#' matches the pooling that produced `pooled`:
#'
#' * `pair_counts` is the sample-weighted mass of subjects observing each
#'   effect pair, \eqn{\sum_s a_s 1[\text{observed}]}. It is a plain integer
#'   count only when the sample weights `a_s` are all 1 (batch fits, pair
#'   bootstrap); the multiplier schemes make it fractional. It drives the
#'   policies on the unweighted branch, where `pooled` is a sum over subjects.
#' * `pair_ess` is Kish's effective sample size of the same pairs under the
#'   precision weights, \eqn{(\sum_s h_s)^2 / \sum_s h_s^2}. It is scale-free,
#'   so it is already on the subject-count scale. It drives the `mask` and
#'   `shrink` policies on the precision-weighted branch. With unit precision and
#'   equal subject weights `pair_ess` equals `pair_counts`, so the two branches
#'   apply identical thresholds in that case.
#'
#' `missingness = "rescale"` turns the pooled sum into a per-pair mean on both
#' branches. Unweighted, that is `numerator / pair_counts`; under precision
#' weighting it is `numerator / pair_weight`, i.e. `pooled` divided by the
#' observed subject mass that scaled it up. Dividing by `pair_ess` instead
#' would apply the per-pair normalization twice.
#'
#' @param pool_cache Optional list from `.dkge_pool_cache()` holding the
#'   sample-independent `structural`/`reliability` pair matrices.
#' @noRd
.dkge_pool_effect_moments <- function(moments,
                                      subject_weights,
                                      obs_masks,
                                      effect_precision,
                                      effect_method = "none",
                                      missingness = c("none", "rescale", "mask", "shrink"),
                                      miss_args = list(),
                                      sample_weights = NULL,
                                      pool_cache = NULL) {
  missingness <- match.arg(missingness)
  S <- length(moments)
  q <- nrow(moments[[1]])
  if (is.null(sample_weights)) sample_weights <- rep(1, S)
  if (length(subject_weights) != S || length(sample_weights) != S) {
    stop("Subject and sample weights must match the number of moments.",
         call. = FALSE)
  }
  if (!is.null(pool_cache) &&
      (!identical(pool_cache$q, q) || length(pool_cache$structural) != S)) {
    pool_cache <- NULL
  }
  if (is.null(pool_cache)) {
    pool_cache <- .dkge_pool_cache(obs_masks, effect_precision, q)
  }

  numerator <- matrix(0, q, q)
  pair_counts <- matrix(0, q, q)
  pair_weight <- matrix(0, q, q)
  pair_weight_sq <- matrix(0, q, q)
  pair_mass <- matrix(0, q, q)

  for (s in seq_len(S)) {
    pair_counts <- pair_counts + sample_weights[s] * pool_cache$structural[[s]]
    sw <- sample_weights[s] * subject_weights[s]
    pair_mass <- pair_mass + sw * pool_cache$structural[[s]]
    h <- sw * pool_cache$reliability[[s]]
    numerator <- numerator + h * moments[[s]]
    pair_weight <- pair_weight + h
    pair_weight_sq <- pair_weight_sq + h^2
  }

  active_weighting <- !identical(effect_method, "none")
  weight_scale <- mean(sample_weights)
  if (!is.finite(weight_scale) || weight_scale <= 0) weight_scale <- 1
  if (active_weighting) {
    valid <- pair_weight > 0
    # `per_pair` is the precision-weighted mean over the subjects contributing
    # to each pair. `pooled` scales that mean by observed subject mass
    # (equals cohort mass under full coverage; avoids S/k inflation when
    # only k of S subjects observe the pair).
    per_pair <- matrix(0, q, q)
    per_pair[valid] <- numerator[valid] / pair_weight[valid]
    pooled <- pair_mass * per_pair
    pair_ess <- matrix(0, q, q)
    pair_ess[valid] <- pair_weight[valid]^2 / pair_weight_sq[valid]

    if (identical(missingness, "rescale")) {
      # Turn the observed-mass sum back into a per-pair mean. The divisor is
      # the mass that was multiplied in above -- NOT `pair_ess`. With unit
      # precision and equal subject weights this reproduces the unweighted
      # branch exactly (`numerator / pair_counts`).
      pooled <- per_pair
    } else if (identical(missingness, "mask")) {
      min_pairs <- .dkge_miss_scalar(miss_args$min_pairs, "min_pairs", 1)
      pooled[pair_ess < min_pairs] <- 0
    } else if (identical(missingness, "shrink")) {
      gamma <- .dkge_miss_scalar(miss_args$gamma, "gamma", 1)
      max_ess <- max(pair_ess)
      rel <- if (max_ess > 0) (pair_ess / max_ess)^gamma else pair_ess
      diagonal_target <- diag(diag(per_pair), q, q)
      pooled <- rel * per_pair + (1 - rel) * diagonal_target
    }
  } else {
    pooled <- numerator
    pair_ess <- pair_counts
    pooled <- .dkge_apply_missingness(
      pooled,
      pair_counts = pair_counts,
      missingness = missingness,
      miss_args = miss_args,
      weight_scale = weight_scale,
      pair_weight = pair_weight
    )
  }

  pooled <- (pooled + t(pooled)) / 2
  effect_names <- rownames(moments[[1]])
  if (!is.null(effect_names)) {
    dn <- list(effect_names, effect_names)
    dimnames(pooled) <- dn
    dimnames(pair_counts) <- dn
    dimnames(pair_weight) <- dn
    dimnames(pair_ess) <- dn
  }
  list(
    pooled = pooled,
    pair_counts = pair_counts,
    pair_weight = pair_weight,
    pair_ess = pair_ess
  )
}

#' @noRd
.dkge_spectral_diagnostics <- function(M, tolerance = 1e-10) {
  values <- eigen((M + t(M)) / 2, symmetric = TRUE, only.values = TRUE)$values
  scale <- max(1, max(abs(values)))
  negative <- values[values < -tolerance * scale]
  list(
    min_eigenvalue = min(values),
    negative_count = length(negative),
    negative_mass = sum(abs(negative)),
    total_spectral_mass = sum(abs(values))
  )
}

#' @noRd
.dkge_build_moment_pool <- function(subjects,
                                    B_list,
                                    Omega_list,
                                    voxel_weights,
                                    obs_masks,
                                    subject_weights,
                                    effect_precision,
                                    effect_method,
                                    R,
                                    Khalf,
                                    missingness,
                                    miss_args,
                                    debias = c("none", "analytic", "split_half"),
                                    sample_weights = NULL,
                                    contribs = TRUE,
                                    pool_cache = NULL) {
  debias <- match.arg(debias)
  S <- length(B_list)
  analytic <- identical(debias, "analytic")
  # `moments_raw` and `noise_moments` only differ from `moments` under analytic
  # debiasing; keeping them otherwise would store two more S x q x q copies of
  # the same numbers.
  moments_raw <- if (analytic) vector("list", S) else NULL
  noise_moments <- if (analytic) vector("list", S) else NULL
  moments <- vector("list", S)

  for (s in seq_len(S)) {
    vw <- .dkge_subject_voxel_weights(voxel_weights, s, ncol(B_list[[s]]),
                                      subjects[[s]]$id)
    raw_s <- if (identical(debias, "split_half")) {
      .dkge_split_effect_moment(subjects[[s]]$split_betas,
                                Omega_list[[s]], vw, obs_masks[[s]],
                                subject_id = subjects[[s]]$id)
    } else {
      .dkge_effect_moment(B_list[[s]], Omega_list[[s]], vw, obs_masks[[s]])
    }
    if (analytic) {
      noise_s <- .dkge_effect_noise_moment(subjects[[s]], Omega_list[[s]], vw,
                                           obs_masks[[s]])
      moments_raw[[s]] <- raw_s
      noise_moments[[s]] <- noise_s
      m <- raw_s - noise_s
      moments[[s]] <- (m + t(m)) / 2
    } else {
      moments[[s]] <- raw_s
    }
  }

  pool_cache <- pool_cache %||% .dkge_pool_cache(obs_masks, effect_precision,
                                                 nrow(moments[[1]]))
  pool <- .dkge_pool_effect_moments(
    moments = moments,
    subject_weights = subject_weights,
    obs_masks = obs_masks,
    effect_precision = effect_precision,
    effect_method = effect_method,
    missingness = missingness,
    miss_args = miss_args,
    sample_weights = sample_weights,
    pool_cache = pool_cache
  )
  A <- R %*% Khalf
  Chat <- .dkge_transform_effect_moment(pool$pooled, R, Khalf, A = A)
  contrib_list <- if (isTRUE(contribs)) {
    lapply(moments, .dkge_transform_effect_moment, R = R, Khalf = Khalf, A = A)
  } else {
    NULL
  }

  c(pool, list(
    Chat = Chat,
    moments = moments,
    moments_raw = moments_raw,
    noise_moments = noise_moments,
    contribs = contrib_list,
    pool_cache = pool_cache,
    diagnostics = list(
      effect = .dkge_spectral_diagnostics(pool$pooled),
      transformed = .dkge_spectral_diagnostics(Chat)
    )
  ))
}

#' Re-pool the q-space sufficient statistics stored on a fitted model
#'
#' @noRd
.dkge_repool_fit <- function(fit, sample_weights = NULL, indices = NULL,
                             missingness = NULL, miss_args = NULL) {
  S <- length(fit$Btil)
  indices <- indices %||% seq_len(S)
  moments <- fit[["effect_moments"]]
  if (is.null(moments)) return(NULL)
  subject_ids <- fit$subject_ids %||% as.character(seq_len(S))
  masks <- .dkge_obs_masks_from_provenance(fit$provenance, subject_ids,
                                           nrow(fit$K))
  if (is.null(masks)) {
    masks <- replicate(S, rep(TRUE, nrow(fit$K)), simplify = FALSE)
  }
  precision <- fit$effect_precision
  if (is.null(precision)) precision <- lapply(masks, as.numeric)
  sample_weights <- sample_weights %||% rep(1, length(indices))
  if (length(sample_weights) == S && length(indices) != S) {
    sample_weights <- sample_weights[indices]
  }
  missingness <- missingness %||% fit$missingness %||% "none"
  miss_args <- miss_args %||% fit$miss_args %||% list()

  # `fit$pool_cache` holds the sample-independent pair matrices, so a
  # resampling loop never recomputes tcrossprod(mask) / tcrossprod(sqrt(p)).
  cache <- fit$pool_cache
  if (!is.null(cache) && !identical(indices, seq_len(S))) {
    cache <- list(structural = cache$structural[indices],
                  reliability = cache$reliability[indices],
                  q = cache$q)
  }

  pool <- .dkge_pool_effect_moments(
    moments = moments[indices],
    subject_weights = fit$weights[indices],
    obs_masks = masks[indices],
    effect_precision = precision[indices],
    effect_method = fit$effect_weight_spec$method %||% "none",
    missingness = missingness,
    miss_args = miss_args,
    sample_weights = sample_weights,
    pool_cache = cache
  )
  pool$Chat <- .dkge_transform_effect_moment(pool$pooled, fit$R, fit$Khalf)
  pool
}
