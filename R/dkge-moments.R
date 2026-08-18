# dkge-moments.R
# Central q-space moment construction, debiasing, and reliability pooling.

#' @noRd
.dkge_scale_effect_columns <- function(B, weights) {
  if (is.null(weights) || !length(weights)) return(B)
  if (length(weights) != ncol(B)) {
    weights <- rep(weights, length.out = ncol(B))
  }
  sweep(B, 2L, sqrt(pmax(as.numeric(weights), 0)), "*")
}

#' Apply raw-effect masking and factorizable right-side weights
#'
#' This is the canonical subject block convention shared by moment, subject-
#' energy, and physical-factor construction. Observation masks are applied
#' before any effect-coordinate transform.
#'
#' @noRd
.dkge_right_weighted_beta <- function(B, Omega = NULL, voxel_weights = NULL,
                                      obs_mask = NULL) {
  B <- as.matrix(B)
  q <- nrow(B)
  obs <- .dkge_observed_rows(obs_mask, q)
  if (length(obs) < q) B[-obs, ] <- 0
  B <- .dkge_scale_effect_columns(B, voxel_weights)

  if (is.null(Omega)) return(B)
  if (is.vector(Omega)) {
    if (length(Omega) != ncol(B)) {
      stop("Diagonal Omega must have one entry per beta column.", call. = FALSE)
    }
    return(sweep(B, 2L, sqrt(pmax(as.numeric(Omega), 0)), "*"))
  }
  Omega <- as.matrix(Omega)
  if (!identical(dim(Omega), c(ncol(B), ncol(B)))) {
    stop("Full Omega must be square with one row per beta column.",
         call. = FALSE)
  }
  B %*% sqrtm_sym(Omega)
}

#' @noRd
.dkge_effect_moment <- function(B, Omega = NULL, voxel_weights = NULL,
                                obs_mask = NULL) {
  B <- as.matrix(B)
  Bwork <- .dkge_right_weighted_beta(B, Omega, voxel_weights, obs_mask)
  M <- tcrossprod(Bwork)
  M <- (M + t(M)) / 2
  if (!is.null(rownames(B))) dimnames(M) <- list(rownames(B), rownames(B))
  M
}

#' @noRd
.dkge_split_effect_moment <- function(split_betas, Omega = NULL,
                                      voxel_weights = NULL,
                                      obs_mask = NULL) {
  if (is.null(split_betas) || !is.list(split_betas) || length(split_betas) != 2L) {
    stop("Split-half debiasing requires two stored split beta matrices for every subject.",
         call. = FALSE)
  }
  B1 <- as.matrix(split_betas[[1]])
  B2 <- as.matrix(split_betas[[2]])
  if (!identical(dim(B1), dim(B2))) {
    stop("Stored split beta matrices must have identical dimensions.", call. = FALSE)
  }
  B1 <- .dkge_right_weighted_beta(B1, Omega, voxel_weights, obs_mask)
  B2 <- .dkge_right_weighted_beta(B2, Omega, voxel_weights, obs_mask)
  cross <- B1 %*% t(B2)
  out <- (cross + t(cross)) / 2
  if (!is.null(rownames(B1))) dimnames(out) <- list(rownames(B1), rownames(B1))
  out
}

#' @noRd
.dkge_noise_trace <- function(subject, Omega = NULL, voxel_weights = NULL) {
  if (!is.null(subject$noise_trace)) {
    scope <- subject$noise_trace_scope %||% "preweighted"
    if (!identical(scope, "unweighted")) {
      return(as.numeric(subject$noise_trace))
    }
    weights_are_identity <- function(x, n) {
      if (is.null(x) || !length(x)) return(TRUE)
      if (is.vector(x)) {
        return(length(x) == n && all(abs(as.numeric(x) - 1) <= 1e-12))
      }
      x <- as.matrix(x)
      identical(dim(x), c(n, n)) &&
        max(abs(x - diag(n))) <= 1e-12
    }
    n_spatial <- length(subject$residual_variance %||% numeric(0))
    if (n_spatial > 0L && weights_are_identity(Omega, n_spatial) &&
        weights_are_identity(voxel_weights, n_spatial)) {
      return(as.numeric(subject$noise_trace))
    }
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
  if (!is.null(subject$residual_df) &&
      (!is.finite(subject$residual_df) || subject$residual_df <= 0)) {
    stop(sprintf(
      "Subject '%s' requires positive residual degrees of freedom for analytic debiasing.",
      subject$id %||% "(unnamed)"
    ), call. = FALSE)
  }
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
  if (length(obs) < q) {
    Lambda[-obs, ] <- 0
    Lambda[, -obs] <- 0
  }
  out <- .dkge_noise_trace(subject, Omega, voxel_weights) * Lambda
  (out + t(out)) / 2
}

#' @noRd
.dkge_transform_effect_moment <- function(M, R, Khalf) {
  out <- Khalf %*% t(R) %*% M %*% R %*% Khalf
  out <- (out + t(out)) / 2
  if (!is.null(dimnames(Khalf))) dimnames(out) <- dimnames(Khalf)
  out
}

#' Pool subject effect moments
#'
#' This is the single pooling engine used by batch and resampled fits. With no
#' effect weighting it preserves the historical missingness policies, now in
#' unmixed effect space. With count or explicit precision weighting it computes
#' a pair-normalised precision-weighted mean and rescales it by the total subject
#' mass so that unit precision exactly recovers the legacy full-coverage sum.
#'
#' @noRd
.dkge_pool_effect_moments <- function(moments,
                                      subject_weights,
                                      obs_masks,
                                      effect_precision,
                                      effect_method = "none",
                                      missingness = c("none", "rescale", "mask", "shrink"),
                                      miss_args = list(),
                                      sample_weights = NULL) {
  missingness <- match.arg(missingness)
  S <- length(moments)
  q <- nrow(moments[[1]])
  if (is.null(sample_weights)) sample_weights <- rep(1, S)
  if (length(subject_weights) != S || length(sample_weights) != S) {
    stop("Subject and sample weights must match the number of moments.",
         call. = FALSE)
  }
  if (is.null(obs_masks)) obs_masks <- replicate(S, rep(TRUE, q), simplify = FALSE)

  numerator <- matrix(0, q, q)
  pair_counts <- matrix(0, q, q)
  pair_weight <- matrix(0, q, q)
  pair_weight_sq <- matrix(0, q, q)
  active_weighting <- !identical(effect_method, "none")

  for (s in seq_len(S)) {
    mask <- as.logical(obs_masks[[s]])
    if (length(mask) != q) mask <- rep(TRUE, q)
    p <- effect_precision[[s]]
    if (length(p) != q) stop("Effect precision must have length q.", call. = FALSE)
    p <- pmax(as.numeric(p), 0)
    if (active_weighting) mask <- mask & p > 0
    structural <- tcrossprod(as.numeric(mask))
    pair_counts <- pair_counts + sample_weights[s] * structural

    reliability <- tcrossprod(sqrt(p)) * structural
    h <- sample_weights[s] * subject_weights[s] * reliability
    numerator <- numerator + h * moments[[s]]
    pair_weight <- pair_weight + h
    pair_weight_sq <- pair_weight_sq + h^2
  }

  if (active_weighting) {
    pooled <- matrix(0, q, q)
    valid <- pair_weight > 0
    cohort_mass <- sum(sample_weights * subject_weights)
    pooled[valid] <- cohort_mass * numerator[valid] / pair_weight[valid]
    pair_ess <- matrix(0, q, q)
    pair_ess[valid] <- pair_weight[valid]^2 / pair_weight_sq[valid]

    if (identical(missingness, "mask")) {
      min_pairs <- miss_args$min_pairs %||% 1
      if (!is.numeric(min_pairs) || length(min_pairs) != 1L ||
          !is.finite(min_pairs) || min_pairs < 0) {
        stop("`miss_args$min_pairs` must be a non-negative numeric scalar.",
             call. = FALSE)
      }
      pooled[pair_ess < min_pairs] <- 0
    } else if (identical(missingness, "shrink")) {
      gamma <- miss_args$gamma %||% 1
      if (!is.numeric(gamma) || length(gamma) != 1L ||
          !is.finite(gamma) || gamma < 0) {
        stop("`miss_args$gamma` must be a non-negative numeric scalar.",
             call. = FALSE)
      }
      max_ess <- max(pair_ess)
      rel <- if (max_ess > 0) (pair_ess / max_ess)^gamma else pair_ess
      diagonal_target <- diag(diag(pooled), q, q)
      pooled <- rel * pooled + (1 - rel) * diagonal_target
    }
  } else {
    pooled <- numerator
    pair_ess <- pair_counts
    pooled <- .dkge_apply_missingness(
      pooled,
      pair_counts = pair_counts,
      missingness = missingness,
      miss_args = miss_args
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
                                    sample_weights = NULL) {
  debias <- match.arg(debias)
  S <- length(B_list)
  moments_raw <- vector("list", S)
  noise_moments <- vector("list", S)
  moments <- vector("list", S)

  for (s in seq_len(S)) {
    vw <- if (is.list(voxel_weights)) voxel_weights[[s]] else voxel_weights
    moments_raw[[s]] <- if (identical(debias, "split_half")) {
      .dkge_split_effect_moment(subjects[[s]]$split_betas,
                                Omega_list[[s]], vw, obs_masks[[s]])
    } else {
      .dkge_effect_moment(B_list[[s]], Omega_list[[s]], vw, obs_masks[[s]])
    }
    noise_moments[[s]] <- if (identical(debias, "analytic")) {
      .dkge_effect_noise_moment(subjects[[s]], Omega_list[[s]], vw,
                                obs_masks[[s]])
    } else {
      matrix(0, nrow(B_list[[s]]), nrow(B_list[[s]]),
             dimnames = dimnames(moments_raw[[s]]))
    }
    moments[[s]] <- moments_raw[[s]] - noise_moments[[s]]
    moments[[s]] <- (moments[[s]] + t(moments[[s]])) / 2
  }

  pool <- .dkge_pool_effect_moments(
    moments = moments,
    subject_weights = subject_weights,
    obs_masks = obs_masks,
    effect_precision = effect_precision,
    effect_method = effect_method,
    missingness = missingness,
    miss_args = miss_args,
    sample_weights = sample_weights
  )
  Chat <- .dkge_transform_effect_moment(pool$pooled, R, Khalf)
  contribs <- lapply(moments, .dkge_transform_effect_moment,
                     R = R, Khalf = Khalf)

  c(pool, list(
    Chat = Chat,
    Chat_sym = Chat,
    moments = moments,
    moments_raw = moments_raw,
    noise_moments = noise_moments,
    contribs = contribs,
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
  sample_weights <- sample_weights %||% rep(1, length(indices))
  if (length(sample_weights) == S && length(indices) != S) {
    sample_weights <- sample_weights[indices]
  }
  integer_multiplicity <- length(sample_weights) == length(indices) &&
    all(is.finite(sample_weights)) && all(sample_weights >= 0) &&
    all(abs(sample_weights - round(sample_weights)) < 1e-12)
  if (integer_multiplicity && !is.null(fit$subjects)) {
    expanded <- rep(indices, times = as.integer(round(sample_weights)))
    if (length(expanded) >= 2L) {
      refit <- .dkge_refit_training_subjects(
        fit, expanded,
        missingness = missingness %||% fit$missingness %||% "none",
        miss_args = miss_args %||% fit$miss_args %||% list()
      )
      return(list(
        pooled = refit$effect_moment,
        pair_counts = refit$pair_counts,
        pair_weight = refit$pair_weight,
        pair_ess = refit$pair_ess,
        Chat = refit$Chat,
        R = refit$R,
        effect_precision = refit$effect_precision,
        effect_precision_diagnostics = refit$effect_precision_diagnostics,
        subject_weights = refit$weights,
        expanded_indices = expanded,
        multiplicities = sample_weights,
        refit = "literal_subject_multiset"
      ))
    }
  }
  moments <- fit$effect_moments
  if (is.null(moments)) return(NULL)
  subject_ids <- fit$subject_ids %||% as.character(seq_len(S))
  masks <- .dkge_obs_masks_from_provenance(fit$provenance, subject_ids,
                                           nrow(fit$K))
  if (is.null(masks)) {
    masks <- replicate(S, rep(TRUE, nrow(fit$K)), simplify = FALSE)
  }
  precision <- fit$effect_precision
  if (is.null(precision)) precision <- lapply(masks, as.numeric)
  missingness <- missingness %||% fit$missingness %||% "none"
  miss_args <- miss_args %||% fit$miss_args %||% list()

  pool <- .dkge_pool_effect_moments(
    moments = moments[indices],
    subject_weights = fit$weights[indices],
    obs_masks = masks[indices],
    effect_precision = precision[indices],
    effect_method = fit$effect_weight_spec$method %||% "none",
    missingness = missingness,
    miss_args = miss_args,
    sample_weights = sample_weights
  )
  pool$Chat <- .dkge_transform_effect_moment(pool$pooled, fit$R, fit$Khalf)
  pool
}
