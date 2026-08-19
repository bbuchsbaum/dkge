# dkge-effect-weights.R
# Subject-by-effect reliability specifications for q-space pooling.

#' Specify effect-level reliability weighting
#'
#' Effect weights control how unequally estimated subject-by-effect rows enter
#' the pooled second moment. They are distinct from [dkge_weights()], which
#' weights voxels or parcels, and from `w_method`, which supplies one scalar per
#' subject.
#'
#' @param method Effect precision source. `"none"` gives every observed effect
#'   equal precision, `"count"` uses `effect_n` stored on each subject,
#'   `"precision"` uses explicitly supplied `effect_precision` values, and
#'   `"random_effects"` combines within-subject uncertainty with a
#'   DerSimonian-Laird between-subject variance estimate.
#' @param within Within-subject variance source for `method = "random_effects"`.
#'   `"noise"` uses `effect_noise_cov` and residual variances, `"count"` uses
#'   `1 / effect_n`, and `"auto"` prefers complete noise statistics before
#'   falling back to complete counts.
#' @param tau_method Between-subject variance estimator. Currently only the
#'   non-negative DerSimonian-Laird moment estimator (`"DL"`) is supported.
#' @param max_ratio Maximum random-effects precision relative to the median
#'   positive precision for the same effect. This cap prevents a single subject
#'   from owning an effect geometry; `Inf` disables it.
#' @param variance_floor Positive lower bound for within-subject variances.
#' @param tau_floor Non-negative lower bound for the estimated between-subject
#'   variance.
#'
#' @return An object of class `dkge_effect_weights`, suitable for the
#'   `effect_weights` argument of [dkge_fit()]. A precision of zero excludes that
#'   subject-effect row. Pairwise reliability is the geometric mean of the two
#'   effect precisions.
#' @details For random effects, effect `c` uses
#'   `p_sc = 1 / (v_sc + tau_c^2)`. The DerSimonian-Laird estimate is
#'   `max(tau_floor, (Q - (k - 1)) / C, 0)`, followed by the documented
#'   median-relative cap. Because subject parcel spaces need not be aligned,
#'   heterogeneity is estimated from each effect row's spatial mean.
#'
#'   The package-wide default remains `"none"`. `"random_effects"` is an
#'   opt-in policy for materially unbalanced trialwise designs.
#' @export
#' @examples
#' dkge_effect_weights("none")
#' dkge_effect_weights("count")
#' dkge_effect_weights("random_effects", within = "auto")
dkge_effect_weights <- function(
    method = c("none", "count", "precision", "random_effects"),
    within = c("auto", "noise", "count"),
    tau_method = "DL",
    max_ratio = 10,
    variance_floor = 1e-8,
    tau_floor = 0) {
  method <- match.arg(method)
  within <- match.arg(within)
  if (!identical(tau_method, "DL")) {
    stop("`tau_method` currently must be 'DL'.", call. = FALSE)
  }
  if (!is.numeric(max_ratio) || length(max_ratio) != 1L ||
      is.na(max_ratio) || max_ratio < 1) {
    stop("`max_ratio` must be at least 1 or Inf.", call. = FALSE)
  }
  if (!is.numeric(variance_floor) || length(variance_floor) != 1L ||
      !is.finite(variance_floor) || variance_floor <= 0) {
    stop("`variance_floor` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(tau_floor) || length(tau_floor) != 1L ||
      !is.finite(tau_floor) || tau_floor < 0) {
    stop("`tau_floor` must be a finite non-negative scalar.", call. = FALSE)
  }
  structure(
    list(
      method = method,
      within = within,
      tau_method = tau_method,
      max_ratio = as.numeric(max_ratio),
      variance_floor = as.numeric(variance_floor),
      tau_floor = as.numeric(tau_floor)
    ),
    class = "dkge_effect_weights"
  )
}

#' @noRd
.dkge_effect_obs_mask <- function(obs_masks, s, q, subject_id) {
  mask <- if (is.null(obs_masks)) rep(TRUE, q) else as.logical(obs_masks[[s]])
  if (length(mask) != q) {
    stop(sprintf(
      "Observation mask for subject '%s' must have one entry per effect (got %d, expected %d).",
      subject_id, length(mask), q
    ), call. = FALSE)
  }
  mask
}

#' @noRd
.dkge_resolve_effect_precision <- function(dataset, spec, obs_masks = NULL,
                                           indices = seq_len(dataset$n_subjects)) {
  if (is.null(spec)) spec <- dkge_effect_weights("none")
  if (!inherits(spec, "dkge_effect_weights")) {
    stop("`effect_weights` must be created by `dkge_effect_weights()`.",
         call. = FALSE)
  }

  method <- spec$method %||% "none"
  q <- dataset$q
  effects <- dataset$effects
  source <- switch(
    method,
    none = NULL,
    count = dataset$effect_n,
    precision = dataset$effect_precision,
    random_effects = NULL,
    stop("Unknown effect-weighting method.", call. = FALSE)
  )

  if (identical(method, "random_effects")) {
    return(.dkge_random_effects_precision(
      dataset, spec, obs_masks = obs_masks, indices = indices
    ))
  }

  out <- vector("list", length(indices))
  for (j in seq_along(indices)) {
    s <- indices[[j]]
    mask <- .dkge_effect_obs_mask(
      obs_masks, s, q, dataset$subject_ids[[s]]
    )

    if (identical(method, "none")) {
      p <- as.numeric(mask)
    } else {
      p <- source[[s]]
      if (is.null(p)) {
        stop(sprintf(
          "Subject '%s' has no `%s` values required by effect_weights='%s'.",
          dataset$subject_ids[[s]],
          if (method == "count") "effect_n" else "effect_precision",
          method
        ), call. = FALSE)
      }
      p <- as.numeric(p)
      if (length(p) != q || any(!is.finite(p)) || any(p < 0)) {
        stop(sprintf(
          "Subject '%s' must have %d finite, non-negative effect precision values.",
          dataset$subject_ids[[s]], q
        ), call. = FALSE)
      }
      p[!mask] <- 0
    }
    names(p) <- effects
    out[[j]] <- p
  }
  names(out) <- dataset$subject_ids[indices]
  attr(out, "diagnostics") <- list(method = method)
  out
}

#' Estimate random-effects subject-by-effect precision
#'
#' The scalar study outcome for an effect is the mean of its beta row across
#' spatial units. Noise-based within variance is the corresponding variance of
#' that mean under spatial independence; count-based variance is `1 / n_sc`.
#' The resulting precision weights the whole subject effect row.
#'
#' @keywords internal
#' @noRd
.dkge_random_effects_precision <- function(dataset, spec, obs_masks = NULL,
                                            indices = seq_len(dataset$n_subjects)) {
  q <- dataset$q
  effects <- dataset$effects
  S <- length(indices)
  masks <- lapply(indices, function(s) {
    .dkge_effect_obs_mask(obs_masks, s, q, dataset$subject_ids[[s]])
  })

  noise_available <- vapply(seq_along(indices), function(j) {
    s <- indices[[j]]
    sigma2 <- dataset$residual_variance[[s]]
    Lambda <- dataset$effect_noise_cov[[s]]
    !is.null(sigma2) && length(sigma2) == ncol(dataset$betas[[s]]) &&
      all(is.finite(sigma2)) && all(sigma2 >= 0) &&
      !is.null(Lambda) && identical(dim(as.matrix(Lambda)), c(q, q)) &&
      all(is.finite(Lambda))
  }, logical(1))
  count_available <- vapply(seq_along(indices), function(j) {
    s <- indices[[j]]
    n <- dataset$effect_n[[s]]
    mask <- masks[[j]]
    !is.null(n) && length(n) == q &&
      all(is.finite(n[mask])) && all(n[mask] >= 0)
  }, logical(1))

  within_source <- spec$within %||% "auto"
  if (identical(within_source, "auto")) {
    within_source <- if (all(noise_available)) {
      "noise"
    } else if (all(count_available)) {
      "count"
    } else {
      stop(paste0(
        "Random-effects weighting requires complete residual/effect covariance ",
        "statistics or complete `effect_n` counts for all selected subjects."
      ), call. = FALSE)
    }
  } else if (identical(within_source, "noise") && !all(noise_available)) {
    stop("Random-effects within='noise' requires residual variance and effect covariance for every subject.",
         call. = FALSE)
  } else if (identical(within_source, "count") && !all(count_available)) {
    stop("Random-effects within='count' requires `effect_n` for every subject.",
         call. = FALSE)
  }

  estimate <- matrix(NA_real_, S, q,
                     dimnames = list(dataset$subject_ids[indices], effects))
  within_variance <- matrix(Inf, S, q, dimnames = dimnames(estimate))
  for (j in seq_along(indices)) {
    s <- indices[[j]]
    B <- as.matrix(dataset$betas[[s]])
    mask <- masks[[j]]
    estimate[j, mask] <- rowMeans(B[mask, , drop = FALSE])
    if (identical(within_source, "count")) {
      n <- as.numeric(dataset$effect_n[[s]])
      valid <- mask & n > 0
      within_variance[j, valid] <- 1 / n[valid]
    } else {
      sigma2 <- as.numeric(dataset$residual_variance[[s]])
      Lambda <- as.matrix(dataset$effect_noise_cov[[s]])
      P <- ncol(B)
      variance_of_mean <- diag(Lambda) * sum(sigma2) / (P^2)
      valid <- mask & is.finite(variance_of_mean) & variance_of_mean >= 0
      within_variance[j, valid] <- variance_of_mean[valid]
    }
  }
  finite <- is.finite(within_variance)
  within_variance[finite] <- pmax(
    within_variance[finite], spec$variance_floor %||% 1e-8
  )

  tau2 <- stats::setNames(numeric(q), effects)
  Q_stat <- stats::setNames(numeric(q), effects)
  k_effect <- stats::setNames(integer(q), effects)
  raw_precision <- matrix(0, S, q, dimnames = dimnames(estimate))
  capped_precision <- raw_precision
  capped <- matrix(FALSE, S, q, dimnames = dimnames(estimate))
  for (c in seq_len(q)) {
    valid <- vapply(seq_len(S), function(j) {
      masks[[j]][c] && is.finite(estimate[j, c]) &&
        is.finite(within_variance[j, c])
    }, logical(1))
    k <- sum(valid)
    k_effect[[c]] <- k
    if (!k) next
    v <- within_variance[valid, c]
    y <- estimate[valid, c]
    fixed_weight <- 1 / v
    fixed_mean <- sum(fixed_weight * y) / sum(fixed_weight)
    Q <- sum(fixed_weight * (y - fixed_mean)^2)
    C_dl <- sum(fixed_weight) -
      sum(fixed_weight^2) / sum(fixed_weight)
    tau <- if (k >= 2L && C_dl > 0) {
      max(0, (Q - (k - 1L)) / C_dl)
    } else {
      0
    }
    tau <- max(tau, spec$tau_floor %||% 0)
    tau2[[c]] <- tau
    Q_stat[[c]] <- Q
    p <- 1 / (v + tau)
    raw_precision[valid, c] <- p

    cap_ratio <- spec$max_ratio %||% 10
    cap <- if (is.finite(cap_ratio)) cap_ratio * stats::median(p) else Inf
    p_capped <- pmin(p, cap)
    capped_precision[valid, c] <- p_capped
    capped[valid, c] <- p > cap
  }

  to_subject_list <- function(x) {
    out <- lapply(seq_len(S), function(j) {
      value <- as.numeric(x[j, ])
      names(value) <- effects
      value
    })
    names(out) <- dataset$subject_ids[indices]
    out
  }
  raw_share <- apply(raw_precision, 2, function(p) {
    if (sum(p) > 0) max(p) / sum(p) else 0
  })
  capped_share <- apply(capped_precision, 2, function(p) {
    if (sum(p) > 0) max(p) / sum(p) else 0
  })
  names(raw_share) <- effects
  names(capped_share) <- effects

  out <- to_subject_list(capped_precision)
  attr(out, "diagnostics") <- list(
    method = "random_effects",
    within_source = within_source,
    tau_method = spec$tau_method %||% "DL",
    tau2 = tau2,
    Q = Q_stat,
    subjects_per_effect = k_effect,
    within_variance = within_variance,
    scalar_effect_estimate = estimate,
    raw_precision = raw_precision,
    precision = capped_precision,
    max_ratio = spec$max_ratio %||% 10,
    capped = capped,
    raw_max_share = raw_share,
    max_share = capped_share,
    assumption = "between-subject heterogeneity is estimated from the spatial mean of each effect row"
  )
  out
}
