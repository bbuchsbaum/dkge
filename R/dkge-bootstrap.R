# dkge-bootstrap.R
# Fast bootstrap approximations for DKGE fits.

#' Subject-level projection bootstrap on a reference support
#'
#' Resamples subject vectors already carried by a typed aligned-map object to
#' quantify between-subject variability without recomputing the group basis or
#' correspondence.
#'
#' @param values_medoid An operator-bound [dkge_aligned_maps()] object returned
#'   by a fitted alignment application. Raw lists and objects from the public
#'   descriptive constructor are refused because they carry no certified
#'   operator application.
#' @param contrast Contrast name or index to bootstrap from `values_medoid`.
#' @param B Number of bootstrap replicates.
#' @param aggregate Aggregation function applied to the resampled subjects
#'   (`"mean"` or `"median"`).
#' @param weights Optional assertion of the immutable subject weights stored in
#'   `values_medoid`. If supplied, it must match that receipt exactly. Only used
#'   when `aggregate = "mean"`.
#' @param seed Optional random seed for reproducibility.
#' @param voxel_operator Optional matrix that maps medoid vectors to voxel space
#'   (columns = voxels). When supplied, summaries in voxel space are also
#'   returned.
#' @param return_samples Logical; when `TRUE` the matrix of bootstrap samples is
#'   returned in the output bundle.
#' @param allow_approximate_alignment Permit an explicitly labelled approximate
#'   aligned-map object. Descriptive/ineligible states are always refused.
#'
#' @return A list containing bootstrap summaries (`mean`, `sd`, `z`, confidence
#'   intervals), and optionally the raw bootstrap draws (reference-support and
#'   voxel space).
#' @export
#' @examples
#' \donttest{
#' toy <- dkge_sim_toy(
#'   factors = list(A = list(L = 2), B = list(L = 3)),
#'   active_terms = c("A", "B"), S = 5, P = 20, snr = 5
#' )
#' fit <- dkge_fit(toy$B_list, toy$X_list, toy$K, rank = 2)
#' # Bootstrap requires transport setup - example shows API
#' }
dkge_bootstrap_projected <- function(values_medoid,
                                     contrast = 1L,
                                     B = 1000L,
                                     aggregate = c("mean", "median"),
                                     weights = NULL,
                                     seed = NULL,
                                     voxel_operator = NULL,
                                     return_samples = TRUE,
                                     allow_approximate_alignment = FALSE) {
  if (!inherits(values_medoid, "dkge_aligned_maps")) {
    .dkge_abort(
      paste0(
        "`dkge_bootstrap_projected()` requires typed `dkge_aligned_maps`; ",
        "raw transported lists have unverified correspondence provenance."
      ),
      "dkge_alignment_ineligible_error"
    )
  }
  aligned_maps <- values_medoid
  .dkge_validate_aligned_maps(aligned_maps)
  .dkge_assert_alignment_eligible(
    aligned_maps,
    allow_approximate = allow_approximate_alignment
  )
  contrast_index <- if (is.character(contrast)) {
    match(contrast, aligned_maps$contrast_ids)
  } else {
    as.integer(contrast)
  }
  if (length(contrast_index) != 1L || is.na(contrast_index) ||
      contrast_index < 1L || contrast_index > length(aligned_maps$values)) {
    .dkge_abort("`contrast` does not identify one aligned map.",
                "dkge_aligned_maps_error")
  }
  Y <- aligned_maps$values[[contrast_index]]
  values_medoid <- lapply(seq_len(nrow(Y)), function(i) Y[i, ])
  aggregate <- match.arg(aggregate)
  Q <- length(values_medoid[[1]])
  S <- length(values_medoid)
  if (!all(vapply(values_medoid, length, integer(1)) == Q)) {
    stop("All medoid vectors must have the same length.")
  }
  if (!is.null(seed)) set.seed(seed)

  declared_weights <- as.numeric(aligned_maps$subject_weights)
  if (!is.null(weights) && !isTRUE(all.equal(
    as.numeric(weights), declared_weights, tolerance = 0
  ))) {
    .dkge_abort(
      "`weights` conflicts with the immutable aligned-map weighting receipt.",
      "dkge_alignment_weighting_error"
    )
  }
  weights <- declared_weights
  if (aggregate != "mean") {
    warning("Subject weights are ignored when aggregate != 'mean'.", call. = FALSE)
  }

  boot_medoid <- matrix(NA_real_, nrow = B, ncol = Q)
  if (aggregate == "mean" && is.null(weights)) {
    draws <- matrix(sample.int(S, size = B * S, replace = TRUE), nrow = B, ncol = S)
    boot_medoid <- matrix(0, B, Q)
    for (s in seq_len(S)) {
      boot_medoid <- boot_medoid + Y[draws[, s], , drop = FALSE]
    }
    boot_medoid <- boot_medoid / S
  } else if (aggregate == "mean") {
    draws <- matrix(sample.int(S, size = B * S, replace = TRUE), nrow = B, ncol = S)
    weight_mat <- matrix(weights[draws], nrow = B, ncol = S)
    boot_medoid <- matrix(0, B, Q)
    denom <- rowSums(weight_mat) + 1e-12
    for (s in seq_len(S)) {
      boot_medoid <- boot_medoid + (weight_mat[, s] * Y[draws[, s], , drop = FALSE])
    }
    boot_medoid <- boot_medoid / denom
  } else {
    for (b in seq_len(B)) {
      idx <- sample.int(S, size = S, replace = TRUE)
      boot_medoid[b, ] <- apply(Y[idx, , drop = FALSE], 2, stats::median)
    }
  }

  summary_medoid <- .dkge_bootstrap_summary(boot_medoid)

  voxel_res <- NULL
  if (!is.null(voxel_operator)) {
    voxel_operator <- as.matrix(voxel_operator)
    stopifnot(nrow(voxel_operator) == Q)
    boot_voxel <- boot_medoid %*% voxel_operator
    summary_voxel <- .dkge_bootstrap_summary(boot_voxel)
    voxel_res <- c(summary_voxel, list(boot = if (return_samples) boot_voxel else NULL))
  }

  list(
    B = B,
    medoid = c(summary_medoid, list(boot = if (return_samples) boot_medoid else NULL)),
    voxel = voxel_res,
    metadata = list(
      alignment_status = aligned_maps$eligibility$status,
      alignment_reason = aligned_maps$eligibility$reason,
      approximate_override = isTRUE(allow_approximate_alignment),
      support_hash = aligned_maps$support_hash,
      contrast = aligned_maps$contrast_ids[[contrast_index]],
      subject_weighting = aligned_maps$subject_weighting
    )
  )
}

#' Multiplier bootstrap in the design space (q-space)
#'
#' Reweights subject contributions with i.i.d. multiplier weights, recomputes
#' the tiny qxq eigendecomposition, and propagates contrasts to an identified
#' reference (and optionally voxel) support using a validated typed fitted
#' alignment. The cohort-trained truncated span makes this route approximate;
#' callers must opt in with `allow_approximate_alignment = TRUE`.
#' Fit-level MFA/energy weights act only on the reweighted pooled moment used to
#' estimate each latent basis. The returned group map has an equal-subject base
#' estimand: bootstrap multipliers resample subjects, but fit-level moment
#' weights are not reused as second-level aggregation weights.
#'
#' @inheritParams dkge_bootstrap_projected
#' @param fit A fitted `dkge` object.
#' @param contrasts Contrast specification accepted by [dkge_contrast()].
#' @param scheme Multiplier distribution (`"poisson"`, `"exp"`, or `"bayes"`).
#'   Poisson draws are conditioned on at least one positive subject multiplier,
#'   so every bootstrap replicate contains a non-empty resampled cohort.
#' @param ridge Optional ridge added to the reweighted compressed covariance.
#' @param align Logical; when `TRUE` the resampled bases are aligned to the
#'   baseline basis via K-Procrustes before contrasts are evaluated.
#' @param allow_reflection Passed to [dkge_procrustes_K()] when aligning bases.
#' @param transport_cache Required typed fitted alignment object produced by
#'   [dkge_prepare_alignment()]. Untyped/legacy caches and automatic full-fit
#'   correspondence are refused.
#' @param mapper Deprecated; automatic mapper fitting is no longer permitted.
#' @param centroids Deprecated; correspondence must be supplied through
#'   `transport_cache`.
#' @param sizes Deprecated compatibility argument; correspondence and masses
#'   must already be fixed in `transport_cache`.
#' @param medoid Deprecated compatibility argument.
#' @param voxel_operator Optional matrix mapping medoid values to voxels.
#' @param allow_approximate_alignment Permit a typed alignment labelled
#'   approximate. Descriptive/ineligible alignment is always refused.
#' @param ... Deprecated compatibility arguments; ignored.
#'
#' @return List containing per-contrast bootstrap summaries and the transport
#'   cache employed during resampling.
#' @export
dkge_bootstrap_qspace <- function(fit,
                                   contrasts,
                                   B = 1000L,
                                   scheme = c("poisson", "exp", "bayes"),
                                   ridge = 0,
                                   align = TRUE,
                                   allow_reflection = FALSE,
                                   seed = NULL,
                                   transport_cache = NULL,
                                   mapper = "sinkhorn",
                                   centroids = NULL,
                                   sizes = NULL,
                                   medoid = 1L,
                                   voxel_operator = NULL,
                                   allow_approximate_alignment = FALSE,
                                   ...) {
  stopifnot(inherits(fit, "dkge"))
  scheme <- match.arg(scheme)
  estimator_eligibility <- dkge_alignment_eligibility(
    feature_source = "geometry_only",
    estimator_source = "same_data_rank_truncated",
    source_verified = TRUE
  )
  estimator_eligibility$reason <- paste0(
    estimator_eligibility$reason,
    "; the multiplier bootstrap re-estimates a cohort-trained truncated span"
  )
  .dkge_assert_alignment_eligible(
    estimator_eligibility,
    allow_approximate = allow_approximate_alignment
  )
  if (!is.null(seed)) set.seed(seed)

  contrast_list <- .normalize_contrasts(contrasts, fit)
  n_contrasts <- length(contrast_list)
  if (!n_contrasts) stop("At least one contrast required.")

  S <- length(fit$Btil)
  q <- nrow(fit$U)
  r <- ncol(fit$U)

  Bmodel <- lapply(seq_along(fit$Btil), function(s) {
    .dkge_apply_fit_spatial(fit, fit$Btil[[s]], subject = s)
  })
  cache <- .dkge_bootstrap_prepare_cache(
    fit, transport_cache, allow_approximate_alignment,
    contrast_list = contrast_list
  )
  operators <- cache$operators
  Q <- cache$reference_support$n_locations

  voxel_operator <- if (is.null(voxel_operator)) NULL else as.matrix(voxel_operator)
  if (!is.null(voxel_operator) && nrow(voxel_operator) != Q) {
    stop("voxel_operator must have as many rows as medoid clusters.")
  }

  KBtil_t <- lapply(Bmodel, function(Bts) t(fit$K %*% Bts))
  Kctil_list <- lapply(contrast_list, function(c) {
    ctil <- backsolve(fit$R, c, transpose = FALSE)
    fit$K %*% ctil
  })

  boot_medoid <- lapply(seq_len(n_contrasts), function(i) matrix(NA_real_, B, Q))
  if (!is.null(voxel_operator)) {
    boot_voxel <- lapply(seq_len(n_contrasts), function(i) matrix(NA_real_, B, ncol(voxel_operator)))
  } else {
    boot_voxel <- NULL
  }

  weights_base <- as.numeric(fit$weights)
  contribs <- fit$contribs
  contrib_matrix <- vapply(contribs, function(M) as.numeric(M), numeric(q * q))
  kernel_support <- fit$kernel_support_projector %||%
    .dkge_kernel_geometry(fit$K)$support_projector

  for (b in seq_len(B)) {
    xi <- .dkge_bootstrap_multipliers(scheme, S)
    moment_coeff <- weights_base * xi
    Chat_vec <- contrib_matrix %*% moment_coeff
    Chat_b <- matrix(Chat_vec, q, q)
    if (ridge > 0) {
      Chat_b <- Chat_b + ridge * kernel_support
    }
    Chat_b <- (Chat_b + t(Chat_b)) / 2

    eig <- eigen(Chat_b, symmetric = TRUE)
    Vb <- eig$vectors[, seq_len(r), drop = FALSE]
    Ub <- fit$Kihalf %*% Vb
    Ub <- dkge_k_orthonormalize(Ub, fit$K)
    if (align) {
      pr <- dkge_procrustes_K(fit$U, Ub, fit$K, allow_reflection = allow_reflection)
      Ub <- pr$U_aligned
    }
    corr_diag <- diag(t(fit$U) %*% fit$K %*% Ub)
    corr_sign <- ifelse(corr_diag < 0, -1, 1)
    Ub <- sweep(Ub, 2, corr_sign, `*`)

    A_list <- lapply(KBtil_t, function(mat) mat %*% Ub)
    for (idx_con in seq_len(n_contrasts)) {
      alpha_b <- as.numeric(crossprod(Ub, Kctil_list[[idx_con]]))
      subject_maps <- matrix(0, S, Q)
      for (s in seq_len(S)) {
        v_s <- as.numeric(A_list[[s]] %*% alpha_b)
        subject_maps[s, ] <- as.numeric(t(operators[[s]]) %*% v_s)
      }
      boot_medoid[[idx_con]][b, ] <-
        .dkge_equal_subject_multiplier_mean(subject_maps, xi)
      if (!is.null(boot_voxel)) {
        boot_voxel[[idx_con]][b, ] <- boot_medoid[[idx_con]][b, ] %*% voxel_operator
      }
    }
  }

  summary_list <- vector("list", n_contrasts)
  names(summary_list) <- names(contrast_list)
  for (i in seq_len(n_contrasts)) {
    medoid_sum <- .dkge_bootstrap_summary(boot_medoid[[i]])
    if (!is.null(boot_voxel)) {
      voxel_sum <- .dkge_bootstrap_summary(boot_voxel[[i]])
      summary_list[[i]] <- c(medoid = list(medoid_sum),
                             voxel = list(voxel_sum),
                             list(boot_medoid = boot_medoid[[i]],
                                  boot_voxel = boot_voxel[[i]]))
    } else {
      summary_list[[i]] <- c(medoid_sum, list(boot = boot_medoid[[i]]))
    }
  }

  list(
    method = "qspace_multiplier",
    scheme = scheme,
    B = B,
    contrasts = names(contrast_list),
    summary = summary_list,
    cache = cache,
    metadata = list(
      alignment_status = cache$eligibility$status,
      alignment_reason = cache$eligibility$reason,
      estimator_status = estimator_eligibility$status,
      estimator_reason = estimator_eligibility$reason,
      approximate_override = isTRUE(allow_approximate_alignment),
      moment_weighting = list(
        method = fit$w_method %||% "stored_fit_weights",
        weights = unname(weights_base)
      ),
      group_weighting = "equal_subject",
      support_hash = cache$reference_support$structural_hash
    )
  )
}

#' Analytic first-order bootstrap in the design space
#'
#' Uses the stored full eigendecomposition to apply first-order perturbations
#' for each bootstrap draw. When the perturbation exceeds the validity region,
#' the method falls back to the exact multiplier bootstrap for that replicate.
#'
#' @inheritParams dkge_bootstrap_qspace
#' @param perturb_tol Maximum absolute coefficient tolerated in the eigenvector
#'   perturbation; larger changes trigger a fallback to the full eigensolve.
#' @param gap_tol Minimum eigen-gap tolerated (in absolute value) before
#'   triggering a fallback to the full eigensolve.
#' @details As in [dkge_bootstrap_qspace()], stored fit weights affect only the
#'   reweighted pooled moment. Subject maps are aggregated with an equal-subject
#'   base estimand using the bootstrap multipliers alone.
#'
#' @return Same structure as [dkge_bootstrap_qspace()] with additional metadata
#'   on the number of fallbacks used.
#' @export
dkge_bootstrap_analytic <- function(fit,
                                    contrasts,
                                    B = 1000L,
                                    scheme = c("poisson", "exp", "bayes"),
                                    ridge = 0,
                                    align = TRUE,
                                    allow_reflection = FALSE,
                                    seed = NULL,
                                    transport_cache = NULL,
                                    mapper = "sinkhorn",
                                    centroids = NULL,
                                    sizes = NULL,
                                    medoid = 1L,
                                    voxel_operator = NULL,
                                    perturb_tol = 0.2,
                                    gap_tol = 1e-6,
                                    allow_approximate_alignment = FALSE,
                                    ...) {
  stopifnot(inherits(fit, "dkge"))
  scheme <- match.arg(scheme)
  estimator_eligibility <- .dkge_rank_truncated_estimator_eligibility("analytic")
  .dkge_assert_alignment_eligible(
    estimator_eligibility,
    allow_approximate = allow_approximate_alignment
  )
  if (!is.null(seed)) set.seed(seed)

  solver_type <- fit$solver
  if (is.null(solver_type)) solver_type <- "pooled"
  if (!identical(solver_type, "pooled")) {
    warning("Analytic bootstrap requires solver = 'pooled'; falling back to q-space bootstrap.",
            call. = FALSE)
    return(dkge_bootstrap_qspace(fit, contrasts, B = B, scheme = scheme, ridge = ridge,
                                 align = align, allow_reflection = allow_reflection,
                                 seed = seed, transport_cache = transport_cache,
                                 mapper = mapper, centroids = centroids, sizes = sizes,
                                 medoid = medoid, voxel_operator = voxel_operator,
                                 allow_approximate_alignment = allow_approximate_alignment,
                                 ...))
  }

  if (is.null(fit$eig_vectors_full) || is.null(fit$eig_values_full)) {
    warning("Full eigendecomposition not stored on fit; falling back to q-space bootstrap.",
            call. = FALSE)
    return(dkge_bootstrap_qspace(fit, contrasts, B = B, scheme = scheme, ridge = ridge,
                                 align = align, allow_reflection = allow_reflection,
                                 seed = NULL, transport_cache = transport_cache,
                                 mapper = mapper, centroids = centroids, sizes = sizes,
                                 medoid = medoid, voxel_operator = voxel_operator,
                                 allow_approximate_alignment = allow_approximate_alignment,
                                 ...))
  }

  contrast_list <- .normalize_contrasts(contrasts, fit)
  n_contrasts <- length(contrast_list)
  if (!n_contrasts) stop("At least one contrast required.")

  S <- length(fit$Btil)
  q <- nrow(fit$U)
  r <- ncol(fit$U)

  Bmodel <- lapply(seq_along(fit$Btil), function(s) {
    .dkge_apply_fit_spatial(fit, fit$Btil[[s]], subject = s)
  })
  cache <- .dkge_bootstrap_prepare_cache(
    fit, transport_cache, allow_approximate_alignment,
    contrast_list = contrast_list
  )
  operators <- cache$operators
  Q <- cache$reference_support$n_locations

  voxel_operator <- if (is.null(voxel_operator)) NULL else as.matrix(voxel_operator)
  if (!is.null(voxel_operator) && nrow(voxel_operator) != Q) {
    stop("voxel_operator must have as many rows as medoid clusters.")
  }

  KBtil_t <- lapply(Bmodel, function(Bts) t(fit$K %*% Bts))
  Kctil_list <- lapply(contrast_list, function(c) {
    ctil <- backsolve(fit$R, c, transpose = FALSE)
    fit$K %*% ctil
  })

  boot_medoid <- lapply(seq_len(n_contrasts), function(i) matrix(NA_real_, B, Q))
  if (!is.null(voxel_operator)) {
    boot_voxel <- lapply(seq_len(n_contrasts), function(i) matrix(NA_real_, B, ncol(voxel_operator)))
  } else {
    boot_voxel <- NULL
  }

  weights_base <- as.numeric(fit$weights)
  contribs <- fit$contribs
  V_full <- fit$eig_vectors_full
  lambda_full <- fit$eig_values_full

  fallback_count <- 0L
  kernel_support <- fit$kernel_support_projector %||%
    .dkge_kernel_geometry(fit$K)$support_projector

  for (b in seq_len(B)) {
    xi <- .dkge_bootstrap_multipliers(scheme, S)

    delta_chat <- matrix(0, q, q)
    for (s in seq_len(S)) {
      delta_chat <- delta_chat + (xi[s] - 1) * weights_base[s] * contribs[[s]]
    }
    if (ridge > 0) {
      delta_chat <- delta_chat + ridge * kernel_support
    }
    delta_chat <- (delta_chat + t(delta_chat)) / 2

    Ub <- .dkge_bootstrap_analytic_basis(fit, V_full, lambda_full, delta_chat,
                                         gap_tol = gap_tol, perturb_tol = perturb_tol)
    if (is.null(Ub)) {
      fallback_count <- fallback_count + 1L
      Chat_b <- fit$Chat + delta_chat
      eig <- eigen(Chat_b, symmetric = TRUE)
      Ub <- fit$Kihalf %*% eig$vectors[, seq_len(r), drop = FALSE]
      Ub <- dkge_k_orthonormalize(Ub, fit$K)
    }

    if (align) {
      pr <- dkge_procrustes_K(fit$U, Ub, fit$K, allow_reflection = allow_reflection)
      Ub <- pr$U_aligned
    }
    corr_diag <- diag(t(fit$U) %*% fit$K %*% Ub)
    corr_sign <- ifelse(corr_diag < 0, -1, 1)
    Ub <- sweep(Ub, 2, corr_sign, `*`)

    subject_maps <- matrix(0, S, Q)
    # A_s depends only on the bootstrap basis, not the contrast; hoist it.
    A_list <- lapply(seq_len(S), function(s) KBtil_t[[s]] %*% Ub)
    for (idx_con in seq_len(n_contrasts)) {
      alpha_b <- as.numeric(crossprod(Ub, Kctil_list[[idx_con]]))
      for (s in seq_len(S)) {
        v_s <- as.numeric(A_list[[s]] %*% alpha_b)
        subject_maps[s, ] <- as.numeric(t(operators[[s]]) %*% v_s)
      }
      boot_medoid[[idx_con]][b, ] <-
        .dkge_equal_subject_multiplier_mean(subject_maps, xi)
      if (!is.null(boot_voxel)) {
        boot_voxel[[idx_con]][b, ] <- boot_medoid[[idx_con]][b, ] %*% voxel_operator
      }
    }
  }

  summary_list <- vector("list", n_contrasts)
  names(summary_list) <- names(contrast_list)
  for (i in seq_len(n_contrasts)) {
    medoid_sum <- .dkge_bootstrap_summary(boot_medoid[[i]])
    if (!is.null(boot_voxel)) {
      voxel_sum <- .dkge_bootstrap_summary(boot_voxel[[i]])
      summary_list[[i]] <- c(medoid = list(medoid_sum),
                             voxel = list(voxel_sum),
                             list(boot_medoid = boot_medoid[[i]],
                                  boot_voxel = boot_voxel[[i]]))
    } else {
      summary_list[[i]] <- c(medoid_sum, list(boot = boot_medoid[[i]]))
    }
  }

  list(
    method = "qspace_analytic",
    scheme = scheme,
    B = B,
    contrasts = names(contrast_list),
    summary = summary_list,
    cache = cache,
    fallbacks = fallback_count,
    metadata = list(
      alignment_status = cache$eligibility$status,
      alignment_reason = cache$eligibility$reason,
      estimator_status = estimator_eligibility$status,
      estimator_reason = estimator_eligibility$reason,
      approximate_override = isTRUE(allow_approximate_alignment),
      moment_weighting = list(
        method = fit$w_method %||% "stored_fit_weights",
        weights = unname(weights_base)
      ),
      group_weighting = "equal_subject",
      support_hash = cache$reference_support$structural_hash
    )
  )
}

# ---------------------------------------------------------------------------
# Internal helpers ----------------------------------------------------------
# ---------------------------------------------------------------------------

.dkge_bootstrap_summary <- function(mat) {
  mean_map <- colMeans(mat)
  sd_map <- apply(mat, 2, stats::sd)
  z_map <- mean_map / (sd_map + 1e-6)
  ci <- t(apply(mat, 2, stats::quantile, probs = c(0.025, 0.975)))
  list(mean = mean_map, sd = sd_map, z = z_map, ci = ci)
}

.dkge_bootstrap_multipliers <- function(scheme, S,
                                        poisson_draw = stats::rpois,
                                        max_redraws = 1000L) {
  switch(scheme,
         poisson = {
           for (attempt in seq_len(max_redraws)) {
             draw <- as.numeric(poisson_draw(S, lambda = 1))
             if (length(draw) != S || any(!is.finite(draw)) ||
                 any(draw < 0)) {
               .dkge_abort(
                 "Poisson multiplier generator returned an invalid draw.",
                 "dkge_bootstrap_weighting_error"
               )
             }
             if (sum(draw) > 0) return(draw)
           }
           .dkge_abort(
             "Could not generate a non-empty Poisson multiplier cohort.",
             "dkge_bootstrap_weighting_error"
           )
         },
         exp = stats::rexp(S, rate = 1),
         bayes = {
           u <- stats::runif(S)
           u / mean(u)
         })
}

.dkge_equal_subject_multiplier_mean <- function(subject_maps, multipliers) {
  subject_maps <- as.matrix(subject_maps)
  multipliers <- as.numeric(multipliers)
  if (nrow(subject_maps) != length(multipliers) ||
      any(!is.finite(subject_maps)) || any(!is.finite(multipliers)) ||
      any(multipliers < 0)) {
    .dkge_abort(
      "Bootstrap subject maps and multipliers are incompatible.",
      "dkge_bootstrap_weighting_error"
    )
  }
  total <- sum(multipliers)
  if (!is.finite(total) || total <= 0) {
    .dkge_abort(
      "Bootstrap multipliers must define a non-empty subject cohort.",
      "dkge_bootstrap_weighting_error"
    )
  }
  colSums(subject_maps * multipliers) / (total + 1e-12)
}

.dkge_bootstrap_prepare_cache <- function(fit, cache,
                                          allow_approximate_alignment = FALSE,
                                          contrast_list = NULL) {
  if (is.null(cache) || !inherits(cache, "dkge_fitted_alignment")) {
    .dkge_abort(
      paste0(
        "Q-space bootstrap requires a typed `dkge_fitted_alignment`; ",
        "automatic or legacy full-fit correspondence is not inferentially ",
        "eligible."
      ),
      "dkge_alignment_ineligible_error"
    )
  }
  .dkge_validate_fitted_alignment_object(cache)
  .dkge_assert_alignment_eligible(
    cache,
    allow_approximate = allow_approximate_alignment
  )
  expected_fit_binding <- .dkge_alignment_fit_binding(fit)
  cached_fit_binding <- cache$preprocessing$fit_binding %||% NULL
  if (is.null(cached_fit_binding) ||
      !identical(cached_fit_binding, expected_fit_binding)) {
    .dkge_abort(
      paste0(
        "Bootstrap fit does not match the fit bound into the typed fitted ",
        "alignment. Same-sized subject cohorts are not interchangeable."
      ),
      "dkge_alignment_cache_mismatch"
    )
  }
  cached_family <- cache$preprocessing$contrast_family_binding %||% NULL
  if (!is.null(cached_family)) {
    expected_family <- .dkge_alignment_contrast_family_binding(contrast_list)
    if (!identical(cached_family, expected_family)) {
      .dkge_abort(
        paste0(
          "Typed correspondence was fitted for a different contrast ",
          "family; its fixed plans cannot be reused for this bootstrap."
        ),
        "dkge_alignment_cache_mismatch"
      )
    }
  }
  fit_ids <- as.character(fit$subject_ids %||% seq_along(fit$Btil))
  cache_ids <- as.character(cache$subject_ids %||% seq_along(cache$operators))
  if (length(cache$operators) != length(fit_ids) ||
      !identical(cache_ids, fit_ids)) {
    .dkge_abort(
      "Bootstrap fit subjects do not match the fitted-alignment provenance.",
      "dkge_alignment_cache_mismatch"
    )
  }
  cache
}

.dkge_bootstrap_analytic_basis <- function(fit, V_full, lambda_full, delta_chat,
                                           gap_tol = 1e-6, perturb_tol = 0.2) {
  q <- nrow(V_full)
  r <- ncol(fit$U)
  H <- t(V_full) %*% delta_chat %*% V_full
  V_new <- matrix(0, q, r)
  for (j in seq_len(r)) {
    gaps <- lambda_full[j] - lambda_full
    gaps[j] <- NA
    if (any(abs(gaps) < gap_tol, na.rm = TRUE)) {
      return(NULL)
    }
    coeffs <- rep(0, q)
    coeffs[-j] <- H[-j, j] / gaps[-j]
    if (any(abs(coeffs[-j]) > perturb_tol, na.rm = TRUE)) {
      return(NULL)
    }
    V_new[, j] <- V_full[, j] + V_full %*% coeffs
  }
  V_ortho <- qr.Q(qr(V_new))
  fit$Kihalf %*% V_ortho
}
