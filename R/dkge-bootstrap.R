# dkge-bootstrap.R
# Fast bootstrap approximations for DKGE fits.

#' Subject-level projection bootstrap in medoid space
#'
#' Resamples transported subject vectors (already aligned in the medoid
#' parcellation) to quantify between-subject variability without recomputing the
#' group basis.
#'
#' @param values_medoid List of length `S` where each element is a numeric vector
#'   defined on the medoid parcellation (e.g. LOSO contrast values).
#' @param B Number of bootstrap replicates.
#' @param aggregate Aggregation function applied to the resampled subjects
#'   (`"mean"` or `"median"`).
#' @param weights Optional subject weights applied when computing the resampled
#'   mean. Only used when `aggregate = "mean"`.
#' @param seed Optional random seed for reproducibility.
#' @param voxel_operator Optional matrix that maps medoid vectors to voxel space
#'   (columns = voxels). When supplied, summaries in voxel space are also
#'   returned.
#' @param return_samples Logical; when `TRUE` the matrix of bootstrap samples is
#'   returned in the output bundle.
#'
#' @return A list containing bootstrap summaries (`mean`, `sd`, `z`, confidence
#'   intervals), and optionally the raw bootstrap draws (medoid and voxel space).
#' @export
dkge_bootstrap_projected <- function(values_medoid,
                                     B = 1000L,
                                     aggregate = c("mean", "median"),
                                     weights = NULL,
                                     seed = NULL,
                                     voxel_operator = NULL,
                                     return_samples = TRUE) {
  stopifnot(is.list(values_medoid), length(values_medoid) > 0)
  aggregate <- match.arg(aggregate)
  Q <- length(values_medoid[[1]])
  S <- length(values_medoid)
  if (!all(vapply(values_medoid, length, integer(1)) == Q)) {
    stop("All medoid vectors must have the same length.")
  }
  if (!is.null(seed)) set.seed(seed)

  Y <- do.call(rbind, values_medoid)
  if (!is.null(weights)) {
    stopifnot(length(weights) == S)
    weights <- as.numeric(weights)
    if (aggregate != "mean") {
      warning("Subject weights are ignored when aggregate != 'mean'.", call. = FALSE)
    }
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
    voxel = voxel_res
  )
}

#' Multiplier bootstrap in the design space (q-space)
#'
#' Reweights subject contributions with i.i.d. multiplier weights, recomputes
#' the tiny qxq eigendecomposition, and propagates contrasts to the medoid (and
#' optionally voxel) space using cached transport operators.
#'
#' @inheritParams dkge_bootstrap_projected
#' @param fit A fitted `dkge` object.
#' @param contrasts Contrast specification accepted by [dkge_contrast()].
#' @param scheme Multiplier distribution (`"poisson"`, `"exp"`, or `"bayes"`).
#' @param ridge Optional ridge added to the reweighted compressed covariance.
#' @param align Logical; when `TRUE` the resampled bases are aligned to the
#'   baseline basis via K-Procrustes before contrasts are evaluated.
#' @param allow_reflection Passed to [dkge_procrustes_K()] when aligning bases.
#' @param transport_cache Optional cache from [dkge_prepare_transport()].
#' @param mapper Mapper specification used when a cache is not supplied.
#' @param centroids Optional centroids overriding those stored on the fit.
#' @param sizes Optional list of cluster weights passed to
#'   [dkge_prepare_transport()].
#' @param medoid Medoid index used during transport.
#' @param voxel_operator Optional matrix mapping medoid values to voxels.
#' @param ... Additional arguments forwarded to [dkge_prepare_transport()] when
#'   the transport cache needs to be built.
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
                                   ...) {
  stopifnot(inherits(fit, "dkge"))
  scheme <- match.arg(scheme)
  if (!is.null(seed)) set.seed(seed)

  contrast_list <- .normalize_contrasts(contrasts, fit)
  n_contrasts <- length(contrast_list)
  if (!n_contrasts) stop("At least one contrast required.")

  S <- length(fit$Btil)
  q <- nrow(fit$U)
  r <- ncol(fit$U)

  loadings_ref <- lapply(fit$Btil, function(Bts) t(Bts) %*% fit$K %*% fit$U)
  cache <- .dkge_bootstrap_prepare_cache(fit, transport_cache, mapper, centroids,
                                         loadings_ref, sizes, medoid, ...)
  operators <- cache$operators
  Q <- ncol(operators[[medoid]])

  voxel_operator <- if (is.null(voxel_operator)) NULL else as.matrix(voxel_operator)
  if (!is.null(voxel_operator) && nrow(voxel_operator) != Q) {
    stop("voxel_operator must have as many rows as medoid clusters.")
  }

  KBtil_t <- lapply(fit$Btil, function(Bts) t(fit$K %*% Bts))
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

  for (b in seq_len(B)) {
    xi <- .dkge_bootstrap_multipliers(scheme, S)
    coeff <- weights_base * xi
    Chat_vec <- contrib_matrix %*% coeff
    Chat_b <- matrix(Chat_vec, q, q)
    if (ridge > 0) {
      diag(Chat_b) <- diag(Chat_b) + ridge
    }
    Chat_b <- (Chat_b + t(Chat_b)) / 2

    eig <- eigen(Chat_b, symmetric = TRUE)
    Vb <- eig$vectors[, seq_len(r), drop = FALSE]
    Ub <- fit$Kihalf %*% Vb
    if (align) {
      pr <- dkge_procrustes_K(fit$U, Ub, fit$K, allow_reflection = allow_reflection)
      Ub <- pr$U_aligned
    }
    corr_diag <- diag(t(fit$U) %*% fit$K %*% Ub)
    corr_sign <- ifelse(corr_diag < 0, -1, 1)
    Ub <- sweep(Ub, 2, corr_sign, `*`)

    A_list <- lapply(KBtil_t, function(mat) mat %*% Ub)
    weights_boot <- coeff
    w_sum <- sum(weights_boot)
    if (!is.finite(w_sum) || w_sum <= 0) w_sum <- 1

    for (idx_con in seq_len(n_contrasts)) {
      alpha_b <- as.numeric(crossprod(Ub, Kctil_list[[idx_con]]))
      subject_maps <- matrix(0, S, Q)
      for (s in seq_len(S)) {
        v_s <- as.numeric(A_list[[s]] %*% alpha_b)
        subject_maps[s, ] <- as.numeric(t(operators[[s]]) %*% v_s)
      }
      boot_medoid[[idx_con]][b, ] <- colSums(subject_maps * weights_boot) / (w_sum + 1e-12)
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
    cache = cache
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
                                    ...) {
  stopifnot(inherits(fit, "dkge"))
  scheme <- match.arg(scheme)
  if (!is.null(seed)) set.seed(seed)

  if (is.null(fit$eig_vectors_full) || is.null(fit$eig_values_full)) {
    warning("Full eigendecomposition not stored on fit; falling back to q-space bootstrap.",
            call. = FALSE)
    return(dkge_bootstrap_qspace(fit, contrasts, B = B, scheme = scheme, ridge = ridge,
                                 align = align, allow_reflection = allow_reflection,
                                 seed = NULL, transport_cache = transport_cache,
                                 mapper = mapper, centroids = centroids, sizes = sizes,
                                 medoid = medoid, voxel_operator = voxel_operator, ...))
  }

  contrast_list <- .normalize_contrasts(contrasts, fit)
  n_contrasts <- length(contrast_list)
  if (!n_contrasts) stop("At least one contrast required.")

  S <- length(fit$Btil)
  q <- nrow(fit$U)
  r <- ncol(fit$U)

  loadings_ref <- lapply(fit$Btil, function(Bts) t(Bts) %*% fit$K %*% fit$U)
  cache <- .dkge_bootstrap_prepare_cache(fit, transport_cache, mapper, centroids,
                                         loadings_ref, sizes, medoid, ...)
  operators <- cache$operators
  Q <- ncol(operators[[medoid]])

  voxel_operator <- if (is.null(voxel_operator)) NULL else as.matrix(voxel_operator)
  if (!is.null(voxel_operator) && nrow(voxel_operator) != Q) {
    stop("voxel_operator must have as many rows as medoid clusters.")
  }

  KBtil_t <- lapply(fit$Btil, function(Bts) t(fit$K %*% Bts))
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

  for (b in seq_len(B)) {
    xi <- .dkge_bootstrap_multipliers(scheme, S)

    delta_chat <- matrix(0, q, q)
    for (s in seq_len(S)) {
      delta_chat <- delta_chat + (xi[s] - 1) * weights_base[s] * contribs[[s]]
    }
    if (ridge > 0) {
      diag(delta_chat) <- diag(delta_chat) + ridge
    }
    delta_chat <- (delta_chat + t(delta_chat)) / 2

    Ub <- .dkge_bootstrap_analytic_basis(fit, V_full, lambda_full, delta_chat,
                                         gap_tol = gap_tol, perturb_tol = perturb_tol)
    if (is.null(Ub)) {
      fallback_count <- fallback_count + 1L
      Chat_b <- fit$Chat + delta_chat
      eig <- eigen(Chat_b, symmetric = TRUE)
      Ub <- fit$Kihalf %*% eig$vectors[, seq_len(r), drop = FALSE]
    }

    if (align) {
      pr <- dkge_procrustes_K(fit$U, Ub, fit$K, allow_reflection = allow_reflection)
      Ub <- pr$U_aligned
    }
    corr_diag <- diag(t(fit$U) %*% fit$K %*% Ub)
    corr_sign <- ifelse(corr_diag < 0, -1, 1)
    Ub <- sweep(Ub, 2, corr_sign, `*`)

    subject_maps <- matrix(0, S, Q)
    weights_boot <- weights_base * xi
    w_sum <- sum(weights_boot)
    if (!is.finite(w_sum) || w_sum <= 0) w_sum <- 1

    for (idx_con in seq_len(n_contrasts)) {
      alpha_b <- as.numeric(crossprod(Ub, Kctil_list[[idx_con]]))
      for (s in seq_len(S)) {
        A_s <- KBtil_t[[s]] %*% Ub
        v_s <- as.numeric(A_s %*% alpha_b)
        subject_maps[s, ] <- as.numeric(t(operators[[s]]) %*% v_s)
      }
      boot_medoid[[idx_con]][b, ] <- colSums(subject_maps * weights_boot) / (w_sum + 1e-12)
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
    fallbacks = fallback_count
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

.dkge_bootstrap_multipliers <- function(scheme, S) {
  switch(scheme,
         poisson = stats::rpois(S, lambda = 1),
         exp = stats::rexp(S, rate = 1),
         bayes = {
           u <- stats::runif(S)
           u / mean(u)
         })
}

.dkge_bootstrap_prepare_cache <- function(fit, cache, mapper, centroids,
                                          loadings, sizes, medoid, ...) {
  if (!is.null(cache)) {
    return(cache)
  }
  mapper <- if (is.null(mapper)) "sinkhorn" else mapper
  dkge_prepare_transport(fit,
                         centroids = centroids,
                         loadings = loadings,
                         sizes = sizes,
                         mapper = mapper,
                         medoid = medoid,
                         ...)
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

