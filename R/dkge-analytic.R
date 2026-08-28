# dkge-analytic.R
# Analytic (first-order) leave-one-out approximation for DKGE contrasts

.dkge_analytic_reason_levels <- c(
  "analytic",
  "pair_normalized_pooling",
  "covariance_aware_moment",
  "nonuniform_voxel_weights",
  "missing_full_decomposition",
  "dimension_mismatch",
  "eigengap",
  "perturbation_magnitude",
  "fallback"
)

.dkge_analytic_diagnostic <- function(reason,
                                      min_eigengap = NA_real_,
                                      max_abs_coeff = NA_real_,
                                      threshold_eigengap = NA_real_,
                                      threshold_coeff = NA_real_) {
  if (length(reason) != 1L || !reason %in% .dkge_analytic_reason_levels) {
    stop("Internal DKGE error: unknown analytic fallback reason.",
         call. = FALSE)
  }
  list(
    reason = reason,
    min_eigengap = as.numeric(min_eigengap),
    max_abs_coeff = as.numeric(max_abs_coeff),
    threshold_eigengap = as.numeric(threshold_eigengap),
    threshold_coeff = as.numeric(threshold_coeff)
  )
}

#' Analytic LOSO contrast using eigenvalue perturbation
#'
#' Approximates leave-one-subject-out contrast values using first-order
#' eigenvalue perturbation theory, avoiding full eigen-decomposition per subject.
#'
#' @param fit A `dkge` object from [dkge_fit()]
#' @param s Subject index (1-based) to leave out
#' @param contrasts Contrast vector in the original design basis (length q)
#' @param tol Numerical tolerance for determining when to fall back to full eigen
#' @param fallback If TRUE, fall back to full eigen when perturbation may be unstable
#' @param ridge Ridge regularization parameter (default 0)
#'
#' @return List with fields:
#'   - `v`: Cluster contrast values for subject s
#'   - `alpha`: Contrast coordinates in latent space
#'   - `basis`: Approximated held-out basis U^(-s)
#'   - `method`: "analytic" or "fallback" if full eigen was used
#'
#' @details
#' This function implements the first-order eigenvalue perturbation approximation
#' described in the paper. For the held-out compressed covariance:
#'
#' Delta_s = Chat^(-s) - Chat
#'
#' where `Chat^(-s)` is formed with the same fold-local subject-weight
#' normalization and shrinkage as exact LOSO. The eigenvalues and eigenvectors
#' are updated using:
#' - deltalambda_j = v_j^T Delta_s v_j (eigenvalue shift)
#' - deltav_j = Sum over k!=j of (v_k^T Delta_s v_j)/(lambda_j - lambda_k) v_k (eigenvector rotation)
#'
#' After the exact fold moment has been constructed, this replaces its O(q^3)
#' eigen-decomposition with an O(q^2 r) first-order eigensystem update, where r
#' is the rank. The approximation is accurate when:
#' 1. The fold perturbation norm is small relative to the fitted eigensystem
#' 2. Eigenvalue gaps are large (well-separated components)
#' 3. The resulting first-order rotation coefficients remain below the gate
#'
#' When these conditions are violated (detected via condition number or
#' eigenvalue gaps), the function can fall back to full eigen-decomposition.
#'
#' Fallback diagnostics use a closed reason vocabulary with this precedence:
#' `pair_normalized_pooling`, `covariance_aware_moment`, `nonuniform_voxel_weights`,
#' `missing_full_decomposition`, `dimension_mismatch`, `eigengap`, and
#' `perturbation_magnitude`. Structural reasons are checked before numerical
#' perturbation thresholds, so a large perturbation cannot mask the more basic
#' fact that the fitted moment does not support this approximation. The
#' q-by-q perturbation itself is exact for the fold pooling contract; only its
#' eigensystem update is first-order. A successful approximation reports
#' `analytic`.
#'
#' @examples
#' \donttest{
#' toy <- dkge_sim_toy(
#'   factors = list(A = list(L = 2), B = list(L = 3)),
#'   active_terms = c("A", "B"), S = 4, P = 20, snr = 5
#' )
#' fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 2)
#' c_vec <- c(1, -1, 0, 0, 0)
#' result <- dkge_analytic_loso(fit, s = 1, contrasts = c_vec)
#' result$method
#' }
#'
#' @references
#' Golub, G. H., & Van Loan, C. F. (2013). Matrix computations (4th ed.).
#' Stewart, G. W., & Sun, J. (1990). Matrix perturbation theory.
#'
#' @keywords internal
#' @export
dkge_analytic_loso <- function(fit, s, contrasts, tol = 1e-6, fallback = TRUE, ridge = 0) {
  stopifnot(inherits(fit, "dkge"), s >= 1L, s <= length(fit$Btil))
  .dkge_assert_crossfit_estimator_supported(
    fit, "Analytic LOSO contrast estimation"
  )
  S <- length(fit$Btil)
  q <- nrow(fit$U)
  r <- ncol(fit$U)
  stopifnot(length(contrasts) == q)

  nonlinear_pool <-
    !identical(fit$effect_weight_spec$method %||% "none", "none") ||
    !identical(fit$missingness %||% "none", "none")
  if (nonlinear_pool) {
    diag_info <- .dkge_analytic_diagnostic("pair_normalized_pooling")
    return(.dkge_analytic_fallback(fit, s, contrasts, ridge,
                                   reason = "pair_normalized_pooling",
                                   diagnostic = diag_info))
  }

  if (!identical(fit$debias %||% "none", "none")) {
    diag_info <- .dkge_analytic_diagnostic("covariance_aware_moment")
    return(.dkge_analytic_fallback(fit, s, contrasts, ridge,
                                   reason = "covariance_aware_moment",
                                   diagnostic = diag_info))
  }

  voxel_weight_payload <- fit$voxel_weights_subject %||% fit$voxel_weights
  if (!is.null(voxel_weight_payload)) {
    is_uniform_one <- function(weights) {
      weights <- as.numeric(weights)
      length(weights) > 0L && all(is.finite(weights)) &&
        isTRUE(all.equal(weights, rep(1, length(weights)), tolerance = 1e-6))
    }
    uniform <- if (is.list(voxel_weight_payload)) {
      length(voxel_weight_payload) == S &&
        all(vapply(voxel_weight_payload, is_uniform_one, logical(1)))
    } else {
      is_uniform_one(voxel_weight_payload)
    }
    if (!uniform) {
      diag_info <- .dkge_analytic_diagnostic("nonuniform_voxel_weights")
      return(.dkge_analytic_fallback(fit, s, contrasts, ridge,
                                     reason = "nonuniform_voxel_weights",
                                     diagnostic = diag_info))
    }
  }

  perturbation <- .dkge_analytic_fold_perturbation(fit, s, ridge = ridge)
  Delta_s <- perturbation$delta

  V_full <- fit$eig_vectors_full
  lambda_full <- fit$eig_values_full
  if (is.null(V_full) || is.null(lambda_full)) {
    diag_info <- .dkge_analytic_diagnostic("missing_full_decomposition")
    return(.dkge_analytic_fallback(fit, s, contrasts, ridge,
                                   reason = "missing_full_decomposition",
                                   diagnostic = diag_info))
  }

  if (ncol(V_full) != q || length(lambda_full) != q) {
    diag_info <- .dkge_analytic_diagnostic("dimension_mismatch")
    return(.dkge_analytic_fallback(fit, s, contrasts, ridge,
                                   reason = "dimension_mismatch",
                                   diagnostic = diag_info))
  }

  # Precompute exact fold-perturbation couplings H = V^T Delta_s V.
  H <- t(V_full) %*% Delta_s %*% V_full
  H <- (H + t(H)) / 2

  lambda_new <- lambda_full[seq_len(r)]
  V_new <- matrix(0, q, r)
  gap_tol <- max(tol, 1e-8)
  perturb_tol <- 0.1

  min_gap_observed <- Inf
  max_coeff_observed <- -Inf

  for (j in seq_len(r)) {
    v_j <- V_full[, j]
    delta_lambda_j <- H[j, j]
    lambda_new[j] <- lambda_full[j] + delta_lambda_j

    gaps <- lambda_full[j] - lambda_full
    gaps[j] <- NA
    min_gap_j <- suppressWarnings(min(abs(gaps), na.rm = TRUE))
    if (!is.finite(min_gap_j)) {
      min_gap_j <- NA_real_
    }
    if (is.finite(min_gap_j)) {
      min_gap_observed <- min(min_gap_observed, min_gap_j)
    }
    if (any(abs(gaps) < gap_tol, na.rm = TRUE)) {
      diag_info <- .dkge_analytic_diagnostic(
        "eigengap",
        min_eigengap = min_gap_j,
        max_abs_coeff = if (is.finite(max_coeff_observed)) {
          max_coeff_observed
        } else {
          NA_real_
        },
        threshold_eigengap = gap_tol,
        threshold_coeff = perturb_tol
      )
      return(.dkge_analytic_fallback(fit, s, contrasts, ridge,
                                     reason = "eigengap",
                                     diagnostic = diag_info))
    }
    coeffs <- rep(0, q)
    coeffs[-j] <- H[-j, j] / gaps[-j]
    max_coeff_j <- suppressWarnings(max(abs(coeffs[-j]), na.rm = TRUE))
    if (!is.finite(max_coeff_j)) {
      max_coeff_j <- NA_real_
    }
    if (!is.na(max_coeff_j)) {
      max_coeff_observed <- max(max_coeff_observed, max_coeff_j, na.rm = TRUE)
    }
    if (any(abs(coeffs[-j]) > perturb_tol, na.rm = TRUE)) {
      diag_info <- .dkge_analytic_diagnostic(
        "perturbation_magnitude",
        min_eigengap = if (is.finite(min_gap_observed)) {
          min_gap_observed
        } else {
          NA_real_
        },
        max_abs_coeff = max_coeff_j,
        threshold_eigengap = gap_tol,
        threshold_coeff = perturb_tol
      )
      return(.dkge_analytic_fallback(fit, s, contrasts, ridge,
                                     reason = "perturbation_magnitude",
                                     diagnostic = diag_info))
    }
    delta_v <- V_full %*% coeffs
    V_new[, j] <- v_j + delta_v
  }

  # Orthonormalize in Euclidean metric
  V_ortho <- qr.Q(qr(V_new))

  U_minus <- fit$Kihalf %*% V_ortho
  KU_minus <- fit$K %*% U_minus

  c_tilde <- backsolve(fit$R, contrasts, transpose = FALSE)
  alpha <- t(U_minus) %*% fit$K %*% c_tilde

  Bts <- fit$Btil[[s]]
  Bmodel <- .dkge_apply_fit_spatial(fit, Bts, subject = s)
  A_s <- t(Bmodel) %*% KU_minus
  v_s <- as.numeric(A_s %*% alpha)

  diag_info <- .dkge_analytic_diagnostic(
    "analytic",
    min_eigengap = if (is.finite(min_gap_observed)) {
      min_gap_observed
    } else {
      NA_real_
    },
    max_abs_coeff = if (is.finite(max_coeff_observed)) {
      max_coeff_observed
    } else {
      NA_real_
    },
    threshold_eigengap = gap_tol,
    threshold_coeff = perturb_tol
  )

  train_ids <- perturbation$train_ids
  preprocessing <- .dkge_alignment_preprocessing_receipt(
    fit, s, Bts, voxel_weights = NULL
  )
  receipt <- .dkge_make_alignment_receipt(
    fit = fit,
    subject = s,
    train_ids = train_ids,
    basis = U_minus,
    evals = lambda_new,
    loadings = A_s,
    preprocessing = preprocessing,
    alphas = list(contrast1 = as.numeric(alpha)),
    fold_index = s,
    holdout = s,
    subject_weights = perturbation$subject_weights,
    subject_weight_source = perturbation$subject_weight_source,
    method = "analytic",
    eligible = FALSE,
    eligibility_reason = "analytic_basis_is_approximate"
  )

  list(
    v = v_s,
    alpha = alpha,
    basis = U_minus,
    evals = lambda_new,
    loadings = A_s,
    method = "analytic",
    diagnostic = diag_info,
    alignment_receipt = receipt
  )
}

#' Exact fold-local q-space perturbation for analytic LOSO
#'
#' The subject weights in a held-out fold are normalized and shrunken over the
#' training cohort. Consequently, the perturbation is generally not
#' `-fit$weights[s] * fit$contribs[[s]]`. This helper forms the exact target
#' used by LOSO before the analytic eigensystem approximation is applied.
#'
#' @keywords internal
#' @noRd
.dkge_analytic_fold_perturbation <- function(fit, s, ridge = 0) {
  S <- length(fit$Btil)
  train_ids <- setdiff(seq_len(S), s)
  ctx <- .dkge_fold_weight_context(
    fit,
    train_ids,
    ridge = ridge,
    missingness = fit$missingness %||% "none",
    miss_args = fit$miss_args %||% list()
  )
  fold_chat <- as.matrix(ctx$Chat)
  full_chat <- as.matrix(fit$Chat)
  if (!identical(dim(fold_chat), dim(full_chat)) ||
      any(!is.finite(fold_chat)) || any(!is.finite(full_chat))) {
    stop("Analytic LOSO requires finite conformable full and fold moments.",
         call. = FALSE)
  }
  delta <- fold_chat - full_chat
  list(
    delta = (delta + t(delta)) / 2,
    fold_chat = (fold_chat + t(fold_chat)) / 2,
    train_ids = train_ids,
    subject_weights = ctx$subject_weights,
    subject_weight_source = ctx$subject_weight_source,
    weight_evaluation = ctx$weights
  )
}


#' Fallback to full eigen-decomposition
#'
#' @param fit dkge object
#' @param s subject index
#' @param contrasts contrast vector
#' @param ridge ridge parameter
#' @return Same as dkge_loso_contrast but with method="fallback"
#' @keywords internal
#' @noRd
.dkge_analytic_fallback <- function(fit, s, contrasts, ridge = 0,
                                    reason = "fallback",
                                    diagnostic = NULL) {
  result <- dkge_loso_contrast(fit, s, contrasts, ridge)
  result$method <- "fallback"
  result$alignment_receipt$estimation_method <- "analytic_fallback_exact"
  result$alignment_receipt$inference$eligible <- TRUE
  result$alignment_receipt$inference$reason <- "exact_loso_fallback"
  diag_out <- diagnostic %||% .dkge_analytic_diagnostic(reason)
  if (!reason %in% .dkge_analytic_reason_levels) {
    stop("Internal DKGE error: unknown analytic fallback reason.",
         call. = FALSE)
  }
  diag_out$reason <- reason
  result$diagnostic <- diag_out
  result
}

#' Gram-Schmidt orthonormalization in K-metric
#'
#' @param V Matrix with columns to orthonormalize
#' @param K Metric matrix
#' @return Orthonormalized matrix where V^T K V = I
#' @keywords internal
#' @noRd
.gram_schmidt_k_metric <- function(V, K) {
  n <- nrow(V)
  r <- ncol(V)
  Q <- matrix(0, n, r)

  for (j in seq_len(r)) {
    v <- V[, j]
    # Remove projections onto previous vectors
    if (j > 1) {
      for (i in seq_len(j - 1)) {
        q_i <- Q[, i]
        # Projection in K-metric: <v, q_i>_K / <q_i, q_i>_K
        # Since Q is already K-orthonormal, <q_i, q_i>_K = 1
        proj <- as.numeric(t(v) %*% K %*% q_i)
        v <- v - proj * q_i
      }
    }
    # Normalize in K-metric
    norm_k <- sqrt(as.numeric(t(v) %*% K %*% v))
    if (norm_k > 1e-10) {
      Q[, j] <- v / norm_k
    } else {
      # Vector became zero - shouldn't happen with stable input
      warning("Gram-Schmidt: vector ", j, " became null")
      Q[, j] <- 0
    }
  }
  Q
}

#' Analytic cross-fitting for multiple contrasts
#'
#' Internal implementation called by dkge_contrast() for method="analytic".
#' Uses first-order perturbation theory to approximate LOSO contrasts.
#'
#' @param fit dkge object
#' @param contrast_list List of normalized contrasts
#' @param ridge Ridge parameter added to the exact fold moment before the
#'   first-order eigensystem update (and used by an exact fallback).
#' @param parallel Logical; enables future.apply-based parallelism for
#'   per-subject computations (requires future.apply)
#' @param verbose Print progress
#' @param tol Tolerance for perturbation stability
#' @param fallback Allow fallback to full eigen when unstable
#' @param ... Additional arguments
#' @return List with values, metadata, etc.
#' @keywords internal
.dkge_contrast_analytic_impl <- function(fit, contrast_list, ridge,
                                   parallel, verbose,
                                   tol = 1e-6, fallback = TRUE, align = TRUE, ...) {
  S <- length(fit$Btil)
  n_contrasts <- length(contrast_list)

  verbose_flag <- .dkge_verbose(verbose)

  if (verbose_flag) {
    message(sprintf("Computing %d contrast(s) via analytic LOSO for %d subjects",
                   n_contrasts, S))
  }

  values <- vector("list", n_contrasts)
  names(values) <- names(contrast_list)
  bases <- vector("list", S)
  alphas <- vector("list", n_contrasts)
  methods_list <- vector("list", n_contrasts)
  diagnostics_list <- vector("list", n_contrasts)
  receipt_sets <- vector("list", n_contrasts)
  subject_indices <- seq_len(S)

  for (i in seq_along(contrast_list)) {
    subject_results <- .dkge_apply(
      subject_indices,
      function(s) {
        if (!parallel && verbose_flag && s %% 10 == 0) {
          message(sprintf("  Subject %d/%d", s, S))
        }
        res <- dkge_analytic_loso(fit, s, contrast_list[[i]],
                                  tol = tol, fallback = fallback, ridge = ridge)
        list(
          value = res$v,
          alpha = as.numeric(res$alpha),
          method = res$method,
          basis = res$basis,
          diagnostic = res$diagnostic %||% list(),
          receipt = res$alignment_receipt
        )
      },
      parallel = parallel
    )

    values[[i]] <- lapply(subject_results, `[[`, "value")
    alpha_mat <- do.call(rbind, lapply(subject_results, `[[`, "alpha"))
    rownames(alpha_mat) <- NULL
    alphas[[i]] <- alpha_mat
    methods_list[[i]] <- vapply(subject_results, function(x) x$method, character(1))
    diagnostics_list[[i]] <- lapply(subject_results, `[[`, "diagnostic")
    receipt_sets[[i]] <- lapply(subject_results, `[[`, "receipt")
    names(receipt_sets[[i]]) <- as.character(fit$subject_ids %||% subject_indices)

    if (i == 1) {
      bases <- lapply(subject_results, `[[`, "basis")
    }
  }

  aligned_bases <- bases
  rotations <- vector("list", length(bases))
  consensus <- NULL
  procrustes <- NULL
  if (align && length(bases) > 0) {
    align_obj <- dkge_align_bases_K(bases, fit$K, allow_reflection = TRUE)
    aligned_bases <- align_obj$U_aligned
    rotations <- align_obj$R
    weights <- fit$weights %||% rep(1, length(bases))
    consensus <- dkge_consensus_basis_K(bases, fit$K,
                                        weights = weights,
                                        allow_reflection = TRUE)
    procrustes <- list(alignment = align_obj, consensus = consensus)
  }

  if (verbose_flag) {
    fallback_counts <- vapply(methods_list, function(m) sum(m == "fallback"), integer(1))
    if (any(fallback_counts > 0)) {
      message(sprintf("Fallback used for %d/%d subjects (max %.1f%% across contrasts)",
                      sum(fallback_counts), S * n_contrasts,
                      max(100 * fallback_counts / S)))
    }
  }

  fallback_rates <- vapply(methods_list, function(m) mean(m == "fallback"), numeric(1))

  subject_ids <- fit$subject_ids %||% seq_len(S)
  fallback_rows <- list()
  for (i in seq_along(diagnostics_list)) {
    contrast_name <- names(contrast_list)[i]
    if (is.null(contrast_name) || !nzchar(contrast_name)) {
      contrast_name <- paste0("contrast", i)
    }
    diag_entries <- diagnostics_list[[i]]
    for (s in seq_along(diag_entries)) {
      diag <- diag_entries[[s]] %||% list()
      fallback_rows[[length(fallback_rows) + 1L]] <- data.frame(
        contrast = contrast_name,
        subject = subject_ids[s],
        reason = diag$reason %||% NA_character_,
        min_eigengap = diag$min_eigengap %||% NA_real_,
        max_abs_coeff = diag$max_abs_coeff %||% NA_real_,
        threshold_eigengap = diag$threshold_eigengap %||% NA_real_,
        threshold_coeff = diag$threshold_coeff %||% NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  fallback_detail <- if (length(fallback_rows)) {
    do.call(rbind, fallback_rows)
  } else {
    data.frame(
      contrast = character(0),
      subject = character(0),
      reason = character(0),
      min_eigengap = numeric(0),
      max_abs_coeff = numeric(0),
      threshold_eigengap = numeric(0),
      threshold_coeff = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  list(
    values = values,
    method = "analytic",
    contrasts = contrast_list,
    metadata = list(
      bases = bases,
      aligned_bases = aligned_bases,
      rotations = rotations,
      alphas = alphas,
      methods = methods_list,
      diagnostics = diagnostics_list,
      alignment_receipts = .dkge_merge_alignment_receipt_sets(
        receipt_sets,
        names(contrast_list) %||% paste0("contrast", seq_along(contrast_list)),
        method = "analytic"
      ),
      tol = tol,
      fallback = fallback,
      fallback_rates = fallback_rates,
      fallback_detail = fallback_detail,
      procrustes = procrustes
    )
  )
}
