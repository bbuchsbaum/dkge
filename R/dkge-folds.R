# dkge-folds.R
# Shared fold-building helpers for LOSO/K-fold cross-fitting.

.dkge_subject_loader_weights <- function(loader_weights, Bts) {
  if (is.null(loader_weights) || length(loader_weights) == 0L) {
    return(NULL)
  }
  w_s <- as.numeric(loader_weights)
  if (length(w_s) != ncol(Bts)) {
    w_s <- rep(w_s, length.out = ncol(Bts))
  }
  w_s
}

#' Build held-out fold bases and loaders
#'
#' Internal utility that re-computes DKGE bases for a collection of held-out
#' subject sets (folds) and optionally caches subject-specific projection
#' loaders. The resulting structure can be reused by contrast, classification,
#' or other cross-fitting modules to avoid duplicating eigen-solves.
#'
#' @param fit dkge object returned by [dkge()] or [dkge_fit()].
#' @param assignments List of integer vectors identifying the subjects held out
#'   in each fold. Subject indices are 1-based.
#' @param ridge Optional ridge added to the held-out Chat matrix before the
#'   eigen decomposition.
#' @param align Logical; when `TRUE`, compute K-orthogonal Procrustes alignment
#'   and a consensus basis across folds for interpretability.
#' @param loader_scope Either `"heldout"` (default) to cache loaders only for the
#'   subjects in each fold, or `"all"` to cache loaders for every subject under
#'   every fold (useful when the caller needs training-set projections).
#' @param verbose Logical; emit progress messages when `TRUE`.
#'
#' @return List with fields used by downstream consumers:
#'   - `folds`: list per fold containing the held-out subjects, raw and aligned
#'     bases, eigenvalues, rotation, and cached loaders.
#'   - `assignments`: original assignment list used to build the folds.
#'   - `align`: whether alignment was requested.
#'   - `consensus`: consensus K-orthogonal basis (when `align = TRUE`).
#'   - `loader_scope`: scope used for caching loaders.
#' @keywords internal
#' @noRd
.dkge_build_fold_bases <- function(fit,
                                   assignments,
                                   ridge = 0,
                                   align = TRUE,
                                   loader_scope = c("heldout", "all"),
                                   verbose = FALSE,
                                   weights = NULL,
                                   missingness = NULL,
                                   miss_args = NULL) {
  stopifnot(inherits(fit, "dkge"))
  stopifnot(is.list(assignments), length(assignments) >= 1)
  loader_scope <- match.arg(loader_scope)
  missingness <- missingness %||% fit$missingness %||% "none"
  missingness <- match.arg(missingness, c("none", "rescale", "mask", "shrink"))
  miss_args <- miss_args %||% fit$miss_args %||% list()

  S <- length(fit$Btil)
  q <- nrow(fit$U)
  r <- ncol(fit$U)
  subject_ids <- fit$subject_ids %||% seq_len(S)

  assignments <- lapply(assignments, function(idx) {
    idx <- sort(unique(as.integer(idx)))
    if (length(idx) == 0L) {
      stop("Each fold must hold out at least one subject.")
    }
    if (any(idx < 1L) || any(idx > S)) {
      stop("Fold assignments contain invalid subject indices.")
    }
    idx
  })

  n_folds <- length(assignments)
  verbose_flag <- .dkge_verbose(verbose)

  weight_spec <- weights %||% fit$weight_spec %||% dkge_weights(adapt = "none")
  stopifnot(inherits(weight_spec, "dkge_weights"))
  fold_bases <- vector("list", n_folds)
  fold_evals <- vector("list", n_folds)
  fold_loaders <- vector("list", n_folds)
  fold_weight_info <- vector("list", n_folds)
  fold_pair_counts <- vector("list", n_folds)
  fold_contexts <- vector("list", n_folds)
  recycled_subjects <- character(0)

  for (fold_idx in seq_len(n_folds)) {
    holdout <- assignments[[fold_idx]]
    train_ids <- setdiff(seq_len(S), holdout)

    if (verbose_flag) {
      message(sprintf("Building fold %d/%d (holdout: %s)",
                      fold_idx, n_folds,
                      paste(subject_ids[holdout], collapse = ", ")))
    }

    ctx <- .dkge_fold_weight_context(fit,
                                     train_ids,
                                     weight_spec,
                                     ridge = ridge,
                                     missingness = missingness,
                                     miss_args = miss_args)
    Chat_minus <- ctx$Chat
    weight_eval <- ctx$weights

    if (!is.null(ctx$fit)) {
      eig_values <- ctx$fit$eig_values_full
      U_fold <- ctx$fit$U
      if (ncol(U_fold) != r) {
        stop(sprintf(
          "Fold %d retained %d components but the reference fit has %d.",
          fold_idx, ncol(U_fold), r
        ), call. = FALSE)
      }
    } else {
      eig_fold <- eigen(Chat_minus, symmetric = TRUE)
      eig_values <- eig_fold$values
      U_fold <- fit$Kihalf %*% eig_fold$vectors[, seq_len(r), drop = FALSE]
    }

    fold_bases[[fold_idx]] <- U_fold
    fold_evals[[fold_idx]] <- eig_values

    loader_weights <- weight_eval$total

    subject_scope <- if (loader_scope == "all") seq_len(S) else holdout
    loader_list <- vector("list", length(subject_scope))
    names(loader_list) <- as.character(subject_scope)

    for (j in seq_along(subject_scope)) {
      s <- subject_scope[[j]]
      Bts <- if (!is.null(ctx$R)) {
        if (identical(fit$effect_scaling %||% "pooled_design", "none")) {
          fit$Braw[[s]]
        } else {
          t(ctx$R) %*% fit$Braw[[s]]
        }
      } else {
        fit$Btil[[s]]
      }
      w_s <- .dkge_subject_loader_weights(loader_weights, Bts)
      if (!is.null(loader_weights) && length(loader_weights) != ncol(Bts)) {
        if (length(loader_weights) > 1L) {
          recycled_subjects <- unique(c(recycled_subjects, subject_ids[s]))
        }
      }
      Bw <- if (is.null(w_s) || length(w_s) == 0L) Bts else {
        sweep(Bts, 2L, sqrt(pmax(w_s, 0)), "*")
      }
      A_s <- .dkge_basis_loadings(
        Bts, U_fold, fit$K, input_scale = "standardized",
        voxel_weights = w_s
      )
      Y_s <- Bw %*% A_s
      loader_list[[j]] <- list(
        subject = s,
        A = A_s,
        Y = Y_s,
        n_cluster = ncol(Bts)
      )
    }

    fold_loaders[[fold_idx]] <- loader_list
    fold_weight_info[[fold_idx]] <- list(
      prior = weight_eval$prior,
      adapt = weight_eval$adapt,
      total = weight_eval$total,
      total_subject = {
        ws <- weight_eval$total_subject
        if (!is.null(ws)) {
          names(ws) <- as.character(train_ids)
        }
        ws
      },
      spec = weight_spec,
      w_prior = weight_eval$prior,
      w_adapt = weight_eval$adapt,
      w_total = weight_eval$total,
      w_total_subject = weight_eval$total_subject,
      weight_spec = weight_spec
    )
    fold_pair_counts[fold_idx] <- list(ctx$pair_counts)
    fold_contexts[[fold_idx]] <- ctx
  }

  aligned_bases <- fold_bases
  rotations <- vector("list", n_folds)
  consensus <- NULL
  alignment <- NULL

  if (align && n_folds > 0) {
    align_obj <- dkge_align_bases_K(fold_bases, fit$K, allow_reflection = FALSE)
    aligned_bases <- align_obj$U_aligned
    rotations <- align_obj$R
    alignment <- align_obj

    fold_weights <- vapply(assignments, length, numeric(1))
    if (any(is.na(fold_weights)) || sum(fold_weights) <= 0) {
      fold_weights <- rep(1, n_folds)
    }
    consensus <- dkge_consensus_basis_K(fold_bases, fit$K,
                                        weights = fold_weights,
                                        allow_reflection = FALSE)
  }

  folds <- vector("list", n_folds)
  for (fold_idx in seq_len(n_folds)) {
    ctx <- fold_contexts[[fold_idx]]
    folds[[fold_idx]] <- list(
      index = fold_idx,
      subjects = assignments[[fold_idx]],
      basis = fold_bases[[fold_idx]],
      basis_aligned = aligned_bases[[fold_idx]],
      rotation = rotations[[fold_idx]],
      evals = fold_evals[[fold_idx]],
      loaders = fold_loaders[[fold_idx]],
      weights = fold_weight_info[[fold_idx]],
      U_minus = fold_bases[[fold_idx]],
      D_minus = fold_evals[[fold_idx]],
      pair_counts = fold_pair_counts[[fold_idx]],
      pair_weight = ctx$pair_weight,
      pair_ess = ctx$pair_ess,
      effect_precision = ctx$effect_precision,
      effect_precision_diagnostics = ctx$effect_precision_diagnostics,
      R = ctx$R %||% fit$R,
      ruler_source = ctx$ruler_source %||% fit$ruler_source,
      estimator = ctx$estimator %||% list(
        effect_weights = fit$effect_weight_spec$method %||% "none",
        debias = fit$debias %||% "none",
        missingness = missingness,
        effect_scaling = fit$effect_scaling %||% "pooled_design",
        train_subjects = subject_ids[setdiff(seq_len(S), assignments[[fold_idx]])]
      ),
      missingness = missingness,
      miss_args = miss_args
    )
  }

  list(
    folds = folds,
    assignments = assignments,
    align = align,
    consensus = consensus,
    alignment = alignment,
    loader_scope = loader_scope,
    weight_spec = weight_spec,
    missingness = missingness,
    miss_args = miss_args
  ) -> result

  if (length(recycled_subjects)) {
    warning(sprintf("Per-subject voxel weights recycled for: %s",
                    paste(recycled_subjects, collapse = ", ")))
    attr(result, "recycled_weights_subjects") <- recycled_subjects
  }
  result
}

#' Build fold loaders using the global (full-data) basis
#'
#' For `mode = "cell"` classification, every subject is projected onto the
#' global `fit$U` rather than fold-specific leave-one-out bases.  This helper
#' produces the same list structure as [.dkge_build_fold_bases()] so that the
#' downstream CV loop can use it without modification.
#'
#' @param fit dkge object.
#' @param assignments List of integer vectors (one per fold, subjects held out).
#' @keywords internal
#' @noRd
.dkge_build_global_fold_loaders <- function(fit, assignments) {
  stopifnot(inherits(fit, "dkge"))
  S <- length(fit$Btil)
  U_global <- fit$U
  n_folds <- length(assignments)

  # Build one loader per subject — same for every fold
  loader_template <- vector("list", S)
  names(loader_template) <- as.character(seq_len(S))
  for (s in seq_len(S)) {
    Bts <- fit$Btil[[s]]
    A_s <- .dkge_basis_loadings(Bts, U_global, fit$K,
                                input_scale = "standardized")
    Y_s <- Bts %*% A_s
    loader_template[[s]] <- list(
      subject = s,
      A = A_s,
      Y = Y_s,
      n_cluster = ncol(Bts)
    )
  }

  folds <- vector("list", n_folds)
  for (fold_idx in seq_len(n_folds)) {
    holdout <- sort(unique(as.integer(assignments[[fold_idx]])))
    folds[[fold_idx]] <- list(
      index = fold_idx,
      subjects = holdout,
      basis = U_global,
      basis_aligned = U_global,
      rotation = NULL,
      evals = fit$evals,
      loaders = loader_template,
      weights = NULL,
      U_minus = U_global,
      D_minus = fit$evals,
      pair_counts = NULL,
      missingness = "none",
      miss_args = list()
    )
  }

  list(
    folds = folds,
    assignments = assignments,
    align = FALSE,
    consensus = NULL,
    alignment = NULL,
    loader_scope = "all",
    weight_spec = fit$weight_spec %||% dkge_weights(adapt = "none"),
    missingness = "none",
    miss_args = list()
  )
}

#' @noRd
.dkge_normalize_folds <- function(folds, fit) {
  S <- length(fit$Btil)
  if (is.null(folds)) {
    return(list(assignments = lapply(seq_len(S), function(s) s), folds = NULL))
  }
  if (is.numeric(folds) && length(folds) == 1) {
    fold_obj <- dkge_define_folds(fit, type = "subject", k = folds)
  } else if (inherits(folds, "dkge_folds")) {
    fold_obj <- folds
  } else {
    fold_obj <- as_dkge_folds(folds, fit_or_data = fit)
    if (!inherits(fold_obj, "dkge_folds")) {
      stop("folds must be an integer k or convertible via as_dkge_folds().", call. = FALSE)
    }
  }
  list(assignments = fold_obj$assignments, folds = fold_obj)
}

#' @noRd
.dkge_subset_refit_spec <- function(fit, train_ids, weight_spec = NULL,
                                    ridge = 0, missingness = NULL,
                                    miss_args = NULL) {
  S <- length(fit$Btil)
  spec <- fit$refit_spec %||% list()
  weights <- weight_spec %||% spec$weights %||% fit$weight_spec %||%
    dkge_weights(adapt = "none")
  if (!is.null(weights$B_list2) && is.list(weights$B_list2) &&
      length(weights$B_list2) == S) {
    weights$B_list2 <- weights$B_list2[train_ids]
  }
  jd_mask <- spec$jd_mask
  if (is.list(jd_mask) && length(jd_mask) == S) {
    jd_mask <- jd_mask[train_ids]
  }
  list(
    w_method = spec$w_method %||% fit$w_method %||% "none",
    w_tau = spec$w_tau %||% fit$w_tau %||% 0.3,
    ridge = (spec$ridge %||% fit$ridge_input %||% 0) + ridge,
    rank = min(spec$rank %||% ncol(fit$U), nrow(fit$K)),
    effect_scaling = spec$effect_scaling %||% fit$effect_scaling %||%
      "pooled_design",
    cpca_blocks = spec$cpca_blocks,
    cpca_T = spec$cpca_T,
    cpca_part = spec$cpca_part %||% fit$cpca$part %||% "none",
    cpca_ridge = spec$cpca_ridge %||% fit$cpca$ridge %||% 0,
    weights = weights,
    effect_weights = spec$effect_weights %||% fit$effect_weight_spec %||%
      dkge_effect_weights("none"),
    debias = spec$debias %||% fit$debias %||% "none",
    solver = spec$solver %||% fit$solver %||% "pooled",
    missingness = missingness %||% spec$missingness %||%
      fit$missingness %||% "none",
    miss_args = miss_args %||% spec$miss_args %||% fit$miss_args %||% list(),
    jd_control = spec$jd_control %||% dkge_jd_control(),
    jd_mask = jd_mask,
    jd_init = spec$jd_init
  )
}

#' Refit the complete estimator on a training-subject multiset
#'
#' Subject indices may repeat, which gives bootstrap multiplicities exactly the
#' same semantics as literal subject duplication. No held-out subject-derived
#' ruler, reliability, or adaptive-weight quantity enters this fit.
#'
#' @keywords internal
#' @noRd
.dkge_refit_training_subjects <- function(fit,
                                           train_ids,
                                           weight_spec = NULL,
                                           ridge = 0,
                                           missingness = NULL,
                                           miss_args = NULL) {
  stopifnot(inherits(fit, "dkge"))
  S <- length(fit$Btil)
  train_ids <- as.integer(train_ids)
  if (length(train_ids) < 1L || anyNA(train_ids) ||
      any(train_ids < 1L | train_ids > S)) {
    stop("A training refit requires at least one valid subject index.",
         call. = FALSE)
  }
  if (is.null(fit$subjects) || length(fit$subjects) != S) {
    stop("This fit does not retain the subject records required for leakage-free refitting.",
         call. = FALSE)
  }

  subjects <- lapply(seq_along(train_ids), function(j) {
    subject <- fit$subjects[[train_ids[[j]]]]
    subject$id <- paste0(subject$id %||% paste0("subject", train_ids[[j]]),
                         "__refit", j)
    subject
  })
  spec <- .dkge_subset_refit_spec(
    fit, train_ids, weight_spec = weight_spec, ridge = ridge,
    missingness = missingness, miss_args = miss_args
  )
  data <- if (length(subjects) >= 2L) {
    dkge_data(subjects)
  } else {
    subject <- subjects[[1L]]
    provenance <- .dkge_full_coverage_provenance(subjects)
    structure(list(
      subjects = subjects,
      betas = list(subject$beta),
      designs = list(subject$design),
      omega = list(subject$omega),
      subject_ids = subject$id,
      effects = subject$effects,
      q = length(subject$effects),
      n_subjects = 1L,
      cluster_ids = list(subject$cluster_ids),
      observed_rows = list(subject$observed_rows),
      effect_n = list(subject$effect_n),
      effect_precision = list(subject$effect_precision),
      effect_noise_cov = list(subject$effect_noise_cov),
      residual_variance = list(subject$residual_variance),
      residual_df = list(subject$residual_df),
      noise_trace = list(subject$noise_trace),
      noise_trace_scope = list(subject$noise_trace_scope),
      residual_sum_squares = list(subject$residual_sum_squares),
      effect_information = list(subject$effect_information),
      effect_score = list(subject$effect_score),
      error_model = list(subject$error_model),
      split_betas = list(subject$split_betas),
      split_provenance = list(subject$split_provenance),
      split_reliability = list(subject$split_reliability),
      provenance = provenance
    ), class = "dkge_data")
  }
  dkge_fit(
    data = data,
    K = fit$K,
    w_method = spec$w_method,
    w_tau = spec$w_tau,
    ridge = spec$ridge,
    rank = spec$rank,
    effect_scaling = spec$effect_scaling,
    keep_X = FALSE,
    force_qspace = TRUE,
    cpca_blocks = spec$cpca_blocks,
    cpca_T = spec$cpca_T,
    cpca_part = spec$cpca_part,
    cpca_ridge = spec$cpca_ridge,
    weights = spec$weights,
    effect_weights = spec$effect_weights,
    debias = spec$debias,
    solver = spec$solver,
    missingness = spec$missingness,
    miss_args = spec$miss_args,
    jd_control = spec$jd_control,
    jd_mask = spec$jd_mask,
    jd_init = spec$jd_init
  )
}

#' @noRd
.dkge_fold_weight_context <- function(fit,
                                      train_ids,
                                      weight_spec = NULL,
                                      ridge = 0,
                                      missingness = NULL,
                                      miss_args = NULL) {
  if (is.null(fit$subjects) || length(fit$subjects) != length(fit$Btil)) {
    return(.dkge_fold_weight_context_legacy(
      fit, train_ids, weight_spec = weight_spec, ridge = ridge,
      missingness = missingness, miss_args = miss_args
    ))
  }
  fold_fit <- .dkge_refit_training_subjects(
    fit, train_ids, weight_spec = weight_spec, ridge = ridge,
    missingness = missingness, miss_args = miss_args
  )
  weight_eval <- list(
    prior = fold_fit$voxel_weights_prior,
    adapt = fold_fit$voxel_weights_adapt,
    total = fold_fit$voxel_weights,
    total_subject = fold_fit$voxel_weights_subject
  )
  list(
    Chat = fold_fit$Chat,
    weights = weight_eval,
    train_ids = train_ids,
    weight_spec = fold_fit$weight_spec,
    subject_weights = fold_fit$weights,
    pair_counts = fold_fit$pair_counts,
    pair_weight = fold_fit$pair_weight,
    pair_ess = fold_fit$pair_ess,
    effect_precision = fold_fit$effect_precision,
    effect_precision_diagnostics = fold_fit$effect_precision_diagnostics,
    missingness = fold_fit$missingness,
    miss_args = fold_fit$miss_args,
    R = fold_fit$R,
    ruler_source = fold_fit$ruler_source,
    fit = fold_fit,
    estimator = list(
      effect_weights = fold_fit$effect_weight_spec$method %||% "none",
      debias = fold_fit$debias %||% "none",
      missingness = fold_fit$missingness %||% "none",
      effect_scaling = fold_fit$effect_scaling %||% "pooled_design",
      train_subjects = fit$subject_ids[train_ids]
    )
  )
}

#' Legacy fallback for fit objects created before subject records were retained
#' @noRd
.dkge_fold_weight_context_legacy <- function(fit,
                                             train_ids,
                                             weight_spec = NULL,
                                             ridge = 0,
                                             missingness = NULL,
                                             miss_args = NULL) {
  stopifnot(inherits(fit, "dkge"))
  weight_spec <- weight_spec %||% fit$weight_spec %||% dkge_weights(adapt = "none")
  stopifnot(inherits(weight_spec, "dkge_weights"))
  missingness <- missingness %||% fit$missingness %||% "none"
  missingness <- match.arg(missingness, c("none", "rescale", "mask", "shrink"))
  miss_args <- miss_args %||% fit$miss_args %||% list()

  kernel_payload <- .dkge_weight_kernel_payload(fit$K, fit$kernel_info)
  B_train <- fit$Btil[train_ids]
  Omega_train <- fit$Omega[train_ids]
  if (!is.null(weight_spec$B_list2)) {
    weight_spec$B_list2 <- weight_spec$B_list2[train_ids]
  }
  subject_ids <- fit$subject_ids %||% seq_along(fit$Btil)
  obs_masks_all <- .dkge_obs_masks_from_provenance(fit$provenance,
                                                   subject_ids,
                                                   nrow(fit$K))
  if (is.null(obs_masks_all)) {
    obs_masks_all <- replicate(length(fit$Btil), rep(TRUE, nrow(fit$K)),
                               simplify = FALSE)
  }
  obs_masks_train <- if (is.null(obs_masks_all)) NULL else obs_masks_all[train_ids]

  weight_eval <- .dkge_resolve_voxel_weights(weight_spec, B_train, kernel_payload)
  voxel_weights_train <- weight_eval$total_subject %||% weight_eval$total

  Braw_all <- fit$Braw
  if (is.null(Braw_all)) {
    Braw_all <- if (identical(fit$effect_scaling, "none") || is.null(fit$R)) {
      fit$Btil
    } else {
      lapply(fit$Btil, function(B) forwardsolve(t(fit$R), B))
    }
  }
  subjects_all <- fit$subjects
  if (is.null(subjects_all)) {
    subjects_all <- lapply(seq_along(Braw_all), function(s) list(
      id = subject_ids[[s]], effect_noise_cov = NULL,
      residual_variance = NULL, noise_trace = NULL
    ))
  }
  effect_precision_all <- fit$effect_precision
  if (is.null(effect_precision_all)) {
    effect_precision_all <- lapply(obs_masks_all, as.numeric)
  }

  R_fit <- fit$R %||% diag(nrow(fit$K))
  subject_weights <- if (is.null(fit$w_method)) {
    # Backward-compatible path for legacy/minimal fit objects whose stored
    # weights have no reproducible weighting specification.
    fit$weights[train_ids]
  } else {
    .dkge_subject_weights(
      Braw_all[train_ids],
      Omega_train,
      fit$Khalf,
      fit$w_method,
      fit$w_tau %||% 0.3,
      obs_masks = obs_masks_train,
      R = R_fit,
      voxel_weights = voxel_weights_train
    )
  }
  accum <- .dkge_build_moment_pool(
    subjects = subjects_all[train_ids],
    B_list = Braw_all[train_ids],
    Omega_list = Omega_train,
    voxel_weights = voxel_weights_train,
    obs_masks = obs_masks_train,
    subject_weights = subject_weights,
    effect_precision = effect_precision_all[train_ids],
    effect_method = fit$effect_weight_spec$method %||% "none",
    R = R_fit,
    Khalf = fit$Khalf,
    missingness = missingness,
    miss_args = miss_args,
    debias = fit$debias %||% "none"
  )
  Chat <- accum$Chat

  if (ridge > 0) Chat <- Chat + ridge * diag(nrow(Chat))
  Chat <- (Chat + t(Chat)) / 2

  list(
    Chat = Chat,
    weights = weight_eval,
    train_ids = train_ids,
    weight_spec = weight_spec,
    subject_weights = subject_weights,
    pair_counts = accum$pair_counts,
    pair_weight = accum$pair_weight,
    pair_ess = accum$pair_ess,
    missingness = missingness,
    miss_args = miss_args
  )
}
