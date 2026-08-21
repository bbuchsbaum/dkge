# dkge-folds.R
# Shared fold-building helpers for LOSO/K-fold cross-fitting.

.dkge_subject_loader_weights <- function(loader_weights, Bts) {
  if (is.null(loader_weights) || length(loader_weights) == 0L) {
    return(NULL)
  }
  if (!is.numeric(loader_weights) || any(!is.finite(loader_weights)) ||
      any(loader_weights < 0)) {
    .dkge_abort(
      "Subject loader weights must be finite non-negative numeric values.",
      "dkge_weight_domain_error"
    )
  }
  w_s <- as.numeric(loader_weights)
  if (length(w_s) == ncol(Bts)) {
    return(w_s)
  }
  if (length(w_s) >= 1L && all(is.finite(w_s)) &&
      all(w_s == w_s[[1L]])) {
    return(rep(w_s[[1L]], ncol(Bts)))
  }
  .dkge_abort(
    sprintf(
      "Subject loader weights have length %d but the subject block has %d columns (parcels); only an exactly matching vector or a constant profile may be re-expanded.",
      length(w_s), ncol(Bts)
    ),
    "dkge_weight_dimension_error"
  )
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
                                   missingness = c("none", "rescale", "mask", "shrink"),
                                   miss_args = list()) {
  stopifnot(inherits(fit, "dkge"))
  stopifnot(is.list(assignments), length(assignments) >= 1)
  loader_scope <- match.arg(loader_scope)
  missingness <- match.arg(missingness)

  S <- length(fit$Btil)
  q <- nrow(fit$U)
  r <- ncol(fit$U)
  subject_ids <- fit$subject_ids %||% seq_len(S)

  assignments <- lapply(assignments, function(idx) {
    if (!is.numeric(idx) || !length(idx) || any(!is.finite(idx)) ||
        any(idx != trunc(idx))) {
      .dkge_abort(
        "Each fold must contain finite integer subject indices before coercion.",
        "dkge_fold_partition_error"
      )
    }
    idx <- sort(unique(as.integer(idx)))
    if (any(idx < 1L) || any(idx > S)) {
      .dkge_abort(
        sprintf("Fold assignments must contain subject indices from 1 through %d.", S),
        "dkge_fold_partition_error"
      )
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
  fold_training_rank <- integer(n_folds)
  fold_moment_tolerance <- numeric(n_folds)
  roots <- .dkge_psd_roots(fit$K)

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

    basis <- .dkge_select_k_basis(
      Chat_minus,
      K = fit$K,
      roots = roots,
      requested_rank = r,
      label = sprintf("Fold %d training", fold_idx)
    )
    fold_bases[[fold_idx]] <- basis$U
    fold_evals[[fold_idx]] <- basis$values
    fold_training_rank[[fold_idx]] <- basis$effective_rank
    fold_moment_tolerance[[fold_idx]] <- basis$moment_tolerance
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
  }

  # Cross-fitted consumers need one component dimension. Use the rank shared by
  # every training fold, never the rank of the full-data fit.
  common_rank <- if (n_folds) {
    min(vapply(fold_bases, ncol, integer(1)))
  } else {
    0L
  }
  if (common_rank < r) {
    message(sprintf(
      "Fold training rank caps the cross-fitted basis at %d component(s) (full-data rank %d).",
      common_rank, r
    ))
  }

  for (fold_idx in seq_len(n_folds)) {
    U_fold <- fold_bases[[fold_idx]][, seq_len(common_rank), drop = FALSE]
    fold_bases[[fold_idx]] <- U_fold
    if (common_rank > 0L) {
      gram_error <- max(abs(crossprod(U_fold, fit$K %*% U_fold) -
                              diag(common_rank)))
      if (!is.finite(gram_error) || gram_error > 1e-7) {
        .dkge_abort(
          sprintf(
            "Fold %d basis failed the K-orthonormal postcondition (maximum Gram error %.3e).",
            fold_idx, gram_error
          ),
          "dkge_kernel_geometry_error"
        )
      }
    }

    holdout <- assignments[[fold_idx]]
    loader_weights <- fold_weight_info[[fold_idx]]$total
    subject_scope <- if (loader_scope == "all") seq_len(S) else holdout
    loader_list <- vector("list", length(subject_scope))
    names(loader_list) <- as.character(subject_scope)

    for (j in seq_along(subject_scope)) {
      s <- subject_scope[[j]]
      Bts <- fit$Btil[[s]]
      w_s <- .dkge_subject_loader_weights(loader_weights, Bts)
      Bw <- if (is.null(w_s) || length(w_s) == 0L) {
        Bts
      } else {
        sweep(Bts, 2L, sqrt(pmax(w_s, 0)), "*")
      }
      A_s <- t(Bw) %*% fit$K %*% U_fold
      Y_s <- Bw %*% A_s
      loader_list[[j]] <- list(
        subject = s,
        A = A_s,
        Y = Y_s,
        n_cluster = ncol(Bts)
      )
    }
    fold_loaders[[fold_idx]] <- loader_list
  }

  aligned_bases <- fold_bases
  rotations <- vector("list", n_folds)
  consensus <- NULL
  alignment <- NULL

  if (align && n_folds > 0 && common_rank > 0L) {
    align_obj <- dkge_align_bases_K(fold_bases, fit$K, allow_reflection = TRUE)
    aligned_bases <- align_obj$U_aligned
    rotations <- align_obj$R
    alignment <- align_obj

    fold_weights <- vapply(assignments, length, numeric(1))
    if (any(is.na(fold_weights)) || sum(fold_weights) <= 0) {
      fold_weights <- rep(1, n_folds)
    }
    consensus <- dkge_consensus_basis_K(fold_bases, fit$K,
                                        weights = fold_weights,
                                         allow_reflection = TRUE)
  }

  folds <- vector("list", n_folds)
  for (fold_idx in seq_len(n_folds)) {
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
      training_rank = fold_training_rank[[fold_idx]],
      effective_rank = common_rank,
      moment_tolerance = fold_moment_tolerance[[fold_idx]],
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
    effective_rank = common_rank,
    loader_scope = loader_scope,
    weight_spec = weight_spec,
    missingness = missingness,
    miss_args = miss_args
  ) -> result

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
    A_s <- t(Bts) %*% fit$K %*% U_global
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
    holdout <- assignments[[fold_idx]]
    if (!is.numeric(holdout) || !length(holdout) ||
        any(!is.finite(holdout)) || any(holdout != trunc(holdout))) {
      .dkge_abort(
        "Each fold must contain finite integer subject indices before coercion.",
        "dkge_fold_partition_error"
      )
    }
    holdout <- sort(unique(as.integer(holdout)))
    if (any(holdout < 1L) || any(holdout > S)) {
      .dkge_abort(
        sprintf("Fold assignments must contain subject indices from 1 through %d.", S),
        "dkge_fold_partition_error"
      )
    }
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

#' Require one assessment for every subject in subject-level consumers
#'
#' @keywords internal
#' @noRd
.dkge_require_unique_assessments <- function(fold_obj, consumer,
                                             n_subjects = NULL) {
  assignments <- fold_obj$assignments %||% list()
  assessment_ids <- unlist(assignments, use.names = FALSE)
  if (anyDuplicated(assessment_ids)) {
    .dkge_abort(
      sprintf(
        "%s does not support repeated assessment sets; use nonoverlapping folds until an explicit aggregation policy is available.",
        consumer
      ),
      "dkge_fold_partition_error"
    )
  }
  n_subjects <- n_subjects %||% fold_obj$metadata$n_subjects
  if (!is.null(n_subjects)) {
    covered <- sort(unique(as.integer(assessment_ids)))
    if (!identical(covered, seq_len(n_subjects))) {
      .dkge_abort(
        sprintf(
          "%s does not support incomplete assessment sets; the supplied folds cover %d of %d subjects. Use a complete nonoverlapping partition until an explicit subset-labeling policy is available.",
          consumer, length(covered), n_subjects
        ),
        "dkge_fold_partition_error"
      )
    }
  }
  invisible(fold_obj)
}

#' @noRd
.dkge_normalize_folds <- function(folds, fit, consumer = "This operation") {
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
  .dkge_require_unique_assessments(fold_obj, consumer, n_subjects = S)
  list(assignments = fold_obj$assignments, folds = fold_obj)
}

#' Do a fold's voxel weights reproduce the ones used at fit time?
#'
#' When they do, `fit$effect_moments` are exactly the per-subject moments this
#' fold would recompute from the raw betas, so `.dkge_repool_fit()` can be used
#' instead (O(S q^2) instead of O(S q^2 P)).
#'
#' @keywords internal
#' @noRd
.dkge_voxel_weights_match <- function(fit, voxel_weights_train, train_ids) {
  if (!length(train_ids)) return(FALSE)
  moments <- fit[["effect_moments"]]
  if (is.null(moments) || length(moments) < max(train_ids)) return(FALSE)
  if (is.null(fit$Khalf) || is.null(fit$R)) return(FALSE)
  fit_weights <- fit$voxel_weights_subject %||% fit$voxel_weights
  same <- function(a, b) {
    if (is.null(a) && is.null(b)) return(TRUE)
    if (is.null(a) || is.null(b)) return(FALSE)
    isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 0))
  }
  if (is.list(voxel_weights_train)) {
    if (!is.list(fit_weights) || length(fit_weights) < max(train_ids)) return(FALSE)
    if (length(voxel_weights_train) != length(train_ids)) return(FALSE)
    return(all(vapply(seq_along(train_ids), function(j) {
      same(fit_weights[[train_ids[[j]]]], voxel_weights_train[[j]])
    }, logical(1))))
  }
  if (is.list(fit_weights)) return(FALSE)
  same(fit_weights, voxel_weights_train)
}

#' @noRd
.dkge_fold_weight_context <- function(fit,
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
  subject_weights <- fit$weights[train_ids]

  # Reliability weighting cross-references a second run (weight_spec$B_list2),
  # one entry per subject. It must be subset to the training subjects so its
  # length matches B_train; otherwise .dkge_adapt_weights() errors on every
  # fold (and, if lengths happened to align, would leak held-out run-2 data).
  if (!is.null(weight_spec$B_list2)) {
    weight_spec$B_list2 <- weight_spec$B_list2[train_ids]
  }

  weight_eval <- .dkge_resolve_voxel_weights(weight_spec, B_train, kernel_payload)
  voxel_weights_train <- weight_eval$total_subject %||% weight_eval$total

  # When the fold's voxel weights match the ones the fit used, the per-subject
  # raw-effect moments are unchanged and only the pooling has to be redone.
  if (.dkge_voxel_weights_match(fit, voxel_weights_train, train_ids)) {
    pool <- .dkge_repool_fit(fit, indices = train_ids,
                             missingness = missingness, miss_args = miss_args)
    if (!is.null(pool)) {
      Chat <- pool$Chat
      if (ridge > 0) Chat <- Chat + ridge * diag(nrow(Chat))
      Chat <- (Chat + t(Chat)) / 2
      return(list(
        Chat = Chat,
        weights = weight_eval,
        train_ids = train_ids,
        weight_spec = weight_spec,
        pair_counts = pool$pair_counts,
        pair_weight = pool$pair_weight,
        pair_ess = pool$pair_ess,
        missingness = missingness,
        miss_args = miss_args
      ))
    }
  }

  subject_ids <- fit$subject_ids %||% seq_along(fit$Btil)
  obs_masks_all <- .dkge_obs_masks_from_provenance(fit$provenance,
                                                   subject_ids,
                                                   nrow(fit$K))
  if (is.null(obs_masks_all)) {
    obs_masks_all <- replicate(length(fit$Btil), rep(TRUE, nrow(fit$K)),
                               simplify = FALSE)
  }
  obs_masks_train <- obs_masks_all[train_ids]

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
    debias = fit$debias %||% "none",
    contribs = FALSE
  )
  Chat <- accum$Chat

  if (ridge > 0) Chat <- Chat + ridge * diag(nrow(Chat))
  Chat <- (Chat + t(Chat)) / 2

  list(
    Chat = Chat,
    weights = weight_eval,
    train_ids = train_ids,
    weight_spec = weight_spec,
    pair_counts = accum$pair_counts,
    pair_weight = accum$pair_weight,
    pair_ess = accum$pair_ess,
    missingness = missingness,
    miss_args = miss_args
  )
}
