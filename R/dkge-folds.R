# dkge-folds.R
# Shared fold-building helpers for LOSO/K-fold cross-fitting.

.dkge_crossfit_estimator_status <- function(fit) {
  solver <- fit$solver %||% "pooled"
  cpca_part <- fit$cpca$part %||% "none"
  supported <- identical(solver, "pooled") && is.null(fit$cpca)
  list(
    supported = supported,
    solver = solver,
    cpca_part = cpca_part,
    reason = if (supported) {
      "pooled_non_cpca_estimator"
    } else {
      paste0("solver='", solver, "', cpca_part='", cpca_part, "'")
    }
  )
}

.dkge_assert_crossfit_estimator_supported <- function(
    fit, operation = "DKGE cross-fitting") {
  status <- .dkge_crossfit_estimator_status(fit)
  if (!isTRUE(status$supported)) {
    .dkge_abort(
      paste0(
        operation, " currently supports only `solver = 'pooled'` with ",
        "`cpca_part = 'none'`; this fit uses ", status$reason, ". ",
        "Replaying CPCA or joint diagonalization inside every training fold ",
        "is not implemented, so an ordinary pooled eigensolve cannot be ",
        "labelled an exact held-out estimator."
      ),
      "dkge_crossfit_estimator_error"
    )
  }
  invisible(status)
}

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
  if (all(w_s == w_s[[1L]])) {
    return(rep(w_s[[1L]], ncol(Bts)))
  }
  .dkge_abort(
    sprintf(
      paste0(
        "Subject loader weights have length %d but the subject block has ",
        "%d columns (parcels); only an exactly matching vector or a ",
        "constant profile may be re-expanded."
      ),
      length(w_s), ncol(Bts)
    ),
    "dkge_weight_dimension_error"
  )
}

#' Stable identity for inferential alignment payloads
#'
#' Alignment receipts deliberately hash the serialized numerical object, rather
#' than a rounded summary, because a cache or transport plan is only reusable
#' when it was fitted from exactly the same structural input.
#'
#' @keywords internal
#' @noRd
.dkge_object_hash <- function(x) {
  digest::digest(x, algo = "xxhash64", serialize = TRUE)
}

#' Diagnose a fold basis at the fitted truncation boundary
#'
#' @keywords internal
#' @noRd
.dkge_fold_basis_diagnostics <- function(U, evals, K, reference = NULL) {
  r <- ncol(U)
  scale <- max(abs(evals), 0, na.rm = TRUE)
  gap <- if (length(evals) > r) {
    as.numeric(evals[[r]] - evals[[r + 1L]])
  } else {
    NA_real_
  }
  out <- list(
    rank = r,
    eigengap = gap,
    relative_eigengap = if (is.finite(gap) && scale > 0) gap / scale else NA_real_,
    k_orthonormality_error = norm(crossprod(U, K %*% U) - diag(r), "F"),
    procrustes_residual = NA_real_,
    principal_cosines = NULL
  )
  if (!is.null(reference)) {
    pr <- dkge_procrustes_K(reference, U, K, allow_reflection = TRUE)
    delta <- pr$U_aligned - reference
    out$procrustes_residual <- sqrt(max(0, sum(delta * (K %*% delta)))) /
      sqrt(max(1, r))
    out$principal_cosines <- pr$cosines
  }
  out
}

#' Record the exact transformations used to construct one held-out loading
#'
#' @keywords internal
#' @noRd
.dkge_alignment_preprocessing_receipt <- function(fit, s, Btil,
                                                   voxel_weights = NULL) {
  P <- ncol(Btil)
  cluster_order <- colnames(Btil)
  if (is.null(cluster_order)) cluster_order <- seq_len(P)
  op <- .dkge_fit_spatial_operator(fit, subject = s, n_cols = P)
  spatial_payload <- if (is.null(op)) {
    list(active = FALSE, lambda = 0, fingerprint = NULL,
         domain_mode = "identity", operator_binding = NULL)
  } else {
    list(
      active = isTRUE(op$lambda > 0),
      lambda = as.numeric(op$lambda),
      fingerprint = op$fingerprint %||% .dkge_object_hash(
        list(L = op$L, lambda = op$lambda, domain = op$domain)
      ),
      operator_binding = .dkge_spatial_operator_binding(op),
      domain_mode = op$domain_mode %||% NA_character_,
      subject_id = op$subject_id %||% NA_character_
    )
  }
  list(
    schema_version = "2.1.0",
    standardized_beta_source = "fit$Btil",
    effect_scaling = fit$effect_scaling %||% "legacy_or_unspecified",
    ruler_hash = .dkge_object_hash(fit$R %||% diag(nrow(Btil))),
    beta_hash = .dkge_object_hash(Btil),
    voxel_weights_applied = !is.null(voxel_weights),
    voxel_weights = if (is.null(voxel_weights)) NULL else as.numeric(voxel_weights),
    voxel_weights_hash = .dkge_object_hash(
      if (is.null(voxel_weights)) NULL else as.numeric(voxel_weights)
    ),
    spatial = spatial_payload,
    effect_separability = list(
      transform_verified = (fit$effect_scaling %||% "legacy_or_unspecified") %in%
        c("pooled_design", "none"),
      assumption_status = "declared_model_assumption_not_empirically_verified",
      effect_transform = if (identical(fit$effect_scaling, "pooled_design")) {
        "R_transpose"
      } else if (identical(fit$effect_scaling, "none")) {
        "identity"
      } else {
        "unverified"
      },
      parcel_transform = paste0(
        "effect-invariant column scaling",
        if (isTRUE(spatial_payload$active)) {
          " followed by an effect-invariant linear spatial solve"
        } else {
          ""
        }
      ),
      covariance_contract = paste0(
        "Cov(Bmodel[,p], Bmodel[,q]) = rho[p,q] * ",
        "R' Lambda R after the recorded parcel-only transforms"
      )
    ),
    cluster_order = cluster_order,
    cluster_order_hash = .dkge_object_hash(cluster_order),
    n_clusters = P
  )
}

#' Construct one typed fold-safe alignment receipt
#'
#' @keywords internal
#' @noRd
.dkge_make_alignment_receipt <- function(fit, subject, train_ids, basis,
                                         evals, loadings, preprocessing,
                                         alphas, fold_index = subject,
                                         holdout = subject,
                                         subject_weights = NULL,
                                         subject_weight_source = NULL,
                                         method = "loso",
                                         reference_basis = NULL,
                                         eligible = TRUE,
                                         eligibility_reason = "exact_heldout_basis") {
  S <- length(fit$Btil)
  subject_ids <- as.character(fit$subject_ids %||% seq_len(S))
  basis_hash <- .dkge_object_hash(basis)
  # The projector is invariant to an admissible right rotation of the basis.
  subspace_projector <- basis %*% crossprod(basis, fit$K)
  train_ids <- as.integer(train_ids)
  sw <- if (is.null(subject_weights)) NULL else as.numeric(subject_weights)
  if (!is.null(sw)) names(sw) <- subject_ids[train_ids]
  structure(
    list(
      schema_version = "1.0.0",
      estimation_method = method,
      subject_index = as.integer(subject),
      subject_id = subject_ids[[subject]],
      fold_index = as.integer(fold_index),
      holdout_subject_indices = as.integer(holdout),
      holdout_subject_ids = subject_ids[holdout],
      training_subject_indices = train_ids,
      training_subject_ids = subject_ids[train_ids],
      basis = basis,
      basis_id = paste0("dkge-basis-", basis_hash),
      basis_hash = basis_hash,
      subspace_hash = .dkge_object_hash(subspace_projector),
      eigenvalues = as.numeric(evals),
      basis_diagnostics = .dkge_fold_basis_diagnostics(
        basis, evals, fit$K, reference = reference_basis
      ),
      loadings = loadings,
      loadings_hash = .dkge_object_hash(loadings),
      alphas = alphas,
      preprocessing = preprocessing,
      training_subject_weights = sw,
      subject_weight_source = subject_weight_source,
      inference = list(
        eligible = isTRUE(eligible),
        reason = eligibility_reason,
        own_subject_excluded_from_basis = subject %in% holdout &&
          !subject %in% train_ids,
        adaptive_group_alignment_exact = FALSE,
        statement = paste0(
          "Cross-fitting excludes the subject from its estimation basis; ",
          "it does not by itself make adaptive group alignment finite-sample exact."
        )
      )
    ),
    class = c("dkge_alignment_receipt", "list")
  )
}

#' Assemble one receipt per held-out subject from cached fold loaders
#'
#' @keywords internal
#' @noRd
.dkge_alignment_receipts_from_folds <- function(fit, fold_info, alphas,
                                                method = c("loso", "kfold")) {
  method <- match.arg(method)
  S <- length(fit$Btil)
  subject_ids <- as.character(fit$subject_ids %||% seq_len(S))
  contrast_names <- names(alphas) %||% paste0("contrast", seq_along(alphas))
  receipts <- vector("list", S)
  names(receipts) <- subject_ids
  reference_basis <- fold_info$consensus$U %||% NULL

  for (fold in fold_info$folds) {
    train_ids <- fold$training_subjects
    for (s in fold$subjects) {
      if (!is.null(receipts[[s]])) {
        .dkge_abort(
          sprintf("Subject '%s' occurs in more than one held-out fold.",
                  subject_ids[[s]]),
          "dkge_alignment_receipt_error"
        )
      }
      loader <- fold$loaders[[as.character(s)]]
      alpha_list <- lapply(seq_along(alphas), function(i) {
        as.numeric(alphas[[i]][fold$index, , drop = TRUE])
      })
      names(alpha_list) <- contrast_names
      receipts[[s]] <- .dkge_make_alignment_receipt(
        fit = fit,
        subject = s,
        train_ids = train_ids,
        basis = fold$basis,
        evals = fold$evals,
        loadings = loader$A,
        preprocessing = loader$preprocessing,
        alphas = alpha_list,
        fold_index = fold$index,
        holdout = fold$subjects,
        subject_weights = fold$training_subject_weights,
        subject_weight_source = fold$subject_weight_source,
        method = method,
        reference_basis = reference_basis,
        eligible = TRUE,
        eligibility_reason = "exact_heldout_basis"
      )
    }
  }
  if (any(vapply(receipts, is.null, logical(1)))) {
    missing <- subject_ids[vapply(receipts, is.null, logical(1))]
    .dkge_abort(
      sprintf("No fold-safe alignment receipt exists for subject(s): %s.",
              paste(missing, collapse = ", ")),
      "dkge_alignment_receipt_error"
    )
  }
  structure(
    receipts,
    class = c("dkge_alignment_receipts", "list"),
    estimation_method = method,
    inferential_contract = paste0(
      "Each loading uses a basis fitted without its held-out subject. ",
      "Adaptive correspondence remains a separate inferential choice."
    )
  )
}

#' Merge per-contrast alignment receipts when a method evaluates subjects
#' separately for each contrast
#'
#' @keywords internal
#' @noRd
.dkge_merge_alignment_receipt_sets <- function(receipt_sets, contrast_names,
                                               method = "analytic") {
  if (!length(receipt_sets)) {
    return(structure(list(), class = c("dkge_alignment_receipts", "list"),
                     estimation_method = method))
  }
  S <- length(receipt_sets[[1]])
  out <- vector("list", S)
  subject_names <- names(receipt_sets[[1]])
  if (is.null(subject_names)) subject_names <- as.character(seq_len(S))
  names(out) <- subject_names
  for (s in seq_len(S)) {
    candidates <- lapply(receipt_sets, `[[`, s)
    base <- candidates[[1]]
    same_structural_source <- length(unique(vapply(
      candidates,
      function(x) paste(x$basis_hash, x$loadings_hash,
                         x$preprocessing$beta_hash, sep = ":"),
      character(1)
    ))) == 1L
    alpha_list <- lapply(candidates, function(x) as.numeric(x$alphas[[1]]))
    names(alpha_list) <- contrast_names
    base$alphas <- alpha_list
    base$per_contrast_estimation_method <- vapply(
      candidates, `[[`, character(1), "estimation_method"
    )
    names(base$per_contrast_estimation_method) <- contrast_names
    if (!same_structural_source) {
      base$inference$eligible <- FALSE
      base$inference$reason <- "contrast_specific_alignment_sources_disagree"
    } else if (!all(vapply(candidates, function(x) {
      isTRUE(x$inference$eligible)
    }, logical(1)))) {
      base$inference$eligible <- FALSE
      reasons <- unique(vapply(candidates, function(x) {
        x$inference$reason %||% "unknown"
      }, character(1)))
      base$inference$reason <- paste(reasons, collapse = ";")
    }
    out[[s]] <- base
  }
  structure(
    out,
    class = c("dkge_alignment_receipts", "list"),
    estimation_method = method,
    inferential_contract = paste0(
      "Analytic receipts are descriptive for adaptive functional alignment ",
      "unless every subject-contrast path used an exact held-out fallback."
    )
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
  .dkge_assert_crossfit_estimator_supported(fit, "Fold-basis construction")
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

    eig_fold <- eigen(Chat_minus, symmetric = TRUE)
    fold_contract <- .dkge_spectral_contract(eig_fold$values)
    fold_rank <- min(fit$kernel_rank %||% qr(fit$K)$rank,
                     fold_contract$rank)
    if (fold_rank < r) {
      .dkge_abort(
        sprintf(
          paste0(
            "Training fold %d has effective rank %d, below fitted rank %d. ",
            "Refit or cross-validate at rank <= %d."
          ),
          fold_idx, fold_rank, r, fold_rank
        ),
        "dkge_fold_rank_error"
      )
    }
    U_fold <- fit$Kihalf %*% eig_fold$vectors[, seq_len(r), drop = FALSE]
    U_fold <- dkge_k_orthonormalize(U_fold, fit$K)

    fold_bases[[fold_idx]] <- U_fold
    fold_evals[[fold_idx]] <- eig_fold$values

    loader_weights <- weight_eval$total

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
      Bmodel <- .dkge_apply_fit_spatial(fit, Bw, subject = s)
      A_s <- t(Bmodel) %*% fit$K %*% U_fold
      Y_s <- Bmodel %*% A_s
      loader_list[[j]] <- list(
        subject = s,
        A = A_s,
        Y = Y_s,
        n_cluster = ncol(Bts),
        preprocessing = .dkge_alignment_preprocessing_receipt(
          fit, s, Bts, voxel_weights = w_s
        )
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
      weight_spec = weight_spec,
      subject_weights = ctx$subject_weights,
      subject_weight_scores_raw = ctx$subject_weight_scores_raw,
      subject_weight_usable = ctx$subject_weight_usable,
      subject_weight_source = ctx$subject_weight_source
    )
    fold_pair_counts[fold_idx] <- list(ctx$pair_counts)
  }

  aligned_bases <- fold_bases
  rotations <- vector("list", n_folds)
  consensus <- NULL
  alignment <- NULL

  if (align && n_folds > 0) {
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
      training_subjects = setdiff(seq_len(S), assignments[[fold_idx]]),
      basis = fold_bases[[fold_idx]],
      basis_aligned = aligned_bases[[fold_idx]],
      rotation = rotations[[fold_idx]],
      evals = fold_evals[[fold_idx]],
      loaders = fold_loaders[[fold_idx]],
      weights = fold_weight_info[[fold_idx]],
      U_minus = fold_bases[[fold_idx]],
      D_minus = fold_evals[[fold_idx]],
      pair_counts = fold_pair_counts[[fold_idx]],
      training_subject_weights = {
        sw <- fold_weight_info[[fold_idx]]$subject_weights
        if (!is.null(sw)) {
          names(sw) <- as.character(setdiff(seq_len(S), assignments[[fold_idx]]))
        }
        sw
      },
      subject_weight_scores_raw = fold_weight_info[[fold_idx]]$subject_weight_scores_raw,
      subject_weight_usable = fold_weight_info[[fold_idx]]$subject_weight_usable,
      subject_weight_source = fold_weight_info[[fold_idx]]$subject_weight_source,
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

  result
}

#' Build fold loaders using the global (full-data) basis
#'
#' For `mode = "cell"` classification, every subject is projected onto the
#' global `fit$U` rather than fold-specific leave-one-out bases.  This helper
#' produces the same list structure as `.dkge_build_fold_bases()` so that the
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
    Bmodel <- .dkge_apply_fit_spatial(fit, Bts, subject = s)
    A_s <- t(Bmodel) %*% fit$K %*% U_global
    Y_s <- Bmodel %*% A_s
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
        paste0(
          "%s does not support repeated assessment sets; use nonoverlapping ",
          "folds until an explicit aggregation policy is available."
        ),
        consumer
      ),
      "dkge_fold_partition_error"
    )
  }
  n_subjects <- n_subjects %||% fold_obj$metadata$n_subjects
  if (!is.null(n_subjects) && is.finite(n_subjects)) {
    covered <- sort(unique(as.integer(assessment_ids)))
    if (!identical(covered, seq_len(as.integer(n_subjects)))) {
      .dkge_abort(
        sprintf(
          paste0(
            "%s does not support incomplete assessment sets; the supplied ",
            "folds cover %d of %d subjects. Use a complete nonoverlapping ",
            "partition until an explicit subset-labeling policy is available."
          ),
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

#' Recover and normalize subject weights for one training fold
#'
#' New fits retain raw per-subject scores. Older fits retain only the final
#' full-cohort normalized and shrunken weights; because raw score scale cancels
#' on renormalization, those scores can be recovered up to scale whenever
#' `w_tau < 1`. At `w_tau = 1`, equal usable-subject weights are exact.
#'
#' @keywords internal
#' @noRd
.dkge_fold_subject_weights <- function(fit, train_ids, obs_masks_all = NULL) {
  S <- length(fit$Btil)
  if (!is.numeric(train_ids) || anyNA(train_ids) ||
      any(train_ids != as.integer(train_ids)) ||
      any(train_ids < 1L) || any(train_ids > S)) {
    .dkge_abort("`train_ids` must contain valid subject indices.",
                "dkge_fold_subject_index_error")
  }
  train_ids <- as.integer(train_ids)
  tau <- fit$w_tau %||% 0

  raw <- fit$subject_weight_scores_raw
  usable <- fit$subject_weight_usable
  source <- "stored_raw"

  if (is.null(raw) || length(raw) != S ||
      is.null(usable) || length(usable) != S) {
    final_obj <- fit$weights %||% numeric()
    final <- if (is.numeric(final_obj) && is.atomic(final_obj)) {
      as.numeric(final_obj)
    } else {
      numeric()
    }
    if (length(final) == S && all(is.finite(final)) && all(final >= 0)) {
      usable <- final > 0
      raw <- numeric(S)
      if (tau < 1) {
        raw[usable] <- pmax((final[usable] - tau) / (1 - tau), 0)
      } else {
        raw[usable] <- 1
      }
      source <- "legacy_inversion"
    }
  }

  if (!is.null(raw) && length(raw) == S &&
      !is.null(usable) && length(usable) == S && any(usable[train_ids])) {
    weights <- .dkge_normalize_subject_weights(
      raw[train_ids], usable[train_ids], tau
    )
    return(list(
      weights = weights,
      raw = as.numeric(raw[train_ids]),
      usable = as.logical(usable[train_ids]),
      source = source
    ))
  }

  # Last-resort compatibility for malformed/minimal legacy fixtures. Recompute
  # from the training data when the necessary inputs exist; otherwise retain
  # the historical explicit equal-weight fallback.
  can_recompute <- !is.null(fit$Khalf) && length(fit$Omega) == S &&
    !is.null(fit$w_method)
  if (can_recompute) {
    spatial_train <- if (is.null(fit$spatial)) {
      vector("list", length(train_ids))
    } else {
      (fit$spatial$operators %||% vector("list", S))[train_ids]
    }
    score_info <- .dkge_subject_weight_scores(
      fit$Btil[train_ids], fit$Omega[train_ids], fit$Khalf, fit$w_method,
      obs_masks = if (is.null(obs_masks_all)) NULL else obs_masks_all[train_ids],
      spatial_list = spatial_train
    )
    return(list(
      weights = .dkge_normalize_subject_weights(
        score_info$raw, score_info$usable, tau
      ),
      raw = score_info$raw,
      usable = score_info$usable,
      source = "recomputed"
    ))
  }

  list(
    weights = rep(1, length(train_ids)),
    raw = rep(1, length(train_ids)),
    usable = rep(TRUE, length(train_ids)),
    source = "equal_fallback"
  )
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
  subject_ids <- fit$subject_ids %||% seq_along(fit$Btil)
  obs_masks_all <- .dkge_obs_masks_from_provenance(fit$provenance,
                                                   subject_ids,
                                                   nrow(fit$K))
  if (is.null(obs_masks_all)) {
    obs_masks_all <- replicate(length(fit$Btil), rep(TRUE, nrow(fit$K)),
                               simplify = FALSE)
  }
  subject_weight_info <- .dkge_fold_subject_weights(
    fit, train_ids, obs_masks_all = obs_masks_all
  )
  subject_weights <- subject_weight_info$weights
  equal_weight_fallback <- !length(subject_weights) ||
    any(!is.finite(subject_weights)) || sum(subject_weights) <= 0
  if (equal_weight_fallback) {
    # A zero-weight training fold has no moment and previously produced an
    # arbitrary basis from the zero matrix. Match the package's downstream
    # summary policy by falling back to explicit equal weights instead.
    subject_weights <- rep(1, length(train_ids))
  }

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
  if (!equal_weight_fallback &&
      .dkge_voxel_weights_match(fit, voxel_weights_train, train_ids)) {
    pool <- .dkge_repool_fit(fit, indices = train_ids,
                             subject_weights = subject_weights,
                             missingness = missingness, miss_args = miss_args)
    if (!is.null(pool)) {
      Chat <- pool$Chat
      if (ridge > 0) {
        support <- fit$kernel_support_projector %||%
          .dkge_kernel_geometry(fit$K)$support_projector
        Chat <- Chat + ridge * support
      }
      Chat <- (Chat + t(Chat)) / 2
      return(list(
        Chat = Chat,
        weights = weight_eval,
        subject_weights = subject_weights,
        subject_weight_scores_raw = subject_weight_info$raw,
        subject_weight_usable = subject_weight_info$usable,
        subject_weight_source = subject_weight_info$source,
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
    spatial_list = if (is.null(fit$spatial)) {
      vector("list", length(train_ids))
    } else {
      fit$spatial$operators[train_ids]
    },
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

  if (ridge > 0) {
    support <- fit$kernel_support_projector %||%
      .dkge_kernel_geometry(fit$K)$support_projector
    Chat <- Chat + ridge * support
  }
  Chat <- (Chat + t(Chat)) / 2

  list(
    Chat = Chat,
    weights = weight_eval,
    subject_weights = subject_weights,
    subject_weight_scores_raw = subject_weight_info$raw,
    subject_weight_usable = subject_weight_info$usable,
    subject_weight_source = subject_weight_info$source,
    train_ids = train_ids,
    weight_spec = weight_spec,
    pair_counts = accum$pair_counts,
    pair_weight = accum$pair_weight,
    pair_ess = accum$pair_ess,
    missingness = missingness,
    miss_args = miss_args
  )
}
