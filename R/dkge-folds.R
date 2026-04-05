# dkge-folds.R
# Shared fold-building helpers for LOSO/K-fold cross-fitting.

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

    eig_fold <- eigen(Chat_minus, symmetric = TRUE)
    U_fold <- fit$Kihalf %*% eig_fold$vectors[, seq_len(r), drop = FALSE]

    fold_bases[[fold_idx]] <- U_fold
    fold_evals[[fold_idx]] <- eig_fold$values

    loader_weights <- weight_eval$total

    subject_scope <- if (loader_scope == "all") seq_len(S) else holdout
    loader_list <- vector("list", length(subject_scope))
    names(loader_list) <- as.character(subject_scope)

    for (j in seq_along(subject_scope)) {
      s <- subject_scope[[j]]
      Bts <- fit$Btil[[s]]
      w_s <- loader_weights
      if (!is.null(w_s) && length(w_s) != ncol(Bts)) {
        if (length(w_s) > 1L) {
          recycled_subjects <- unique(c(recycled_subjects, subject_ids[s]))
        }
        w_s <- rep(w_s, length.out = ncol(Bts))
      }
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
.dkge_fold_weight_context <- function(fit,
                                      train_ids,
                                      weight_spec = NULL,
                                      ridge = 0,
                                      missingness = c("none", "rescale", "mask", "shrink"),
                                      miss_args = list()) {
  stopifnot(inherits(fit, "dkge"))
  weight_spec <- weight_spec %||% fit$weight_spec %||% dkge_weights(adapt = "none")
  stopifnot(inherits(weight_spec, "dkge_weights"))
  missingness <- match.arg(missingness)

  kernel_payload <- .dkge_weight_kernel_payload(fit$K, fit$kernel_info)
  B_train <- fit$Btil[train_ids]
  Omega_train <- fit$Omega[train_ids]
  subject_weights <- fit$weights[train_ids]

  weight_eval <- .dkge_resolve_voxel_weights(weight_spec, B_train, kernel_payload)
  voxel_weights_train <- weight_eval$total_subject %||% weight_eval$total

  accum <- .dkge_accumulate_chat(B_train, Omega_train, fit$Khalf, subject_weights,
                                 voxel_weights = voxel_weights_train)
  Chat <- accum$Chat

  pair_counts <- NULL
  prov <- fit$provenance
  masks <- prov$obs_mask %||% NULL
  if (!is.null(masks) && length(train_ids) > 0L) {
    subject_ids <- fit$subject_ids %||% seq_along(fit$Btil)
    if (!is.null(names(masks))) {
      order_idx <- match(subject_ids, names(masks))
      if (!anyNA(order_idx)) {
        masks <- masks[order_idx]
      }
    }
    if (length(masks) >= max(train_ids)) {
      freq_mat <- matrix(0L, nrow = nrow(Chat), ncol = ncol(Chat))
      dimnames(freq_mat) <- dimnames(Chat)
      if (!is.null(prov$effect_ids) && length(prov$effect_ids) == nrow(Chat)) {
        dimnames(freq_mat) <- list(prov$effect_ids, prov$effect_ids)
      }
      for (idx in train_ids) {
        mask <- masks[[idx]]
        if (is.null(mask)) next
        mask <- as.logical(mask)
        if (length(mask) != nrow(Chat)) next
        sel <- which(mask)
        if (!length(sel)) next
        freq_mat[sel, sel] <- freq_mat[sel, sel] + 1L
      }
      pair_counts <- freq_mat

      if (!identical(missingness, "none")) {
        pc_safe <- pmax(freq_mat, 1L)
        if (identical(missingness, "rescale")) {
          Chat <- Chat / pc_safe
          Chat[freq_mat == 0L] <- 0
        } else if (identical(missingness, "mask")) {
          threshold <- miss_args$min_pairs %||% 1L
          mask_zero <- freq_mat < threshold
          Chat[mask_zero] <- 0
        } else if (identical(missingness, "shrink")) {
          rescaled <- Chat / pc_safe
          rescaled[freq_mat == 0L] <- 0
          max_pc <- max(freq_mat)
          gamma <- miss_args$gamma %||% 1
          if (max_pc <= 0) {
            weights_mat <- matrix(0, nrow = nrow(Chat), ncol = ncol(Chat))
          } else {
            weights_mat <- (freq_mat / max_pc)^gamma
          }
          if (!is.null(dimnames(Chat))) {
            dimnames(weights_mat) <- dimnames(Chat)
          }
          diag_part <- diag(diag(Chat), nrow = nrow(Chat), ncol = ncol(Chat))
          if (!is.null(dimnames(Chat))) {
            dimnames(diag_part) <- dimnames(Chat)
          }
          Chat <- weights_mat * rescaled + (1 - weights_mat) * diag_part
        }
      }
    }
  }

  if (ridge > 0) Chat <- Chat + ridge * diag(nrow(Chat))
  Chat <- (Chat + t(Chat)) / 2

  list(
    Chat = Chat,
    weights = weight_eval,
    train_ids = train_ids,
    weight_spec = weight_spec,
    pair_counts = pair_counts,
    missingness = missingness,
    miss_args = miss_args
  )
}
