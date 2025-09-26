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
                                   weights = NULL) {
  stopifnot(inherits(fit, "dkge"))
  stopifnot(is.list(assignments), length(assignments) >= 1)
  loader_scope <- match.arg(loader_scope)

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
  recycled_subjects <- character(0)

  for (fold_idx in seq_len(n_folds)) {
    holdout <- assignments[[fold_idx]]
    train_ids <- setdiff(seq_len(S), holdout)

    if (verbose_flag) {
      message(sprintf("Building fold %d/%d (holdout: %s)",
                      fold_idx, n_folds,
                      paste(subject_ids[holdout], collapse = ", ")))
    }

    ctx <- .dkge_fold_weight_context(fit, train_ids, weight_spec, ridge = ridge)
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
      D_minus = fold_evals[[fold_idx]]
    )
  }

  list(
    folds = folds,
    assignments = assignments,
    align = align,
    consensus = consensus,
    alignment = alignment,
    loader_scope = loader_scope,
    weight_spec = weight_spec
  ) -> result

  if (length(recycled_subjects)) {
    warning(sprintf("Per-subject voxel weights recycled for: %s",
                    paste(recycled_subjects, collapse = ", ")))
    attr(result, "recycled_weights_subjects") <- recycled_subjects
  }
  result
}
