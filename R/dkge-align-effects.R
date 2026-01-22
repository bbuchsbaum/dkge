# dkge-align-effects.R
# Fold-aware helpers for aligning subject effect kernels under partial overlap.

#' Align effect kernels across subjects with partial overlap
#'
#' Builds a fold-aware common effect geometry using masked averaging over the
#' training subjects, then completes each subject's kernel on the shared effect
#' index via Nyström, shrinkage, or intersection.
#'
#' @param K_list List of per-subject symmetric kernels; `K_list[[s]]` has
#'   dimensions `|O_s| x |O_s|`.
#' @param effects List of character vectors. `effects[[s]]` gives the effect IDs
#'   associated with the rows and columns of `K_list[[s]]`.
#' @param subject_ids Optional character vector naming subjects; defaults to the
#'   list names or sequential labels.
#' @param folds Optional cross-validation structure defining held-out subjects
#'   per fold (see Details).
#' @param mode Completion mode. One of \code{"nystrom"}, \code{"shrinkage"}, or
#'   \code{"intersection"}.
#' @param weights Optional numeric vector of subject weights used when pooling
#'   the training kernels.
#' @param ridge Ridge factor used when inverting training blocks for Nyström.
#' @param alpha Shrinkage weight applied in \code{mode = "shrinkage"}.
#' @param ensure_psd Logical; when `TRUE` (default) project pooled and completed
#'   matrices to the PSD cone.
#' @param psd_tol Eigenvalue floor expressed as a fraction of the largest
#'   eigenvalue when projecting to PSD.
#' @param min_train_coverage Drop effects observed by fewer than this many
#'   training subjects when forming the union.
#' @param intersection_scope When \code{mode = "intersection"}, restrict the
#'   intersection to \code{"all_subjects"} (default) or \code{"train_only"}.
#' @param effect_prior Optional PSD matrix indexed by effect IDs used to seed
#'   zero-coverage entries of the group kernel.
#' @param prior_weight Blend factor in [0, 1] applied to `effect_prior` when
#'   available.
#' @param verbose Logical; emit messages when `TRUE`.
#'
#' @details
#' The `folds` argument accepts: `NULL` (single context), a `dkge_folds` object,
#' a data frame with columns `subject` and `fold`, or a list whose elements name
#' the held-out subjects. Fold-specific results are returned under
#' `result$folds[[f]]` with training/test indices attached.
#'
#' @return When `folds = NULL`, a list with fields
#'   - `K_aligned`: list of aligned `n x n` kernels per subject
#'   - `effect_ids`: character vector of shared effect IDs
#'   - `G`: pooled training kernel (when applicable)
#'   - `obs_mask`: list of logical vectors indicating observed effects per
#'     subject
#'   - `pair_counts`: integer matrix of training coverage per effect pair
#'   - `coverage`: data frame summarising training coverage per effect
#'   - `mode`: completion mode used.
#'
#'   When folds are supplied, returns `list(folds = list(...))` where each fold
#'   entry includes the same fields along with `train_idx` and `test_idx`.
#'
#' @examples
#' K_list <- list(s1 = diag(5), s2 = diag(4), s3 = diag(5))
#' effects <- list(
#'   s1 = c("a", "b", "c", "d", "e"),
#'   s2 = c("a", "b", "c", "d"),
#'   s3 = c("b", "c", "d", "e", "f")
#' )
#' aligned <- dkge_align_effects(K_list, effects, mode = "intersection")
#' length(aligned$K_aligned)
#'
#' @export
dkge_align_effects <- function(K_list,
                               effects,
                               subject_ids = NULL,
                               folds = NULL,
                               mode = c("nystrom", "shrinkage", "intersection"),
                               weights = NULL,
                               ridge = 1e-6,
                               alpha = 0.25,
                               ensure_psd = TRUE,
                               psd_tol = 1e-10,
                               min_train_coverage = 1L,
                               intersection_scope = c("all_subjects", "train_only"),
                               effect_prior = NULL,
                               prior_weight = 0,
                               verbose = FALSE) {
  mode <- match.arg(mode)
  intersection_scope <- match.arg(intersection_scope)

  stopifnot(is.list(K_list), length(K_list) >= 1L)
  stopifnot(is.list(effects), length(effects) == length(K_list))

  S <- length(K_list)
  if (is.null(subject_ids)) {
    subject_ids <- names(K_list)
    if (is.null(subject_ids)) {
      subject_ids <- as.character(seq_len(S))
    }
  } else {
    stopifnot(length(subject_ids) == S)
  }

  if (is.null(weights)) {
    weights <- rep(1, S)
  } else {
    stopifnot(length(weights) == S)
    weights <- as.numeric(weights)
  }

  for (i in seq_len(S)) {
    Ki <- K_list[[i]]
    Ei <- effects[[i]]
    if (!is.matrix(Ki) || nrow(Ki) != ncol(Ki)) {
      stop(sprintf("K_list[[%d]] must be a square matrix.", i), call. = FALSE)
    }
    if (nrow(Ki) != length(Ei)) {
      stop(sprintf("Length of effects[[%d]] must match nrow(K_list[[%d]]).", i, i),
           call. = FALSE)
    }
    K_list[[i]] <- 0.5 * (Ki + t(Ki))
    effects[[i]] <- as.character(Ei)
  }

  fold_list <- .dkge_align_expand_folds(folds, subject_ids)

  align_context <- function(train_idx, test_idx, tag) {
    if (mode == "intersection") {
      if (intersection_scope == "all_subjects") {
        union_ids <- Reduce(intersect, effects)
      } else {
        union_ids <- Reduce(intersect, effects[train_idx])
      }
      if (length(union_ids) == 0L) {
        stop("Intersection of effects is empty for mode='intersection'.", call. = FALSE)
      }
      cover_counts <- NULL
    } else {
      union_ids <- unique(unlist(effects[train_idx], use.names = FALSE))
      if (!length(union_ids)) {
        stop("Training subjects contribute zero effects; cannot form union.", call. = FALSE)
      }
      cover_counts <- table(unlist(lapply(effects[train_idx], unique), use.names = FALSE))
      keep <- names(cover_counts)[cover_counts >= min_train_coverage]
      union_ids <- intersect(union_ids, keep)
      if (!length(union_ids)) {
        stop("All effects filtered by min_train_coverage; relax threshold.", call. = FALSE)
      }
      if (!is.null(effect_prior) && prior_weight > 0) {
        prior_ids <- rownames(effect_prior)
        if (is.null(prior_ids) || !identical(prior_ids, colnames(effect_prior))) {
          stop("effect_prior must have matching row and column names.", call. = FALSE)
        }
        test_union <- unique(unlist(effects[test_idx], use.names = FALSE))
        add_ids <- intersect(setdiff(test_union, union_ids), prior_ids)
        if (length(add_ids)) {
          base_len <- length(union_ids)
          union_ids <- c(union_ids, add_ids)
          if (length(add_ids) > base_len) {
            warning(sprintf("Prior contributed %d additional effects beyond training union in context %s.",
                            length(add_ids), tag), call. = FALSE)
          }
        }
      }
    }
    union_ids <- unique(union_ids)

    n <- length(union_ids)
    pair_counts <- matrix(0L, n, n, dimnames = list(union_ids, union_ids))
    if (length(train_idx) > 0L) {
      for (s in train_idx) {
        idx <- match(effects[[s]], union_ids)
        idx <- idx[!is.na(idx)]
        if (length(idx)) {
          pair_counts[idx, idx] <- pair_counts[idx, idx] + 1L
        }
      }
    }

    G <- NULL
    if (mode != "intersection") {
      G <- .dkge_align_group_kernel(K_list, effects, train_idx, union_ids,
                                    weights = weights, ensure_psd = ensure_psd,
                                    psd_tol = psd_tol, verbose = verbose)
      if (!is.null(effect_prior) && prior_weight > 0) {
        prior_ids <- rownames(effect_prior)
        Pidx <- match(union_ids, prior_ids)
        use_prior <- !is.na(Pidx)
        if (any(use_prior)) {
          G_prior <- effect_prior[Pidx[use_prior], Pidx[use_prior], drop = FALSE]
          G0 <- G[use_prior, use_prior, drop = FALSE]
          G[use_prior, use_prior] <- (1 - prior_weight) * G0 + prior_weight * G_prior
          if (ensure_psd) {
            G[use_prior, use_prior] <- .dkge_align_psd_project(G[use_prior, use_prior], tol = psd_tol)
          }
        }
      }
    }

    K_aligned <- vector("list", S)
    obs_mask <- vector("list", S)
    names(K_aligned) <- subject_ids
    names(obs_mask) <- subject_ids

    for (s in seq_len(S)) {
      idx <- match(effects[[s]], union_ids)
      keep <- !is.na(idx)
      obs <- idx[keep]
      order_obs <- order(obs)
      obs_ord <- obs[order_obs]
      mask <- rep(FALSE, n)
      if (length(obs)) mask[obs] <- TRUE
      obs_mask[[s]] <- mask

      if (!length(obs)) {
        if (mode == "intersection") {
          warning(sprintf("Subject %s has no overlap with intersection (context %s); returning zeros.",
                          subject_ids[[s]], tag), call. = FALSE)
        }
        Khat <- matrix(0, n, n)
        dimnames(Khat) <- list(union_ids, union_ids)
        K_aligned[[s]] <- Khat
        next
      }

      if (mode == "intersection") {
        Ki <- K_list[[s]][keep, keep, drop = FALSE]
        Ki <- Ki[order_obs, order_obs, drop = FALSE]
        Khat <- matrix(0, n, n)
        Khat[obs_ord, obs_ord] <- Ki
      } else if (mode == "shrinkage") {
        Ki <- K_list[[s]][keep, keep, drop = FALSE]
        Ki <- Ki[order_obs, order_obs, drop = FALSE]
        Khat <- G
        if (length(obs_ord)) {
          Gobs <- G[obs_ord, obs_ord, drop = FALSE]
          Khat[obs_ord, obs_ord] <- (1 - alpha) * Ki + alpha * Gobs
        }
      } else {
        Ki <- K_list[[s]][keep, keep, drop = FALSE]
        Ki <- Ki[order_obs, order_obs, drop = FALSE]
        Khat <- .dkge_align_complete_subject(Ki, obs_ord, G, ridge = ridge)
      }

      if (ensure_psd) {
        Khat <- .dkge_align_psd_project(Khat, tol = psd_tol)
      } else {
        Khat <- 0.5 * (Khat + t(Khat))
      }
      dimnames(Khat) <- list(union_ids, union_ids)
      K_aligned[[s]] <- Khat
    }

    subj_counts <- if (is.null(cover_counts)) NA_integer_ else as.integer(cover_counts[union_ids])
    coverage <- data.frame(
      effect = union_ids,
      train_subjects = subj_counts,
      subjects = subj_counts
    )

    list(
      K_aligned = K_aligned,
      effect_ids = union_ids,
      G = G,
      obs_mask = obs_mask,
      pair_counts = pair_counts,
      coverage = coverage,
      train_idx = train_idx,
      test_idx = test_idx,
      mode = mode
    )
  }

  if (is.null(fold_list)) {
    ctx <- align_context(seq_len(S), integer(0), "all")
    return(ctx)
  }

  fold_out <- vector("list", length(fold_list))
  names(fold_out) <- paste0("fold_", seq_along(fold_list))
  for (f in seq_along(fold_list)) {
    test_idx <- fold_list[[f]]
    train_idx <- setdiff(seq_len(S), test_idx)
    if (!length(train_idx)) {
      stop("Each fold must retain at least one training subject.", call. = FALSE)
    }
    fold_out[[f]] <- align_context(train_idx, test_idx, names(fold_out)[[f]])
  }
  list(folds = fold_out)
}

.dkge_align_expand_folds <- function(folds, subject_ids) {
  if (is.null(folds)) {
    return(NULL)
  }
  S <- length(subject_ids)
  coerce_idx <- function(x) {
    if (is.null(x)) {
      integer(0)
    } else if (is.character(x)) {
      idx <- match(x, subject_ids)
      if (anyNA(idx)) {
        stop("Fold specification references unknown subject IDs.", call. = FALSE)
      }
      sort(unique(idx))
    } else {
      idx <- as.integer(x)
      if (any(idx < 1L) || any(idx > S)) {
        stop("Fold specification contains invalid subject indices.", call. = FALSE)
      }
      sort(unique(idx))
    }
  }

  if (inherits(folds, "dkge_folds") && !is.null(folds$assignments)) {
    return(lapply(folds$assignments, coerce_idx))
  }

  if (is.data.frame(folds) && all(c("subject", "fold") %in% names(folds))) {
    return(lapply(split(folds$subject, folds$fold), coerce_idx))
  }

  if (is.list(folds)) {
    return(lapply(folds, coerce_idx))
  }

  stop("Unsupported fold specification. Provide a dkge_folds object, data frame, or list.",
       call. = FALSE)
}

.dkge_align_group_kernel <- function(K_list,
                                     effects,
                                     train_idx,
                                     union_ids,
                                     weights,
                                     ensure_psd,
                                     psd_tol,
                                     verbose) {
  n <- length(union_ids)
  Sacc <- matrix(0, n, n)
  Cacc <- matrix(0, n, n)
  for (s in train_idx) {
    idx <- match(effects[[s]], union_ids)
    keep <- !is.na(idx)
    idx <- idx[keep]
    if (!length(idx)) next
    Ki <- K_list[[s]][keep, keep, drop = FALSE]
    Ki <- 0.5 * (Ki + t(Ki))
    w <- weights[s]
    Sacc[idx, idx] <- Sacc[idx, idx] + w * Ki
    Cacc[idx, idx] <- Cacc[idx, idx] + w
  }

  nz <- Cacc > 0
  G <- matrix(0, n, n)
  G[nz] <- Sacc[nz] / Cacc[nz]
  G <- 0.5 * (G + t(G))
  if (ensure_psd) {
    G <- .dkge_align_psd_project(G, tol = psd_tol)
  }
  dimnames(G) <- list(union_ids, union_ids)
  if (.dkge_verbose(verbose)) {
    message(sprintf("Aligned union contains %d effects; %d training subjects.",
                    n, length(train_idx)))
  }
  G
}

.dkge_align_complete_subject <- function(Ks,
                                         obs_idx,
                                         G,
                                         ridge) {
  n <- nrow(G)
  Khat <- matrix(0, n, n)
  Ks <- 0.5 * (Ks + t(Ks))
  obs_idx <- sort(obs_idx)
  Khat[obs_idx, obs_idx] <- Ks
  if (length(obs_idx) == n) {
    return(0.5 * (Khat + t(Khat)))
  }
  miss_idx <- setdiff(seq_len(n), obs_idx)
  GOO <- G[obs_idx, obs_idx, drop = FALSE]
  GUO <- G[miss_idx, obs_idx, drop = FALSE]
  scale0 <- mean(diag(GOO))
  if (!is.finite(scale0) || scale0 <= 0) scale0 <- 1
  eps <- ridge * scale0 + 1e-12
  A <- GOO + diag(eps, nrow(GOO))
  invA <- tryCatch(chol2inv(chol(A)), error = function(...) solve(A, diag(nrow(A))))
  K_UO <- GUO %*% invA %*% Ks
  Khat[miss_idx, obs_idx] <- K_UO
  Khat[obs_idx, miss_idx] <- t(K_UO)
  Khat[miss_idx, miss_idx] <- GUO %*% invA %*% Ks %*% invA %*% t(GUO)
  0.5 * (Khat + t(Khat))
}

.dkge_align_psd_project <- function(M, tol = 1e-10) {
  M <- 0.5 * (M + t(M))
  eig <- eigen(M, symmetric = TRUE)
  vals <- eig$values
  vmax <- max(vals, 0)
  floor <- tol * if (vmax > 0) vmax else 1
  vals2 <- pmax(vals, floor)
  Mpsd <- eig$vectors %*% (vals2 * t(eig$vectors))
  0.5 * (Mpsd + t(Mpsd))
}
