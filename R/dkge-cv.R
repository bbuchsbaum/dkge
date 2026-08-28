# dkge-cv.R
# Cross-validation helpers and diagnostics for DKGE fits.

#' Compute per-component variance explained
#'
#' Returns the standard deviation, variance, and cumulative variance explained
#' by the DKGE components extracted in [dkge_fit()].
#'
#' @param fit A `dkge` object.
#' @param relative_to Compute variance proportions relative to "kept" (default,
#'   only the retained components) or "total" (all possible components).
#' @return Data frame with columns `component`, `sdev`, `variance`,
#'   `prop_var`, and `cum_prop_var`.
#' @examples
#' toy <- dkge_sim_toy(
#'   factors = list(A = list(L = 2), B = list(L = 3)),
#'   active_terms = c("A", "B"), S = 3, P = 15, snr = 5
#' )
#' fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 2)
#' dkge_variance_explained(fit)
#' @export
dkge_variance_explained <- function(fit, relative_to = c("kept", "total")) {
  stopifnot(inherits(fit, "dkge"))
  relative_to <- match.arg(relative_to)
  evals_all <- fit$evals %||% rep(0, nrow(fit$U))
  sdev <- fit$sdev %||% sqrt(pmax(evals_all[seq_len(fit$rank)], 0))
  variance <- sdev^2
  total_kept <- sum(variance)
  total_all <- sum(pmax(evals_all, 0))
  denom <- if (relative_to == "kept") total_kept else total_all
  prop <- if (denom > 0) variance / denom else rep(NA_real_, length(variance))
  cum_prop <- cumsum(prop)
  data.frame(
    component = seq_along(sdev),
    sdev = sdev,
    variance = variance,
    prop_var = prop,
    cum_prop_var = cum_prop
  )
}

#' Summarize DKGE diagnostics
#'
#' Provides a compact list of variance explained, subject weights, rank, and
#' kernel-support metadata for quick inspection.
#'
#' @param fit A `dkge` object.
#' @return List with variance table, subject weights, rank info, and scalar
#'   kernel diagnostics, including numerical rank/nullity, condition, status,
#'   participation-ratio effective rank, leading-eigenvalue share, and any
#'   model-level spatial-regularization provenance.
#' @examples
#' toy <- dkge_sim_toy(
#'   factors = list(A = list(L = 2), B = list(L = 3)),
#'   active_terms = c("A", "B"), S = 3, P = 15, snr = 5
#' )
#' fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 2)
#' diag <- dkge_diagnostics(fit)
#' names(diag)
#' @export
dkge_diagnostics <- function(fit) {
  stopifnot(inherits(fit, "dkge"))
  if (!is.null(fit$spatial)) {
    .dkge_spatial_fit_payload(fit$spatial, validate = TRUE)
  }
  voxel_stats <- if (!is.null(fit$voxel_weights)) {
    vw <- fit$voxel_weights
    list(mean = mean(vw), sd = stats::sd(vw), min = min(vw), max = max(vw))
  } else NULL
  list(
    variance = dkge_variance_explained(fit),
    weights = fit$weights,
    rank = fit$rank,
    q = nrow(fit$U),
    kernel = fit$kernel_diagnostics %||%
      .dkge_kernel_diagnostics(.dkge_kernel_geometry(fit$K)),
    n_subjects = length(fit$Btil),
    voxel_weights = voxel_stats,
    weight_spec = fit$weight_spec,
    spatial = if (is.null(fit$spatial)) NULL else list(
      active = fit$spatial$active,
      requested = fit$spatial$requested %||% (fit$spatial$lambda > 0),
      effective = fit$spatial$effective %||% fit$spatial$active,
      fully_effective = fit$spatial$fully_effective %||% fit$spatial$active,
      status = fit$spatial$status %||%
        if (isTRUE(fit$spatial$active)) "active" else "inactive",
      lambda = fit$spatial$lambda,
      shared = fit$spatial$shared,
      diagnostics = fit$spatial$diagnostics,
      provenance = fit$spatial$provenance
    )
  )
}

#' One standard-error rule selection helper
#'
#' Aggregates cross-validation scores by parameter setting and returns both the
#' best-performing parameter and the one within one standard error of the best.
#'
#' @param scores Data frame containing per-fold scores.
#' @param param_col Column name identifying the tuning parameter.
#' @param metric_col Column name holding the metric (larger is better).
#' @return List with `best`, `pick`, and `summary` table of mean/se by parameter.
#' @examples
#' scores <- data.frame(param = rep(1:3, each = 2), score = c(0.2, 0.1, 0.25, 0.2, 0.24, 0.23))
#' dkge_one_se(scores, param_col = "param", metric_col = "score")$pick
#' @export
dkge_one_se <- function(scores, param_col = "param", metric_col = "score") {
  stopifnot(is.data.frame(scores))
  agg <- aggregate(scores[[metric_col]], by = list(scores[[param_col]]),
                   FUN = function(x) {
                     n <- length(x)
                     m <- mean(x)
                     se <- if (n > 1) stats::sd(x) / sqrt(n) else 0
                     c(mean = m, se = se)
                   })
  param <- agg[[1]]
  stats_mat <- agg[[2]]
  if (is.list(stats_mat)) {
    stats_mat <- do.call(rbind, stats_mat)
  } else {
    stats_mat <- as.matrix(stats_mat)
  }
  if (is.null(colnames(stats_mat))) {
    colnames(stats_mat) <- c("mean", "se")
  }
  means <- stats_mat[, "mean"]
  ses <- stats_mat[, "se"]
  best_idx <- which.max(means)
  threshold <- means[best_idx] - ses[best_idx]
  pick_idx <- min(which(means >= threshold))
  list(
    best = param[best_idx],
    pick = param[pick_idx],
    summary = data.frame(param = param, mean = means, se = ses)
  )
}

#' Resolve one candidate-independent validation geometry for cross-validation
#'
#' @keywords internal
#' @noRd
.dkge_cv_validation_geometry <- function(validation_K, q) {
  kind <- if (is.null(validation_K)) "effect_space_identity" else "custom"
  if (is.null(validation_K)) validation_K <- diag(q)
  if (!is.matrix(validation_K) || any(dim(validation_K) != c(q, q))) {
    stop(sprintf("`validation_K` must be NULL or a %d by %d matrix.", q, q),
         call. = FALSE)
  }
  geometry <- .dkge_kernel_geometry(validation_K)
  if (geometry$rank == 0L) {
    .dkge_abort(
      "`validation_K` has numerical rank zero and cannot define a held-out score.",
      "dkge_cv_validation_error"
    )
  }
  list(
    roots = geometry,
    diagnostics = c(list(kind = kind), .dkge_kernel_diagnostics(geometry))
  )
}

#' Validate one CV candidate's effect-space dimensions
#'
#' @keywords internal
#' @noRd
.dkge_cv_candidate_geometry <- function(K, q, label = "candidate kernel") {
  if (!is.matrix(K) || !is.numeric(K) || any(dim(K) != c(q, q))) {
    stop(sprintf("%s must be a finite numeric %d by %d matrix.", label, q, q),
         call. = FALSE)
  }
  .dkge_kernel_geometry(K)
}

#' Euclidean orthonormal basis for a matrix column span
#'
#' @keywords internal
#' @noRd
.dkge_cv_span_basis <- function(W, tol = 1e-10) {
  W <- as.matrix(W)
  if (!ncol(W) || !nrow(W)) return(matrix(0, nrow(W), 0L))
  sv <- svd(W, nu = min(dim(W)), nv = 0L)
  if (!length(sv$d) || max(sv$d) <= 0) return(matrix(0, nrow(W), 0L))
  keep <- sv$d > tol * max(sv$d)
  sv$u[, keep, drop = FALSE]
}

#' Score a learned effect-space span in a fixed validation geometry
#'
#' @keywords internal
#' @noRd
.dkge_cv_score_fixed <- function(Bw, U, validation_roots) {
  Xs <- validation_roots$Khalf %*% Bw
  total <- sum(Xs^2)
  if (!is.finite(total) || total <= 0) return(NA_real_)
  Q <- .dkge_cv_span_basis(validation_roots$Khalf %*% U)
  captured <- if (ncol(Q)) sum(crossprod(Q, Xs)^2) else 0
  pmin(1, pmax(0, captured / total))
}

#' Apply a fixed held-out spatial metric without candidate smoothing
#'
#' @keywords internal
#' @noRd
.dkge_cv_apply_spatial_metric <- function(B, Omega = NULL) {
  if (is.null(Omega)) return(B)
  P <- ncol(B)
  if (is.vector(Omega)) {
    if (length(Omega) != P) {
      stop("Held-out diagonal Omega must match the beta-block width.",
           call. = FALSE)
    }
    return(sweep(B, 2L, sqrt(pmax(as.numeric(Omega), 0)), "*"))
  }
  Omega <- as.matrix(Omega)
  if (!identical(dim(Omega), c(P, P))) {
    stop("Held-out full Omega must be square and match the beta-block width.",
         call. = FALSE)
  }
  B %*% sqrtm_sym(Omega)
}

#' Numerical rank and candidate basis for one fold moment
#'
#' @keywords internal
#' @noRd
.dkge_cv_fold_basis <- function(Chat, fit, rank) {
  eg <- eigen((Chat + t(Chat)) / 2, symmetric = TRUE)
  scale <- max(eg$values, 0)
  eig_tol <- if (scale > 0) 1e-10 * scale else 0
  available <- min(fit$kernel_rank %||% .dkge_kernel_geometry(fit$K)$rank,
                   sum(eg$values > eig_tol))
  if (rank > available) {
    return(list(basis = NULL, eigen = eg, available_rank = available))
  }
  U <- fit$Kihalf %*% eg$vectors[, seq_len(rank), drop = FALSE]
  U <- dkge_k_orthonormalize(U, fit$K)
  list(basis = U, eigen = eg, available_rank = available)
}

#' Validate candidate ranks
#'
#' @keywords internal
#' @noRd
.dkge_cv_ranks <- function(ranks) {
  if (!is.numeric(ranks) || !length(ranks) || any(!is.finite(ranks)) ||
      any(ranks < 1L) || any(ranks != as.integer(ranks))) {
    stop("`ranks` must contain one or more positive integers.", call. = FALSE)
  }
  sort(unique(as.integer(ranks)))
}

#' Warn when a held-out criterion has saturated
#'
#' @keywords internal
#' @noRd
.dkge_cv_saturation <- function(scores, threshold = 0.999) {
  scores <- scores[is.finite(scores)]
  saturated <- length(scores) >= 2L && all(scores >= threshold)
  if (saturated) {
    .dkge_warn(
      sprintf(
        paste0(
          "All comparable held-out scores are >= %.3f in the fixed validation ",
          "geometry; this criterion provides little discrimination among candidates."
        ),
        threshold
      ),
      "dkge_cv_saturation_warning"
    )
  }
  saturated
}

#' Diagnose spectral concentration relative to a selected latent rank
#'
#' Predictive CV can legitimately prefer a kernel that concentrates most of its
#' mass in fewer directions than the selected latent rank. That is not a rank
#' failure, but it is a separability warning: distinct effect queries may become
#' nearly proportional even though the kernel is algebraically full rank.
#'
#' @keywords internal
#' @noRd
.dkge_cv_validate_kernel_concentration_threshold <- function(threshold) {
  if (!is.null(threshold) &&
      (!is.numeric(threshold) || length(threshold) != 1L ||
       !is.finite(threshold) || threshold <= 0 || threshold > 1)) {
    stop("`kernel_concentration_threshold` must be NULL or one finite scalar in (0, 1].",
         call. = FALSE)
  }
  invisible(NULL)
}

.dkge_cv_kernel_concentration <- function(geometry, rank, threshold = 0.75,
                                          kernel = NULL) {
  .dkge_cv_validate_kernel_concentration_threshold(threshold)
  rank <- as.integer(rank)[[1L]]
  ratio <- if (rank > 0L) geometry$effective_rank_pr / rank else NA_real_
  concentrated <- !is.null(threshold) && rank > 1L && is.finite(ratio) &&
    ratio < threshold
  if (concentrated) {
    label <- if (is.null(kernel)) "Selected kernel" else
      sprintf("Selected kernel %s", shQuote(kernel))
    .dkge_warn(
      sprintf(
        paste0(
          "%s has participation-ratio effective rank %.3f for selected latent ",
          "rank %d (ratio %.3f < warning threshold %.3f). Predictive CV can ",
          "favor stable common directions without preserving contrast ",
          "separability; inspect the kernel spectral diagnostics and run ",
          "`dkge_contrast_diagnostics()` on the planned contrasts."
        ),
        label, geometry$effective_rank_pr, rank, ratio, threshold
      ),
      "dkge_cv_kernel_concentration_warning"
    )
  }
  list(
    concentrated = concentrated,
    effective_rank_pr = geometry$effective_rank_pr,
    selected_rank = rank,
    effective_rank_to_selected_rank = ratio,
    warning_threshold = threshold
  )
}

#' LOSO cross-validation for rank selection
#'
#' Evaluates candidate ranks by recomputing LOSO bases and measuring explained
#' energy on the held-out subject in one candidate-independent validation
#' geometry. By default this is ordinary effect space after pooled-design
#' scaling, not the candidate kernel's own metric.
#'
#' @param B_list List of qxP subject beta matrices.
#' @param X_list List of Txq subject design matrices.
#' @param K qxq design kernel.
#' @param ranks Integer vector of ranks to evaluate.
#' @param Omega_list Optional list of spatial weights.
#' @param ridge Optional ridge parameter passed to [dkge_fit()].
#' @param w_method Subject-level weighting scheme passed to [dkge_fit()].
#' @param w_tau Shrinkage parameter toward equal weights passed to [dkge_fit()].
#' @param spatial Optional fixed [dkge_spatial_regularizer()] applied inside
#'   every candidate fit and training fold.
#' @param validation_K Optional fixed q by q PSD matrix used to score every
#'   candidate. `NULL` (default) uses the identity. It affects scoring only;
#'   candidate `K` still determines the learned training subspace.
#' @param kernel_rank_policy Kernel-support policy for selection. The default,
#'   `"full"`, requires `rank(K) = q`. Use `"allow_singular"` only when every
#'   candidate is intentionally allowed to define a quotient effect space.
#' @param kernel_concentration_threshold Optional warning threshold for spectral
#'   concentration. The default `0.75` warns when the selected kernel's
#'   participation-ratio effective rank divided by the selected latent rank is
#'   below `0.75`. This does not change the CV score or selection; use `NULL` to
#'   disable the warning.
#' @return List containing the one-SE selection (`pick`), the best rank,
#'   aggregated and per-fold score tables, kernel/validation diagnostics,
#'   excluded ranks, and a saturation flag.
#' @examples
#' \donttest{
#' toy <- dkge_sim_toy(
#'   factors = list(cond = list(L = 3)),
#'   active_terms = "cond", S = 4, P = 15, snr = 5
#' )
#' cv <- dkge_cv_rank_loso(toy$B_list, toy$X_list, toy$K, ranks = 1:2)
#' cv$pick
#' }
#' @export
dkge_cv_rank_loso <- function(B_list, X_list, K, ranks,
                              Omega_list = NULL, ridge = 0,
                              w_method = "mfa_sigma1", w_tau = 0.3,
                              validation_K = NULL,
                              kernel_rank_policy = c("full", "allow_singular"),
                              kernel_concentration_threshold = 0.75,
                              spatial = NULL) {
  stopifnot(length(B_list) == length(X_list))
  kernel_rank_policy <- match.arg(kernel_rank_policy)
  .dkge_cv_validate_kernel_concentration_threshold(kernel_concentration_threshold)
  S <- length(B_list)
  q <- nrow(B_list[[1]])
  ranks <- .dkge_cv_ranks(ranks)
  kernel_geometry <- .dkge_cv_candidate_geometry(K, q)
  if (kernel_geometry$rank == 0L) {
    .dkge_abort("Candidate kernel has numerical rank zero.",
                "dkge_kernel_rank_error")
  }
  if (kernel_rank_policy == "full" && kernel_geometry$nullity > 0L) {
    .dkge_abort(
      sprintf(
        paste0(
          "Candidate kernel has numerical rank %d of %d. Full-rank CV is the ",
          "default; use `kernel_rank_policy = \"allow_singular\"` only for an ",
          "intentional quotient effect space."
        ),
        kernel_geometry$rank, q
      ),
      "dkge_cv_kernel_rank_error"
    )
  }
  validation <- .dkge_cv_validation_geometry(validation_K, q)

  kernel_invalid <- ranks > kernel_geometry$rank
  if (any(kernel_invalid)) {
    .dkge_warn(
      sprintf(
        paste0(
          "Dropping rank(s) %s: candidate kernel rank is %d of %d, so those ",
          "latent dimensions do not exist."
        ),
        paste(ranks[kernel_invalid], collapse = ", "), kernel_geometry$rank, q
      ),
      "dkge_cv_rank_warning"
    )
  }
  ranks_fit <- ranks[!kernel_invalid]
  if (!length(ranks_fit)) {
    .dkge_abort(
      sprintf("No requested rank is admissible for kernel rank %d.", kernel_geometry$rank),
      "dkge_cv_rank_error"
    )
  }

  base <- dkge_fit(B_list, X_list, K, Omega_list = Omega_list,
                  w_method = w_method, w_tau = w_tau,
                  ridge = ridge, rank = max(ranks_fit), spatial = spatial)

  rows <- vector("list", S * length(ranks_fit))
  row_id <- 1L
  for (s in seq_len(S)) {
    train_ids <- setdiff(seq_len(S), s)
    ctx <- .dkge_fold_weight_context(base, train_ids, ridge = ridge)
    Bts <- base$Btil[[s]]
    loader_weights <- .dkge_subject_loader_weights(ctx$weights$total, Bts)
    Bw <- if (is.null(loader_weights)) {
      Bts
    } else {
      sweep(Bts, 2L, sqrt(pmax(loader_weights, 0)), "*")
    }
    Bw <- .dkge_apply_fit_spatial(base, Bw, subject = s)
    for (r in ranks_fit) {
      fold <- .dkge_cv_fold_basis(ctx$Chat, base, r)
      score <- if (is.null(fold$basis)) {
        NA_real_
      } else {
        .dkge_cv_score_fixed(Bw, fold$basis, validation$roots)
      }
      rows[[row_id]] <- data.frame(
        subject = s, rank = r, rank_used = if (is.null(fold$basis)) fold$available_rank else r,
        score = score, admissible = !is.null(fold$basis) && is.finite(score)
      )
      row_id <- row_id + 1L
    }
  }
  tab <- do.call(rbind, rows)
  fold_ok <- vapply(ranks_fit, function(r) {
    rows_r <- tab[tab$rank == r, , drop = FALSE]
    nrow(rows_r) == S && all(rows_r$admissible)
  }, logical(1))
  fold_invalid <- ranks_fit[!fold_ok]
  if (length(fold_invalid)) {
    .dkge_warn(
      sprintf(
        "Dropping rank(s) %s because at least one training fold has lower effective rank.",
        paste(fold_invalid, collapse = ", ")
      ),
      "dkge_cv_fold_rank_warning"
    )
  }
  usable_ranks <- ranks_fit[fold_ok]
  if (!length(usable_ranks)) {
    .dkge_abort("No requested rank is estimable in every training fold.",
                "dkge_cv_fold_rank_error")
  }
  usable <- tab[tab$rank %in% usable_ranks, , drop = FALSE]
  sel <- dkge_one_se(usable, param_col = "rank", metric_col = "score")
  sel$summary$rank_used <- sel$summary$param
  saturated <- .dkge_cv_saturation(sel$summary$mean)
  concentration <- .dkge_cv_kernel_concentration(
    kernel_geometry, sel$pick,
    threshold = kernel_concentration_threshold
  )
  list(
    pick = sel$pick,
    best = sel$best,
    table = sel$summary,
    raw = tab,
    kernel = .dkge_kernel_diagnostics(kernel_geometry),
    validation = validation$diagnostics,
    kernel_rank_policy = kernel_rank_policy,
    inadmissible_ranks = sort(unique(c(ranks[kernel_invalid], fold_invalid))),
    saturated = saturated,
    kernel_concentration = concentration
  )
}

#' LOSO kernel grid search
#' 
#' Evaluates a named list of candidate design kernels using LOSO explained
#' energy at a fixed rank. Every candidate is scored in the same fixed
#' validation geometry, so a kernel cannot inflate its score by collapsing the
#' target metric it is judged against.
#'
#' @inheritParams dkge_cv_rank_loso
#' @param K_grid Named list of candidate kernels.
#' @param rank Rank used for evaluation.
#' @return List with the pick, best kernel, candidate audit table (including
#'   kernel rank, nullity, condition, spectral-concentration metrics,
#'   admissibility, and exclusion reason), raw fold scores, fixed validation
#'   diagnostics, saturation flag, and selected-kernel concentration summary.
#' @examples
#' \donttest{
#' toy <- dkge_sim_toy(
#'   factors = list(cond = list(L = 3)),
#'   active_terms = "cond", S = 4, P = 15, snr = 5
#' )
#' q <- nrow(toy$K)
#' K_grid <- list(base = toy$K, identity = diag(q))
#' cv <- dkge_cv_kernel_grid(toy$B_list, toy$X_list, K_grid, rank = 1)
#' cv$pick
#' }
#' @export
dkge_cv_kernel_grid <- function(B_list, X_list, K_grid, rank,
                                Omega_list = NULL, ridge = 0,
                                w_method = "mfa_sigma1", w_tau = 0.3,
                                validation_K = NULL,
                                kernel_rank_policy = c("full", "allow_singular"),
                                kernel_concentration_threshold = 0.75,
                                spatial = NULL) {
  stopifnot(is.list(K_grid), length(K_grid) >= 1)
  kernel_rank_policy <- match.arg(kernel_rank_policy)
  .dkge_cv_validate_kernel_concentration_threshold(kernel_concentration_threshold)
  if (is.null(names(K_grid)) || any(!nzchar(names(K_grid))) || anyDuplicated(names(K_grid))) {
    stop("`K_grid` must have unique, non-empty names.", call. = FALSE)
  }
  rank <- .dkge_cv_ranks(rank)
  if (length(rank) != 1L) stop("`rank` must be one positive integer.", call. = FALSE)
  q <- nrow(B_list[[1]])
  validation <- .dkge_cv_validation_geometry(validation_K, q)
  rows <- list()
  summaries <- list()
  geometries <- list()

  for (nm in names(K_grid)) {
    Kc <- K_grid[[nm]]
    geometry <- .dkge_cv_candidate_geometry(Kc, q, sprintf("Kernel %s", shQuote(nm)))
    geometries[[nm]] <- geometry
    reason <- NA_character_
    if (geometry$rank == 0L) {
      reason <- "kernel rank is zero"
    } else if (kernel_rank_policy == "full" && geometry$nullity > 0L) {
      reason <- sprintf(
        "kernel rank %d is below q = %d under the full-rank policy",
        geometry$rank, q
      )
    } else if (rank > geometry$rank) {
      reason <- sprintf("requested rank %d exceeds kernel rank %d", rank, geometry$rank)
    }
    if (!is.na(reason)) {
      summaries[[nm]] <- data.frame(
        kernel = nm, mean = NA_real_, se = NA_real_,
        kernel_rank = geometry$rank, kernel_nullity = geometry$nullity,
        kernel_condition = geometry$condition,
        kernel_effective_rank_pr = geometry$effective_rank_pr,
        kernel_effective_rank_fraction = geometry$effective_rank_fraction,
        kernel_leading_eigenvalue_share = geometry$leading_eigenvalue_share,
        kernel_effective_rank_to_rank = geometry$effective_rank_pr / rank,
        rank_requested = rank,
        rank_used = min(rank, geometry$rank), admissible = FALSE,
        reason = reason, stringsAsFactors = FALSE
      )
      next
    }

    base <- dkge_fit(B_list, X_list, Kc, Omega_list = Omega_list,
                    w_method = w_method, w_tau = w_tau,
                    ridge = ridge, rank = rank, spatial = spatial)
    S <- length(B_list)
    candidate_rows <- list()

    for (s in seq_len(S)) {
      train_ids <- setdiff(seq_len(S), s)
      ctx <- .dkge_fold_weight_context(base, train_ids, ridge = ridge)
      fold <- .dkge_cv_fold_basis(ctx$Chat, base, rank)

      Bts <- base$Btil[[s]]
      loader_weights <- .dkge_subject_loader_weights(ctx$weights$total, Bts)
      Bw <- if (is.null(loader_weights)) Bts else sweep(Bts, 2L, sqrt(pmax(loader_weights, 0)), "*")
      Bw <- .dkge_apply_fit_spatial(base, Bw, subject = s)
      ev <- if (is.null(fold$basis)) {
        NA_real_
      } else {
        .dkge_cv_score_fixed(Bw, fold$basis, validation$roots)
      }
      candidate_rows[[s]] <- data.frame(
        kernel = nm, subject = s, rank = rank,
        rank_used = if (is.null(fold$basis)) fold$available_rank else rank,
        score = ev, admissible = !is.null(fold$basis) && is.finite(ev)
      )
    }
    candidate_tab <- do.call(rbind, candidate_rows)
    rows[[nm]] <- candidate_tab
    admissible <- nrow(candidate_tab) == S && all(candidate_tab$admissible)
    scores <- candidate_tab$score[candidate_tab$admissible]
    summaries[[nm]] <- data.frame(
      kernel = nm,
      mean = if (admissible) mean(scores) else NA_real_,
      se = if (admissible && length(scores) > 1L) stats::sd(scores) / sqrt(length(scores)) else if (admissible) 0 else NA_real_,
      kernel_rank = geometry$rank,
      kernel_nullity = geometry$nullity,
      kernel_condition = geometry$condition,
      kernel_effective_rank_pr = geometry$effective_rank_pr,
      kernel_effective_rank_fraction = geometry$effective_rank_fraction,
      kernel_leading_eigenvalue_share = geometry$leading_eigenvalue_share,
      kernel_effective_rank_to_rank = geometry$effective_rank_pr / rank,
      rank_requested = rank,
      rank_used = if (admissible) rank else min(candidate_tab$rank_used),
      admissible = admissible,
      reason = if (admissible) NA_character_ else "at least one training fold has lower effective rank",
      stringsAsFactors = FALSE
    )
  }

  tab <- if (length(rows)) do.call(rbind, rows) else data.frame(
    kernel = character(0), subject = integer(0), rank = integer(0),
    rank_used = integer(0), score = numeric(0), admissible = logical(0)
  )
  table <- do.call(rbind, summaries)
  rownames(table) <- NULL
  excluded <- table$kernel[!table$admissible]
  if (length(excluded)) {
    details <- paste0(excluded, " (", table$reason[!table$admissible], ")")
    .dkge_warn(
      sprintf("Excluded inadmissible kernel candidate(s): %s.", paste(details, collapse = "; ")),
      "dkge_cv_kernel_rank_warning"
    )
  }
  eligible <- table[table$admissible & is.finite(table$mean), , drop = FALSE]
  if (!nrow(eligible)) {
    .dkge_abort("No kernel candidate is estimable at the requested rank in every fold.",
                "dkge_cv_kernel_rank_error")
  }
  best_idx <- which.max(eligible$mean)
  # Kernels are nominal: the one-SE "first within tolerance" rule reduces to an
  # arbitrary alphabetical tie-break, so select the best-scoring kernel directly.
  pick_idx <- best_idx
  saturated <- .dkge_cv_saturation(eligible$mean)
  selected_kernel <- eligible$kernel[[pick_idx]]
  concentration <- .dkge_cv_kernel_concentration(
    geometries[[selected_kernel]], rank,
    threshold = kernel_concentration_threshold,
    kernel = selected_kernel
  )

  list(
    pick = selected_kernel,
    best = eligible$kernel[best_idx],
    table = table,
    raw = tab,
    validation = validation$diagnostics,
    kernel_rank_policy = kernel_rank_policy,
    saturated = saturated,
    kernel_concentration = concentration
  )
}

#' LOSO selection of spatial regularization strength
#'
#' Selects the Laplacian penalty `lambda` while keeping the spatial graph,
#' design kernel, latent rank, and validation geometry fixed. Each candidate is
#' fitted inside each LOSO training fold, but scored against the **unsmoothed**
#' held-out beta block in a candidate-independent effect-space geometry. Fixed
#' effect scaling, adaptive location weights, and `Omega_list` remain part of
#' that held-out geometry; only the candidate Laplacian solve is omitted. Thus a
#' larger `lambda` cannot improve its score merely by smoothing or shrinking the
#' same field used in the denominator.
#'
#' The returned `pick` is the largest (smoothest) candidate whose mean score is
#' within one standard error of the best candidate. Include `0` in `lambdas` to
#' compare against the unsmoothed model. A saturated score is reported and
#' warned about in the same way as other DKGE CV helpers.
#'
#' Spatial CV fails closed when any positive candidate is requested but every
#' supplied graph is edgeless: all candidate resolvents would be the identity,
#' so `lambda` is not identifiable from the score. If only some
#' subject-specific graphs are edgeless, CV warns once with their identifiers
#' and continues with status `partial`.
#'
#' @inheritParams dkge_cv_rank_loso
#' @param spatial A [dkge_spatial_regularizer()] supplying the fixed graph. Its
#'   stored `lambda` is ignored while evaluating `lambdas`.
#' @param lambdas Non-negative finite candidate penalties.
#' @param rank One positive latent rank used for every candidate.
#' @param effect_scaling Effect-space scaling passed to every candidate
#'   [dkge_fit()]. Use the same setting planned for the final refit.
#' @return A list with selected `pick`, unconstrained `best`, summary `table`,
#'   per-fold `raw` scores, the selected `spatial` specification, validation and
#'   kernel diagnostics, the one-SE threshold, recorded `fit_settings`, and a
#'   saturation flag.
#' @export
#' @examples
#' \donttest{
#' toy <- dkge_sim_toy(
#'   factors = list(cond = list(L = 3)), active_terms = "cond",
#'   S = 4, P = 15, snr = 4
#' )
#' coords <- cbind(x = seq_len(ncol(toy$B_list[[1]])), y = 0, z = 0)
#' spatial <- dkge_spatial_regularizer(
#'   coords, lambda = 1, dthresh = 1.01, nnk = 3,
#'   weight_mode = "binary"
#' )
#' cv <- dkge_cv_spatial_grid(
#'   toy$B_list, toy$X_list, toy$K, spatial,
#'   lambdas = c(0, 0.25, 1), rank = 1
#' )
#' cv$pick
#' }
dkge_cv_spatial_grid <- function(B_list, X_list, K, spatial, lambdas, rank,
                                 Omega_list = NULL, ridge = 0,
                                 w_method = "mfa_sigma1", w_tau = 0.3,
                                 effect_scaling = c("pooled_design", "none"),
                                 validation_K = NULL,
                                 kernel_rank_policy = c("full", "allow_singular"),
                                 kernel_concentration_threshold = 0.75) {
  stopifnot(is.list(B_list), is.list(X_list),
            length(B_list) == length(X_list), length(B_list) >= 2L)
  if (!inherits(spatial, "dkge_spatial_regularizer")) {
    .dkge_abort("`spatial` must be created by `dkge_spatial_regularizer()`.",
                "dkge_spatial_spec_error")
  }
  .dkge_spatial_spec_payload(spatial, validate = TRUE)
  if (!is.numeric(lambdas) || !length(lambdas) ||
      any(!is.finite(lambdas)) || any(lambdas < 0)) {
    .dkge_abort("`lambdas` must contain non-negative finite numbers.",
                "dkge_spatial_spec_error")
  }
  lambdas <- sort(unique(as.numeric(lambdas)))
  cv_topology <- .dkge_spatial_topology(
    spatial$laplacians,
    if (any(lambdas > 0)) max(lambdas) else 0
  )
  if (identical(cv_topology$status, "inert")) {
    .dkge_abort(
      paste0(
        "Spatial CV cannot distinguish positive `lambda` candidates because ",
        "all supplied graphs have no edges: every resolvent is the identity. ",
        "Increase `dthresh` to the scale of `coords` or supply a connected ",
        "Laplacian before tuning."
      ),
      "dkge_cv_spatial_inert_error"
    )
  }
  if (identical(cv_topology$status, "partial")) {
    .dkge_spatial_warn_topology(
      cv_topology,
      spatial$source,
      spatial$construction$dthresh %||% NULL
    )
  }
  topology_checked <- identical(cv_topology$status, "partial")
  rank <- .dkge_cv_ranks(rank)
  if (length(rank) != 1L) stop("`rank` must be one positive integer.", call. = FALSE)
  effect_scaling <- match.arg(effect_scaling)

  kernel_rank_policy <- match.arg(kernel_rank_policy)
  .dkge_cv_validate_kernel_concentration_threshold(kernel_concentration_threshold)
  q <- nrow(B_list[[1]])
  geometry <- .dkge_cv_candidate_geometry(K, q)
  if (geometry$rank == 0L) {
    .dkge_abort("Candidate kernel has numerical rank zero.",
                "dkge_kernel_rank_error")
  }
  if (kernel_rank_policy == "full" && geometry$nullity > 0L) {
    .dkge_abort(
      sprintf(
        paste0(
          "Candidate kernel has numerical rank %d of %d. Full-rank CV is the ",
          "default; use `kernel_rank_policy = \"allow_singular\"` only for an ",
          "intentional quotient effect space."
        ),
        geometry$rank, q
      ),
      "dkge_cv_kernel_rank_error"
    )
  }
  if (rank > geometry$rank) {
    .dkge_abort(
      sprintf("Requested rank %d exceeds kernel rank %d.", rank, geometry$rank),
      "dkge_cv_rank_error"
    )
  }
  validation <- .dkge_cv_validation_geometry(validation_K, q)
  S <- length(B_list)
  candidate_rows <- vector("list", length(lambdas))

  for (i in seq_along(lambdas)) {
    lambda <- lambdas[[i]]
    candidate <- .dkge_spatial_with_lambda(
      spatial, lambda, topology_checked = topology_checked
    )
    base <- dkge_fit(
      B_list, X_list, K,
      Omega_list = Omega_list,
      w_method = w_method,
      w_tau = w_tau,
      ridge = ridge,
      rank = rank,
      effect_scaling = effect_scaling,
      spatial = candidate
    )
    fold_rows <- vector("list", S)
    for (s in seq_len(S)) {
      train_ids <- setdiff(seq_len(S), s)
      ctx <- .dkge_fold_weight_context(base, train_ids, ridge = ridge)
      fold <- .dkge_cv_fold_basis(ctx$Chat, base, rank)

      # Deliberately do not call `.dkge_apply_fit_spatial()` here. Candidate
      # lambda changes the training basis, while every candidate is judged on
      # the same raw held-out field and fixed validation metric.
      Bheld <- base$Btil[[s]]
      loader_weights <- .dkge_subject_loader_weights(ctx$weights$total, Bheld)
      if (!is.null(loader_weights)) {
        Bheld <- sweep(Bheld, 2L, sqrt(pmax(loader_weights, 0)), "*")
      }
      Bheld <- .dkge_cv_apply_spatial_metric(Bheld, Omega_list[[s]])
      score <- if (is.null(fold$basis)) {
        NA_real_
      } else {
        .dkge_cv_score_fixed(Bheld, fold$basis, validation$roots)
      }
      fold_rows[[s]] <- data.frame(
        lambda = lambda,
        subject = s,
        rank = rank,
        rank_used = if (is.null(fold$basis)) fold$available_rank else rank,
        score = score,
        admissible = !is.null(fold$basis) && is.finite(score),
        stringsAsFactors = FALSE
      )
    }
    candidate_rows[[i]] <- do.call(rbind, fold_rows)
  }

  raw <- do.call(rbind, candidate_rows)
  table <- do.call(rbind, lapply(lambdas, function(lambda) {
    rows <- raw[raw$lambda == lambda, , drop = FALSE]
    admissible <- nrow(rows) == S && all(rows$admissible)
    scores <- rows$score[rows$admissible]
    data.frame(
      lambda = lambda,
      mean = if (admissible) mean(scores) else NA_real_,
      se = if (admissible && length(scores) > 1L) {
        stats::sd(scores) / sqrt(length(scores))
      } else if (admissible) {
        0
      } else {
        NA_real_
      },
      rank = rank,
      rank_used = if (admissible) rank else min(rows$rank_used),
      admissible = admissible,
      reason = if (admissible) NA_character_ else
        "at least one training fold has lower effective rank",
      stringsAsFactors = FALSE
    )
  }))
  eligible <- table[table$admissible & is.finite(table$mean), , drop = FALSE]
  if (!nrow(eligible)) {
    .dkge_abort("No spatial candidate is estimable in every training fold.",
                "dkge_cv_spatial_rank_error")
  }
  best_idx <- which.max(eligible$mean)
  best <- eligible$lambda[[best_idx]]
  threshold <- eligible$mean[[best_idx]] - eligible$se[[best_idx]]
  one_se <- eligible[eligible$mean >= threshold, , drop = FALSE]
  pick <- max(one_se$lambda)
  saturated <- .dkge_cv_saturation(eligible$mean)
  concentration <- .dkge_cv_kernel_concentration(
    geometry, rank, threshold = kernel_concentration_threshold
  )

  list(
    pick = pick,
    best = best,
    table = table,
    raw = raw,
    spatial = .dkge_spatial_with_lambda(
      spatial, pick, topology_checked = topology_checked
    ),
    validation = validation$diagnostics,
    heldout_geometry = "raw_beta_block",
    heldout_spatial_metric = if (is.null(Omega_list)) "identity" else
      "Omega_list",
    selection_rule = "largest_lambda_within_one_se",
    one_se_threshold = threshold,
    fit_settings = list(
      rank = rank,
      ridge = ridge,
      w_method = w_method,
      w_tau = w_tau,
      effect_scaling = effect_scaling
    ),
    kernel = .dkge_kernel_diagnostics(geometry),
    kernel_rank_policy = kernel_rank_policy,
    saturated = saturated,
    kernel_concentration = concentration
  )
}

#' Pooled design-space covariance and Cholesky factor
#'
#' Computes the qxq pooled design covariance and the corresponding Cholesky factor
#' needed for fast kernel alignment screening.
#'
#' @inheritParams dkge_cv_rank_loso
#' @return List containing `C` (pooled covariance in the ruler metric), `R`
#'   (upper-triangular Cholesky factor), and `G` (pooled design Gram matrix).
#' @examples
#' toy <- dkge_sim_toy(
#'   factors = list(cond = list(L = 3)),
#'   active_terms = "cond", S = 3, P = 10, snr = 5
#' )
#' pooled <- dkge_pooled_cov_q(toy$B_list, toy$X_list)
#' dim(pooled$C)
#' @keywords internal
#' @export
dkge_pooled_cov_q <- function(B_list, X_list, Omega_list = NULL,
                              spatial = NULL) {
  stopifnot(length(B_list) == length(X_list))
  S <- length(B_list)
  q <- nrow(B_list[[1]])

  G_list <- lapply(X_list, crossprod)
  G <- Reduce(`+`, G_list)
  diag(G) <- diag(G) + 1e-10
  R <- chol(G)

  if (is.null(Omega_list)) Omega_list <- vector("list", S)
  spatial_fit <- .dkge_resolve_spatial(
    spatial, B_list,
    subject_ids = names(B_list) %||% paste0("subject", seq_len(S))
  )

  C <- matrix(0, q, q)
  for (s in seq_len(S)) {
    Bt <- t(R) %*% B_list[[s]]
    Bt <- .dkge_spatial_apply_betas(
      Bt, spatial_fit$operators[[s]] %||% NULL
    )
    Omega <- Omega_list[[s]]
    if (is.null(Omega)) {
      C <- C + Bt %*% t(Bt)
    } else if (is.vector(Omega)) {
      C <- C + (Bt * rep(Omega, each = q)) %*% t(Bt)
    } else {
      C <- C + Bt %*% Omega %*% t(Bt)
    }
  }
  C <- (C + t(C)) / 2
  list(C = C, R = R, G = G)
}

#' Kernel alignment pre-screening
#'
#' Ranks candidate kernels by their alignment with the pooled design-space
#' covariance produced by [dkge_pooled_cov_q()].
#'
#' @param K_grid Named list of qxq kernels.
#' @param C Pooled covariance matrix.
#' @param normalize_k Logical; if `TRUE`, kernels are scaled to unit trace before
#'   alignment.
#' @param top_k Number of kernels to retain.
#' @return Data frame sorted by decreasing alignment; the `top` attribute carries
#'   the names of the retained kernels.
#' @examples
#' toy <- dkge_sim_toy(
#'   factors = list(cond = list(L = 3)),
#'   active_terms = "cond", S = 3, P = 10, snr = 5
#' )
#' pooled <- dkge_pooled_cov_q(toy$B_list, toy$X_list)
#' q <- nrow(toy$K)
#' K_grid <- list(base = toy$K, identity = diag(q))
#' dkge_kernel_prescreen(K_grid, pooled$C, top_k = 1)
#' @export
dkge_kernel_prescreen <- function(K_grid, C, normalize_k = TRUE, top_k = 3) {
  stopifnot(is.list(K_grid), length(K_grid) >= 1)
  scores <- lapply(names(K_grid), function(nm) {
    K <- K_grid[[nm]]
    if (normalize_k) {
      tr <- sum(diag(K))
      if (tr > 0) K <- K / tr
    }
    data.frame(kernel = nm, align = dkge_kernel_alignment(K, C))
  })
  tab <- do.call(rbind, scores)
  tab <- tab[order(-tab$align), , drop = FALSE]
  attr(tab, "top") <- head(tab$kernel, min(top_k, nrow(tab)))
  tab
}

#' Combined kernel and rank selection via pre-screening and LOSO CV
#'
#' Runs kernel alignment pre-screening followed by LOSO explained-variance
#' cross-validation, applying the one-standard-error rule to pick a kernel/rank
#' pair.
#'
#' @inheritParams dkge_cv_rank_loso
#' @param K_grid Named list of candidate kernels.
#' @param ranks Integer vector of ranks to evaluate.
#' @param top_k Number of kernels to keep after pre-screening.
#' @return List with the selected `kernel` and `rank`, alignment and CV tables
#'   carrying kernel-support and spectral-concentration diagnostics, per-kernel
#'   selections, exclusions, fixed validation diagnostics, saturation flag, and
#'   selected-kernel concentration summary.
#' @examples
#' \donttest{
#' toy <- dkge_sim_toy(
#'   factors = list(cond = list(L = 3)),
#'   active_terms = "cond", S = 4, P = 15, snr = 5
#' )
#' q <- nrow(toy$K)
#' K_grid <- list(base = toy$K, identity = diag(q))
#' sel <- dkge_cv_kernel_rank(toy$B_list, toy$X_list, K_grid, ranks = 1:2)
#' sel$pick
#' }
#' @export
dkge_cv_kernel_rank <- function(B_list, X_list, K_grid, ranks,
                                Omega_list = NULL, ridge = 0,
                                w_method = "mfa_sigma1", w_tau = 0.3,
                                top_k = 3, validation_K = NULL,
                                kernel_rank_policy = c("full", "allow_singular"),
                                kernel_concentration_threshold = 0.75,
                                spatial = NULL) {
  stopifnot(is.list(K_grid), length(K_grid) >= 1)
  kernel_rank_policy <- match.arg(kernel_rank_policy)
  .dkge_cv_validate_kernel_concentration_threshold(kernel_concentration_threshold)
  if (is.null(names(K_grid)) || any(!nzchar(names(K_grid))) || anyDuplicated(names(K_grid))) {
    stop("`K_grid` must have unique, non-empty names.", call. = FALSE)
  }
  if (!is.numeric(top_k) || length(top_k) != 1L || !is.finite(top_k) ||
      top_k < 1L || top_k != as.integer(top_k)) {
    stop("`top_k` must be one positive integer.", call. = FALSE)
  }
  top_k <- as.integer(top_k)
  ranks <- .dkge_cv_ranks(ranks)
  q <- nrow(B_list[[1]])
  validation <- .dkge_cv_validation_geometry(validation_K, q)
  kernel_geometries <- lapply(names(K_grid), function(nm) {
    .dkge_cv_candidate_geometry(K_grid[[nm]], q, sprintf("Kernel %s", shQuote(nm)))
  })
  names(kernel_geometries) <- names(K_grid)

  exclusion_reason <- vapply(names(K_grid), function(nm) {
    geometry <- kernel_geometries[[nm]]
    if (geometry$rank == 0L) {
      return("kernel rank is zero")
    }
    if (kernel_rank_policy == "full" && geometry$nullity > 0L) {
      return(sprintf(
        "kernel rank %d is below q = %d under the full-rank policy",
        geometry$rank, q
      ))
    }
    if (!any(ranks <= geometry$rank)) {
      return(sprintf("all requested ranks exceed kernel rank %d", geometry$rank))
    }
    NA_character_
  }, character(1))

  pooled <- dkge_pooled_cov_q(B_list, X_list, Omega_list,
                              spatial = spatial)
  alignment <- dkge_kernel_prescreen(
    K_grid, pooled$C, normalize_k = TRUE, top_k = length(K_grid)
  )
  alignment$kernel_rank <- vapply(alignment$kernel, function(nm) {
    kernel_geometries[[nm]]$rank
  }, integer(1))
  alignment$kernel_nullity <- vapply(alignment$kernel, function(nm) {
    kernel_geometries[[nm]]$nullity
  }, integer(1))
  alignment$kernel_condition <- vapply(alignment$kernel, function(nm) {
    kernel_geometries[[nm]]$condition
  }, numeric(1))
  alignment$kernel_effective_rank_pr <- vapply(alignment$kernel, function(nm) {
    kernel_geometries[[nm]]$effective_rank_pr
  }, numeric(1))
  alignment$kernel_effective_rank_fraction <- vapply(alignment$kernel, function(nm) {
    kernel_geometries[[nm]]$effective_rank_fraction
  }, numeric(1))
  alignment$kernel_leading_eigenvalue_share <- vapply(alignment$kernel, function(nm) {
    kernel_geometries[[nm]]$leading_eigenvalue_share
  }, numeric(1))
  alignment$rank_policy_admissible <- is.na(exclusion_reason[alignment$kernel])
  eligible_screen <- alignment$kernel[alignment$rank_policy_admissible]
  keep <- head(eligible_screen, min(top_k, length(eligible_screen)))
  attr(alignment, "top") <- keep

  cv_rows <- list()
  picks <- list()
  excluded <- as.list(exclusion_reason[!is.na(exclusion_reason)])

  for (i in seq_along(keep)) {
    nm <- keep[[i]]
    geometry <- kernel_geometries[[nm]]
    candidate_ranks <- ranks[ranks <= geometry$rank]
    cv <- dkge_cv_rank_loso(
      B_list, X_list, K_grid[[nm]], candidate_ranks,
      Omega_list = Omega_list, ridge = ridge,
      w_method = w_method, w_tau = w_tau,
      validation_K = validation_K,
      kernel_rank_policy = kernel_rank_policy,
      kernel_concentration_threshold = NULL,
      spatial = spatial
    )
    summary_tbl <- cv$table %||% cv$summary
    stopifnot(!is.null(summary_tbl), all(c("param", "mean", "se") %in% names(summary_tbl)))
    tmp <- data.frame(kernel = nm,
                      rank = summary_tbl$param,
                      mean = summary_tbl$mean,
                      se = summary_tbl$se,
                      kernel_rank = geometry$rank,
                      kernel_nullity = geometry$nullity,
                      kernel_condition = geometry$condition,
                      kernel_effective_rank_pr = geometry$effective_rank_pr,
                      kernel_effective_rank_fraction = geometry$effective_rank_fraction,
                      kernel_leading_eigenvalue_share = geometry$leading_eigenvalue_share,
                      kernel_effective_rank_to_rank = geometry$effective_rank_pr /
                        summary_tbl$param)
    cv_rows[[nm]] <- tmp
    idx <- tmp$rank == cv$pick
    picks[[nm]] <- data.frame(kernel = nm,
                              rank = cv$pick,
                              score = tmp$mean[idx],
                              se = tmp$se[idx],
                              kernel_rank = geometry$rank,
                              kernel_nullity = geometry$nullity,
                              kernel_condition = geometry$condition,
                              kernel_effective_rank_pr = geometry$effective_rank_pr,
                              kernel_effective_rank_fraction = geometry$effective_rank_fraction,
                              kernel_leading_eigenvalue_share = geometry$leading_eigenvalue_share,
                              kernel_effective_rank_to_rank = geometry$effective_rank_pr /
                                cv$pick)
  }

  if (length(excluded)) {
    detail <- paste0(names(excluded), " (", unlist(excluded, use.names = FALSE), ")")
    .dkge_warn(
      sprintf("Excluded inadmissible kernel candidate(s): %s.", paste(detail, collapse = "; ")),
      "dkge_cv_kernel_rank_warning"
    )
  }
  if (!length(cv_rows) || !length(picks)) {
    .dkge_abort("No kernel-rank candidate is admissible.",
                "dkge_cv_kernel_rank_error")
  }
  cv_table <- do.call(rbind, cv_rows)
  pick_df <- do.call(rbind, picks)
  rownames(cv_table) <- NULL
  rownames(pick_df) <- NULL

  best_idx <- which.max(pick_df$score)
  threshold <- pick_df$score[best_idx] - pick_df$se[best_idx]
  candidates <- pick_df[pick_df$score >= threshold, , drop = FALSE]
  selected <- candidates[order(candidates$rank, -candidates$score), ][1, ]
  saturated <- .dkge_cv_saturation(pick_df$score)
  concentration <- .dkge_cv_kernel_concentration(
    kernel_geometries[[selected$kernel]], selected$rank,
    threshold = kernel_concentration_threshold,
    kernel = selected$kernel
  )

  list(
    pick = list(kernel = selected$kernel, rank = selected$rank),
    tables = list(alignment = alignment, cv = cv_table),
    picks_per_kernel = pick_df,
    excluded = excluded,
    validation = validation$diagnostics,
    kernel_rank_policy = kernel_rank_policy,
    saturated = saturated,
    kernel_concentration = concentration
  )
}
