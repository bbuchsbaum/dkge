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

#' Summarise DKGE diagnostics
#'
#' Provides a compact list of variance explained, subject weights, and rank
#' metadata for quick inspection.
#'
#' @param fit A `dkge` object.
#' @return List with variance table, subject weights, and rank info.
#' @export
dkge_diagnostics <- function(fit) {
  stopifnot(inherits(fit, "dkge"))
  voxel_stats <- if (!is.null(fit$voxel_weights)) {
    vw <- fit$voxel_weights
    list(mean = mean(vw), sd = stats::sd(vw), min = min(vw), max = max(vw))
  } else NULL
  list(
    variance = dkge_variance_explained(fit),
    weights = fit$weights,
    rank = fit$rank,
    q = nrow(fit$U),
    n_subjects = length(fit$Btil),
    voxel_weights = voxel_stats,
    weight_spec = fit$weight_spec
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

#' LOSO cross-validation for rank selection
#'
#' Evaluates candidate ranks by recomputing LOSO bases and measuring explained
#' variance on the held-out subject in the K^{1/2} metric.
#'
#' @param B_list List of qxP subject beta matrices.
#' @param X_list List of Txq subject design matrices.
#' @param K qxq design kernel.
#' @param ranks Integer vector of ranks to evaluate.
#' @param Omega_list Optional list of spatial weights.
#' @param ridge Optional ridge parameter passed to [dkge_fit()].
#' @param w_method Subject-level weighting scheme passed to [dkge_fit()].
#' @param w_tau Shrinkage parameter toward equal weights passed to [dkge_fit()].
#' @return List containing the one-SE selection (`pick`), the best rank, and the
#'   aggregated score table.
#' @export
dkge_cv_rank_loso <- function(B_list, X_list, K, ranks,
                              Omega_list = NULL, ridge = 0,
                              w_method = "mfa_sigma1", w_tau = 0.3) {
  stopifnot(length(B_list) == length(X_list), length(ranks) >= 1)
  S <- length(B_list)
  q <- nrow(B_list[[1]])

  fit_fun <- if (exists("dkge_fit_fast")) get("dkge_fit_fast") else dkge_fit
  base <- fit_fun(B_list, X_list, K, Omega_list = Omega_list,
                  w_method = w_method, w_tau = w_tau,
                  ridge = ridge, rank = max(ranks))
  Khalf <- base$Khalf

  rows <- vector("list", S * length(ranks))
  row_id <- 1L
  for (s in seq_len(S)) {
    train_ids <- setdiff(seq_len(S), s)
    ctx <- .dkge_fold_weight_context(base, train_ids, ridge = ridge)
    eg <- eigen(ctx$Chat, symmetric = TRUE)
    loader_weights <- ctx$weights$total
    for (r in ranks) {
      Uminus <- base$Kihalf %*% eg$vectors[, seq_len(r), drop = FALSE]
      score <- {
        Bts <- base$Btil[[s]]
        Bw <- if (is.null(loader_weights)) Bts else sweep(Bts, 2L, sqrt(pmax(loader_weights, 0)), "*")
        Xs <- base$Khalf %*% Bw
        V <- base$Khalf %*% Uminus
        Xhat <- V %*% (t(Uminus) %*% base$K %*% Bw)
        sum(Xhat * Xhat) / (sum(Xs * Xs) + 1e-12)
      }
      rows[[row_id]] <- data.frame(subject = s, rank = r, score = score)
      row_id <- row_id + 1L
    }
  }
  tab <- do.call(rbind, rows)
  sel <- dkge_one_se(tab, param_col = "rank", metric_col = "score")
  list(pick = sel$pick, best = sel$best, table = sel$summary, raw = tab)
}

#' LOSO kernel grid search
#' 
#' Evaluates a named list of candidate design kernels using LOSO explained
#' variance at a fixed rank.
#'
#' @inheritParams dkge_cv_rank_loso
#' @param K_grid Named list of candidate kernels.
#' @param rank Rank used for evaluation.
#' @return List with the one-SE pick, best kernel, summary table, and raw scores.
#' @export
dkge_cv_kernel_grid <- function(B_list, X_list, K_grid, rank,
                                Omega_list = NULL, ridge = 0,
                                w_method = "mfa_sigma1", w_tau = 0.3) {
  stopifnot(is.list(K_grid), length(K_grid) >= 1)
  rows <- list()

  for (nm in names(K_grid)) {
    Kc <- K_grid[[nm]]
    base <- dkge_fit(B_list, X_list, Kc, Omega_list = Omega_list,
                     w_method = w_method, w_tau = w_tau,
                     ridge = ridge, rank = rank)
    Khalf <- base$Khalf
    S <- length(B_list)

    for (s in seq_len(S)) {
      train_ids <- setdiff(seq_len(S), s)
      ctx <- .dkge_fold_weight_context(base, train_ids, ridge = ridge)
      eg <- eigen(ctx$Chat, symmetric = TRUE)
      loader_weights <- ctx$weights$total
      Uminus <- base$Kihalf %*% eg$vectors[, seq_len(rank), drop = FALSE]

      Bts <- base$Btil[[s]]
      Bw <- if (is.null(loader_weights)) Bts else sweep(Bts, 2L, sqrt(pmax(loader_weights, 0)), "*")
      V <- Khalf %*% Uminus
      Xs <- Khalf %*% Bw
      Xhat <- V %*% (t(Uminus) %*% base$K %*% Bw)
      ev <- sum(Xhat * Xhat) / (sum(Xs * Xs) + 1e-12)
      rows[[length(rows) + 1L]] <- data.frame(kernel = nm, subject = s, score = ev)
    }
  }

  tab <- do.call(rbind, rows)
  agg <- aggregate(tab$score, by = list(tab$kernel),
                   FUN = function(x) c(mean = mean(x), se = stats::sd(x) / sqrt(length(x))))
  kernels <- agg[[1]]
  stats_raw <- agg[[2]]
  stats_mat <- if (is.list(stats_raw)) do.call(rbind, stats_raw) else as.matrix(stats_raw)
  if (is.null(colnames(stats_mat))) {
    colnames(stats_mat) <- c("mean", "se")
  }
  means <- stats_mat[, "mean"]
  ses <- stats_mat[, "se"]
  best_idx <- which.max(means)
  threshold <- means[best_idx] - ses[best_idx]
  pick_idx <- min(which(means >= threshold))

  list(
    pick = kernels[pick_idx],
    best = kernels[best_idx],
    table = data.frame(kernel = kernels, mean = means, se = ses),
    raw = tab
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
#' @export
dkge_pooled_cov_q <- function(B_list, X_list, Omega_list = NULL) {
  stopifnot(length(B_list) == length(X_list))
  S <- length(B_list)
  q <- nrow(B_list[[1]])

  G_list <- lapply(X_list, crossprod)
  G <- Reduce(`+`, G_list)
  diag(G) <- diag(G) + 1e-10
  R <- chol(G)

  if (is.null(Omega_list)) Omega_list <- vector("list", S)

  C <- matrix(0, q, q)
  for (s in seq_len(S)) {
    Bt <- t(R) %*% B_list[[s]]
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
#' @return List with the selected `kernel` and `rank`, alignment and CV tables,
#'   and a per-kernel summary of scores at the selected rank.
#' @export
dkge_cv_kernel_rank <- function(B_list, X_list, K_grid, ranks,
                                Omega_list = NULL, ridge = 0,
                                w_method = "mfa_sigma1", w_tau = 0.3,
                                top_k = 3) {
  stopifnot(is.list(K_grid), length(K_grid) >= 1)

  pooled <- dkge_pooled_cov_q(B_list, X_list, Omega_list)
  alignment <- dkge_kernel_prescreen(K_grid, pooled$C, normalize_k = TRUE, top_k = top_k)
  keep <- attr(alignment, "top")

  cv_rows <- list()
  picks <- vector("list", length(keep))

  for (i in seq_along(keep)) {
    nm <- keep[[i]]
    cv <- dkge_cv_rank_loso(B_list, X_list, K_grid[[nm]], ranks,
                             Omega_list = Omega_list, ridge = ridge,
                             w_method = w_method, w_tau = w_tau)
    tmp <- data.frame(kernel = nm,
                      rank = cv$summary$param,
                      mean = cv$summary$score,
                      se = cv$summary$se)
    cv_rows[[i]] <- tmp
    idx <- tmp$rank == cv$pick
    picks[[i]] <- data.frame(kernel = nm,
                             rank = cv$pick,
                             score = tmp$mean[idx],
                             se = tmp$se[idx])
  }

  cv_table <- do.call(rbind, cv_rows)
  pick_df <- do.call(rbind, picks)

  best_idx <- which.max(pick_df$score)
  threshold <- pick_df$score[best_idx] - pick_df$se[best_idx]
  candidates <- pick_df[pick_df$score >= threshold, , drop = FALSE]
  selected <- candidates[order(candidates$rank, -candidates$score), ][1, ]

  list(
    pick = list(kernel = selected$kernel, rank = selected$rank),
    tables = list(alignment = alignment, cv = cv_table),
    picks_per_kernel = pick_df
  )
}
