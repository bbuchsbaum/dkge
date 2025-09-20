
# dkge-cv.R
# Cross-validation helpers for rank and kernel selection.

#' One-SE selection
#' @param grid data.frame with columns metric (higher better) and param
#' @return list(best, pick) -- best param and one-SE pick
#' @export
dkge_one_se <- function(scores, param_col = "param", metric_col = "score", se_col = NULL) {
  # scores: data.frame with per-fold metrics; compute mean, se per param
  agg <- aggregate(scores[[metric_col]], by = list(scores[[param_col]]), FUN = function(x) c(mean=mean(x), se=sd(x)/sqrt(length(x))))
  # unpack
  param <- agg[[1]]
  m  <- sapply(agg[[2]], function(v) v[["mean"]])
  se <- sapply(agg[[2]], function(v) v[["se"]])
  best_idx <- which.max(m)
  thr <- m[best_idx] - se[best_idx]
  pick_idx <- min(which(m >= thr))
  list(best = param[best_idx], pick = param[pick_idx], summary = data.frame(param=param, mean=m, se=se))
}

#' LOSO CV for rank: explained variance on held-out subject
#' @param B_list, X_list, K as in dkge_fit
#' @param ranks integer vector to try
#' @param ridge, w_method, w_tau passed to dkge_fit
#' @return list(pick, table)
#' @export
dkge_cv_rank_loso <- function(B_list, X_list, K, ranks,
                              Omega_list = NULL, ridge = 0,
                              w_method = "mfa_sigma1", w_tau = 0.3) {
  S <- length(B_list)
  q <- nrow(B_list[[1]])
  # pre-fit once to get pooled design Cholesky and K roots & per-subject contribs
  base <- dkge_fit(B_list, X_list, K, Omega_list, w_method=w_method, w_tau=w_tau, ridge=ridge, rank=max(ranks))
  Khalf <- base$Khalf

  # helper: explained variance on subject s given U_minus (rank r)
  ev_subject <- function(s, Uminus) {
    Bts <- base$Btil[[s]]             # q×P
    Xs  <- Khalf %*% Bts              # q×P in K^{1/2} metric
    # Projection: P_K = U U' K  ⇒ K^{1/2} P_K = (K^{1/2} U) (U' K)
    V   <- base$Khalf %*% Uminus      # q×r
    Xhat <- V %*% (t(Uminus) %*% base$K %*% Bts)  # q×P in K^{1/2} space
    num <- sum(Xhat * Xhat)
    den <- sum(Xs * Xs) + 1e-12
    num / den
  }

  rows <- list()
  for (s in seq_len(S)) {
    # LOSO: rebuild Chat without s, then cut to each rank
    Chat_minus <- base$Chat - base$weights[s]*base$contribs[[s]]
    if (ridge > 0) Chat_minus <- Chat_minus + ridge * diag(q)
    eg <- eigen((Chat_minus + t(Chat_minus))/2, symmetric = TRUE)
    for (r in ranks) {
      Uminus <- base$Kihalf %*% eg$vectors[, seq_len(r), drop = FALSE]
      ev <- ev_subject(s, Uminus)
      rows[[length(rows)+1L]] <- data.frame(subject=s, rank=r, score=ev)
    }
  }
  tab <- do.call(rbind, rows)
  sel <- dkge_one_se(tab, param_col="rank", metric_col="score")
  list(pick = sel$pick, best = sel$best, table = sel$summary, raw = tab)
}

#' Grid search for kernel hyper-parameters with LOSO CV
#' @param K_grid list of candidate K matrices (named)
#' @param ... forwarded to dkge_cv_rank_loso (use a fixed rank or a picked rank)
#' @export
dkge_cv_kernel_grid <- function(B_list, X_list, K_grid, rank,
                                Omega_list = NULL, ridge = 0,
                                w_method = "mfa_sigma1", w_tau = 0.3) {
  rows <- list()
  for (nm in names(K_grid)) {
    K <- K_grid[[nm]]
    # Fit base to get ruler and K roots
    base <- dkge_fit(B_list, X_list, K, Omega_list, w_method=w_method, w_tau=w_tau, ridge=ridge, rank=rank)
    Khalf <- base$Khalf
    S <- length(B_list); q <- nrow(B_list[[1]])
    for (s in seq_len(S)) {
      Chat_minus <- base$Chat - base$weights[s]*base$contribs[[s]]
      if (ridge > 0) Chat_minus <- Chat_minus + ridge * diag(q)
      eg <- eigen((Chat_minus + t(Chat_minus))/2, symmetric = TRUE)
      Uminus <- base$Kihalf %*% eg$vectors[, seq_len(rank), drop = FALSE]
      # Explained variance on s
      Bts <- base$Btil[[s]]
      V   <- Khalf %*% Uminus
      Xs  <- Khalf %*% Bts
      Xhat <- V %*% (t(Uminus) %*% base$K %*% Bts)
      ev <- sum(Xhat*Xhat) / (sum(Xs*Xs) + 1e-12)
      rows[[length(rows)+1L]] <- data.frame(kernel=nm, subject=s, score=ev)
    }
  }
  tab <- do.call(rbind, rows)
  agg <- aggregate(tab$score, by=list(tab$kernel), FUN=function(x) c(mean=mean(x), se=sd(x)/sqrt(length(x))))
  kernel <- agg[[1]]
  m <- sapply(agg[[2]], function(v) v[["mean"]]); se <- sapply(agg[[2]], function(v) v[["se"]])
  best <- kernel[which.max(m)]
  pick <- kernel[min(which(m >= max(m) - se[which.max(m)]))]
  list(pick = pick, best = best, table = data.frame(kernel=kernel, mean=m, se=se), raw = tab)
}
