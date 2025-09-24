# dkge-metrics.R
# Metric utilities for DKGE classification.

.dkge_metric_accuracy <- function(truth, pred, weights = NULL) {
  stopifnot(length(truth) == length(pred))
  w <- weights %||% rep(1, length(truth))
  w <- as.numeric(w)
  correct <- as.numeric(truth == pred)
  sum(w * correct) / sum(w)
}

.dkge_metric_balanced_accuracy <- function(truth, pred, class_levels = NULL) {
  if (is.null(class_levels)) class_levels <- sort(unique(truth))
  recalls <- numeric(length(class_levels))
  for (i in seq_along(class_levels)) {
    cls <- class_levels[[i]]
    idx <- which(truth == cls)
    if (!length(idx)) {
      recalls[[i]] <- NA_real_
    } else {
      recalls[[i]] <- mean(pred[idx] == cls)
    }
  }
  mean(recalls, na.rm = TRUE)
}

.dkge_metric_logloss <- function(truth, prob, class_levels = NULL, eps = 1e-12, weights = NULL) {
  prob <- as.matrix(prob)
  n <- length(truth)
  if (nrow(prob) != n) stop("Probability matrix has incompatible number of rows.")
  if (is.null(class_levels)) class_levels <- colnames(prob)
  if (is.null(class_levels)) class_levels <- sort(unique(truth))

  w <- weights %||% rep(1, n)
  w <- as.numeric(w)
  total <- sum(w)
  if (total <= 0) return(NA_real_)

  idx <- match(truth, class_levels)
  if (any(is.na(idx))) {
    stop("Some truth labels not found in probability columns.")
  }
  p <- prob[cbind(seq_len(n), idx)]
  p <- pmin(pmax(p, eps), 1 - eps)
  - sum(w * log(p)) / total
}

.dkge_metric_brier <- function(truth, prob, class_levels = NULL, weights = NULL) {
  prob <- as.matrix(prob)
  n <- length(truth)
  if (nrow(prob) != n) stop("Probability matrix has incompatible number of rows.")
  if (is.null(class_levels)) class_levels <- colnames(prob)
  if (is.null(class_levels)) class_levels <- sort(unique(truth))

  w <- weights %||% rep(1, n)
  w <- as.numeric(w)

  target <- matrix(0, n, length(class_levels))
  colnames(target) <- class_levels
  target[cbind(seq_len(n), match(truth, class_levels))] <- 1

  sqerr <- rowSums((prob - target)^2)
  sum(w * sqerr) / sum(w)
}

.dkge_compute_metric <- function(metric, truth, pred, prob, class_levels, weights = NULL) {
  metric <- match.arg(metric, c("accuracy", "balanced_accuracy", "logloss", "brier"))
  switch(metric,
    accuracy = .dkge_metric_accuracy(truth, pred, weights = weights),
    balanced_accuracy = .dkge_metric_balanced_accuracy(truth, pred, class_levels = class_levels),
    logloss = .dkge_metric_logloss(truth, prob, class_levels = class_levels, weights = weights),
    brier = .dkge_metric_brier(truth, prob, class_levels = class_levels, weights = weights)
  )
}

.dkge_empirical_pval <- function(observed, null_dist, tail = c("greater", "less")) {
  tail <- match.arg(tail)
  if (!length(null_dist)) return(NA_real_)
  if (tail == "greater") {
    (1 + sum(null_dist >= observed)) / (length(null_dist) + 1)
  } else {
    (1 + sum(null_dist <= observed)) / (length(null_dist) + 1)
  }
}
