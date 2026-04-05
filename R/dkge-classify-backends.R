# dkge-classify-backends.R
# Classifier backends (LDA, logistic regression) and prediction helpers.

.dkge_resolve_lambda <- function(lambda_candidate, lambda_fun, target_name, fold_idx, method, fallback = NULL) {
  base <- lambda_candidate %||% fallback
  if (!is.null(lambda_fun)) {
    val <- lambda_fun(target = target_name,
                      fold = fold_idx,
                      method = method,
                      default = base)
    if (is.null(val)) {
      val <- base
    }
  } else {
    val <- base
  }
  if (is.null(val)) {
    return(NULL)
  }
  if (!is.numeric(val) || length(val) != 1 || val <= 0) {
    stop("lambda returned by control must be a positive numeric scalar")
  }
  as.numeric(val)
}

.dkge_fit_lda_classifier <- function(X, y, lambda = NULL, class_weights = c("none", "balanced", "inverse")) {
  class_weights <- match.arg(class_weights)
  y <- droplevels(y)
  classes <- levels(y)
  r <- ncol(X)
  k <- length(classes)
  lambda <- if (is.null(lambda)) 1e-4 else as.numeric(lambda)

  counts <- table(y)
  means <- matrix(0, k, r)
  rownames(means) <- classes
  cov_mat <- matrix(0, r, r)

  for (i in seq_along(classes)) {
    cls <- classes[[i]]
    idx <- which(y == cls)
    Xi <- X[idx, , drop = FALSE]
    mu <- colMeans(Xi)
    means[i, ] <- mu
    centered <- sweep(Xi, 2, mu, FUN = "-")
    cov_mat <- cov_mat + t(centered) %*% centered
  }
  denom <- max(nrow(X) - k, 1)
  cov_mat <- cov_mat / denom
  cov_mat <- cov_mat + lambda * diag(r)
  inv_cov <- tryCatch(solve(cov_mat), error = function(...) {
    solve(cov_mat + 1e-6 * diag(r))
  })

  priors <- switch(class_weights,
    none = counts / sum(counts),
    balanced = rep(1 / k, k),
    inverse = {
      w <- 1 / as.numeric(counts)
      w / sum(w)
    }
  )
  names(priors) <- classes

  list(type = "lda",
       classes = classes,
       means = means,
       inv_cov = inv_cov,
       priors = priors)
}

.dkge_sample_weights <- function(y, scheme = c("none", "balanced", "inverse")) {
  scheme <- match.arg(scheme)
  n <- length(y)
  if (n == 0) return(numeric(0))
  counts <- table(y)
  weights <- switch(scheme,
    none = rep(1, n),
    balanced = {
      w <- length(counts) / as.numeric(counts[as.character(y)])
      w
    },
    inverse = {
      w <- 1 / as.numeric(counts[as.character(y)])
      w
    }
  )
  weights / mean(weights)
}

.dkge_fit_logit_classifier <- function(X, y, lambda = NULL,
                                       class_weights = c("none", "balanced", "inverse"),
                                       max_iter = 50,
                                       tol = 1e-6) {
  class_weights <- match.arg(class_weights)
  y <- droplevels(y)
  classes <- levels(y)
  k <- length(classes)
  lambda <- if (is.null(lambda)) 1e-2 else as.numeric(lambda)
  weights <- .dkge_sample_weights(y, class_weights)
  X_aug <- cbind(`(Intercept)` = 1, X)
  r <- ncol(X_aug)
  betas <- matrix(0, r, k)
  penalty <- diag(c(0, rep(1, r - 1)), nrow = r)

  for (j in seq_len(k)) {
    y_indicator <- y == classes[[j]]
    y_bin <- as.numeric(y_indicator)
    beta <- rep(0, r)
    if (!any(y_indicator) || all(y_indicator)) {
      next
    }
    for (iter in seq_len(max_iter)) {
      eta <- drop(X_aug %*% beta)
      p <- stats::plogis(eta)
      w <- weights * pmax(p * (1 - p), 1e-6)
      z <- eta + (y_bin - p) / pmax(p * (1 - p), 1e-6)
      sqrt_w <- sqrt(w)
      Xw <- X_aug * sqrt_w
      zw <- sqrt_w * z
      G <- crossprod(Xw) + lambda * penalty
      b <- crossprod(Xw, zw)
      beta_new <- tryCatch({
        qr.solve(G, b)
      }, error = function(e) {
        qr.solve(G + 1e-6 * diag(r), b)
      })
      if (max(abs(beta_new - beta)) < tol) {
        beta <- beta_new
        break
      }
      beta <- beta_new
    }
    betas[, j] <- beta
  }

  list(type = "logit",
       classes = classes,
       beta = betas,
       lambda = lambda)
}

.dkge_predict_classifier <- function(model, X, class_levels) {
  if (model$type == "lda") {
    probs <- .dkge_predict_lda(model, X, class_levels)
  } else if (model$type == "logit") {
    probs <- .dkge_predict_logit(model, X, class_levels)
  } else {
    stop("Unsupported classifier type")
  }
  pred <- class_levels[max.col(probs, ties.method = "first")]
  list(class = pred, prob = probs)
}

.dkge_predict_lda <- function(model, X, class_levels) {
  classes <- model$classes
  inv_cov <- model$inv_cov
  means <- model$means
  priors <- model$priors
  scores <- matrix(0, nrow(X), length(classes))
  for (i in seq_along(classes)) {
    mu <- means[i, ]
    linear <- drop(X %*% (inv_cov %*% mu))
    const <- as.numeric(-0.5 * drop(mu %*% inv_cov %*% mu) + log(priors[[i]]))
    scores[, i] <- linear + const
  }
  max_score <- apply(scores, 1, max)
  exp_scores <- exp(scores - max_score)
  probs_partial <- exp_scores / rowSums(exp_scores)
  colnames(probs_partial) <- classes
  # expand to full class_levels
  probs <- matrix(0, nrow(X), length(class_levels))
  colnames(probs) <- class_levels
  match_idx <- match(classes, class_levels)
  probs[, match_idx] <- probs_partial
  probs[!is.finite(probs)] <- 0
  row_totals <- rowSums(probs)
  zero_rows <- which(row_totals <= .Machine$double.eps)
  if (length(zero_rows)) {
    probs[zero_rows, ] <- 1 / length(class_levels)
  }
  row_totals <- rowSums(probs)
  probs <- sweep(probs, 1, pmax(row_totals, .Machine$double.eps), "/")
  probs
}

.dkge_predict_logit <- function(model, X, class_levels) {
  X_aug <- cbind(`(Intercept)` = 1, X)
  scores <- X_aug %*% model$beta
  classes <- model$classes
  if (!length(classes)) {
    probs <- matrix(1 / length(class_levels), nrow(X), length(class_levels))
    colnames(probs) <- class_levels
    return(probs)
  }
  colnames(scores) <- classes
  max_score <- apply(scores, 1, max)
  scores_centered <- scores - max_score
  exp_scores <- exp(scores_centered)
  row_sums <- rowSums(exp_scores)
  zero_idx <- which(row_sums == 0)
  if (length(zero_idx)) {
    exp_scores[zero_idx, ] <- 1
    row_sums <- rowSums(exp_scores)
  }
  probs_partial <- exp_scores / row_sums
  colnames(probs_partial) <- classes
  probs <- matrix(0, nrow(X), length(class_levels))
  colnames(probs) <- class_levels
  match_idx <- match(classes, class_levels)
  probs[, match_idx] <- probs_partial
  probs
}
