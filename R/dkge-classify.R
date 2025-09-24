# dkge-classify.R
# Classification add-on for DKGE fits.

#' Cross-validated classification on DKGE effect patterns
#'
#' @param fit dkge object.
#' @param targets Target specification consumed by [dkge_targets()] or a list of
#'   `dkge_target` objects.
#' @param method Classifier backend: "lda" (default) or "logit" (one-vs-rest
#'   ridge logistic regression).
#' @param folds Cross-fitting specification. `NULL` (default) performs LOSO.
#'   Integer values request subject-level K-fold. A `dkge_folds` object is also
#'   accepted.
#' @param lambda Optional ridge parameter passed to the classifier backend.
#' @param metric Evaluation metric(s); defaults to c("accuracy", "logloss").
#' @param mode Decoding mode: "auto" (default), "cell", or "delta".
#' @param residualize Forwarded to [dkge_targets()].
#' @param collapse Forwarded to [dkge_targets()] for factor collapsing.
#' @param restrict_factors Optional factor subset for `spec = "fullcell"`.
#' @param n_perm Number of permutations for empirical p-values (default 0).
#' @param scope Override permutation scope (otherwise target scope used).
#' @param class_weights Class weighting scheme ("none", "balanced", "inverse").
#' @param ridge Optional ridge added when recomputing held-out bases (default 0).
#' @param control Optional list of advanced controls (power users). Recognised
#'   entries: `lambda_grid` (numeric vector of candidate penalties) and
#'   `lambda_fun` (function returning a lambda per target/fold with signature
#'   `function(target, fold, method, default)`). Defaults to `NULL`, leaving the
#'   standard `lambda` behaviour unchanged.
#' @param blocks Optional vector identifying within-subject blocks (e.g., runs or
#'   sessions) used to constrain permutations. Length must match the number of
#'   subjects in `fit` when supplied.
#' @param parallel Logical; reserved for future parallelism hooks.
#' @param verbose Logical; print progress messages.
#' @param seed Optional random seed applied before permutations.
#' @export
dkge_classify <- function(fit,
                          targets,
                          method = c("lda", "logit"),
                          folds = NULL,
                          lambda = NULL,
                          metric = c("accuracy", "logloss"),
                          mode = c("auto", "cell", "delta"),
                          residualize = TRUE,
                          collapse = NULL,
                          restrict_factors = NULL,
                          n_perm = 0L,
                          scope = NULL,
                          class_weights = c("none", "balanced", "inverse"),
                          ridge = 0,
                          control = NULL,
                          blocks = NULL,
                          parallel = FALSE,
                          verbose = FALSE,
                          seed = NULL) {
  stopifnot(inherits(fit, "dkge"))
  method <- match.arg(method)
  metric <- unique(metric)
  mode <- match.arg(mode)
  class_weights <- match.arg(class_weights)
  if (!is.null(seed)) set.seed(seed)

  control <- control %||% list()
  if (!is.list(control)) {
    stop("`control` must be a list or NULL")
  }
  lambda_grid <- control$lambda_grid %||% NULL
  if (!is.null(lambda_grid)) {
    if (!is.numeric(lambda_grid) || any(lambda_grid <= 0)) {
      stop("`control$lambda_grid` must contain positive numeric values")
    }
    lambda_grid <- sort(unique(lambda_grid))
  }
  lambda_fun <- control$lambda_fun %||% NULL
  if (!is.null(lambda_grid) && !is.null(lambda_fun)) {
    stop("Specify at most one of `control$lambda_grid` or `control$lambda_fun`.")
  }

  if (!is.null(blocks)) {
    if (length(blocks) != length(fit$Btil)) {
      stop("`blocks` must have length equal to the number of subjects in `fit`.")
    }
    blocks <- as.character(blocks)
  }

  target_list <- .dkge_prepare_target_list(fit, targets,
                                           residualize = residualize,
                                           collapse = collapse,
                                           restrict_factors = restrict_factors,
                                           scope_override = scope)
  if (!length(target_list)) {
    stop("No valid targets supplied")
  }

  fold_bundle <- .dkge_prepare_folds(fit, folds)
  fold_info <- .dkge_build_fold_bases(fit,
                                      assignments = fold_bundle$assignments,
                                      ridge = ridge,
                                      align = FALSE,
                                      loader_scope = "all",
                                      verbose = verbose)

  results <- vector("list", length(target_list))
  names(results) <- vapply(target_list, `[[`, character(1), "name")

  for (i in seq_along(target_list)) {
    target <- target_list[[i]]
    target_mode <- if (mode == "auto") .dkge_choose_target_mode(target) else mode
    target_scope <- scope %||% target$scope
    if (target_mode == "cell") {
      if (nrow(target$weight_matrix) < 2) {
        warning(sprintf("Target '%s' has fewer than two classes; skipping.", target$name))
        next
      }
      res <- .dkge_classify_cell_target(fit, target, fold_info,
                                        method = method,
                                        lambda = lambda,
                                        metric = metric,
                                        class_weights = class_weights,
                                        n_perm = n_perm,
                                        scope = target_scope,
                                        control = list(lambda_grid = lambda_grid, lambda_fun = lambda_fun),
                                        blocks = blocks,
                                        parallel = parallel,
                                        verbose = verbose)
    } else if (target_mode == "delta") {
      res <- .dkge_classify_delta_target(fit, target, fold_info,
                                         lambda = lambda,
                                         metric = metric,
                                         n_perm = n_perm,
                                         scope = target_scope,
                                         control = list(lambda_fun = lambda_fun),
                                         verbose = verbose)
    } else {
      warning(sprintf("Target '%s': unsupported mode '%s'; skipping.", target$name, target_mode))
      next
    }
    results[[i]] <- res
  }

  control_out <- if (length(control) == 0) NULL else control

  structure(list(
    fit = fit,
    targets = target_list,
    results = results,
    method = method,
    metric = metric,
    folds = fold_bundle$folds,
    n_perm = n_perm,
    class_weights = class_weights,
    control = control_out,
    blocks = blocks
  ), class = "dkge_classification")
}

.dkge_choose_target_mode <- function(target) {
  if (!is.null(target$class_labels) && length(target$class_labels) >= 2) {
    "cell"
  } else {
    "delta"
  }
}

.dkge_prepare_target_list <- function(fit, targets,
                                      residualize, collapse,
                                      restrict_factors,
                                      scope_override = NULL) {
  if (inherits(targets, "dkge_target")) {
    target_list <- list(targets)
  } else if (is.list(targets) && all(vapply(targets, inherits, logical(1), "dkge_target"))) {
    target_list <- targets
  } else if (is.matrix(targets)) {
    target_list <- list(.dkge_wrap_direct_target(targets))
  } else {
    target_list <- dkge_targets(fit, targets,
                                residualize = residualize,
                                collapse = collapse,
                                restrict_factors = restrict_factors)
  }

  for (i in seq_along(target_list)) {
    if (!inherits(target_list[[i]], "dkge_target")) {
      stop("All targets must inherit from class 'dkge_target'.")
    }
    if (!is.null(scope_override)) {
      target_list[[i]]$scope <- scope_override
    }
    if (is.null(target_list[[i]]$class_labels)) {
      target_list[[i]]$class_labels <- paste0("class", seq_len(nrow(target_list[[i]]$weight_matrix)))
    }
  }
  target_list
}

.dkge_wrap_direct_target <- function(weight_matrix) {
  stopifnot(is.matrix(weight_matrix))
  nm <- paste0("target", seq_len(1))
  target <- list(
    name = nm,
    factors = character(0),
    labels = data.frame(),
    class_labels = if (!is.null(rownames(weight_matrix))) rownames(weight_matrix) else paste0("class", seq_len(nrow(weight_matrix))),
    weight_matrix = weight_matrix,
    indicator = NULL,
    residualized = NA,
    collapse = NULL,
    scope = "within_subject"
  )
  class(target) <- c("dkge_target", "list")
  target
}

.dkge_prepare_folds <- function(fit, folds) {
  S <- length(fit$Btil)
  if (is.null(folds)) {
    assignments <- lapply(seq_len(S), function(s) s)
    fold_obj <- NULL
  } else if (is.numeric(folds) && length(folds) == 1) {
    fold_obj <- dkge_define_folds(fit, type = "subject", k = folds)
    assignments <- fold_obj$assignments
  } else if (inherits(folds, "dkge_folds")) {
    fold_obj <- folds
    assignments <- folds$assignments
  } else {
    stop("folds must be NULL, an integer K, or a dkge_folds object.")
  }
  list(assignments = assignments, folds = fold_obj)
}

.dkge_classify_cell_target <- function(fit, target, fold_info,
                                       method, lambda, metric,
                                       class_weights, n_perm, scope,
                                       control, blocks,
                                       parallel, verbose) {
  S <- length(fit$Btil)
  n_classes <- nrow(target$weight_matrix)
  rank_dim <- ncol(fit$U)
  if (rank_dim == 0) {
    stop("DKGE rank is zero; cannot run classification.")
  }

  # Map subjects to folds
  subject_fold <- integer(S)
  for (fold in fold_info$folds) {
    subject_fold[fold$subjects] <- fold$index
  }

  subject_labels <- fit$subject_ids %||% paste0("subject", seq_len(S))
  class_labels <- target$class_labels
  n_rows <- S * n_classes

  control <- control %||% list()
  lambda_grid <- control$lambda_grid %||% NULL
  lambda_fun <- control$lambda_fun %||% NULL

  row_data <- data.frame(
    row_id = seq_len(n_rows),
    subject_idx = rep(seq_len(S), each = n_classes),
    subject_label = rep(subject_labels, each = n_classes),
    class_idx = rep(seq_len(n_classes), times = S),
    class_label = rep(class_labels, times = S),
    fold = rep(subject_fold, each = n_classes),
    stringsAsFactors = FALSE
  )
  if (!is.null(blocks)) {
    row_data$block <- rep(blocks, each = n_classes)
  } else {
    row_data$block <- rep("1", n_rows)
  }

  feature_cache <- lapply(seq_along(fold_info$folds), function(..) vector("list", S))

  compute_features <- function(fold_idx, subject_idx) {
    cached <- feature_cache[[fold_idx]][[subject_idx]]
    if (!is.null(cached)) return(cached)
    loader <- fold_info$folds[[fold_idx]]$loaders[[as.character(subject_idx)]]
    if (is.null(loader)) {
      stop(sprintf("Loader missing for subject %d in fold %d", subject_idx, fold_idx))
    }
    Z <- target$weight_matrix %*% loader$Y
    feature_cache[[fold_idx]][[subject_idx]] <<- Z
    Z
  }

  run_cv <- function(labels_vec, record = FALSE, lambda_value = NULL) {
    pred_class <- if (record) rep(NA_character_, n_rows) else NULL
    prob_mat <- if (record) matrix(NA_real_, n_rows, n_classes) else NULL
    if (record) colnames(prob_mat) <- class_labels

    obs_truth <- character(0)
    obs_pred <- character(0)
    obs_prob <- NULL

    folds <- fold_info$folds
    for (fold in folds) {
      holdout <- fold$subjects
      train_subjects <- setdiff(seq_len(S), holdout)
      train_rows <- which(row_data$subject_idx %in% train_subjects)
      test_rows <- which(row_data$subject_idx %in% holdout)
      if (!length(train_rows) || !length(test_rows)) next

      X_train <- matrix(0, length(train_rows), rank_dim)
      for (i in seq_along(train_rows)) {
        r <- train_rows[[i]]
        s <- row_data$subject_idx[[r]]
        c_idx <- row_data$class_idx[[r]]
        Z <- compute_features(fold$index, s)
        X_train[i, ] <- Z[c_idx, ]
      }
      y_train <- factor(labels_vec[train_rows], levels = class_labels)
      available_classes <- levels(droplevels(y_train))
      if (length(available_classes) < length(class_labels)) {
        next
      }

      X_test <- matrix(0, length(test_rows), rank_dim)
      for (i in seq_along(test_rows)) {
        r <- test_rows[[i]]
        s <- row_data$subject_idx[[r]]
        c_idx <- row_data$class_idx[[r]]
        Z <- compute_features(fold$index, s)
        X_test[i, ] <- Z[c_idx, ]
      }
      y_test <- factor(labels_vec[test_rows], levels = class_labels)

      lambda_fold <- lambda_value %||% lambda
      lambda_fold <- .dkge_resolve_lambda(lambda_fold, lambda_fun, target$name, fold$index, method, lambda_fold)
      classifier <- switch(method,
        lda = .dkge_fit_lda_classifier(X_train, y_train, lambda = lambda_fold, class_weights = class_weights),
        logit = .dkge_fit_logit_classifier(X_train, y_train, lambda = lambda_fold, class_weights = class_weights)
      )
      pred <- .dkge_predict_classifier(classifier, X_test, class_labels)

      if (record) {
        pred_class[test_rows] <- as.character(pred$class)
        prob_mat[test_rows, ] <- pred$prob
      }
      obs_truth <- c(obs_truth, as.character(y_test))
      obs_pred <- c(obs_pred, as.character(pred$class))
      obs_prob <- rbind(obs_prob, pred$prob)
    }

    if (!length(obs_pred)) {
      metrics <- setNames(rep(NA_real_, length(metric)), metric)
    } else {
      metrics <- setNames(numeric(length(metric)), metric)
      for (j in seq_along(metric)) {
        metrics[[j]] <- .dkge_compute_metric(metric[[j]], obs_truth, obs_pred,
                                             prob = obs_prob,
                                             class_levels = class_labels,
                                             weights = NULL)
      }
    }

    if (record) {
      list(metrics = metrics,
           pred = pred_class,
           prob = prob_mat)
    } else {
      list(metrics = metrics)
    }
  }

  lambda_selected <- lambda
  if (!is.null(lambda_grid)) {
    metric_primary <- metric[[1]]
    metric_direction <- if (metric_primary %in% c("logloss", "brier")) "min" else "max"
    scores <- numeric(length(lambda_grid))
    for (idx_lambda in seq_along(lambda_grid)) {
      res <- run_cv(row_data$class_label, record = FALSE, lambda_value = lambda_grid[[idx_lambda]])
      scores[idx_lambda] <- res$metrics[[metric_primary]]
    }
    if (metric_direction == "max") {
      scores_adj <- replace(scores, is.na(scores), -Inf)
      lambda_selected <- lambda_grid[[which.max(scores_adj)]]
    } else {
      scores_adj <- replace(scores, is.na(scores), Inf)
      lambda_selected <- lambda_grid[[which.min(scores_adj)]]
    }
  }

  observed <- run_cv(row_data$class_label, record = TRUE, lambda_value = lambda_selected)
  perm_matrix <- matrix(NA_real_, n_perm, length(metric))
  if (n_perm > 0) {
    scope_use <- scope %||% "within_subject"
    for (b in seq_len(n_perm)) {
      perm_labels <- .dkge_permute_labels(row_data$class_label,
                                          scope = scope_use,
                                          subjects = row_data$subject_idx,
                                          blocks = if (!is.null(blocks)) row_data$block else NULL)
      perm_res <- run_cv(perm_labels, record = FALSE, lambda_value = lambda_selected)
      perm_matrix[b, ] <- perm_res$metrics
    }
    colnames(perm_matrix) <- metric
  }

  p_values <- if (n_perm > 0) {
    setNames(numeric(length(metric)), metric)
  } else NULL
  if (!is.null(p_values)) {
    for (j in seq_along(metric)) {
      tail_dir <- if (metric[[j]] %in% c("logloss", "brier")) "less" else "greater"
      p_values[[j]] <- .dkge_empirical_pval(observed$metrics[[metric[[j]]]],
                                            perm_matrix[, j],
                                            tail = tail_dir)
    }
  }

  list(
    target = target,
    mode = "cell",
    metrics = observed$metrics,
    permutations = perm_matrix,
    p_values = p_values,
    predictions = observed$pred,
    probabilities = observed$prob,
    row_data = row_data,
    lambda = lambda_selected
  )
}

.dkge_classify_delta_target <- function(fit, target, fold_info,
                                        lambda, metric, n_perm,
                                        scope, control, verbose) {
  S <- length(fit$Btil)
  n_classes <- nrow(target$weight_matrix)
  rank_dim <- ncol(fit$U)
  if (rank_dim == 0) {
    stop("DKGE rank is zero; cannot run classification.")
  }
  if (n_classes != 2) {
    warning(sprintf("Target '%s': delta mode currently supports exactly two classes; skipping.", target$name))
    return(NULL)
  }

  control <- control %||% list()
  lambda_fun <- control$lambda_fun %||% NULL
  lambda_dir <- lambda %||% 1e-3
  lambda_dir <- .dkge_resolve_lambda(lambda_dir, lambda_fun, target$name, fold_idx = NA_integer_, method = "delta", fallback = lambda_dir)
  if (is.null(lambda_dir)) lambda_dir <- 1e-3
  folds <- fold_info$folds

  delta_cache <- lapply(folds, function(fold) {
    out <- vector("list", S)
    for (s in seq_len(S)) {
      loader <- fold$loaders[[as.character(s)]]
      if (is.null(loader)) {
        stop(sprintf("Loader missing for subject %d in fold %d", s, fold$index))
      }
      Z <- target$weight_matrix %*% loader$Y
      out[[s]] <- Z[1, ] - Z[2, ]
    }
    out
  })

  run_delta <- function(signs) {
    scores <- rep(NA_real_, S)
    probs <- rep(NA_real_, S)
    for (fold in folds) {
      idx <- fold$index
      holdout <- fold$subjects
      train_idx <- setdiff(seq_len(S), holdout)
      train_delta <- do.call(rbind, lapply(train_idx, function(s) signs[s] * delta_cache[[idx]][[s]]))
      if (nrow(train_delta) == 0) next
      mu <- colMeans(train_delta)
      if (nrow(train_delta) <= 1) {
        cov_mat <- diag(rank_dim)
      } else {
        cov_mat <- stats::cov(train_delta)
      }
      cov_mat <- cov_mat + lambda_dir * diag(rank_dim)
      direction <- tryCatch({
        qr.solve(cov_mat, mu)
      }, error = function(e) {
        qr.solve(cov_mat + 1e-6 * diag(rank_dim), mu)
      })
      for (s in holdout) {
        delta_vec <- signs[s] * delta_cache[[idx]][[s]]
        score <- sum(delta_vec * direction)
        scores[s] <- score
        probs[s] <- stats::plogis(score)
      }
    }
    list(scores = scores, probs = probs)
  }

  observed <- run_delta(rep(1, S))

  supported_metrics <- c("accuracy", "logloss", "brier")
  unsupported <- setdiff(metric, supported_metrics)
  if (length(unsupported)) {
    warning(sprintf("Metric(s) %s not supported in delta mode; returning NA.",
                    paste(unsupported, collapse = ", ")))
  }

  metrics_out <- setNames(rep(NA_real_, length(metric)), metric)
  valid_idx <- which(!is.na(observed$scores))
  if (length(valid_idx)) {
    if ("accuracy" %in% metric) {
      metrics_out["accuracy"] <- mean(observed$scores[valid_idx] > 0)
    }
    if ("logloss" %in% metric) {
      metrics_out["logloss"] <- mean(-log(pmax(observed$probs[valid_idx], 1e-8)))
    }
    if ("brier" %in% metric) {
      metrics_out["brier"] <- mean((observed$probs[valid_idx] - 1)^2)
    }
  }

  perm_matrix <- matrix(NA_real_, n_perm, length(metric))
  if (n_perm > 0) {
    scope_use <- scope %||% "signflip"
    if (!identical(scope_use, "signflip") && verbose) {
      message("Delta mode permutations use subject-wise sign flips regardless of scope setting.")
    }
    for (b in seq_len(n_perm)) {
      signs <- sample(c(-1, 1), S, replace = TRUE)
      perm_res <- run_delta(signs)
      perm_scores <- setNames(rep(NA_real_, length(metric)), metric)
      if (length(valid_idx)) {
        if ("accuracy" %in% metric) {
          perm_scores["accuracy"] <- mean(perm_res$scores[valid_idx] > 0)
        }
        if ("logloss" %in% metric) {
          perm_scores["logloss"] <- mean(-log(pmax(perm_res$probs[valid_idx], 1e-8)))
        }
        if ("brier" %in% metric) {
          perm_scores["brier"] <- mean((perm_res$probs[valid_idx] - 1)^2)
        }
      }
      perm_matrix[b, ] <- perm_scores
    }
    colnames(perm_matrix) <- metric
  }

  p_values <- if (n_perm > 0) setNames(numeric(length(metric)), metric) else NULL
  if (!is.null(p_values)) {
    for (j in seq_along(metric)) {
      tail_dir <- if (metric[[j]] %in% c("logloss", "brier")) "less" else "greater"
      p_values[[j]] <- .dkge_empirical_pval(metrics_out[[metric[[j]]]],
                                            perm_matrix[, j],
                                            tail = tail_dir)
    }
  }

  class_labels <- target$class_labels
  predictions <- ifelse(observed$scores > 0, class_labels[1], class_labels[2])

  list(
    target = target,
    mode = "delta",
    metrics = metrics_out,
    permutations = perm_matrix,
    p_values = p_values,
    subject_scores = observed$scores,
    subject_probabilities = observed$probs,
    predictions = predictions,
    lambda = lambda_dir
  )
}

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
    y_bin <- as.numeric(y == classes[[j]])
    beta <- rep(0, r)
    if (!any(y_bin) || all(y_bin == 1)) {
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

#' @export
print.dkge_classification <- function(x, ...) {
  cat("DKGE Classification\n")
  cat("--------------------\n")
  cat(sprintf("Targets: %d\n", length(x$results)))
  cat(sprintf("Classifier: %s\n", x$method))
  cat(sprintf("Metrics: %s\n", paste(x$metric, collapse = ", ")))
  cat(sprintf("Permutations: %d\n", x$n_perm))
  for (nm in names(x$results)) {
    res <- x$results[[nm]]
    if (is.null(res)) {
      cat(sprintf("  %s: <skipped>\n", nm))
      next
    }
    metric_str <- paste(sprintf("%s=%.3f", names(res$metrics), res$metrics), collapse = ", ")
    cat(sprintf("  %s: %s\n", nm, metric_str))
    if (!is.null(res$p_values)) {
      pv_str <- paste(sprintf("%s p=%.3f", names(res$p_values), res$p_values), collapse = ", ")
      cat(sprintf("    %s\n", pv_str))
    }
  }
  invisible(x)
}

#' @export
as.data.frame.dkge_classification <- function(x, row.names = NULL, optional = FALSE, ...) {
  rows <- list()
  for (nm in names(x$results)) {
    res <- x$results[[nm]]
    if (is.null(res)) next
    for (metric_nm in names(res$metrics)) {
      rows[[length(rows) + 1L]] <- data.frame(
        target = nm,
        metric = metric_nm,
        value = res$metrics[[metric_nm]],
        p_value = if (!is.null(res$p_values)) res$p_values[[metric_nm]] else NA_real_,
        n_perm = x$n_perm,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) {
    df <- data.frame(target = character(0), metric = character(0), value = numeric(0),
                     p_value = numeric(0), n_perm = integer(0), stringsAsFactors = FALSE)
  } else {
    df <- do.call(rbind, rows)
  }
  if (!is.null(row.names)) rownames(df) <- row.names else rownames(df) <- NULL
  df
}
