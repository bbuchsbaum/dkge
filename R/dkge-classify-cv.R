# dkge-classify-cv.R
# Cross-validation loops for DKGE classification (cell and delta targets).

.dkge_classify_cell_target <- function(fit, target, fold_info,
                                       mode, method, lambda, metric,
                                       class_weights, n_perm, scope,
                                       control, blocks,
                                       parallel, verbose,
                                       standardize_within_fold) {
  if (identical(mode, "cell_cross")) {
    warning(
      "mode='cell_cross' is not yet fully implemented: it currently behaves as ",
      "mode='cell' with standardize_within_fold=TRUE. ",
      "A distinct cross-subject projection protocol is planned for a future release.",
      call. = FALSE
    )
  }
  S <- length(fit$Btil)
  n_classes <- nrow(target$weight_matrix)
  rank_dim <- ncol(fit$U)
  if (rank_dim == 0) {
    stop("DKGE rank is zero; cannot run classification.")
  }

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

  supported_metrics <- c("accuracy", "balanced_accuracy", "logloss", "brier")
  unsupported <- setdiff(metric, supported_metrics)
  if (length(unsupported)) {
    warning(sprintf("Target '%s': metric(s) %s not supported in cell mode; returning NA.",
                    target$name, paste(unsupported, collapse = ", ")),
            call. = FALSE)
  }
  metric_eval <- intersect(metric, supported_metrics)

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

  make_scaler <- function(X) {
    mu <- colMeans(X)
    sigma <- apply(X, 2, stats::sd)
    sigma[!is.finite(sigma) | sigma < 1e-8] <- 1
    list(center = mu, scale = sigma)
  }

  apply_scaler <- function(X, scaler) {
    centered <- sweep(X, 2, scaler$center, "-")
    sweep(centered, 2, scaler$scale, "/")
  }

  run_cv <- function(labels_vec, record = FALSE, lambda_value = NULL) {
    pred_class <- if (record) rep(NA_character_, n_rows) else NULL
    prob_mat <- if (record) matrix(NA_real_, n_rows, n_classes) else NULL
    if (record && length(prob_mat)) colnames(prob_mat) <- class_labels

    obs_truth <- character(0)
    obs_pred <- character(0)
    obs_prob <- NULL

    folds <- fold_info$folds
    fold_diag <- if (record) vector("list", length(folds)) else NULL

    for (fold in folds) {
      idx <- fold$index
      holdout <- fold$subjects
      train_subjects <- setdiff(seq_len(S), holdout)
      train_rows <- which(row_data$subject_idx %in% train_subjects)
      test_rows <- which(row_data$subject_idx %in% holdout)

      diag_entry <- NULL
      if (record) {
        diag_entry <- list(
          fold = idx,
          mode = mode,
          standardized = standardize_within_fold,
          holdout_subjects = subject_labels[holdout],
          lambda = NA_real_,
          class_counts_train = setNames(rep(0L, n_classes), class_labels),
          class_counts_test = setNames(rep(0L, n_classes), class_labels),
          confusion = NULL,
          skipped = FALSE,
          reason = NULL
        )
      }

      if (!length(train_rows) || !length(test_rows)) {
        if (record) {
          diag_entry$skipped <- TRUE
          diag_entry$reason <- "empty_split"
          fold_diag[[idx]] <- diag_entry
        }
        next
      }

      X_train <- matrix(0, length(train_rows), rank_dim)
      for (i in seq_along(train_rows)) {
        r <- train_rows[[i]]
        s <- row_data$subject_idx[[r]]
        c_idx <- row_data$class_idx[[r]]
        Z <- compute_features(fold$index, s)
        X_train[i, ] <- Z[c_idx, ]
      }
      y_train <- factor(labels_vec[train_rows], levels = class_labels)
      train_counts <- table(factor(y_train, levels = class_labels))
      if (record) {
        diag_entry$class_counts_train <- as.integer(train_counts)
        names(diag_entry$class_counts_train) <- class_labels
      }

      available_classes <- levels(droplevels(y_train))
      if (length(available_classes) < length(class_labels)) {
        if (record) {
          diag_entry$skipped <- TRUE
          diag_entry$reason <- "missing_training_classes"
          test_counts <- table(factor(labels_vec[test_rows], levels = class_labels))
          diag_entry$class_counts_test <- as.integer(test_counts)
          names(diag_entry$class_counts_test) <- class_labels
          fold_diag[[idx]] <- diag_entry
        }
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
      if (record) {
        test_counts <- table(factor(y_test, levels = class_labels))
        diag_entry$class_counts_test <- as.integer(test_counts)
        names(diag_entry$class_counts_test) <- class_labels
      }

      if (standardize_within_fold) {
        scaler <- make_scaler(X_train)
        X_train <- apply_scaler(X_train, scaler)
        X_test <- apply_scaler(X_test, scaler)
      }

      lambda_fold <- lambda_value %||% lambda
      lambda_fold <- .dkge_resolve_lambda(lambda_fold, lambda_fun, target$name,
                                          fold$index, method, lambda_fold)
      if (record) {
        diag_entry$lambda <- lambda_fold
      }
      classifier <- switch(method,
        lda = .dkge_fit_lda_classifier(X_train, y_train, lambda = lambda_fold, class_weights = class_weights),
        logit = .dkge_fit_logit_classifier(X_train, y_train, lambda = lambda_fold, class_weights = class_weights)
      )
      pred <- .dkge_predict_classifier(classifier, X_test, class_labels)

      if (record) {
        pred_class[test_rows] <- as.character(pred$class)
        prob_mat[test_rows, ] <- pred$prob
        confusion_tbl <- table(factor(y_test, levels = class_labels),
                               factor(pred$class, levels = class_labels))
        diag_entry$confusion <- matrix(as.integer(confusion_tbl),
                                       nrow = nrow(confusion_tbl),
                                       dimnames = dimnames(confusion_tbl))
        fold_diag[[idx]] <- diag_entry
      }

      obs_truth <- c(obs_truth, as.character(y_test))
      obs_pred <- c(obs_pred, as.character(pred$class))
      obs_prob <- if (is.null(obs_prob)) pred$prob else rbind(obs_prob, pred$prob)
    }

    if (!length(obs_pred) || !length(metric_eval)) {
      metrics <- setNames(rep(NA_real_, length(metric_eval)), metric_eval)
    } else {
      metrics <- setNames(numeric(length(metric_eval)), metric_eval)
      for (j in seq_along(metric_eval)) {
        metrics[[j]] <- .dkge_compute_metric(metric_eval[[j]], obs_truth, obs_pred,
                                             prob = obs_prob,
                                             class_levels = class_labels,
                                             weights = NULL)
      }
    }

    if (record) {
      list(metrics = metrics,
           pred = pred_class,
           prob = prob_mat,
           fold_diag = fold_diag)
    } else {
      list(metrics = metrics)
    }
  }

  lambda_selected <- lambda
  if (!is.null(lambda_grid) && length(metric_eval)) {
    metric_primary <- metric[[1]]
    if (!metric_primary %in% metric_eval) {
      metric_primary <- metric_eval[[1]]
    }
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
  metrics_out <- setNames(rep(NA_real_, length(metric)), metric)
  if (length(metric_eval)) {
    metrics_out[metric_eval] <- observed$metrics
  }

  perm_matrix <- matrix(NA_real_, n_perm, length(metric), dimnames = list(NULL, metric))
  if (n_perm > 0) {
    scope_use <- scope %||% "within_subject"
    for (b in seq_len(n_perm)) {
      perm_labels <- .dkge_permute_labels(row_data$class_label,
                                          scope = scope_use,
                                          subjects = row_data$subject_idx,
                                          blocks = if (!is.null(blocks)) row_data$block else NULL)
      perm_res <- run_cv(perm_labels, record = FALSE, lambda_value = lambda_selected)
      if (length(metric_eval)) {
        perm_matrix[b, metric_eval] <- perm_res$metrics
      }
    }
  }

  p_values <- if (n_perm > 0) setNames(rep(NA_real_, length(metric)), metric) else NULL
  if (!is.null(p_values) && length(metric_eval)) {
    for (nm in metric_eval) {
      obs_val <- metrics_out[[nm]]
      if (is.na(obs_val)) {
        p_values[[nm]] <- NA_real_
      } else {
        tail_dir <- if (nm %in% c("logloss", "brier")) "less" else "greater"
        col_idx <- match(nm, metric)
        p_values[[nm]] <- .dkge_empirical_pval(obs_val, perm_matrix[, col_idx], tail = tail_dir)
      }
    }
  }

  diagnostics <- list(
    folds = observed$fold_diag,
    standardize_within_fold = standardize_within_fold,
    metric_eval = metric_eval
  )

  list(
    target = target,
    mode = mode,
    metrics = metrics_out,
    permutations = perm_matrix,
    p_values = p_values,
    predictions = observed$pred,
    probabilities = observed$prob,
    row_data = row_data,
    lambda = lambda_selected,
    diagnostics = diagnostics
  )
}

.dkge_classify_delta_target <- function(fit, target, fold_info,
                                        lambda, metric, n_perm,
                                        scope, control, verbose,
                                        y = NULL, subject_ids = NULL) {
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

  class_labels <- target$class_labels
  positive_class <- class_labels[[1]]

  label_vec <- y %||% target$y
  y_factor <- NULL
  if (!is.null(label_vec)) {
    if (is.data.frame(label_vec)) {
      if (ncol(label_vec) != 1) {
        stop(sprintf("Target '%s': subject labels must be a single column.", target$name))
      }
      label_vec <- label_vec[[1]]
    }
    if (is.matrix(label_vec)) {
      if (ncol(label_vec) != 1) {
        stop(sprintf("Target '%s': subject labels must be a single column.", target$name))
      }
      label_vec <- label_vec[, 1]
    }
    if (is.list(label_vec) && !is.factor(label_vec)) {
      stop(sprintf("Target '%s': subject labels must be an atomic vector or factor.", target$name))
    }
    if (length(label_vec) != S) {
      stop(sprintf("Target '%s': subject labels must have length %d (got %d).",
                   target$name, S, length(label_vec)))
    }
    y_factor <- factor(label_vec, levels = class_labels)
    if (anyNA(y_factor)) {
      bad_idx <- which(is.na(y_factor))
      offending <- unique(label_vec[bad_idx])
      stop(sprintf("Target '%s': subject labels contain unknown levels: %s",
                   target$name, paste(offending, collapse = ", ")))
    }
  }
  if (!is.null(subject_ids) && !is.null(y_factor)) {
    names(y_factor) <- subject_ids
  }

  control <- control %||% list()
  lambda_fun <- control$lambda_fun %||% NULL
  lambda_dir <- lambda %||% 1e-3
  lambda_dir <- .dkge_resolve_lambda(lambda_dir, lambda_fun, target$name,
                                     fold_idx = NA_integer_, method = "delta",
                                     fallback = lambda_dir)
  if (is.null(lambda_dir)) lambda_dir <- 1e-3
  folds <- fold_info$folds

  supported_metrics <- c("accuracy", "logloss", "brier", "auroc", "ece")
  unsupported <- setdiff(metric, supported_metrics)
  if (length(unsupported)) {
    warning(sprintf("Target '%s': metric(s) %s not supported in delta mode; returning NA.",
                    target$name, paste(unsupported, collapse = ", ")),
            call. = FALSE)
  }
  metric_eval <- intersect(metric, supported_metrics)

  compute_delta_metrics <- function(scores, probs, labels, metrics) {
    out <- setNames(rep(NA_real_, length(metrics)), metrics)
    if (!length(metrics) || is.null(labels)) {
      return(out)
    }
    valid_idx <- which(!is.na(scores) & !is.na(probs) & !is.na(labels))
    if (!length(valid_idx)) {
      return(out)
    }
    scores <- scores[valid_idx]
    probs <- probs[valid_idx]
    labels <- labels[valid_idx]
    truth_pos <- as.numeric(labels == positive_class)

    if ("accuracy" %in% metrics) {
      pred_pos <- as.numeric(scores >= 0)
      out["accuracy"] <- mean(pred_pos == truth_pos)
    }
    if ("logloss" %in% metrics) {
      p <- ifelse(truth_pos == 1, probs, 1 - probs)
      p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
      out["logloss"] <- mean(-log(p))
    }
    if ("brier" %in% metrics) {
      out["brier"] <- mean((probs - truth_pos)^2)
    }
    if ("auroc" %in% metrics) {
      out["auroc"] <- .dkge_metric_auc_binary(probs, labels, positive_class)
    }
    if ("ece" %in% metrics) {
      out["ece"] <- .dkge_metric_ece(probs, labels, positive_class)
    }
    out
  }

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
      cov_mat <- if (nrow(train_delta) <= 1) diag(rank_dim) else stats::cov(train_delta)
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
  if (!is.null(subject_ids)) {
    names(observed$scores) <- subject_ids
    names(observed$probs) <- subject_ids
  }

  metrics_out <- setNames(rep(NA_real_, length(metric)), metric)
  if (!is.null(y_factor)) {
    metrics_eval <- compute_delta_metrics(observed$scores, observed$probs, y_factor, metric_eval)
    metrics_out[names(metrics_eval)] <- metrics_eval
  } else if (length(metric_eval)) {
    warning(sprintf("Target '%s': delta mode requires subject labels; metrics set to NA.",
                    target$name), call. = FALSE)
  }

  perm_matrix <- matrix(NA_real_, n_perm, length(metric), dimnames = list(NULL, metric))
  if (n_perm > 0) {
    if (!is.null(scope) && !identical(scope, "signflip")) {
      stop("Delta mode permutations support only scope = 'signflip'.")
    }
    perm_columns <- match(metric_eval, metric)
    for (b in seq_len(n_perm)) {
      signs <- sample(c(-1, 1), S, replace = TRUE)
      perm_res <- run_delta(signs)
      if (length(metric_eval) && !is.null(y_factor)) {
        perm_metrics <- compute_delta_metrics(perm_res$scores, perm_res$probs, y_factor, metric_eval)
        perm_matrix[b, perm_columns] <- perm_metrics
      }
    }
  }

  p_values <- if (n_perm > 0) setNames(rep(NA_real_, length(metric)), metric) else NULL
  if (!is.null(p_values) && length(metric_eval) && !is.null(y_factor)) {
    tails <- setNames(rep("greater", length(metric_eval)), metric_eval)
    tails[c("logloss", "brier", "ece")] <- "less"
    for (nm in metric_eval) {
      obs_val <- metrics_out[[nm]]
      if (is.na(obs_val)) {
        p_values[[nm]] <- NA_real_
      } else {
        col_idx <- match(nm, metric)
        p_values[[nm]] <- .dkge_empirical_pval(obs_val, perm_matrix[, col_idx], tail = tails[[nm]])
      }
    }
  }

  predictions <- rep(class_labels[2], S)
  pos_idx <- !is.na(observed$scores) & observed$scores >= 0
  predictions[pos_idx] <- positive_class
  predictions <- factor(predictions, levels = class_labels)
  if (!is.null(subject_ids)) {
    names(predictions) <- subject_ids
  }

  total_signflip <- if (S <= 30) 2^S else Inf
  rng_state <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else NULL

  list(
    target = target,
    mode = "delta",
    metrics = metrics_out,
    permutations = perm_matrix,
    p_values = p_values,
    subject_scores = observed$scores,
    subject_probabilities = observed$probs,
    predictions = predictions,
    positive_class = positive_class,
    subject_labels = y_factor,
    lambda = lambda_dir,
    permutation_scope = "signflip",
    n_perm_requested = n_perm,
    n_perm_performed = if (n_perm > 0 && !is.null(y_factor)) n_perm else 0L,
    total_signflip_configurations = total_signflip,
    rng_state = rng_state
  )
}
