# dkge-regress.R
# Multi-output regression wrapper operating on DKGE effect embeddings.

#' Multi-output regression on DKGE effects
#'
#' Fits cross-validated regressors that map DKGE effect embeddings to
#' continuous targets supplied per effect (e.g., semantic vectors, ratings,
#' reference maps). The default engine uses ridge-regularised multi-response
#' regression via \pkg{glmnet}.
#'
#' @param fit A `dkge` object containing an effect-space basis (`fit$U`).
#' @param Y Numeric matrix of shape `(#effects) x (#outputs)` with row names that
#'   identify the effects/items to be modelled.
#' @param effect_ids Optional character vector of effect identifiers to include
#'   (defaults to `rownames(Y)`).
#' @param folds Optional list giving held-out effect IDs for each fold. When
#'   `NULL`, random K-fold splits are generated.
#' @param k Number of folds to create when `folds` is `NULL` (default `5`).
#' @param features Either the string "effects" (use the fitted DKGE effect
#'   embeddings) or a custom function `function(fit, ids, ctx)` returning a
#'   matrix of predictors for the supplied effect IDs. The function may use the
#'   optional `ctx` argument to receive fold-specific context objects.
#' @param engine Regression engine. One of "glmnet" (default ridge multi-output
#'   regression), "lm" (ordinary least squares solved jointly across outputs),
#'   or a list with `fit`/`predict` closures: `list(fit = function(X, Y) {...},
#'   predict = function(model, X) {...})`.
#' @param alpha Elastic-net mixing parameter passed to \pkg{glmnet} (0 for ridge,
#'   1 for lasso).
#' @param standardize Logical; when `TRUE` (default) the \pkg{glmnet} engine
#'   standardises predictors internally.
#' @param seed Random seed used when generating folds.
#' @param return_models Logical; when `TRUE`, attach the fitted model for each
#'   fold to the returned object.
#'
#' @return An object of class `"dkge_regress"` containing predictions (`pred`),
#'   observed targets (`truth`), fold assignments (`folds`), summary metrics
#'   (`metrics`), and optionally the per-fold models (`models`).
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' Q <- 4L
#' effects <- paste0("e", seq_len(Q))
#' betas <- replicate(3, matrix(rnorm(Q * 6), Q, 6, dimnames = list(effects, NULL)), simplify = FALSE)
#' designs <- replicate(3, diag(Q), simplify = FALSE)
#' designs <- lapply(designs, function(X) { colnames(X) <- effects; X })
#' fit <- dkge_fit(betas, designs = designs, K = diag(Q))
#' W <- matrix(rnorm(ncol(fit$U) * 2L), ncol = 2L)
#' rownames(W) <- colnames(fit$U)
#' Y <- fit$U %*% W
#' rownames(Y) <- effects
#' res <- dkge_regress(fit, Y, k = 2, engine = "lm")
#' print(res)
#' }
dkge_regress <- function(fit,
                         Y,
                         effect_ids = rownames(Y),
                         folds = NULL,
                         k = 5,
                         features = "effects",
                         engine = "glmnet",
                         alpha = 0,
                         standardize = TRUE,
                         seed = 1,
                         return_models = FALSE) {
  stopifnot(inherits(fit, "dkge"))
  if (!(is.matrix(Y) || is.data.frame(Y))) {
    stop("`Y` must be a matrix or data.frame of continuous targets.")
  }
  Y <- as.matrix(Y)
  if (is.null(rownames(Y))) {
    stop("`Y` must have row names corresponding to effect IDs.")
  }
  if (is.null(effect_ids)) {
    effect_ids <- rownames(Y)
  }
  missing_in_Y <- setdiff(effect_ids, rownames(Y))
  if (length(missing_in_Y)) {
    stop("These effect IDs were not found in `Y`: ",
         paste(missing_in_Y, collapse = ", "))
  }
  Y <- Y[effect_ids, , drop = FALSE]

  # Default feature extractor -------------------------------------------------
  features_fun_default <- function(fit, ids, ctx = NULL) {
    U <- fit$U
    if (is.null(U) || ncol(U) == 0) {
      stop("`fit$U` is unavailable or rank=0; supply a custom `features` function.")
    }
    row_ids <- rownames(U)
    if (is.null(row_ids)) {
      row_ids <- fit$effects
    }
    if (is.null(row_ids)) {
      stop("Unable to align effects because neither `rownames(fit$U)` nor `fit$effects` are set. Provide `features` manually.")
    }
    idx <- match(ids, row_ids)
    if (anyNA(idx)) {
      missing <- ids[is.na(idx)]
      stop("Some requested effects are missing from the DKGE basis: ",
           paste(missing, collapse = ", "))
    }
    out <- U[idx, , drop = FALSE]
    rownames(out) <- ids
    out
  }

  features_fun <- if (is.character(features)) {
    if (!identical(features, "effects")) {
      stop("Unsupported `features` string. Use 'effects' or provide a custom function.")
    }
    features_fun_default
  } else if (is.function(features)) {
    features
  } else {
    stop("`features` must be 'effects' or a function(fit, ids, ctx).")
  }

  # Fold allocation -----------------------------------------------------------
  if (is.null(folds)) {
    set.seed(seed)
    n <- length(effect_ids)
    if (k < 2 || k > n) {
      stop("`k` must be between 2 and the number of effects.")
    }
    fold_ids <- sample(rep(seq_len(k), length.out = n))
    folds <- lapply(seq_len(k), function(i) effect_ids[fold_ids == i])
  }
  if (!is.list(folds) || !length(folds)) {
    stop("`folds` must be a non-empty list of effect ID vectors.")
  }

  # Optionally provide fold-aware context objects -----------------------------
  get_ctx_fun <- NULL
  if (exists(".dkge_fold_weight_context", mode = "function")) {
    get_ctx_fun <- get(".dkge_fold_weight_context", mode = "function")
  }

  P <- matrix(NA_real_, nrow = nrow(Y), ncol = ncol(Y),
              dimnames = list(effect_ids, colnames(Y)))
  models <- vector("list", length(folds))

  for (f in seq_along(folds)) {
    test_ids <- folds[[f]]
    train_ids <- setdiff(effect_ids, test_ids)
    if (!length(train_ids)) {
      stop("Each fold must leave at least one effect for training.")
    }

    ctx <- if (!is.null(get_ctx_fun)) {
      tryCatch(get_ctx_fun(fit, train_ids), error = function(e) NULL)
    } else {
      NULL
    }

    X_train <- features_fun(fit, train_ids, ctx)
    X_test  <- features_fun(fit, test_ids, ctx)
    Y_train <- Y[train_ids, , drop = FALSE]

    if (!is.matrix(X_train)) X_train <- as.matrix(X_train)
    if (!is.matrix(X_test)) X_test <- as.matrix(X_test)
    rownames(X_train) <- train_ids
    rownames(X_test) <- test_ids
    if (nrow(X_train) != nrow(Y_train)) {
      stop("Row mismatch between training predictors and targets on fold ", f)
    }
    if (!is.null(rownames(X_train)) && !all(rownames(X_train) == train_ids)) {
      X_train <- X_train[match(train_ids, rownames(X_train)), , drop = FALSE]
    }
    if (!is.null(rownames(X_test)) && !all(rownames(X_test) == test_ids)) {
      X_test <- X_test[match(test_ids, rownames(X_test)), , drop = FALSE]
    }

    if (identical(engine, "glmnet")) {
      if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Engine 'glmnet' requires the glmnet package. Install it or choose another engine.")
      }
      fit_obj <- glmnet::cv.glmnet(x = X_train,
                                   y = Y_train,
                                   family = "mgaussian",
                                   alpha = alpha,
                                   standardize = standardize)
      pred <- stats::predict(fit_obj, newx = X_test, s = "lambda.min")
      pred <- pred[, , 1, drop = FALSE]
      P[test_ids, ] <- pred
      if (return_models) models[[f]] <- fit_obj
    } else if (identical(engine, "lm")) {
      coef <- qr.solve(X_train, Y_train)
      P[test_ids, ] <- X_test %*% coef
      if (return_models) {
        models[[f]] <- list(coef = coef, engine = "lm")
      }
    } else if (is.list(engine) && is.function(engine$fit) && is.function(engine$predict)) {
      mod <- engine$fit(X_train, Y_train)
      P[test_ids, ] <- engine$predict(mod, X_test)
      if (return_models) models[[f]] <- mod
    } else {
      stop("Unsupported engine. Use 'glmnet', 'lm', or provide fit/predict closures.")
    }
  }

  # Metrics -------------------------------------------------------------------
  diff <- P - Y
  rmse_micro <- sqrt(mean(diff^2, na.rm = TRUE))
  sse <- sum(diff^2, na.rm = TRUE)
  grand_mean <- matrix(colMeans(Y, na.rm = TRUE), nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE)
  tss <- sum((Y - grand_mean)^2, na.rm = TRUE)
  r2_micro <- if (tss > 0) 1 - sse / tss else NA_real_

  col_sse <- colSums(diff^2, na.rm = TRUE)
  col_means <- colMeans(Y, na.rm = TRUE)
  col_tss <- colSums((Y - matrix(col_means, nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE))^2, na.rm = TRUE)
  r2_macro <- mean(ifelse(col_tss > 0, 1 - col_sse / col_tss, NA_real_), na.rm = TRUE)
  rmse_by_dim <- sqrt(colMeans(diff^2, na.rm = TRUE))

  row_norm <- function(M) sqrt(rowSums(M^2))
  numer <- rowSums(P * Y)
  denom <- row_norm(P) * row_norm(Y)
  cos_scores <- numer / denom
  cos_scores[!is.finite(cos_scores)] <- NA_real_
  cosine_mean <- mean(cos_scores, na.rm = TRUE)

  out <- list(
    pred = P,
    truth = Y,
    folds = folds,
    metrics = list(
      rmse_micro = rmse_micro,
      r2_micro = r2_micro,
      r2_macro = r2_macro,
      rmse_by_dim = rmse_by_dim,
      cosine_mean = cosine_mean
    ),
    models = if (return_models) models else NULL
  )
  class(out) <- "dkge_regress"
  out
}

#' @export
print.dkge_regress <- function(x, ...) {
  cat("dkge_regress results\n")
  cat(sprintf("  RMSE (micro): %.4f\n", x$metrics$rmse_micro))
  cat(sprintf("  R^2  (micro): %.4f\n", x$metrics$r2_micro))
  cat(sprintf("  R^2  (macro): %.4f\n", x$metrics$r2_macro))
  cat(sprintf("  Cosine mean : %.4f\n", x$metrics$cosine_mean))
  invisible(x)
}
