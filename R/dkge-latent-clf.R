#' Cross-fitted linear classifiers in the DKGE latent space
#'
#' Trains a binary linear classifier on DKGE latent features using subject-level
#' folds (LOSO or K-fold). The result stores one weight vector \eqn{\beta^{(-s)}}
#' per held-out subject so downstream information-map routines can construct
#' bias-aware decoder or Haufe maps.
#'
#' @param fit Fitted `dkge` object.
#' @param y Either a factor of length `S` (subject-level labels) or a list of
#'   length `S` supplying per-sample labels for each subject.
#' @param Z_by_subject Optional list of latent feature matrices
#'   (`P_s x r`). When `NULL`, [dkge_project_clusters_to_latent()] is used.
#' @param folds Either an integer `K` or a `dkge_folds` object created with
#'   [dkge_define_folds()].
#' @param model Binary classifier: `"lda"` (pooled covariance linear
#'   discriminant, default), `"ridge_logit"` (ridge-penalised logistic
#'   regression via glmnet when available), or `"lsvm"` (placeholder falling
#'   back to LDA).
#' @param ridge Ridge penalty added to the pooled covariance matrix when using
#'   LDA. Ensures numerical stability when the latent dimension is high.
#' @param level Training granularity: `"subject"` averages each subject's
#'   clusters, while `"sample"` stacks all cluster samples from training
#'   subjects.
#' @param standardize Logical; when `TRUE` (default) latent features are
#'   standardised within each training fold before fitting. Stored weight
#'   vectors are converted back to the original latent scale.
#' @return An object of class `dkge_clf` containing fold models, per-subject
#'   weight vectors, and metadata required by downstream mapping utilities.
#' @importFrom glmnet cv.glmnet
#' @export
#' @examples
#' \dontrun{
#' clf <- dkge_cv_train_latent_classifier(fit, y, folds = 5)
#' length(clf$beta_by_subject)
#' }
dkge_cv_train_latent_classifier <- function(fit, y,
                                            Z_by_subject = NULL,
                                            folds = 5,
                                            model = c("lda", "ridge_logit", "lsvm"),
                                            ridge = 1e-3,
                                            level = c("subject", "sample"),
                                            standardize = TRUE) {
  stopifnot(inherits(fit, "dkge"))
  model <- match.arg(model)
  level <- match.arg(level)

  S <- length(fit$Btil)
  r <- ncol(fit$U)
  stopifnot(r >= 1L)

  if (is.null(Z_by_subject)) {
    Z_by_subject <- dkge_project_clusters_to_latent(fit)
  }
  stopifnot(is.list(Z_by_subject), length(Z_by_subject) == S)
  Z_by_subject <- lapply(Z_by_subject, function(Zs) {
    stopifnot(ncol(Zs) == r)
    as.matrix(Zs)
  })

  if (is.numeric(folds) && length(folds) == 1) {
    folds <- dkge_define_folds(fit, type = "subject", k = folds)
  }
  stopifnot(inherits(folds, "dkge_folds"), length(folds$assignments) == folds$k)

  beta_by_subject <- vector("list", S)
  models_by_fold <- vector("list", folds$k)
  fold_assignment <- integer(S)

  if (level == "subject") {
    stopifnot(length(y) == S)
    y <- as.factor(y)
    if (length(levels(y)) != 2L) {
      stop("Binary classification required; supply a two-level factor for 'y'.")
    }
    Z_bar <- do.call(rbind, lapply(Z_by_subject, colMeans))

    for (k in seq_len(folds$k)) {
      hold <- folds$assignments[[k]]
      train <- setdiff(seq_len(S), hold)
      fold_assignment[hold] <- k

      Z_train_raw <- Z_bar[train, , drop = FALSE]
      y_train <- droplevels(y[train])
      if (length(levels(y_train)) != 2L) {
        stop("Each training fold must contain both classes.")
      }

      std <- if (standardize) .dkge_standardize_matrix(Z_train_raw) else list(X = Z_train_raw, mean = rep(0, r), sd = rep(1, r))
      Z_train_proc <- std$X

      mdl <- .dkge_fit_latent_binary(Z_train_proc, y_train, model = model, ridge = ridge)
      mdl$standardize <- if (standardize) list(mu = std$mean, sd = std$sd) else NULL

      beta_proc <- mdl$beta
      beta_raw <- beta_proc / std$sd
      mdl$beta_proc <- beta_proc
      mdl$beta_raw <- beta_raw
      mdl$beta <- beta_raw

      if (!is.null(mdl$intercept)) {
        mdl$intercept_proc <- mdl$intercept
        shift <- sum(beta_proc * std$mean / std$sd)
        mdl$intercept <- mdl$intercept - shift
      }

      if (identical(mdl$model, "lda")) {
        cls <- mdl$classes
        mu_pos_raw <- colMeans(Z_train_raw[y_train == cls[1], , drop = FALSE])
        mu_neg_raw <- colMeans(Z_train_raw[y_train == cls[2], , drop = FALSE])
        mdl$mu_pos_raw <- mu_pos_raw
        mdl$mu_neg_raw <- mu_neg_raw
        mdl$intercept <- -0.5 * sum(beta_raw * (mu_pos_raw + mu_neg_raw))
      }

      models_by_fold[[k]] <- mdl
      beta_by_subject[hold] <- replicate(length(hold), beta_raw, simplify = FALSE)
    }
  } else {
    stopifnot(is.list(y), length(y) == S)
    for (s in seq_len(S)) {
      stopifnot(length(y[[s]]) == nrow(Z_by_subject[[s]]))
    }

    for (k in seq_len(folds$k)) {
      hold <- folds$assignments[[k]]
      train <- setdiff(seq_len(S), hold)
      fold_assignment[hold] <- k

      Z_train_raw <- do.call(rbind, Z_by_subject[train])
      y_train <- droplevels(as.factor(unlist(y[train])))
      if (length(levels(y_train)) != 2L) {
        stop("Each training fold must contain both classes.")
      }

      std <- if (standardize) .dkge_standardize_matrix(Z_train_raw) else list(X = Z_train_raw, mean = rep(0, r), sd = rep(1, r))
      Z_train_proc <- std$X

      mdl <- .dkge_fit_latent_binary(Z_train_proc, y_train, model = model, ridge = ridge)
      mdl$standardize <- if (standardize) list(mu = std$mean, sd = std$sd) else NULL

      beta_proc <- mdl$beta
      beta_raw <- beta_proc / std$sd
      mdl$beta_proc <- beta_proc
      mdl$beta_raw <- beta_raw
      mdl$beta <- beta_raw

      if (!is.null(mdl$intercept)) {
        mdl$intercept_proc <- mdl$intercept
        shift <- sum(beta_proc * std$mean / std$sd)
        mdl$intercept <- mdl$intercept - shift
      }

      if (identical(mdl$model, "lda")) {
        cls <- mdl$classes
        mu_pos_raw <- colMeans(Z_train_raw[y_train == cls[1], , drop = FALSE])
        mu_neg_raw <- colMeans(Z_train_raw[y_train == cls[2], , drop = FALSE])
        mdl$mu_pos_raw <- mu_pos_raw
        mdl$mu_neg_raw <- mu_neg_raw
        mdl$intercept <- -0.5 * sum(beta_raw * (mu_pos_raw + mu_neg_raw))
      }

      models_by_fold[[k]] <- mdl
      beta_by_subject[hold] <- replicate(length(hold), beta_raw, simplify = FALSE)
    }
  }

  if (any(vapply(beta_by_subject, is.null, logical(1)))) {
    stop("Some subjects were not assigned a weight vector; check fold definitions.")
  }

  structure(list(
    model = model,
    level = level,
    folds = folds,
    models_by_fold = models_by_fold,
    beta_by_subject = beta_by_subject,
    fold_assignment = fold_assignment,
    r = r,
    call = match.call()
  ), class = "dkge_clf")
}

#' @keywords internal
#' @noRd
.dkge_fit_latent_binary <- function(Z, y, model = c("lda", "ridge_logit", "lsvm"), ridge = 1e-3) {
  model <- match.arg(model)
  stopifnot(is.matrix(Z), length(y) == nrow(Z))
  y <- droplevels(as.factor(y))
  if (length(levels(y)) != 2L) {
    stop("Binary classification required for latent classifiers.")
  }

  res <- switch(model,
    lda = .dkge_fit_lda(Z, y, ridge = ridge),
    ridge_logit = .dkge_fit_ridge_logit(Z, y),
    lsvm = {
      warning("'lsvm' placeholder falls back to LDA until a native solver is provided.")
      .dkge_fit_lda(Z, y, ridge = ridge)
    }
  )
  res$model <- model
  res
}

#' @keywords internal
#' @noRd
.dkge_fit_lda <- function(Z, y, ridge = 1e-3) {
  cls <- levels(y)
  Z1 <- Z[y == cls[1], , drop = FALSE]
  Z2 <- Z[y == cls[2], , drop = FALSE]
  mu1 <- colMeans(Z1)
  mu2 <- colMeans(Z2)

  n1 <- nrow(Z1)
  n2 <- nrow(Z2)
  S1 <- if (n1 > 1) stats::cov(Z1) else matrix(0, ncol(Z), ncol(Z))
  S2 <- if (n2 > 1) stats::cov(Z2) else matrix(0, ncol(Z), ncol(Z))
  pooled <- ((n1 - 1) * S1 + (n2 - 1) * S2)
  denom <- max(n1 + n2 - 2, 1)
  Sigma <- pooled / denom
  Sigma <- (Sigma + t(Sigma)) / 2
  diag(Sigma) <- diag(Sigma) + ridge

  beta <- tryCatch(solve(Sigma, mu1 - mu2), error = function(e) {
    eig <- eigen(Sigma, symmetric = TRUE)
    vals <- pmax(eig$values, 1e-8)
    eig$vectors %*% ((t(eig$vectors) %*% (mu1 - mu2)) / vals)
  })
  beta <- as.numeric(beta)
  if (sum(beta * (mu1 - mu2)) < 0) {
    beta <- -beta
  }

  list(
    beta = beta,
    intercept = -0.5 * sum(beta * (mu1 + mu2)),
    Sigma = Sigma,
    mu_pos = mu1,
    mu_neg = mu2,
    classes = cls
  )
}

#' @keywords internal
#' @noRd
.dkge_fit_ridge_logit <- function(Z, y) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    warning("glmnet not available; falling back to LDA.")
    return(.dkge_fit_lda(Z, y, ridge = 1e-3))
  }
  cls <- levels(y)
  ybin <- as.integer(y == cls[1])
  cvfit <- glmnet::cv.glmnet(Z, ybin, family = "binomial", alpha = 0, standardize = FALSE)
  coef_vec <- as.numeric(stats::coef(cvfit, s = "lambda.min"))
  list(
    beta = coef_vec[-1],
    intercept = coef_vec[1],
    Sigma = NULL,
    mu_pos = NULL,
    mu_neg = NULL,
    classes = cls
  )
}
