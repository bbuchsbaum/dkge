# dkge-mapper.R
# Pluggable mapping strategies from subject parcels to a consensus reference.

#' Specify a DKGE mapper strategy
#'
#' @param strategy Mapping strategy identifier ("sinkhorn", "ridge", "ols").
#' @param ... Strategy-specific hyperparameters stored in the specification.
#' @param name Optional user-facing name for diagnostics.
#' @return A mapper specification object.
#' @export
dkge_mapper_spec <- function(strategy = c("sinkhorn", "ridge", "ols"),
                             ..., name = NULL) {
  strategy <- match.arg(strategy)
  spec <- list(strategy = strategy, params = list(...),
               name = name %||% strategy)
class(spec) <- c(paste0("dkge_mapper_spec_", strategy), "dkge_mapper_spec")
  spec
}

#' Create a pluggable DKGE anchor mapper
#'
#' Constructs lightweight mapper objects used by the dense rendering core. The
#' initial implementation exposes a fast kNN barycentric mapper and keeps the
#' interface open for richer transports.
#'
#' @param method Mapper backend identifier. Currently only `"knn"` is
#'   implemented in the core.
#' @param ... Backend-specific parameters stored within the mapper object.
#' @return An S3 mapper descriptor consumed by [fit_mapper()].
#' @export
dkge_mapper <- function(method = c("knn", "sinkhorn", "ridge", "gw"), ...) {
  method <- match.arg(method)
  structure(list(method = method, pars = list(...)),
            class = c(paste0("dkge_mapper_", method), "dkge_mapper"))
}

# Internal registry (reserved for future extensions)
.dkge_mapper_registry <- new.env(parent = emptyenv())

register_mapper <- function(key, constructor) {
  assign(key, constructor, envir = .dkge_mapper_registry)
}

get_mapper_ctor <- function(key) {
  get0(key, envir = .dkge_mapper_registry, ifnotfound = NULL)
}

#' Fit a mapper on subject/reference features
#'
#' @param spec Mapper specification created with [dkge_mapper_spec()] or
#'   [dkge_mapper()].
#' @param ... Strategy-specific arguments (e.g. feature matrices, spatial
#'   coordinates, weights).
#' @return A fitted mapping object (see strategy-specific classes).
#' @export
fit_mapper <- function(spec, ...) {
  UseMethod("fit_mapper", spec)
}

#' Apply a fitted mapper to new source values
#'
#' @param mapping Mapping object returned by [fit_mapper()].
#' @param new_source_vals Matrix or vector of new source values (P_s x K).
#' @param ... Optional arguments used by specific strategies.
#' @return Matrix (Q x K) or vector of mapped values.
#' @export
predict_mapper <- function(mapping, new_source_vals, ...) {
  UseMethod("predict_mapper", mapping)
}

# ---------------------------------------------------------------------------
# kNN barycentric mapper (dense rendering core) -----------------------------
# ---------------------------------------------------------------------------

#' @export
fit_mapper.dkge_mapper_knn <- function(spec,
                                       subj_points,
                                       anchor_points,
                                       subj_feats = NULL,
                                       anchor_feats = NULL,
                                       reliab = NULL,
                                       ...) {
  stopifnot(is.matrix(subj_points), ncol(subj_points) == 3,
            is.matrix(anchor_points), ncol(anchor_points) == 3)

  pars <- spec$pars
  k <- pars$k %||% 8L
  sigx <- pars$sigx %||% 3.0
  sigz <- pars$sigz

  nn <- FNN::get.knnx(anchor_points, subj_points, k = k)
  dist2 <- nn$nn.dist^2
  W <- exp(-dist2 / (2 * sigx^2))

  use_feat <- !is.null(subj_feats) && !is.null(anchor_feats) && !is.null(sigz)
  if (use_feat) {
    stopifnot(nrow(subj_feats) == nrow(subj_points))
    stopifnot(nrow(anchor_feats) == nrow(anchor_points))
    feature_weights <- matrix(1, nrow = nrow(subj_points), ncol = k)
    for (i in seq_len(nrow(subj_points))) {
      js <- nn$nn.index[i, ]
      diffs <- anchor_feats[js, , drop = FALSE] -
        matrix(subj_feats[i, ], nrow = k, ncol = ncol(anchor_feats), byrow = TRUE)
      d2_feat <- rowSums(diffs * diffs)
      feature_weights[i, ] <- exp(-d2_feat / (2 * sigz^2))
    }
    W <- W * feature_weights
  }

  W <- W / (rowSums(W) + 1e-12)

  structure(list(
    type = "knn",
    idx = nn$nn.index,
    weights = W,
    Q = nrow(anchor_points),
    P = nrow(subj_points),
    reliab = reliab
  ), class = "dkge_mapper_fit_knn")
}

#' @export
apply_mapper.dkge_mapper_fit_knn <- function(fitted_mapper,
                                             values,
                                             reliab = NULL,
                                             normalize_by_reliab = TRUE,
                                             ...) {
  vals <- if (is.null(dim(values))) matrix(values, ncol = 1) else as.matrix(values)
  stopifnot(nrow(vals) == fitted_mapper$P)

  r <- reliab %||% fitted_mapper$reliab
  idx <- fitted_mapper$idx
  W <- fitted_mapper$weights
  Q <- fitted_mapper$Q

  if (ncol(vals) > 1) {
    out <- matrix(0, nrow = Q, ncol = ncol(vals))
    for (col in seq_len(ncol(vals))) {
      out[, col] <- Recall(fitted_mapper, vals[, col], reliab = r,
                           normalize_by_reliab = normalize_by_reliab)
    }
    return(out)
  }

  values <- vals[, 1]
  if (!normalize_by_reliab) {
    anchor <- numeric(Q)
    for (t in seq_len(ncol(idx))) {
      js <- idx[, t]
      wt <- W[, t]
      anchor[js] <- anchor[js] + wt * values
    }
    return(anchor)
  }

  if (is.null(r)) {
    r <- rep(1, fitted_mapper$P)
  } else {
    stopifnot(length(r) == fitted_mapper$P)
  }

  y_num <- numeric(Q)
  y_den <- numeric(Q)
  rv <- r * values
  for (t in seq_len(ncol(idx))) {
    js <- idx[, t]
    wt <- W[, t]
    y_num[js] <- y_num[js] + wt * rv
    y_den[js] <- y_den[js] + wt * r
  }
  anchor <- numeric(Q)
  ok <- y_den > 1e-12
  anchor[ok] <- y_num[ok] / y_den[ok]
  anchor
}

# Placeholders for alternative backends (implemented elsewhere)
#' @export
fit_mapper.dkge_mapper_sinkhorn <- function(spec,
                                            subj_points,
                                            anchor_points,
                                            subj_feats = NULL,
                                            anchor_feats = NULL,
                                            reliab = NULL,
                                            anchor_weights = NULL,
                                            ...) {
  stopifnot(is.matrix(subj_points), ncol(subj_points) == 3,
            is.matrix(anchor_points), ncol(anchor_points) == 3)

  pars <- spec$pars
  epsilon <- pars$epsilon %||% 0.05
  max_iter <- pars$max_iter %||% 300
  tol <- pars$tol %||% 1e-6
  sigx <- pars$sigx %||% 5.0
  lambda_xyz <- pars$lambda_xyz %||% 1.0
  sigz <- pars$sigz %||% 1.0
  lambda_feat <- pars$lambda_feat
  if (is.null(lambda_feat)) {
    lambda_feat <- if (!is.null(subj_feats) && !is.null(anchor_feats)) 1.0 else 0
  }

  P <- nrow(subj_points)
  Q <- nrow(anchor_points)

  C <- matrix(0, P, Q)
  if (lambda_xyz > 0) {
    C <- C + lambda_xyz * .dkge_pairwise_sqdist(subj_points / sigx,
                                                anchor_points / sigx)
  }
  if (lambda_feat > 0) {
    stopifnot(!is.null(subj_feats), !is.null(anchor_feats))
    C <- C + lambda_feat * .dkge_pairwise_sqdist(subj_feats / sigz,
                                                 anchor_feats / sigz)
  }

  source_weights <- pars$source_weights %||% reliab
  if (is.null(source_weights)) {
    source_weights <- rep(1, P)
  }
  source_weights <- as.numeric(source_weights)
  stopifnot(length(source_weights) == P)
  if (any(source_weights < 0)) {
    stop("Source weights must be non-negative for Sinkhorn mapper.")
  }
  if (sum(source_weights) <= 0) {
    source_weights <- rep(1, P)
  }
  mu <- source_weights / sum(source_weights)

  anchor_weights <- anchor_weights %||% pars$anchor_weights
  if (is.null(anchor_weights)) {
    anchor_weights <- rep(1, Q)
  }
  anchor_weights <- as.numeric(anchor_weights)
  stopifnot(length(anchor_weights) == Q)
  if (any(anchor_weights < 0)) {
    stop("Anchor weights must be non-negative for Sinkhorn mapper.")
  }
  if (sum(anchor_weights) <= 0) {
    anchor_weights <- rep(1, Q)
  }
  nu <- anchor_weights / sum(anchor_weights)

  plan <- .dkge_sinkhorn_plan(C,
                              mu = mu,
                              nu = nu,
                              epsilon = epsilon,
                              max_iter = max_iter,
                              tol = tol)

  plan_sparse <- Matrix::Matrix(plan, sparse = TRUE)
  transport_cost <- sum(plan * C)

  structure(list(
    type = "sinkhorn",
    plan = plan_sparse,
    cost = C,
    mu = mu,
    nu = nu,
    epsilon = epsilon,
    max_iter = max_iter,
    tol = tol,
    reliab = reliab,
    stats = list(transport_cost = transport_cost),
    P = P,
    Q = Q
  ), class = "dkge_mapper_fit_sinkhorn")
}

#' @export
apply_mapper.dkge_mapper_fit_sinkhorn <- function(fitted_mapper,
                                                  values,
                                                  reliab = NULL,
                                                  normalize_by_reliab = TRUE,
                                                  ...) {
  vals <- if (is.null(dim(values))) matrix(values, ncol = 1) else as.matrix(values)
  stopifnot(nrow(vals) == fitted_mapper$P)

  plan <- fitted_mapper$plan
  if (!inherits(plan, "Matrix")) {
    plan <- Matrix::Matrix(plan, sparse = TRUE)
  }

  plan_t <- Matrix::t(plan)

  if (ncol(vals) == 1) {
    if (!normalize_by_reliab) {
      return(as.numeric(plan_t %*% vals))
    }

    r <- reliab %||% fitted_mapper$reliab %||% fitted_mapper$mu
    r <- as.numeric(r)
    stopifnot(length(r) == fitted_mapper$P)
    if (any(r < 0)) {
      stop("Reliability weights must be non-negative for Sinkhorn mapper.")
    }

    numer <- as.numeric(plan_t %*% (vals * r))
    denom_raw <- as.numeric(plan_t %*% r)
    denom_safe <- pmax(denom_raw, 1e-12)
    out <- numer / denom_safe
    out[denom_raw <= 1e-12] <- 0
    return(out)
  }

  if (!normalize_by_reliab) {
    res <- plan_t %*% vals
    return(as.matrix(res))
  }

  r <- reliab %||% fitted_mapper$reliab %||% fitted_mapper$mu
  r <- as.numeric(r)
  stopifnot(length(r) == fitted_mapper$P)
  if (any(r < 0)) {
    stop("Reliability weights must be non-negative for Sinkhorn mapper.")
  }

  numer <- plan_t %*% (vals * r)
  denom_raw <- as.numeric(plan_t %*% r)
  denom_safe <- pmax(denom_raw, 1e-12)
  res <- sweep(as.matrix(numer), 1, denom_safe, "/")
  if (any(denom_raw <= 1e-12)) {
    res[denom_raw <= 1e-12, ] <- 0
  }
  res
}

#' @export
fit_mapper.dkge_mapper_ridge <- function(spec, ...) {
  stop("Ridge mapper not implemented in the dense rendering core. ",
       "Provide a 'dkge_mapper_ridge' plugin that defines fit_mapper() and apply_mapper().")
}

#' @export
fit_mapper.dkge_mapper_gw <- function(spec, ...) {
  stop("Gromov-Wasserstein mapper not implemented in the dense rendering core. ",
       "Provide a 'dkge_mapper_gw' plugin that defines fit_mapper() and apply_mapper().")
}

#' Apply a fitted mapper to values
#'
#' Alias to [predict_mapper()] that mirrors the terminology used by the
#' dense rendering helpers.
#'
#' @param fitted_mapper Mapping object returned by [fit_mapper()].
#' @param values Numeric vector or matrix of source-space values.
#' @param ... Optional backend-specific arguments (e.g. reliabilities).
#' @return Numeric vector or matrix of mapped values.
#' @export
apply_mapper <- function(fitted_mapper, values, ...) {
  UseMethod("apply_mapper", fitted_mapper)
}

# ---------------------------------------------------------------------------
# Sinkhorn optimal transport mapper ----------------------------------------
# ---------------------------------------------------------------------------

#' @export
fit_mapper.dkge_mapper_spec_sinkhorn <- function(spec, source_feat, source_vals = NULL,
                                                 target_feat,
                                                 source_weights = NULL,
                                                 target_weights = NULL,
                                                 source_xyz = NULL,
                                                 target_xyz = NULL,
                                                 ...) {
  params <- spec$params
  epsilon <- params$epsilon %||% 0.05
  max_iter <- params$max_iter %||% 300
  tol <- params$tol %||% 1e-6
  lambda_emb <- params$lambda_emb %||% 1
  lambda_spa <- params$lambda_spa %||% 0.5
  sigma_mm <- params$sigma_mm %||% 15
  lambda_size <- params$lambda_size %||% 0

  n <- nrow(source_feat)
  m <- nrow(target_feat)

  a <- .dkge_normalize_weights(source_weights, n)
  b <- .dkge_normalize_weights(target_weights, m)

  C <- .dkge_cost_matrix(source_feat, target_feat,
                         X_s = source_xyz,
                         X_ref = target_xyz,
                         lambda_emb = lambda_emb,
                         lambda_spa = lambda_spa,
                         sigma_mm = sigma_mm,
                         sizes_s = a,
                         sizes_ref = b,
                         lambda_size = lambda_size)

  plan <- .dkge_sinkhorn_plan(C, mu = a, nu = b,
                              epsilon = epsilon,
                              max_iter = max_iter,
                              tol = tol)

  mapping <- list(operator = plan,
                  strategy = "sinkhorn",
                  fit_info = list(epsilon = epsilon,
                                  max_iter = max_iter,
                                  tol = tol,
                                  cost = C,
                                  weights = list(source = a, target = b)))
  class(mapping) <- c("dkge_mapping_sinkhorn", "dkge_mapping")
  mapping
}

#' @export
predict_mapper.dkge_mapping_sinkhorn <- function(mapping, new_source_vals, ...) {
  .dkge_apply_operator(mapping$operator, new_source_vals)
}

# ---------------------------------------------------------------------------
# Ridge / OLS mapper --------------------------------------------------------
# ---------------------------------------------------------------------------

#' @export
fit_mapper.dkge_mapper_spec_ridge <- function(spec, source_feat, source_vals = NULL,
                                              target_feat,
                                              source_weights = NULL,
                                              target_weights = NULL,
                                              ...) {
  params <- spec$params
  lambda <- params$lambda %||% 1e-2

  A <- as.matrix(source_feat)
  B <- as.matrix(target_feat)
  r <- ncol(A)

  M <- crossprod(A)
  if (lambda > 0) {
    M <- M + diag(lambda, r)
  }

  solve_rhs <- tryCatch(solve(M, t(B)), error = function(e) {
    solve(M + 1e-8 * diag(r), t(B))
  })
  coeff <- A %*% solve_rhs

  mapping <- list(operator = coeff,
                  strategy = "ridge",
                  fit_info = list(lambda = lambda))
  class(mapping) <- c("dkge_mapping_ridge", "dkge_mapping")
  mapping
}

#' @export
predict_mapper.dkge_mapping_ridge <- function(mapping, new_source_vals, ...) {
  .dkge_apply_operator(mapping$operator, new_source_vals)
}

#' @export
fit_mapper.dkge_mapper_spec_ols <- function(spec, source_feat, source_vals = NULL,
                                            target_feat, ...) {
  params <- spec$params
  params$lambda <- params$lambda %||% 0
  spec$params <- params
  mapping <- fit_mapper.dkge_mapper_spec_ridge(spec, source_feat, source_vals,
                                               target_feat, ...)
  mapping$strategy <- "ols"
  class(mapping) <- c("dkge_mapping_ols", "dkge_mapping")
  mapping
}

#' @export
predict_mapper.dkge_mapping_ols <- predict_mapper.dkge_mapping_ridge

# ---------------------------------------------------------------------------
# Default predict method ----------------------------------------------------
# ---------------------------------------------------------------------------

#' @export
predict_mapper.dkge_mapping <- function(mapping, new_source_vals, ...) {
  if (is.null(mapping$operator)) {
    stop("Mapping does not expose an operator for prediction.")
  }
  .dkge_apply_operator(mapping$operator, new_source_vals)
}

#' @export
apply_mapper.dkge_mapping <- function(fitted_mapper, values, ...) {
  predict_mapper(fitted_mapper, values, ...)
}

# Helper to normalize weight vectors
.dkge_normalize_weights <- function(w, n) {
  if (is.null(w)) {
    return(rep(1 / n, n))
  }
  stopifnot(length(w) == n)
  w <- as.numeric(w)
  if (any(w < 0)) {
    stop("Weights must be non-negative")
  }
  total <- sum(w)
  if (!is.finite(total) || total <= 0) {
    stop("Weights must sum to a positive finite value")
  }
  w / total
}

# Apply operator (P x Q) to source values (P x K) returning Q x K
.dkge_apply_operator <- function(operator, new_source_vals) {
  if (is.null(operator)) stop("Operator is NULL")
  vals <- new_source_vals
  if (is.null(dim(vals))) {
    vals <- matrix(vals, ncol = 1)
  } else {
    vals <- as.matrix(vals)
  }
  if (nrow(vals) != nrow(operator)) {
    stop("Source values do not match operator dimensions")
  }
  mapped <- t(operator) %*% vals
  if (ncol(mapped) == 1) {
    drop(mapped)
  } else {
    mapped
  }
}
