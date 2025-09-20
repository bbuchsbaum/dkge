# dkge-transport.R
# Transport DKGE maps to a common reference parcellation via pluggable mappers.

.dkge_pairwise_sqdist <- function(A, B) {
  n <- nrow(A); m <- nrow(B)
  An <- rowSums(A * A); Bn <- rowSums(B * B)
  outer(An, rep(1, m)) + outer(rep(1, n), Bn) - 2 * (A %*% t(B))
}

.dkge_cost_matrix <- function(Aemb_s, Aemb_ref, X_s = NULL, X_ref = NULL,
                              lambda_emb = 1, lambda_spa = 0.5, sigma_mm = 15,
                              sizes_s = NULL, sizes_ref = NULL, lambda_size = 0) {
  C_emb <- .dkge_pairwise_sqdist(Aemb_s, Aemb_ref)
  C_spa <- 0
  if (!is.null(X_s) && !is.null(X_ref)) {
    C_spa <- .dkge_pairwise_sqdist(X_s / sigma_mm, X_ref / sigma_mm)
  }
  C <- lambda_emb * C_emb + lambda_spa * C_spa
  if (!is.null(sizes_s) && !is.null(sizes_ref) && lambda_size > 0) {
    la <- log(pmax(sizes_s, 1))
    lb <- log(pmax(sizes_ref, 1))
    C <- C + lambda_size * outer(la, lb, function(x, y) (x - y)^2)
  }
  C
}

.dkge_logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) {
    return(m)
  }
  m + log(sum(exp(x - m)))
}

.dkge_sinkhorn_cache <- new.env(parent = emptyenv())

.dkge_sinkhorn_plan <- function(C, mu, nu, epsilon = 0.05, max_iter = 200, tol = 1e-6) {
  stopifnot(is.matrix(C), length(mu) == nrow(C), length(nu) == ncol(C))
  stopifnot(all(mu > 0), all(nu > 0))
  if (abs(sum(mu) - sum(nu)) > 1e-6) {
    stop("mu and nu must sum to the same total mass")
  }

  key <- paste(nrow(C), ncol(C), signif(epsilon, 8), sep = "|")
  state <- .dkge_sinkhorn_cache[[key]]
  log_u_init <- if (!is.null(state)) state$log_u else NULL
  log_v_init <- if (!is.null(state)) state$log_v else NULL

  res <- sinkhorn_plan_cpp(C, mu, nu, epsilon, as.integer(max_iter), tol,
                           log_u_init = log_u_init,
                           log_v_init = log_v_init,
                           keep_duals = TRUE)
  plan <- res$plan
  if (!is.null(res$log_u) && !is.null(res$log_v)) {
    .dkge_sinkhorn_cache[[key]] <- list(log_u = res$log_u,
                                        log_v = res$log_v,
                                        iterations = res$iterations)
  }
  plan
}

# ---------------------------------------------------------------------------
# Helper for mapping lists --------------------------------------------------
# ---------------------------------------------------------------------------

.dkge_prepare_mapping_inputs <- function(loadings, centroids, sizes) {
  S <- length(loadings)
  stopifnot(length(centroids) == S)
  if (is.null(sizes)) {
    sizes <- lapply(loadings, function(A) rep(1, nrow(A)))
  } else {
    stopifnot(length(sizes) == S)
    sizes <- Map(function(sz, A) if (is.null(sz)) rep(1, nrow(A)) else sz,
                 sizes, loadings)
  }
  feat <- lapply(loadings, function(A) {
    A <- as.matrix(A)
    norms <- pmax(sqrt(rowSums(A^2)), 1e-8)
    A / norms
  })
  list(features = feat, sizes = sizes)
}

.dkge_run_mapper <- function(mapper_spec, value_list, feature_list, feature_ref,
                             centroids_list, centroid_ref,
                             sizes_list, size_ref,
                             medoid,
                             operators = NULL) {
  S <- length(feature_list)
  Q <- nrow(feature_ref)
  compute_values <- !is.null(value_list)
  if (compute_values) {
    stopifnot(length(value_list) == S)
  }

  mapped <- if (compute_values) matrix(NA_real_, S, Q) else NULL
  operators_out <- vector("list", S)

  for (s in seq_len(S)) {
    use_cached <- !is.null(operators) && length(operators) >= s &&
      !is.null(operators[[s]])

    if (use_cached) {
      operator <- operators[[s]]
    } else if (s == medoid) {
      operator <- diag(1, Q)
    } else {
      map_fit <- fit_mapper(mapper_spec,
                            source_feat = feature_list[[s]],
                            target_feat = feature_ref,
                            source_weights = sizes_list[[s]],
                            target_weights = size_ref,
                            source_xyz = centroids_list[[s]],
                            target_xyz = centroid_ref)
      operator <- map_fit$operator
    }

    operators_out[[s]] <- operator

    if (compute_values) {
      vals <- value_list[[s]]
      if (s == medoid && length(vals) != Q) {
        stop("Medoid subject values must match reference dimension")
      }
      mapped[s, ] <- as.numeric(.dkge_apply_operator(operator, vals))
    }
  }

  list(
    value = if (compute_values) apply(mapped, 2, stats::median, na.rm = TRUE) else NULL,
    subj_values = mapped,
    plans = operators_out,
    mapper = list(strategy = mapper_spec$strategy,
                  params = mapper_spec$params),
    feature_ref = feature_ref,
    size_ref = size_ref
  )
}

.dkge_resolve_mapper_spec <- function(mapper, method, dots) {
  if (inherits(mapper, "dkge_mapper_spec")) {
    return(mapper)
  }

  if (is.character(mapper)) {
    mapper <- dkge_mapper_spec(mapper)
  }

  if (inherits(mapper, "dkge_mapper_spec")) {
    if (length(dots)) {
      mapper$params[names(dots)] <- dots
    }
    return(mapper)
  }

  strategy <- method %||% "sinkhorn"
  allowed <- c("epsilon", "max_iter", "tol", "lambda_emb", "lambda_spa",
               "sigma_mm", "lambda_size")
  params <- dots[intersect(names(dots), allowed)]
  if (length(params)) {
    do.call(dkge_mapper_spec, c(list(strategy = strategy), params))
  } else {
    dkge_mapper_spec(strategy)
  }
}

.dkge_transport_to_medoid <- function(mapper_spec, value_list, loadings,
                                      centroids, sizes, medoid,
                                      transport_cache = NULL) {
  if (is.null(transport_cache)) {
    prep <- .dkge_prepare_mapping_inputs(loadings, centroids, sizes)
    feature_list <- prep$features
    size_list <- prep$sizes
    operators <- NULL
    feature_ref <- feature_list[[medoid]]
    size_ref <- size_list[[medoid]]
  } else {
    feature_list <- transport_cache$feature_list
    size_list <- transport_cache$size_list
    operators <- transport_cache$operators
    feature_ref <- transport_cache$feature_ref
    size_ref <- transport_cache$size_ref
  }

  centroid_ref <- centroids[[medoid]]
  res <- .dkge_run_mapper(mapper_spec, value_list, feature_list, feature_ref,
                          centroids, centroid_ref, size_list, size_ref,
                          medoid,
                          operators = operators)

  res$feature_list <- feature_list
  res$size_list <- size_list
  res$centroids <- centroids
  res$medoid <- medoid
  res
}

#' Prepare subject-to-medoid transport operators for reuse
#'
#' Computes and caches subject-to-medoid transport matrices so downstream
#' routines (e.g. bootstraps) can reuse a fixed consensus mapping without
#' re-solving the transport problem on every call.
#'
#' @param fit A `dkge` object.
#' @param centroids List of subject centroid matrices. Defaults to the centroids
#'   stored on `fit` or `fit$input`.
#' @param loadings Optional list of subject loadings (`P_s x r`). When omitted,
#'   they are recomputed from `fit$Btil` or the supplied `betas`.
#' @param betas Optional list of subject betas used to recompute loadings when
#'   `loadings` is `NULL`.
#' @param sizes Optional list of cluster masses.
#' @param mapper Mapper specification or shorthand passed to
#'   [dkge_mapper_spec()].
#' @param medoid Index (1-based) of the reference subject.
#' @param ... Additional mapper arguments such as `epsilon` or `lambda_spa`.
#'
#' @return A list containing cached transport objects: `operators`,
#'   `mapper_spec`, `feature_list`, `size_list`, `feature_ref`, `size_ref`,
#'   `centroids`, and `medoid`.
#' @export
dkge_prepare_transport <- function(fit,
                                   centroids = NULL,
                                   loadings = NULL,
                                   betas = NULL,
                                   sizes = NULL,
                                   mapper = "sinkhorn",
                                   medoid = 1L,
                                   ...) {
  stopifnot(inherits(fit, "dkge"))

  if (is.null(loadings)) {
    if (!is.null(betas)) {
      loadings <- dkge_predict_loadings(fit, betas)
    } else if (!is.null(fit$Btil)) {
      loadings <- lapply(fit$Btil, function(Bts) t(Bts) %*% fit$K %*% fit$U)
    } else {
      stop("Provide betas or pre-computed loadings.")
    }
  }

  if (is.null(centroids)) {
    centroids <- fit$centroids %||% fit$input$centroids %||%
      stop("Centroids required for transport preparation; none found in fit or arguments.")
  }

  stopifnot(length(centroids) == length(loadings))

  mapper_spec <- .dkge_resolve_mapper_spec(mapper, method = NULL, dots = list(...))

  prep <- .dkge_prepare_mapping_inputs(loadings, centroids, sizes)
  feature_list <- prep$features
  size_list <- prep$sizes
  feature_ref <- feature_list[[medoid]]
  size_ref <- size_list[[medoid]]
  centroid_ref <- centroids[[medoid]]

  mapper_run <- .dkge_run_mapper(mapper_spec,
                                 value_list = NULL,
                                 feature_list = feature_list,
                                 feature_ref = feature_ref,
                                 centroids_list = centroids,
                                 centroid_ref = centroid_ref,
                                 sizes_list = size_list,
                                 size_ref = size_ref,
                                 medoid = medoid,
                                 operators = NULL)

  list(
    operators = mapper_run$plans,
    mapper_spec = mapper_spec,
    feature_list = feature_list,
    size_list = size_list,
    feature_ref = mapper_run$feature_ref,
    size_ref = mapper_run$size_ref,
    centroids = centroids,
    medoid = medoid
  )
}

# ---------------------------------------------------------------------------
# Public wrappers -----------------------------------------------------------
# ---------------------------------------------------------------------------

#' Transport cluster values to a medoid via entropic Sinkhorn OT
#'
#' @param v_list List of subject-level cluster values (length P_s each).
#' @param A_list List of subject loadings (P_s × r).
#' @param centroids List of subject cluster centroids (each P_s × 3 matrix).
#' @param sizes Optional list of cluster masses (defaults to uniform weights).
#' @param medoid Integer index of the reference subject (1-based).
#' @param lambda_emb,lambda_spa Cost weights for embedding and spatial terms.
#' @param sigma_mm Spatial rescaling (millimetres).
#' @param epsilon,max_iter,tol Sinkhorn parameters.
#' @param transport_cache Optional cache returned by [dkge_prepare_transport()].
#'   When supplied, cached operators are reused and the mapper configuration
#'   stored in the cache takes precedence.
#' @return List containing summary statistics, transported subject maps, and
#'   the per-subject transport operators.
#' @export
dkge_transport_to_medoid_sinkhorn <- function(v_list, A_list, centroids, sizes = NULL,
                                              medoid,
                                              lambda_emb = 1, lambda_spa = 0.5, sigma_mm = 15,
                                              epsilon = 0.05, max_iter = 200, tol = 1e-6,
                                              transport_cache = NULL) {
  if (!is.null(transport_cache)) {
    mapper_spec <- transport_cache$mapper_spec
  } else {
    mapper_spec <- dkge_mapper_spec("sinkhorn",
                                    lambda_emb = lambda_emb,
                                    lambda_spa = lambda_spa,
                                    sigma_mm = sigma_mm,
                                    epsilon = epsilon,
                                    max_iter = max_iter,
                                    tol = tol)
  }
  .dkge_transport_to_medoid(mapper_spec, v_list, A_list, centroids, sizes, medoid,
                            transport_cache = transport_cache)
}

#' @rdname dkge_transport_to_medoid_sinkhorn
#' @param return_plans Logical; if TRUE, include transport plans in the output.
#' @export
dkge_transport_to_medoid_sinkhorn_cpp <- function(v_list, A_list, centroids, sizes = NULL,
                                                  medoid,
                                                  lambda_emb = 1, lambda_spa = 0.5, sigma_mm = 15,
                                                  epsilon = 0.05, max_iter = 300, tol = 1e-6,
                                                  return_plans = FALSE,
                                                  transport_cache = NULL) {
  res <- dkge_transport_to_medoid_sinkhorn(v_list, A_list, centroids, sizes = sizes,
                                           medoid = medoid,
                                           lambda_emb = lambda_emb,
                                           lambda_spa = lambda_spa,
                                           sigma_mm = sigma_mm,
                                           epsilon = epsilon,
                                           max_iter = max_iter,
                                           tol = tol,
                                           transport_cache = transport_cache)
  if (!return_plans) {
    res$plans <- NULL
  }
  res
}

#' Transport component loadings to a medoid parcellation
#'
#' @param fit A `dkge` object used to compute the loadings.
#' @param medoid Integer index of the reference subject (1-based).
#' @param centroids List of subject cluster centroids (each P_s × 3 matrix).
#' @param loadings Optional list of subject loadings (P_s × r). When omitted,
#'   they are recomputed from `betas`.
#' @param betas Optional list of subject betas used to recompute loadings when
#'   `loadings` is `NULL`.
#' @param sizes Optional list of cluster masses (defaults to uniform weights).
#' @param mapper Optional mapper specification created by [dkge_mapper_spec()].
#'   When `NULL`, defaults to Sinkhorn with the supplied parameters.
#' @param method Backwards-compatible alias; currently only "sinkhorn" is
#'   supported.
#' @param transport_cache Optional cache from [dkge_prepare_transport()]. When
#'   supplied, cached operators are reused for all components.
#' @param ... Additional parameters passed when building the default mapper
#'   specification (e.g. `epsilon`, `lambda_emb`).
#' @return List with `group` (medoid cluster vectors per component),
#'   `subjects` (per-subject transported values), and `cache` (transport cache
#'   reused for future calls).
#' @export
dkge_transport_loadings_to_medoid <- function(fit, medoid, centroids,
                                               loadings = NULL,
                                               betas = NULL,
                                               sizes = NULL,
                                               mapper = NULL,
                                               method = c("sinkhorn", "sinkhorn_cpp"),
                                               transport_cache = NULL,
                                               ...) {
  stopifnot(inherits(fit, "dkge"))
  if (is.null(loadings)) {
    if (!is.null(betas)) {
      loadings <- dkge_predict_loadings(fit, betas)
    } else if (!is.null(fit$Btil)) {
      loadings <- lapply(fit$Btil, function(Bts) t(Bts) %*% fit$K %*% fit$U)
    } else {
      stop("Provide betas or pre-computed loadings.")
    }
  }
  S <- length(loadings)
  stopifnot(length(centroids) == S)

  dots <- list(...)
  if (!is.null(transport_cache)) {
    mapper_spec <- transport_cache$mapper_spec
  } else {
    mapper_spec <- .dkge_resolve_mapper_spec(mapper, method = method[1], dots = dots)
  }

  if (is.null(sizes)) {
    sizes <- lapply(loadings, function(A) rep(1, nrow(A)))
  } else {
    sizes <- Map(function(sz, A) if (is.null(sz)) rep(1, nrow(A)) else sz, sizes, loadings)
  }
  rank <- ncol(loadings[[1]])

  subj_vals <- vector("list", rank)
  group_vals <- vector("list", rank)
  cache_local <- transport_cache
  for (j in seq_len(rank)) {
    v_list <- lapply(loadings, function(A) A[, j])
    tr <- .dkge_transport_to_medoid(mapper_spec, v_list, loadings, centroids,
                                    sizes = sizes, medoid = medoid,
                                    transport_cache = cache_local)
    subj_vals[[j]] <- tr$subj_values
    group_vals[[j]] <- tr$value
    if (is.null(cache_local)) {
      cache_local <- list(
        operators = tr$plans,
        mapper_spec = mapper_spec,
        feature_list = tr$feature_list,
        size_list = tr$size_list,
        feature_ref = tr$feature_ref,
        size_ref = tr$size_ref,
        centroids = tr$centroids,
        medoid = tr$medoid
      )
    }
  }
  list(group = group_vals, subjects = subj_vals, cache = cache_local)
}


#' Transport subject contrasts to a medoid parcellation
#'
#' @param fit A `dkge` object used to compute the contrasts.
#' @param contrast_obj A `dkge_contrasts` result.
#' @param medoid Integer index of the reference subject (1-based).
#' @param centroids List of subject cluster centroids (each P_s × 3 matrix).
#' @param loadings Optional list of subject loadings (P_s × r).
#' @param betas Optional list of subject betas used to recompute loadings.
#' @param sizes Optional list of cluster masses.
#' @param mapper Optional mapper specification created by [dkge_mapper_spec()].
#' @param method Backwards-compatible alias; currently only "sinkhorn" is
#'   supported.
#' @param transport_cache Optional cache from [dkge_prepare_transport()]. When
#'   supplied, cached operators are reused for every contrast.
#' @param ... Additional parameters passed when building the default mapper
#'   specification.
#' @return Named list of transport results (one per contrast) with an attached
#'   `cache` element for reuse.
#' @export
dkge_transport_contrasts_to_medoid <- function(fit, contrast_obj, medoid, centroids = NULL,
                                               loadings = NULL,
                                               betas = NULL,
                                               sizes = NULL,
                                               mapper = NULL,
                                               method = c("sinkhorn", "sinkhorn_cpp"),
                                               transport_cache = NULL,
                                               ...) {
  stopifnot(inherits(fit, "dkge"), inherits(contrast_obj, "dkge_contrasts"))
  if (is.null(loadings)) {
    if (!is.null(betas)) {
      loadings <- dkge_predict_loadings(fit, betas)
    } else if (!is.null(fit$Btil)) {
      loadings <- lapply(fit$Btil, function(Bts) t(Bts) %*% fit$K %*% fit$U)
    } else {
      stop("Provide betas or pre-computed loadings.")
    }
  }
  S <- length(loadings)
  if (is.null(centroids)) {
    centroids <- fit$centroids %||% fit$input$centroids %||%
      stop("Centroids required for transport; none found in fit or arguments.")
  }
  stopifnot(length(centroids) == S)

  dots <- list(...)
  if (!is.null(transport_cache)) {
    mapper_spec <- transport_cache$mapper_spec
  } else {
    mapper_spec <- .dkge_resolve_mapper_spec(mapper, method = method[1], dots = dots)
  }

  if (is.null(sizes)) {
    sizes <- lapply(loadings, function(A) rep(1, nrow(A)))
  } else {
    sizes <- Map(function(sz, A) if (is.null(sz)) rep(1, nrow(A)) else sz, sizes, loadings)
  }

  cache_local <- transport_cache
  out <- vector("list", length(contrast_obj$values))
  for (i in seq_along(contrast_obj$values)) {
    tr <- .dkge_transport_to_medoid(mapper_spec,
                                    contrast_obj$values[[i]],
                                    loadings,
                                    centroids,
                                    sizes = sizes,
                                    medoid = medoid,
                                    transport_cache = cache_local)
    out[[i]] <- tr
    if (is.null(cache_local)) {
      cache_local <- list(
        operators = tr$plans,
        mapper_spec = mapper_spec,
        feature_list = tr$feature_list,
        size_list = tr$size_list,
        feature_ref = tr$feature_ref,
        size_ref = tr$size_ref,
        centroids = tr$centroids,
        medoid = tr$medoid
      )
    }
  }
  names(out) <- names(contrast_obj$values)
  attr(out, "cache") <- cache_local
  out
}


#' Paint medoid cluster values back to a label volume
#'
#' @param values Vector or matrix of values defined on the medoid parcellation.
#' @param labels Medoid labels describing the target parcellation.
#' @param out_file Optional path to save the rendered map.
#' @return A `BrainVolume` or file path, depending on `out_file`.
#' @export
dkge_paint_medoid_map <- function(values, labels, out_file = NULL) {
  dkge_write_group_map(values, labels, out_file = out_file)
}
