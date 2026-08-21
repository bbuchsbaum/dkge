# dkge-specs.R
# User-facing helper constructors for orchestration specs.

#' Transport specification helper
#'
#' Builds a validated transport configuration that can be passed to
#' [dkge_pipeline()] or transport utilities. The helper enforces basic argument
#' checks and provides sensible defaults for Sinkhorn-based mapping.
#'
#' @param centroids List of subject-specific centroid matrices (P_s x d).
#' @param sizes Optional list of cluster sizes (one numeric vector per subject).
#' @param medoid Integer index of the medoid subject (default 1).
#' @param method Mapper backend. Default "sinkhorn".
#' @param mapper Optional prefit mapper specification (advanced use).
#' @param epsilon Sinkhorn entropic regularisation parameter.
#' @param max_iter Maximum Sinkhorn iterations.
#' @param tol Convergence tolerance for Sinkhorn scaling.
#' @param lambda_emb Weight on embedding distance in the cost matrix.
#' @param lambda_spa Weight on spatial distance in the cost matrix.
#' @param sigma_mm Spatial scale (in millimetres) used when spatial coordinates
#'   are available.
#' @param lambda_size Weight on size regularisation between clusters.
#' @param value_type Semantics of transported values: `"intensive"` preserves
#'   constant fields; `"extensive"` preserves the sum of source totals.
#' @param warm_start Reuse converged Sinkhorn duals for identical problems.
#' @param provenance Optional [dkge_transport_provenance()] declaration. When
#'   omitted, loading-derived transport is marked descriptive and
#'   [dkge_infer()] will reject it rather than assume inferential validity.
#' @param ... Additional fields stored on the spec (e.g., precomputed loadings
#'   or betas).
#' @return Object with class `dkge_transport_spec`.
#' @export
#' @examples
#' transport <- dkge_transport_spec(centroids = list(matrix(runif(12), 4, 3)))
dkge_transport_spec <- function(centroids,
                                sizes = NULL,
                                medoid = 1L,
                                method = c("sinkhorn", "ridge", "ols", "sinkhorn_cpp"),
                                mapper = NULL,
                                epsilon = 0.05,
                                max_iter = 5000L,
                                tol = 1e-4,
                                lambda_emb = 1,
                                lambda_spa = 0.5,
                                sigma_mm = 15,
                                lambda_size = 0,
                                value_type = c("intensive", "extensive"),
                                warm_start = TRUE,
                                provenance = NULL,
                                ...) {
  method <- match.arg(method)
  if (identical(method, "sinkhorn_cpp")) {
    warning("`sinkhorn_cpp` is deprecated; `sinkhorn` already uses the compiled backend.",
            call. = FALSE)
    method <- "sinkhorn"
  }
  value_type <- match.arg(value_type)
  stopifnot(is.list(centroids), length(centroids) >= 1L)
  if (!is.null(sizes)) {
    stopifnot(is.list(sizes), length(sizes) == length(centroids))
  }
  medoid <- .dkge_validate_positive_integer(medoid, "medoid")
  if (medoid < 1L || medoid > length(centroids)) {
    stop("`medoid` must index one of the provided centroid lists")
  }
  epsilon <- .dkge_validate_positive_scalar(epsilon, "epsilon")
  max_iter <- .dkge_validate_positive_integer(max_iter, "max_iter")
  tol <- .dkge_validate_positive_scalar(tol, "tol")
  lambda_emb <- .dkge_validate_nonnegative_scalar(lambda_emb, "lambda_emb")
  lambda_spa <- .dkge_validate_nonnegative_scalar(lambda_spa, "lambda_spa")
  sigma_mm <- .dkge_validate_positive_scalar(sigma_mm, "sigma_mm")
  lambda_size <- .dkge_validate_nonnegative_scalar(lambda_size, "lambda_size")

  spec <- list(
    centroids = centroids,
    sizes = sizes,
    medoid = medoid,
    method = method,
    mapper = mapper,
    epsilon = epsilon,
    max_iter = max_iter,
    tol = tol,
    lambda_emb = lambda_emb,
    lambda_spa = lambda_spa,
    sigma_mm = sigma_mm,
    lambda_size = lambda_size,
    value_type = value_type,
    warm_start = isTRUE(warm_start),
    provenance = provenance %||% .dkge_data_derived_loading_provenance()
  )
  extra <- list(...)
  if (length(extra)) spec <- c(spec, extra)
  class(spec) <- c("dkge_transport_spec", "list")
  spec
}

#' Inference specification helper
#'
#' @param B Number of permutations for sign-flip inference.
#' @param tail Tail of the test: "two.sided", "greater", or "less".
#' @param center Location statistic for permutations. Only `"mean"` is
#'   implemented by the beta inference service.
#' @return Object with class `dkge_inference_spec`.
#' @export
#' @examples
#' infer <- dkge_inference_spec(B = 1000, tail = "two.sided")
dkge_inference_spec <- function(B = 2000L,
                                tail = c("two.sided", "greater", "less"),
                                center = "mean") {
  B <- .dkge_validate_resample_B(B)
  tail <- match.arg(tail)
  .dkge_validate_inference_compatibility(
    inference = "signflip", correction = "maxT",
    has_transport = FALSE, center = center,
    n_targets = NA_integer_, parallel = FALSE
  )
  structure(list(B = B, tail = tail, center = center),
            class = c("dkge_inference_spec", "list"))
}

#' Classification specification helper
#'
#' @param targets Target specification accepted by [dkge_classify()].
#' @param method Classifier backend ("lda" or "logit").
#' @param folds Optional fold specification.
#' @param lambda Optional ridge parameter.
#' @param metric Classification metrics to report.
#' @param mode Decoding mode passed to [dkge_classify()]: "auto" selects
#'   automatically, "cell" uses a transductive global-basis embedding,
#'   "cell_cross" uses prospective cross-fitted per-cell embeddings, and
#'   "delta" uses pairwise deltas.
#' @param ... Additional fields stored on the spec (e.g., `n_perm`, `scope`).
#' @return Object with class `dkge_classification_spec`.
#' @export
#' @examples
#' cls <- dkge_classification_spec(targets = ~ condition, method = "lda")
dkge_classification_spec <- function(targets,
                                     method = c("lda", "logit"),
                                     folds = NULL,
                                     lambda = NULL,
                                     metric = c("accuracy", "logloss"),
                                     mode = c("auto", "cell", "cell_cross", "delta"),
                                     ...) {
  method <- match.arg(method)
  metric <- match.arg(metric, several.ok = TRUE)
  mode <- match.arg(mode)
  spec <- list(
    targets = targets,
    method = method,
    folds = folds,
    lambda = lambda,
    metric = unique(metric),
    mode = mode
  )
  extra <- list(...)
  if (length(extra)) spec <- c(spec, extra)
  class(spec) <- c("dkge_classification_spec", "list")
  spec
}

#' Print a transport specification
#'
#' @param x A `dkge_transport_spec` object.
#' @param ... Unused; present for S3 method compatibility.
#' @return `x`, invisibly.
#' @method print dkge_transport_spec
#' @export
print.dkge_transport_spec <- function(x, ...) {
  cat("<dkge_transport_spec>", "\n", sep = "")
  cat("  method   :", x$method, "\n")
  cat("  subjects :", length(x$centroids), "\n")
  cat("  medoid   :", x$medoid, "\n")
  if (!is.null(x$epsilon)) cat("  epsilon  :", x$epsilon, "\n")
  invisible(x)
}

#' Print an inference specification
#'
#' @param x A `dkge_inference_spec` object.
#' @param ... Unused; present for S3 method compatibility.
#' @return `x`, invisibly.
#' @method print dkge_inference_spec
#' @export
print.dkge_inference_spec <- function(x, ...) {
  cat("<dkge_inference_spec>", "\n", sep = "")
  cat("  permutations :", x$B, "\n")
  cat("  tail         :", x$tail, "\n")
  cat("  center       :", x$center, "\n")
  invisible(x)
}

#' Print a classification specification
#'
#' @param x A `dkge_classification_spec` object.
#' @param ... Unused; present for S3 method compatibility.
#' @return `x`, invisibly.
#' @method print dkge_classification_spec
#' @export
print.dkge_classification_spec <- function(x, ...) {
  cat("<dkge_classification_spec>", "\n", sep = "")
  cat("  method :", x$method, "\n")
  cat("  targets:", if (is.null(x$targets)) "<unspecified>" else "<supplied>", "\n")
  if (!is.null(x$folds)) cat("  folds  : <custom>\n")
  invisible(x)
}
