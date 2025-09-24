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
#' @param ... Additional fields stored on the spec (e.g., precomputed loadings
#'   or betas).
#' @return Object with class `dkge_transport_spec`.
#' @export
#' @examples
#' transport <- dkge_transport_spec(centroids = list(matrix(runif(12), 4, 3)))
dkge_transport_spec <- function(centroids,
                                sizes = NULL,
                                medoid = 1L,
                                method = c("sinkhorn", "sinkhorn_cpp", "knn"),
                                mapper = NULL,
                                epsilon = 0.05,
                                max_iter = 200L,
                                tol = 1e-6,
                                lambda_emb = 1,
                                lambda_spa = 0.5,
                                sigma_mm = 15,
                                lambda_size = 0,
                                ...) {
  method <- match.arg(method)
  stopifnot(is.list(centroids), length(centroids) >= 1L)
  if (!is.null(sizes)) {
    stopifnot(is.list(sizes), length(sizes) == length(centroids))
  }
  medoid <- as.integer(medoid)
  if (medoid < 1L || medoid > length(centroids)) {
    stop("`medoid` must index one of the provided centroid lists")
  }
  stopifnot(epsilon > 0, max_iter > 0, tol > 0, lambda_emb >= 0,
            lambda_spa >= 0, sigma_mm > 0, lambda_size >= 0)

  spec <- list(
    centroids = centroids,
    sizes = sizes,
    medoid = medoid,
    method = method,
    mapper = mapper,
    epsilon = epsilon,
    max_iter = as.integer(max_iter),
    tol = tol,
    lambda_emb = lambda_emb,
    lambda_spa = lambda_spa,
    sigma_mm = sigma_mm,
    lambda_size = lambda_size
  )
  extra <- list(...)
  if (length(extra)) spec <- c(spec, extra)
  class(spec) <- c("dkge_transport_spec", "list")
  spec
}

#' Inference specification helper
#'
#' @param B Number of permutations for sign-flip inference.
#' @param tail Tail of the test: "two.sided", "left", or "right".
#' @param center Centering method for permutations: "mean", "median", or "none".
#' @return Object with class `dkge_inference_spec`.
#' @export
#' @examples
#' infer <- dkge_inference_spec(B = 1000, tail = "two.sided")
dkge_inference_spec <- function(B = 2000L,
                                tail = c("two.sided", "left", "right"),
                                center = c("mean", "median", "none")) {
  stopifnot(B > 0)
  tail <- match.arg(tail)
  center <- match.arg(center)
  structure(list(B = as.integer(B), tail = tail, center = center),
            class = c("dkge_inference_spec", "list"))
}

#' Classification specification helper
#'
#' @param targets Target specification accepted by [dkge_classify()].
#' @param method Classifier backend ("lda" or "logit").
#' @param folds Optional fold specification.
#' @param lambda Optional ridge parameter.
#' @param metric Classification metrics to report.
#' @param mode Decoding mode passed to [dkge_classify()].
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
                                     mode = c("auto", "cell", "delta"),
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

#' @export
print.dkge_transport_spec <- function(x, ...) {
  cat("<dkge_transport_spec>", "\n", sep = "")
  cat("  method   :", x$method, "\n")
  cat("  subjects :", length(x$centroids), "\n")
  cat("  medoid   :", x$medoid, "\n")
  if (!is.null(x$epsilon)) cat("  epsilon  :", x$epsilon, "\n")
  invisible(x)
}

#' @export
print.dkge_inference_spec <- function(x, ...) {
  cat("<dkge_inference_spec>", "\n", sep = "")
  cat("  permutations :", x$B, "\n")
  cat("  tail         :", x$tail, "\n")
  cat("  center       :", x$center, "\n")
  invisible(x)
}

#' @export
print.dkge_classification_spec <- function(x, ...) {
  cat("<dkge_classification_spec>", "\n", sep = "")
  cat("  method :", x$method, "\n")
  cat("  targets:", if (is.null(x$targets)) "<unspecified>" else "<supplied>", "\n")
  if (!is.null(x$folds)) cat("  folds  : <custom>\n")
  invisible(x)
}
