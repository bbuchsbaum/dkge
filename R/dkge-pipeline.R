# dkge-pipeline.R
# High-level orchestration for DKGE analyses.

#' End-to-end DKGE workflow
#'
#' Fits DKGE (if needed), computes cross-fitted contrasts, optionally transports
#' them to a medoid parcellation, and performs sign-flip inference.
#'
#' @param fit Optional pre-computed `dkge` object. If `NULL`, provide `betas`,
#'   `designs`, and `kernel` to fit inside the pipeline.
#' @param betas,designs,kernel Inputs passed to [dkge()] when `fit` is `NULL`.
#' @param omega Optional spatial weights forwarded to [dkge()].
#' @param contrasts Contrast specification as accepted by [dkge_contrast()].
#' @param transport Optional list configuring transport. Recognised fields:
#'   `method` ("sinkhorn" or "sinkhorn_cpp"), `centroids`, `sizes`, `medoid`,
#'   `lambda_emb`, `lambda_spa`, `sigma_mm`, `epsilon`, `max_iter`, `tol`.
#'   Alternatively, supply precomputed transport matrices via
#'   `transforms`/`matrices` to bypass the medoid transport helper.
#' @param inference Optional list configuring sign-flip inference with entries
#'   `B` (permutations), `tail`, and `center`.
#' @param method Cross-fitting strategy for contrasts (default "loso").
#' @param ridge Optional ridge added during held-out decompositions.
#' @param ... Additional arguments passed to [dkge()] when fitting inside the
#'   pipeline, or to [dkge_contrast()].
#' @return List containing the fit, diagnostics, raw contrast values, transported
#'   maps (if requested), and inference results.
#' @export
dkge_pipeline <- function(fit = NULL,
                          betas = NULL, designs = NULL, kernel = NULL, omega = NULL,
                          contrasts,
                          transport = NULL,
                          inference = list(),
                          method = c("loso", "kfold", "analytic"),
                          ridge = 0,
                          ...) {
  method <- match.arg(method)

  if (is.null(fit)) {
    stopifnot(!is.null(betas), !is.null(designs), !is.null(kernel))
    fit <- dkge(betas, designs = designs, kernel = kernel, omega = omega, ...)
  }
  stopifnot(inherits(fit, "dkge"))

  contrast_results <- dkge_contrast(fit, contrasts, method = method, ridge = ridge, ...)

  transport_results <- NULL
  if (!is.null(transport) && !is.null(transport$centroids)) {
    medoid <- transport$medoid %||% 1L
    mapper_spec <- transport$mapper %||% NULL
    method_arg <- transport$method %||% "sinkhorn"
    mapper_args <- transport[intersect(names(transport),
                                       c("epsilon", "max_iter", "tol",
                                         "lambda_emb", "lambda_spa",
                                         "sigma_mm", "lambda_size"))]
    args <- list(
      fit = fit,
      contrast_obj = contrast_results,
      medoid = medoid,
      centroids = transport$centroids,
      loadings = transport$loadings,
      betas = transport$betas,
      sizes = transport$sizes,
      mapper = mapper_spec,
      method = method_arg
    )
    args <- c(args, mapper_args)
    args <- args[!vapply(args, is.null, logical(1))]
    transport_results <- do.call(dkge_transport_contrasts_to_medoid, args)
  }

  inference_results <- NULL
  if (!is.null(inference)) {
    B <- inference$B %||% 2000
    tail <- inference$tail %||% "two.sided"
    center <- inference$center %||% "mean"

    inference_results <- lapply(seq_along(contrast_results$values), function(i) {
      subj_mat <- if (!is.null(transport_results)) {
        transport_results[[i]]$subj_values
      } else {
        as.matrix(contrast_results, contrast = i)
      }
      dkge_signflip_maxT(subj_mat, B = B, tail = tail, center = center)
    })
    names(inference_results) <- names(contrast_results$values)
  }

  list(
    fit = fit,
    diagnostics = dkge_diagnostics(fit),
    contrasts = contrast_results,
    transport = transport_results,
    inference = inference_results
  )
}
