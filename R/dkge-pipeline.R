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
#' @param transport Either a transport specification/service or `NULL`.
#' @param inference Either an inference specification/service or `NULL`.
#' @param classification Optional specification passed to [dkge_classify()].
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
                          classification = NULL,
                          method = c("loso", "kfold", "analytic"),
                          ridge = 0,
                          ...) {
  method <- match.arg(method)
  extra_args <- list(...)

  if (inherits(transport, "dkge_transport_spec")) {
    transport <- unclass(transport)
  }
  if (inherits(inference, "dkge_inference_spec")) {
    inference <- unclass(inference)
  }
  if (inherits(classification, "dkge_classification_spec")) {
    classification <- unclass(classification)
  }

  if (is.null(fit)) {
    stopifnot(!is.null(betas), !is.null(designs), !is.null(kernel))
    fit_args <- c(list(betas, designs = designs, kernel = kernel, omega = omega),
                  extra_args)
    fit <- do.call(dkge, fit_args)
  }
  stopifnot(inherits(fit, "dkge"))

  contrast_service <- dkge_contrast_service(method = method, ridge = ridge)
  contrast_results <- .dkge_run_contrast_service(contrast_service, fit, contrasts, extra_args)

  transport_service <- if (inherits(transport, "dkge_transport_service")) {
    transport
  } else {
    dkge_transport_service(transport)
  }
  transport_results <- .dkge_run_transport_service(transport_service, fit, contrast_results)

  classification_result <- NULL
  if (!is.null(classification)) {
    if (inherits(classification, "dkge_classification")) {
      classification_result <- classification
    } else if (is.list(classification) && !is.null(classification$targets)) {
      args <- utils::modifyList(list(fit = fit), classification)
      classification_result <- do.call(dkge_classify, args)
    } else {
      classification_result <- dkge_classify(fit, classification)
    }
  }

  inference_service <- if (inherits(inference, "dkge_inference_service")) {
    inference
  } else {
    dkge_inference_service(inference)
  }
  inference_results <- .dkge_run_inference_service(inference_service,
                                                   contrast_results,
                                                   transport_results)

  list(
    fit = fit,
    diagnostics = dkge_diagnostics(fit),
    contrasts = contrast_results,
    transport = transport_results,
    inference = inference_results,
    classification = classification_result
  )
}
