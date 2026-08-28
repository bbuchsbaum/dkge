# dkge-pipeline.R
# High-level orchestration for DKGE analyses.

#' End-to-end DKGE workflow
#'
#' Fits DKGE (if needed), computes cross-fitted contrasts, optionally produces
#' legacy descriptive transport output, and performs native-support sign-flip
#' inference. Pipeline transport and inference cannot be composed. Functional
#' alignment inference uses [dkge_transport_contrasts_to_reference()] followed
#' by [dkge_infer_aligned()].
#'
#' @param fit Optional pre-computed `dkge` object. If `NULL`, provide `betas`,
#'   `designs`, and `kernel` to fit inside the pipeline.
#' @param input Optional DKGE input descriptor created with
#'   [dkge_input_anchor()] or future helpers. When supplied (and `fit` is
#'   `NULL`), `dkge_pipeline()` will build the fit via [dkge_fit_from_input()].
#' @param betas,designs,kernel Inputs passed to [dkge()] when neither `fit` nor
#'   `input` is supplied.
#' @param omega Optional spatial weights forwarded to [dkge()].
#' @param spatial Optional model-level [dkge_spatial_regularizer()] forwarded
#'   only to the raw-beta fitting stage. Current anchor input descriptors do not
#'   expose a physical beta-column domain and therefore reject this argument.
#' @param contrasts Contrast specification as accepted by [dkge_contrast()].
#' @param transport Either a legacy descriptive transport specification/service
#'   or `NULL`. It cannot be combined with `inference`.
#' @param inference Either an inference specification/service or `NULL` (the
#'   default). Same-data rank-truncated inference is approximate and requires an
#'   explicit `allow_approximate_alignment = TRUE` in the inference spec.
#' @param classification Optional specification passed to [dkge_classify()].
#' @param method Cross-fitting strategy for contrasts (default "loso").
#' @param ridge Optional ridge added during held-out decompositions.
#' @param ... Additional arguments passed to [dkge()] when fitting inside the
#'   pipeline, or to [dkge_contrast()].
#' @return List containing the fit, diagnostics, raw contrast values, optional
#'   legacy descriptive maps, and optional native-support inference results.
#' @examples
#' # Simulate toy data
#' toy <- dkge_sim_toy(
#'   factors = list(A = list(L = 2), B = list(L = 3)),
#'   active_terms = c("A", "B"), S = 5, P = 25, snr = 5
#' )
#'
#' # Run pipeline with LOSO contrasts
#' result <- dkge_pipeline(
#'   betas = toy$B_list,
#'   designs = toy$X_list,
#'   kernel = toy$K,
#'   contrasts = c(1, rep(0, 4)),  # first effect
#'   method = "loso"
#' )
#' names(result)
#' @export
dkge_pipeline <- function(fit = NULL,
                          input = NULL,
                          betas = NULL, designs = NULL, kernel = NULL, omega = NULL,
                          spatial = NULL,
                          contrasts,
                          transport = NULL,
                          inference = NULL,
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
  if (!is.null(transport) && !is.null(inference)) {
    .dkge_abort(
      paste0(
        "`dkge_pipeline()` cannot compose its legacy descriptive transport ",
        "with inference. Use `dkge_transport_contrasts_to_reference()` ",
        "followed by `dkge_infer_aligned()`."
      ),
      "dkge_alignment_ineligible_error"
    )
  }

  if (is.null(fit)) {
    if (!is.null(input)) {
      if (!inherits(input, "dkge_input")) {
        stop("`input` must be constructed with dkge_input_* helpers.", call. = FALSE)
      }
      fit_call <- c(list(input = input), extra_args)
      if (!is.null(spatial)) fit_call$spatial <- spatial
      fit <- do.call(dkge_fit_from_input, fit_call)
    } else {
      stopifnot(!is.null(betas), !is.null(designs), !is.null(kernel))
      fit_args <- c(list(betas, designs = designs, K = kernel,
                         Omega_list = omega, spatial = spatial),
                    extra_args)
      fit <- do.call(dkge, fit_args)
    }
  } else if (!is.null(spatial)) {
    .dkge_abort(
      "`spatial` cannot be added to a pre-computed fit; refit with the regularizer.",
      "dkge_spatial_spec_error"
    )
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
