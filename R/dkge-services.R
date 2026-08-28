# dkge-services.R
# Service abstractions used by dkge_pipeline() to orchestrate major steps.

#' Construct a contrast service
#'
#' Packages the arguments required by [dkge_contrast()] into a reusable service
#' object. The service can be passed to [dkge_pipeline()] (via the `service`
#' argument when it gains support) or executed manually with
#' `.dkge_run_contrast_service()`.
#'
#' @param method Cross-fitting strategy ("loso", "kfold", or "analytic").
#' @param ridge Ridge penalty forwarded to [dkge_contrast()].
#' @param ... Additional arguments stored with the service (e.g. `folds`,
#'   `parallel`).
#' @return Object of class `dkge_contrast_service`.
#' @export
#' @examples
#' contrast_srv <- dkge_contrast_service(method = "loso", ridge = 0)
dkge_contrast_service <- function(method = c("loso", "kfold", "analytic"),
                                  ridge = 0,
                                  ...) {
  method <- match.arg(method)
  structure(list(method = method,
                 ridge = ridge,
                 args = list(...)),
            class = "dkge_contrast_service")
}

#' Construct a legacy descriptive transport service
#'
#' Pipeline transport services are retained for descriptive migration output.
#' They cannot be composed with inference; use the typed reference-oriented
#' alignment workflow for inferential maps.
#'
#' @param spec Transport specification (list or `dkge_transport_spec`).
#' @param ... Additional key-value pairs merged into the specification.
#' @return Object of class `dkge_transport_service`.
#' @export
#' @examples
#' transport_srv <- dkge_transport_service(dkge_transport_spec(centroids = list(matrix(0, 2, 3))))
dkge_transport_service <- function(spec = NULL, ...) {
  if (inherits(spec, "dkge_transport_spec")) {
    spec <- unclass(spec)
  }
  extra <- list(...)
  if (length(extra)) {
    spec <- c(spec, extra)
  }
  structure(list(spec = spec), class = "dkge_transport_service")
}

#' Construct an inference service
#'
#' @param spec Inference specification (list or `dkge_inference_spec`).
#' @param ... Additional parameters merged into the specification.
#' @return Object of class `dkge_inference_service`.
#' @export
#' @examples
#' inference_srv <- dkge_inference_service(dkge_inference_spec(B = 1000))
dkge_inference_service <- function(spec = NULL, ...) {
  if (inherits(spec, "dkge_inference_spec")) {
    spec <- unclass(spec)
  }
  extra <- list(...)
  if (length(extra)) spec <- c(spec, extra)
  structure(list(spec = spec), class = "dkge_inference_service")
}

#' Execute a contrast service
#'
#' @param service Contrast service produced by [dkge_contrast_service()].
#' @param fit dkge object.
#' @param contrasts Contrast specification forwarded to [dkge_contrast()].
#' @param extra_args Additional arguments to append when invoking
#'   [dkge_contrast()].
#' @keywords internal
#' @noRd
.dkge_run_contrast_service <- function(service, fit, contrasts, extra_args = list()) {
  if (is.null(service)) {
    stop("Contrast service is NULL; supply dkge_contrast_service().", call. = FALSE)
  }
  args <- c(list(fit = fit,
                 contrasts = contrasts,
                 method = service$method,
                 ridge = service$ridge),
            service$args,
            extra_args)
  do.call(dkge_contrast, args)
}

#' Execute a transport service
#'
#' @param service Transport service produced by [dkge_transport_service()].
#' @param fit dkge object.
#' @param contrast_results Result from [dkge_contrast()].
#' @keywords internal
#' @noRd
.dkge_run_transport_service <- function(service, fit, contrast_results) {
  if (is.null(service) || is.null(service$spec)) {
    return(NULL)
  }
  spec <- service$spec
  if (is.null(spec$centroids)) {
    return(NULL)
  }
  medoid <- spec$medoid %||% 1L
  mapper_spec <- spec$mapper %||% NULL
  method_arg <- spec$method %||% "sinkhorn"
  alignment_mode <- spec$alignment_mode %||% "fold_safe"
  if (identical(alignment_mode, "fold_safe") &&
      (!is.null(spec$loadings) || !is.null(spec$betas))) {
    warning(
      "Transport-service loadings/betas are ignored in fold-safe mode; contrast receipts own alignment features.",
      call. = FALSE
    )
  }
  mapper_args <- spec[intersect(names(spec),
                                c("epsilon", "max_iter", "tol",
                                  "lambda_emb", "lambda_spa",
                                  "sigma_mm", "lambda_size",
                                  "value_type", "warm_start"))]
  args <- list(
    fit = fit,
    contrast_obj = contrast_results,
    medoid = medoid,
    centroids = spec$centroids,
    loadings = if (identical(alignment_mode, "descriptive")) spec$loadings else NULL,
    betas = if (identical(alignment_mode, "descriptive")) spec$betas else NULL,
    sizes = spec$sizes,
    mapper = mapper_spec,
    method = method_arg,
    alignment_mode = alignment_mode,
    transport_cache = spec$transport_cache %||% spec$fitted_alignment
  )
  args <- c(args, mapper_args)
  args <- args[!vapply(args, is.null, logical(1))]
  do.call(.dkge_transport_contrasts_to_reference_core, args)
}

#' Execute an inference service
#'
#' @param service Inference service produced by [dkge_inference_service()].
#' @param contrast_results Result from [dkge_contrast()].
#' @param transport_results Optional transported contrasts.
#' @keywords internal
#' @noRd
.dkge_run_inference_service <- function(service, contrast_results, transport_results = NULL) {
  if (is.null(service) || is.null(service$spec)) {
    return(NULL)
  }
  if (!is.null(transport_results)) {
    .dkge_abort(
      paste0(
        "Pipeline transport is legacy descriptive output and cannot enter ",
        "inference. Use `dkge_transport_contrasts_to_reference()` followed by ",
        "`dkge_infer_aligned()`."
      ),
      "dkge_alignment_ineligible_error"
    )
  }
  spec <- service$spec
  B <- spec$B %||% 2000L
  tail <- spec$tail %||% "two.sided"
  center <- spec$center %||% "mean"
  alpha <- spec$alpha %||% 0.05
  allow_approximate <- isTRUE(spec$allow_approximate_alignment)
  estimator_eligibility <- .dkge_rank_truncated_estimator_eligibility(
    contrast_results$method
  )
  .dkge_assert_alignment_eligible(
    estimator_eligibility,
    allow_approximate = allow_approximate
  )
  mapped_values <- NULL
  res <- .infer_signflip(
    contrast_results,
    n_perm = B,
    correction = "maxT",
    mapped_values = mapped_values,
    tail = tail,
    center = center
  )
  res$significant <- lapply(res$p_adjusted, function(p) p <= alpha)
  res$alpha <- alpha
  res$metadata$estimator <- list(
    status = estimator_eligibility$status,
    reason = estimator_eligibility$reason,
    method = contrast_results$method,
    approximate_override = allow_approximate
  )
  structure(res, class = "dkge_inference")
}

#' Run a dkge service object
#'
#' Generic dispatcher for executing pre-configured service objects. Dispatches
#' to the appropriate typed runner based on the service class.
#'
#' @param service A service object created by `dkge_contrast_service()`,
#'   `dkge_inference_service()`, or `dkge_transport_service()`
#' @param ... Additional arguments passed to the typed runner
#' @return Result from the typed service runner
#' @export
dkge_run_service <- function(service, ...) {
  UseMethod("dkge_run_service")
}

#' @export
#' @rdname dkge_run_service
dkge_run_service.dkge_contrast_service <- function(service, ...) {
  .dkge_run_contrast_service(service, ...)
}

#' @export
#' @rdname dkge_run_service
dkge_run_service.dkge_inference_service <- function(service, ...) {
  .dkge_run_inference_service(service, ...)
}

#' @export
#' @rdname dkge_run_service
dkge_run_service.dkge_transport_service <- function(service, ...) {
  .dkge_run_transport_service(service, ...)
}

#' @export
#' @rdname dkge_run_service
dkge_run_service.default <- function(service, ...) {
  stop(
    "No dkge_run_service method for class '",
    paste(class(service), collapse = "/"), "'.",
    call. = FALSE
  )
}
