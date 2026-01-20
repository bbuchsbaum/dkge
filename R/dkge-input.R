# dkge-input.R
# Canonical input descriptors for DKGE initialisation.

#' Anchor-based DKGE input descriptor
#'
#' Builds an input specification that projects subject item kernels onto a
#' shared anchor basis before fitting DKGE. Use this object with
#' [dkge_fit_from_input()] or supply it to [dkge_pipeline()] via the `input`
#' argument. The descriptor is immutable; downstream calls may extend its
#' `dkge_args` field with additional DKGE fitting options.
#'
#' @param features_list List of subject feature matrices (`n_s \times d`).
#' @param K_item_list List of subject item kernels (`n_s \times n_s` PSD).
#' @param folds Optional fold structure passed to
#'   [dkge_build_anchor_kernels()].
#' @param anchors Optional list overriding anchor selection defaults (see
#'   [dkge_anchor_fit()]).
#' @param design_kernel Optional design kernel supplied to
#'   [dkge_fit_from_kernels()]. Defaults to the identity in effect space.
#' @param dkge_args Optional list of arguments forwarded to the DKGE fitter
#'   after anchor kernels are constructed (e.g. `w_method`, `cpca_part`).
#' @return Object of class `dkge_input_anchor` (inherits from `dkge_input`).
#' @export
dkge_input_anchor <- function(features_list,
                               K_item_list,
                               folds = NULL,
                               anchors = list(),
                               design_kernel = NULL,
                               dkge_args = list()) {
  stopifnot(is.list(features_list), is.list(K_item_list),
            length(features_list) == length(K_item_list),
            is.list(anchors), is.list(dkge_args))
  structure(list(
    features_list = features_list,
    K_item_list = K_item_list,
    folds = folds,
    anchor_args = anchors,
    design_kernel = design_kernel,
    dkge_args = dkge_args
  ), class = c("dkge_input_anchor", "dkge_input"))
}

#' Fit DKGE from an input descriptor
#'
#' Dispatches to the appropriate preprocessing pipeline (anchors, raw betas,
#' etc.) before invoking the DKGE core. Currently supports anchor descriptors
#' created via [dkge_input_anchor()].
#'
#' @param input Object inheriting from `dkge_input`.
#' @param ... Additional arguments merged into the DKGE fitting call (these are
#'   interpreted as DKGE core options, e.g. `w_method`).
#' @return A fitted `dkge` object.
#' @export
dkge_fit_from_input <- function(input, ...) {
  UseMethod("dkge_fit_from_input")
}

#' @export
#' @rdname dkge_fit_from_input
dkge_fit_from_input.dkge_input_anchor <- function(input, ...) {
  extra <- list(...)
  anchor_args <- utils::modifyList(input$anchor_args, extra$anchors %||% list())
  extra$anchors <- NULL
  if (!is.null(extra$design_kernel)) {
    input$design_kernel <- extra$design_kernel
    extra$design_kernel <- NULL
  }
  if (!is.null(extra$folds)) {
    input$folds <- extra$folds
    extra$folds <- NULL
  }
  dkge_args <- utils::modifyList(input$dkge_args, extra)
  call_args <- list(
    features_list = input$features_list,
    K_item_list = input$K_item_list,
    folds = input$folds,
    anchors = anchor_args,
    design_kernel = input$design_kernel,
    dkge_args = dkge_args
  )
  call_args <- call_args[!vapply(call_args, is.null, logical(1))]
  do.call(dkge_anchor_fit, call_args)
}

#' @export
#' @rdname dkge_fit_from_input
dkge_fit_from_input.default <- function(input, ...) {
  stop("Unsupported DKGE input descriptor.", call. = FALSE)
}
