# dkge-utils.R
# Shared helper utilities for DKGE

#' Null-coalescing helper
#'
#' Returns `b` when `a` is `NULL`, otherwise returns `a`.
#'
#' @name grapes-or-or-grapes
#' @keywords internal
NULL

#' @rdname grapes-or-or-grapes
#' @param a Primary value tested for `NULL`.
#' @param b Fallback value returned when `a` is `NULL`.
#' @usage a \%||\% b
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Apply helper with optional parallelism
#'
#' Wraps `lapply()` with an optional future.apply backend so callers can enable
#' `parallel = TRUE` without repeating boilerplate dependency checks.
#'
#' @param X Vector or list to iterate over.
#' @param FUN Function to apply.
#' @param parallel Logical; if `TRUE`, uses `future.apply::future_lapply()`.
#' @param ... Additional arguments passed to the apply backend.
#' @return List of results matching `lapply()` semantics.
#' @keywords internal
.dkge_apply <- function(X, FUN, parallel = FALSE, ...) {
  if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("parallel=TRUE requires the future.apply package; install it or set parallel=FALSE.",
           call. = FALSE)
    }
    future.apply::future_lapply(X, FUN, ...)
  } else {
    lapply(X, FUN, ...)
  }
}

#' Check whether verbose output should be emitted
#'
#' Uses the per-call `verbose` flag combined with the global
#' `options(dkge.verbose = TRUE)` toggle.
#'
#' @keywords internal
.dkge_verbose <- function(verbose) {
  isTRUE(verbose) && isTRUE(getOption("dkge.verbose", TRUE))
}
