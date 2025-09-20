# dkge-utils.R
# Shared helper utilities for DKGE

#' Null-coalescing helper
#'
#' Returns `b` when `a` is `NULL`, otherwise returns `a`.
#'
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a
