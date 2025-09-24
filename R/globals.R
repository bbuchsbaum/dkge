#' Internal global variable declarations
#'
#' Declares non-standard evaluation symbols used inside ggplot2 calls to keep
#' R CMD check quiet.
#' @name dkge-global-variables
#' @keywords internal
NULL

utils::globalVariables(c(
  "anchor",
  "angle_deg",
  "base",
  "component",
  "cumulative",
  "effect",
  "prop_var",
  "subject",
  "value",
  "weight"
))
