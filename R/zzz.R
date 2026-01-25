# zzz.R

#' @useDynLib dkge, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

.onLoad <- function(libname, pkgname) {
  # Optional interoperability: register dkge aligner into neuralign if present.
  if (requireNamespace("neuralign", quietly = TRUE)) {
    try(.dkge_register_neuralign_aligner(), silent = TRUE)
  }
}
