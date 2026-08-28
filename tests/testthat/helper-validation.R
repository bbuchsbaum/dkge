# Locate a validation module under both layouts: the source tree used by
# devtools::test(), where the file is still under inst/, and the installed
# layout used by R CMD check, where inst/ contents sit at the package root.
# Returns NA_character_ when neither resolves, so callers can skip rather than
# error at file-source time.
dkge_validation_path <- function(...) {
  candidates <- c(
    file.path(testthat::test_path("..", ".."), "inst", "validation", ...),
    system.file("validation", ..., package = "dkge")
  )
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) NA_character_ else hit[[1]]
}
