#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: freeze-public-beta-baseline.R SOURCE_ROOT OUTPUT_CSV",
    call. = FALSE
  )
}

source_root <- normalizePath(args[[1]], mustWork = TRUE)
output_csv <- args[[2]]

git <- function(...) {
  out <- system2(
    "git",
    c("-C", shQuote(source_root), ...),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(out, "status") %||% 0L
  if (!identical(status, 0L)) {
    stop(paste(out, collapse = "\n"), call. = FALSE)
  }
  out
}

`%||%` <- function(x, y) if (is.null(x)) y else x

status <- git("status", "--porcelain=v1", "--untracked-files=all")
if (!length(status)) {
  manifest <- data.frame(
    status = character(), path = character(), bytes = numeric(),
    sha256 = character(), source_head = character(),
    stringsAsFactors = FALSE
  )
} else {
  paths <- substring(status, 4L)
  paths <- sub("^.* -> ", "", paths)
  absolute <- file.path(source_root, paths)
  exists <- file.exists(absolute)
  sha256 <- rep(NA_character_, length(paths))
  sha256[exists] <- vapply(
    absolute[exists],
    digest::digest,
    character(1),
    file = TRUE,
    algo = "sha256",
    serialize = FALSE
  )
  bytes <- rep(NA_real_, length(paths))
  bytes[exists] <- file.info(absolute[exists])$size
  manifest <- data.frame(
    status = substring(status, 1L, 2L),
    path = paths,
    bytes = bytes,
    sha256 = sha256,
    source_head = rep(git("rev-parse", "HEAD")[[1]], length(paths)),
    stringsAsFactors = FALSE
  )
}

dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(manifest, output_csv, row.names = FALSE, na = "")
