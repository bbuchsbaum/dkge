#!/usr/bin/env Rscript

# Build the certification site from an exact package tarball and leave a
# machine-verifiable receipt in the generated site.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "inst/validation/functional-alignment-certification/run-pkgdown.R",
    mustWork = TRUE
  )
}
root <- normalizePath(file.path(dirname(script), "..", "..", ".."),
                      mustWork = TRUE)
utils_script <- file.path(dirname(script), "evidence-utils.R")
provenance_script <- file.path(dirname(script), "package-provenance.R")
court_script <- file.path(
  root, "inst", "validation", "functional-alignment-court", "court.R"
)
script_env <- environment()
source(utils_script, local = TRUE)
builder_paths <- c(script, utils_script, court_script, provenance_script)
input_snapshot <- function() {
  list(
    source_tree_sha256 = dkfa_candidate_source_tree_hash(root),
    docs = dkfa_pkgdown_input_manifest(root),
    builder = dkfa_bundle_manifest(builder_paths, root)
  )
}
loaded_binding <- dkfa_bind_loaded_inputs(
  input_snapshot,
  function() {
    source(utils_script, local = script_env)
    source(court_script, local = script_env)
    source(provenance_script, local = script_env)
    invisible(NULL)
  },
  "Pkgdown builder loader"
)

locale_info <- l10n_info()
if (!isTRUE(locale_info[["UTF-8"]])) {
  stop(
    paste0(
      "Certified pkgdown rendering requires a UTF-8 locale; set LC_ALL and ",
      "LANG to a valid UTF-8 locale before invoking this wrapper."
    ),
    call. = FALSE
  )
}
build_locale <- Sys.getlocale()

parse_args <- function(x) {
  x <- x[grepl("^--[^=]+=", x)]
  values <- sub("^--", "", x)
  stats::setNames(
    as.list(sub("^[^=]+=", "", values)),
    sub("=.*$", "", values)
  )
}
cli <- parse_args(commandArgs(trailingOnly = TRUE))
if (!all(c("tarball", "dest-dir") %in% names(cli))) {
  stop("Usage: run-pkgdown.R --tarball=<path> --dest-dir=<path>", call. = FALSE)
}
tarball <- normalizePath(cli[["tarball"]], mustWork = TRUE)
if (dir.exists(tarball)) stop("`tarball` must be a file.", call. = FALSE)
tarball_sha256 <- dkfa_hash_file(tarball)
dest_dir <- normalizePath(cli[["dest-dir"]], mustWork = FALSE)
if (dir.exists(dest_dir) && length(list.files(
  dest_dir, all.files = TRUE, no.. = TRUE
))) {
  stop("Pkgdown destination must be absent or empty: ", dest_dir, call. = FALSE)
}
dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

extract_dir <- tempfile("dkge-pkgdown-tarball-")
dir.create(extract_dir)
on.exit(unlink(extract_dir, recursive = TRUE), add = TRUE)
utils::untar(tarball, exdir = extract_dir)
package_roots <- list.dirs(extract_dir, recursive = FALSE, full.names = TRUE)
if (length(package_roots) != 1L) {
  stop("Tarball must contain exactly one package root.", call. = FALSE)
}
package_root <- package_roots[[1L]]
source_provenance <- dkfa_verify_built_source(root, package_root)
if (!isTRUE(source_provenance$passed)) {
  stop("Tarball source projection does not match the certified source tree.",
       call. = FALSE)
}
docs <- dkfa_apply_docs_overlay(root, package_root)
source_hash <- source_provenance$certified_source_tree_sha256
if (!identical(source_hash, loaded_binding$snapshot$source_tree_sha256)) {
  stop("Tarball projection is not bound to the loaded source snapshot.",
       call. = FALSE)
}
builder <- loaded_binding$snapshot$builder

required_files <- c(
  "index.html",
  "articles/dkge-functional-alignment.html",
  "reference/index.html",
  "reference/dkge_prepare_alignment.html",
  "reference/dkge_fit_functional_template.html",
  "reference/dkge_infer_aligned.html",
  "reference/dkge_render_aligned.html"
)
pkgdown::build_site(
  package_root,
  override = list(destination = dest_dir),
  new_process = FALSE,
  install = TRUE,
  quiet = FALSE
)
missing <- required_files[!file.exists(file.path(dest_dir, required_files))]
if (length(missing)) {
  stop("Pkgdown omitted required page(s): ", paste(missing, collapse = ", "),
       call. = FALSE)
}
dkfa_assert_input_snapshot(
  loaded_binding$snapshot, input_snapshot(), "Pkgdown build"
)
if (!identical(dkfa_hash_file(tarball), tarball_sha256)) {
  stop("The package tarball changed during the pkgdown build.", call. = FALSE)
}
site <- dkfa_site_manifest(dest_dir)
receipt <- list(
  schema_version = "1.4.0",
  tarball = basename(tarball),
  tarball_sha256 = tarball_sha256,
  tarball_source_tree_sha256 = source_hash,
  tarball_raw_source_tree_sha256 =
    source_provenance$built_raw_source_tree_sha256,
  package_source_projection_sha256 =
    source_provenance$package_projection_sha256,
  tarball_shipped_payload_tree_sha256 =
    source_provenance$shipped_payload_tree_sha256,
  description_semantic_sha256 =
    source_provenance$description$semantic_sha256,
  source_projection_checks = as.list(source_provenance$checks),
  docs_input_sha256 = docs$certified_docs_input_sha256,
  packaged_docs_tree_sha256 = docs$packaged_docs_tree_sha256,
  docs_overlay_bundle_sha256 = docs$overlay_bundle_sha256,
  docs_projection_checks = as.list(docs$checks),
  site_content_sha256 = site$tree_sha256,
  required_files = as.list(required_files),
  builder_bundle_sha256 = builder$bundle_sha256,
  builder_files = as.list(builder$files),
  loaded_input_snapshot_sha256 = loaded_binding$snapshot_sha256,
  utf8_locale = TRUE,
  locale = build_locale,
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  R = R.version.string,
  pkgdown = as.character(utils::packageVersion("pkgdown"))
)
jsonlite::write_json(
  receipt,
  file.path(dest_dir, ".dkge-pkgdown-receipt.json"),
  pretty = TRUE, auto_unbox = TRUE
)
cat(dest_dir, "\n")
