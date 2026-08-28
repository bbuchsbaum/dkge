#!/usr/bin/env Rscript

# Provenance helpers for comparing an R CMD build tarball with its live source.
# R CMD build rewrites DESCRIPTION and normally omits pkgdown-only inputs, so a
# raw tree hash cannot establish this relationship. These helpers verify the
# build projection explicitly and bind the documentation overlay separately.

dkfa_normalize_dcf_value <- function(x) {
  trimws(gsub("[[:space:]]+", " ", as.character(x)))
}

dkfa_raw_court_source_tree_hash <- function(root) {
  root <- normalizePath(root, mustWork = TRUE)
  roots <- c("DESCRIPTION", "NAMESPACE", "R", "src")
  paths <- unlist(lapply(roots, function(path) {
    absolute <- file.path(root, path)
    if (dir.exists(absolute)) {
      list.files(absolute, recursive = TRUE, full.names = TRUE,
                 all.files = TRUE, no.. = TRUE)
    } else if (file.exists(absolute)) {
      absolute
    } else {
      character()
    }
  }), use.names = FALSE)
  paths <- paths[file.info(paths)$isdir %in% FALSE]
  paths <- paths[!grepl("\\.(o|so|dll|dylib)$", paths,
                        ignore.case = TRUE)]
  paths <- sort(paths, method = "radix")
  relative <- substring(paths, nchar(root) + 2L)
  records <- paste(
    relative,
    vapply(paths, dkfa_hash_file, character(1)),
    sep = "="
  )
  digest::digest(paste(records, collapse = "\n"), algo = "sha256",
                 serialize = FALSE)
}

dkfa_certified_source_tree_hash <- function(root) {
  dkfa_raw_court_source_tree_hash(root)
}

dkfa_test_validation_manifest <- function(root) {
  dkfa_tree_manifest(
    root,
    roots = c(
      "DESCRIPTION", "NAMESPACE", "R", "src", "tests", "inst", "man",
      "vignettes", "README.md",
      "data-raw/dkge-between-rotation-plan.md",
      "data-raw/dkge-between-rotation-report.md",
      "dev/calibrate-dkge-between-rotation.R"
    ),
    exclude = c(
      "\\.(o|so|dll|dylib)$", "(^|/)src/symbols\\.rds$",
      "(^|/)\\.DS_Store$"
    )
  )
}

dkfa_description_projection <- function(live_path, built_path) {
  live <- read.dcf(normalizePath(live_path, mustWork = TRUE))
  built <- read.dcf(normalizePath(built_path, mustWork = TRUE))
  if (nrow(live) != 1L || nrow(built) != 1L) {
    stop("DESCRIPTION projection requires one DCF record per file.",
         call. = FALSE)
  }

  live_fields <- colnames(live)
  built_fields <- colnames(built)
  allowed_derived <- c(
    "NeedsCompilation", "Packaged", "Author", "Maintainer", "Built",
    "Repository", "Date/Publication"
  )
  missing <- setdiff(live_fields, built_fields)
  unexpected <- setdiff(setdiff(built_fields, live_fields), allowed_derived)
  common <- intersect(live_fields, built_fields)
  live_values <- stats::setNames(
    vapply(live[1L, common, drop = TRUE], dkfa_normalize_dcf_value,
           character(1)),
    common
  )
  built_values <- stats::setNames(
    vapply(built[1L, common, drop = TRUE], dkfa_normalize_dcf_value,
           character(1)),
    common
  )
  changed <- common[live_values != built_values]
  canonical_fields <- sort(live_fields, method = "radix")
  canonical_values <- stats::setNames(
    vapply(live[1L, canonical_fields, drop = TRUE],
           dkfa_normalize_dcf_value, character(1)),
    canonical_fields
  )
  semantic_hash <- digest::digest(
    paste(names(canonical_values), canonical_values, sep = "=",
          collapse = "\n"),
    algo = "sha256", serialize = FALSE
  )
  checks <- c(
    all_live_fields_retained = !length(missing),
    no_unexpected_built_fields = !length(unexpected),
    live_field_values_unchanged = !length(changed)
  )
  list(
    passed = all(checks),
    checks = checks,
    missing_fields = missing,
    unexpected_fields = unexpected,
    changed_fields = changed,
    semantic_sha256 = semantic_hash
  )
}

dkfa_verify_built_source <- function(live_root, built_root,
                                     certified_source_tree_sha256 = NULL) {
  live_root <- normalizePath(live_root, mustWork = TRUE)
  built_root <- normalizePath(built_root, mustWork = TRUE)
  live_source_tree_sha256 <- dkfa_certified_source_tree_hash(live_root)
  requested_source_tree_sha256 <- certified_source_tree_sha256
  if (is.null(requested_source_tree_sha256)) {
    requested_source_tree_sha256 <- live_source_tree_sha256
  }
  description <- dkfa_description_projection(
    file.path(live_root, "DESCRIPTION"),
    file.path(built_root, "DESCRIPTION")
  )
  excludes <- c(
    "\\.(o|so|dll|dylib)$", "(^|/)src/symbols\\.rds$",
    "(^|/)\\.DS_Store$"
  )
  source_roots <- c("NAMESPACE", "R", "src")
  live_files <- dkfa_tree_manifest(live_root, source_roots, excludes)
  built_files <- dkfa_tree_manifest(built_root, source_roots, excludes)
  file_names_match <- identical(names(live_files$files),
                                names(built_files$files))
  file_bytes_match <- file_names_match &&
    identical(live_files$files, built_files$files)
  shipped_payload_roots <- c(
    "tests", "inst", "LICENSE", "NEWS.md", "CONTRIBUTING.md"
  )
  live_shipped_payload <- dkfa_tree_manifest(
    live_root, shipped_payload_roots, excludes
  )
  built_shipped_payload <- dkfa_tree_manifest(
    built_root, shipped_payload_roots, excludes
  )
  shipped_payload_names_match <- identical(
    names(live_shipped_payload$files), names(built_shipped_payload$files)
  )
  shipped_payload_bytes_match <- shipped_payload_names_match &&
    identical(live_shipped_payload$files, built_shipped_payload$files)
  package_payload_roots <- c(
    "NAMESPACE", "R", "src", "tests", "inst", "LICENSE", "NEWS.md",
    "CONTRIBUTING.md", "man", "vignettes", "README.md"
  )
  live_package_payload <- dkfa_tree_manifest(
    live_root, package_payload_roots, excludes
  )
  built_package_payload <- dkfa_tree_manifest(
    built_root, ".", c(excludes, "(^|/)DESCRIPTION$")
  )
  package_payload_names_match <- identical(
    names(live_package_payload$files), names(built_package_payload$files)
  )
  package_payload_bytes_match <- package_payload_names_match &&
    identical(live_package_payload$files, built_package_payload$files)
  checks <- c(
    certified_source_tree = identical(
      live_source_tree_sha256, requested_source_tree_sha256
    ),
    description_projection = isTRUE(description$passed),
    source_file_names = file_names_match,
    source_file_bytes = file_bytes_match,
    shipped_payload_file_names = shipped_payload_names_match,
    shipped_payload_file_bytes = shipped_payload_bytes_match,
    closed_world_package_file_names = package_payload_names_match,
    closed_world_package_file_bytes = package_payload_bytes_match
  )
  projection_hash <- digest::digest(
    paste(
      c(
        paste0("description=", description$semantic_sha256),
        paste0("package_payload=", live_package_payload$tree_sha256),
        paste0("source_files=", live_files$tree_sha256),
        paste0("shipped_payload=", live_shipped_payload$tree_sha256)
      ),
      collapse = "\n"
    ),
    algo = "sha256", serialize = FALSE
  )
  list(
    passed = all(checks),
    checks = checks,
    certified_source_tree_sha256 = live_source_tree_sha256,
    requested_source_tree_sha256 = requested_source_tree_sha256,
    built_raw_source_tree_sha256 = dkfa_raw_court_source_tree_hash(built_root),
    source_files_tree_sha256 = live_files$tree_sha256,
    shipped_payload_tree_sha256 = live_shipped_payload$tree_sha256,
    built_shipped_payload_tree_sha256 = built_shipped_payload$tree_sha256,
    package_payload_tree_sha256 = live_package_payload$tree_sha256,
    built_package_payload_tree_sha256 = built_package_payload$tree_sha256,
    package_projection_sha256 = projection_hash,
    description = description,
    live_only_files = setdiff(names(live_files$files),
                              names(built_files$files)),
    built_only_files = setdiff(names(built_files$files),
                               names(live_files$files)),
    live_only_shipped_payload_files = setdiff(
      names(live_shipped_payload$files), names(built_shipped_payload$files)
    ),
    built_only_shipped_payload_files = setdiff(
      names(built_shipped_payload$files), names(live_shipped_payload$files)
    ),
    live_only_package_payload_files = setdiff(
      names(live_package_payload$files), names(built_package_payload$files)
    ),
    built_only_package_payload_files = setdiff(
      names(built_package_payload$files), names(live_package_payload$files)
    )
  )
}

dkfa_docs_projection <- function(live_root, built_root) {
  live_root <- normalizePath(live_root, mustWork = TRUE)
  built_root <- normalizePath(built_root, mustWork = TRUE)
  packaged_roots <- c("man", "vignettes", "README.md")
  live_packaged <- dkfa_tree_manifest(
    live_root, packaged_roots, exclude = c("(^|/)\\.DS_Store$")
  )
  built_packaged <- dkfa_tree_manifest(
    built_root, packaged_roots, exclude = c("(^|/)\\.DS_Store$")
  )
  names_match <- identical(names(live_packaged$files),
                           names(built_packaged$files))
  bytes_match <- names_match && identical(live_packaged$files,
                                          built_packaged$files)
  overlay_roots <- c("DESCRIPTION", "_pkgdown.yml", "pkgdown")
  overlay <- dkfa_tree_manifest(
    live_root, overlay_roots, exclude = c("(^|/)\\.DS_Store$")
  )
  live_docs <- dkfa_pkgdown_input_manifest(live_root)
  checks <- c(
    packaged_doc_names = names_match,
    packaged_doc_bytes = bytes_match,
    site_config_present = file.exists(file.path(live_root, "_pkgdown.yml"))
  )
  list(
    passed = all(checks),
    checks = checks,
    packaged_docs_tree_sha256 = live_packaged$tree_sha256,
    overlay_bundle_sha256 = overlay$tree_sha256,
    overlay_files = overlay$files,
    certified_docs_input_sha256 = live_docs$tree_sha256,
    live_only_files = setdiff(names(live_packaged$files),
                              names(built_packaged$files)),
    built_only_files = setdiff(names(built_packaged$files),
                               names(live_packaged$files))
  )
}

dkfa_apply_docs_overlay <- function(live_root, built_root) {
  live_root <- normalizePath(live_root, mustWork = TRUE)
  built_root <- normalizePath(built_root, mustWork = TRUE)
  provenance <- dkfa_docs_projection(live_root, built_root)
  if (!isTRUE(provenance$passed)) {
    stop("Built package documentation does not match the live source.",
         call. = FALSE)
  }
  overlay_paths <- names(provenance$overlay_files)
  for (relative in overlay_paths) {
    source_path <- file.path(live_root, relative)
    destination <- file.path(built_root, relative)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(source_path, destination, overwrite = TRUE,
                   copy.mode = TRUE, copy.date = TRUE)) {
      stop("Could not apply documentation overlay: ", relative,
           call. = FALSE)
    }
  }
  overlaid <- dkfa_pkgdown_input_manifest(built_root)
  if (!identical(overlaid$tree_sha256,
                 provenance$certified_docs_input_sha256)) {
    stop("Documentation overlay did not reconstruct the certified inputs.",
         call. = FALSE)
  }
  provenance$overlaid_docs_input_sha256 <- overlaid$tree_sha256
  provenance
}

dkfa_validate_r_cmd_check_log <- function(lines) {
  option_lines <- grep("^\\* using option(s)? ", lines, value = TRUE)
  option_record <- paste(option_lines, collapse = " ")
  required <- c("--no-manual", "--ignore-vignettes")
  allowed <- c(required, "--no-clean")
  # R may print the option list with UTF-8 typographic quotes even in a C
  # locale. Normalize their byte sequences without asking the regex engine to
  # interpret locale-dependent Unicode characters.
  quote_bytes <- lapply(
    c(0x98L, 0x99L, 0x9cL, 0x9dL),
    function(last) rawToChar(as.raw(c(0xe2L, 0x80L, last)))
  )
  for (quote in quote_bytes) {
    option_record <- gsub(
      quote, "'", option_record, fixed = TRUE, useBytes = TRUE
    )
  }
  matches <- gregexpr(
    "--[^[:space:]'\"`]+", option_record, perl = TRUE, useBytes = TRUE
  )
  options <- unique(regmatches(option_record, matches)[[1L]])
  if (identical(options, character(0)) || identical(options, "")) {
    options <- character()
  }
  contains <- function(option) grepl(option, option_record, fixed = TRUE)
  unexpected <- setdiff(options, allowed)
  checks <- c(
    status_ok = any(grepl("^Status: OK$", trimws(lines))),
    options_recorded = length(option_lines) > 0L,
    no_manual = contains("--no-manual"),
    ignore_vignettes = contains("--ignore-vignettes"),
    no_unexpected_options = length(unexpected) == 0L
  )
  list(
    passed = all(checks), checks = checks,
    option_record = option_record,
    recorded_options = options,
    required_options = required,
    allowed_options = allowed,
    unexpected_options = unexpected
  )
}

dkfa_validate_source_test_evidence <- function(
    manifest_path, results_path, session_info_path) {
  paths <- c(
    manifest = manifest_path,
    results = results_path,
    session_info = session_info_path
  )
  missing <- names(paths)[!file.exists(paths)]
  if (length(missing)) {
    stop(
      "Missing source-test evidence: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  manifest <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
  results <- utils::read.csv(results_path, stringsAsFactors = FALSE)
  required_columns <- c(
    "nb", "failed", "error", "warning", "skipped"
  )
  absent_columns <- setdiff(required_columns, names(results))
  if (length(absent_columns)) {
    stop(
      "Source-test results omit required columns: ",
      paste(absent_columns, collapse = ", "),
      call. = FALSE
    )
  }
  count_column <- function(name) {
    value <- suppressWarnings(as.numeric(results[[name]]))
    if (anyNA(value) || any(!is.finite(value)) || any(value < 0)) {
      stop("Source-test result column `", name,
           "` is not a finite non-negative count.", call. = FALSE)
    }
    sum(value)
  }
  expectation_failures <- count_column("failed")
  errors <- count_column("error")
  counts <- c(
    test_blocks = nrow(results),
    expectations = count_column("nb"),
    expectation_failures = expectation_failures,
    errors = errors,
    failed = expectation_failures + errors,
    warnings = count_column("warning"),
    skipped = count_column("skipped")
  )
  same_count <- function(name) {
    value <- manifest[[name]]
    length(value) == 1L && !is.na(value) && isTRUE(all.equal(
      as.numeric(value), as.numeric(counts[[name]]), tolerance = 0
    ))
  }
  hashes <- c(
    results_sha256 = dkfa_hash_file(results_path),
    session_info_sha256 = dkfa_hash_file(session_info_path)
  )
  checks <- c(
    schema = identical(manifest$schema_version, "1.2.0"),
    gate = identical(manifest$gate, "full_source_test_suite"),
    source_loader = identical(
      manifest$source_loader$method, "pkgload::load_all"
    ) && isTRUE(manifest$source_loader$compile),
    nonempty_test_blocks = counts[["test_blocks"]] > 0,
    nonempty_expectations = counts[["expectations"]] > 0,
    test_blocks = same_count("test_blocks"),
    expectations = same_count("expectations"),
    expectation_failures = same_count("expectation_failures"),
    errors = same_count("errors"),
    failed = same_count("failed"),
    warnings = same_count("warnings"),
    skipped = same_count("skipped"),
    passed = identical(isTRUE(manifest$passed), counts[["failed"]] == 0),
    results_hash = identical(manifest$results_sha256,
                             unname(hashes[["results_sha256"]])),
    session_info_hash = identical(
      manifest$session_info_sha256,
      unname(hashes[["session_info_sha256"]])
    )
  )
  list(
    passed = all(checks),
    checks = checks,
    counts = counts,
    hashes = hashes,
    manifest = manifest,
    results = results
  )
}

dkfa_validate_certified_pkgdown_receipt <- function(
    site_dir, tarball_sha256, source_provenance, docs_provenance,
    required_files) {
  receipt_path <- file.path(site_dir, ".dkge-pkgdown-receipt.json")
  if (!file.exists(receipt_path)) {
    return(list(passed = FALSE, reason = "missing pkgdown receipt"))
  }
  receipt <- jsonlite::read_json(receipt_path, simplifyVector = TRUE)
  missing <- required_files[!file.exists(file.path(site_dir, required_files))]
  current_site <- dkfa_site_manifest(site_dir)
  checks <- c(
    schema = identical(receipt$schema_version, "1.4.0"),
    utf8_locale = isTRUE(receipt$utf8_locale),
    locale_recorded = is.character(receipt$locale) &&
      length(receipt$locale) == 1L && !is.na(receipt$locale) &&
      nzchar(receipt$locale),
    tarball = identical(receipt$tarball_sha256, tarball_sha256),
    certified_source = identical(
      receipt$tarball_source_tree_sha256,
      source_provenance$certified_source_tree_sha256
    ),
    raw_built_source = identical(
      receipt$tarball_raw_source_tree_sha256,
      source_provenance$built_raw_source_tree_sha256
    ),
    source_projection = identical(
      receipt$package_source_projection_sha256,
      source_provenance$package_projection_sha256
    ),
    shipped_payload = identical(
      receipt$tarball_shipped_payload_tree_sha256,
      source_provenance$shipped_payload_tree_sha256
    ),
    description_projection = identical(
      receipt$description_semantic_sha256,
      source_provenance$description$semantic_sha256
    ),
    docs = identical(
      receipt$docs_input_sha256,
      docs_provenance$certified_docs_input_sha256
    ),
    docs_overlay = identical(
      receipt$docs_overlay_bundle_sha256,
      docs_provenance$overlay_bundle_sha256
    ),
    site = identical(receipt$site_content_sha256,
                     current_site$tree_sha256),
    required_page_receipt = identical(
      as.character(unlist(receipt$required_files, use.names = FALSE)),
      as.character(required_files)
    ),
    required_pages = !length(missing)
  )
  list(
    passed = all(checks), checks = checks, receipt = receipt,
    missing_required_files = missing,
    site_content_sha256 = current_site$tree_sha256
  )
}
