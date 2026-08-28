#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
} else {
  normalizePath(
    "inst/validation/functional-alignment-certification/run-source-tests.R",
    mustWork = TRUE
  )
}
root <- normalizePath(file.path(dirname(script), "..", "..", ".."),
                      mustWork = TRUE)
court_script <- file.path(
  root, "inst", "validation", "functional-alignment-court", "court.R"
)
utils_script <- file.path(dirname(script), "evidence-utils.R")
provenance_script <- file.path(dirname(script), "package-provenance.R")
script_env <- environment()
source(utils_script, local = TRUE)
input_snapshot <- function() {
  list(
    source_tree_sha256 = dkfa_candidate_source_tree_hash(root),
    test_validation = dkfa_test_validation_manifest(root),
    protocol_history = dkfa_bundle_manifest(
      dkfa_court_harness_paths(root), root
    )
  )
}
loaded_binding <- dkfa_bind_loaded_inputs(
  input_snapshot,
  function() {
    source(utils_script, local = script_env)
    source(court_script, local = script_env)
    source(provenance_script, local = script_env)
    dkfa_load_source_package(root, quiet = TRUE)
  },
  "Source-test loader"
)
source_loader <- loaded_binding$value
source_loader$bound_input_snapshot_sha256 <-
  loaded_binding$snapshot_sha256
source_hash <- loaded_binding$snapshot$source_tree_sha256
test_validation_hash <-
  loaded_binding$snapshot$test_validation$tree_sha256
protocol_history <- dkfa_validate_current_protocol(
  root, candidate_source_tree_sha256 = source_hash
)
if (!isTRUE(protocol_history$passed)) {
  stop("The v9 protocol history is invalid.", call. = FALSE)
}
run_dir <- file.path(
  root, "data-raw", "functional-alignment-certification",
  paste0("package-", substr(source_hash, 1L, 16L), "-",
         substr(test_validation_hash, 1L, 12L))
)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
results_path <- file.path(run_dir, "source-test-results.csv")
if (file.exists(results_path)) {
  stop("Source-test evidence already exists for this package source tree.",
       call. = FALSE)
}

started <- proc.time()[["elapsed"]]
results <- testthat::test_local(
  root, reporter = "summary", stop_on_failure = FALSE,
  load_package = "none"
)
elapsed <- proc.time()[["elapsed"]] - started
end_source_hash <- dkfa_certified_source_tree_hash(root)
end_test_validation_hash <- dkfa_test_validation_manifest(root)$tree_sha256
dkfa_assert_input_snapshot(
  loaded_binding$snapshot, input_snapshot(), "Source-test run"
)
if (!identical(end_source_hash, source_hash) ||
    !identical(end_test_validation_hash, test_validation_hash)) {
  stop(
    "Package or test-validation inputs changed during the source test run.",
    call. = FALSE
  )
}
table <- as.data.frame(results)
table$result <- NULL
utils::write.csv(table, results_path, row.names = FALSE)
expectation_failures <- sum(table$failed)
errors <- sum(table$error)
failed <- expectation_failures + errors
warnings <- sum(table$warning)
skipped <- sum(table$skipped)
session_info_path <- file.path(run_dir, "session-info.txt")
writeLines(capture.output(sessionInfo()), session_info_path, useBytes = TRUE)

optional <- c(
  future = requireNamespace("future", quietly = TRUE),
  future.apply = requireNamespace("future.apply", quietly = TRUE),
  T4transport = requireNamespace("T4transport", quietly = TRUE),
  neuralign = requireNamespace("neuralign", quietly = TRUE),
  patchwork = requireNamespace("patchwork", quietly = TRUE),
  ragg = requireNamespace("ragg", quietly = TRUE),
  pkgdown = requireNamespace("pkgdown", quietly = TRUE)
)
jsonlite::write_json(
  list(
    schema_version = "1.2.0",
    gate = "full_source_test_suite",
    passed = failed == 0L,
    source_loader = source_loader,
    loaded_input_snapshot_sha256 = loaded_binding$snapshot_sha256,
    source_tree_sha256 = source_hash,
    test_validation_tree_sha256 = test_validation_hash,
    protocol_id = "dkge-functional-alignment-certification-v9",
    protocol_v2_sha256 = protocol_history$protocol_v2_sha256,
    protocol_v3_sha256 = protocol_history$protocol_v3_sha256,
    protocol_v3_supersession_sha256 =
      protocol_history$protocol_v3_supersession_sha256,
    protocol_v4_sha256 = protocol_history$protocol_v4_sha256,
    protocol_v4_invalidation_sha256 =
      protocol_history$protocol_v4_invalidation_sha256,
    v4_history_bundle_sha256 =
      protocol_history$v4_history_bundle$bundle_sha256,
    protocol_v5_sha256 = protocol_history$protocol_v5_sha256,
    protocol_v5_supersession_sha256 =
      protocol_history$protocol_v5_supersession_sha256,
    v5_history_bundle_sha256 =
      protocol_history$v5_history_bundle$bundle_sha256,
    protocol_v6_sha256 = protocol_history$protocol_v6_sha256,
    protocol_v6_supersession_sha256 =
      protocol_history$protocol_v6_supersession_sha256,
    v6_history_bundle_sha256 =
      protocol_history$v6_history_bundle$bundle_sha256,
    protocol_v7_sha256 = protocol_history$protocol_v7_sha256,
    protocol_v7_supersession_sha256 =
      protocol_history$protocol_v7_supersession_sha256,
    v7_history_bundle_sha256 =
      protocol_history$v7_history_bundle$bundle_sha256,
    protocol_v8_sha256 = protocol_history$protocol_v8_sha256,
    protocol_v8_supersession_sha256 =
      protocol_history$protocol_v8_supersession_sha256,
    v8_history_bundle_sha256 =
      protocol_history$v8_history_bundle$bundle_sha256,
    protocol_v9_sha256 = protocol_history$protocol_v9_sha256,
    protocol_history_bundle_sha256 =
      loaded_binding$snapshot$protocol_history$bundle_sha256,
    test_blocks = nrow(table),
    expectations = sum(table$nb),
    expectation_failures = expectation_failures,
    errors = errors,
    failed = failed,
    warnings = warnings,
    skipped = skipped,
    results_sha256 = fa_court_hash_file(results_path),
    session_info_sha256 = fa_court_hash_file(session_info_path),
    elapsed_seconds = elapsed,
    optional_backends = as.list(optional),
    skip_classification = list(
      long_between_subject_null = "opt-in DKGE_LONG_TESTS; unrelated to alignment candidate",
      future_apply_absence_branch = "installed backend makes absence-only branch inapplicable",
      T4transport = if (optional[["T4transport"]]) {
        "installed"
      } else {
        "not installed; optional comparison backend skip"
      },
      long_weight_null = "opt-in DKGE_LONG_TESTS; unrelated to alignment candidate"
    ),
    warning_classification = list(
      expected_model_diagnostics = paste(
        "rank reduction, degenerate classification target, reduced-rank input,",
        "and intentionally strict Sinkhorn convergence comparator warnings"
      ),
      optional_plot_backend = if (optional[["patchwork"]]) {
        "installed"
      } else {
        "patchwork absent; tested placeholder fallback"
      },
      environment = "locale and package-built-under-version messages are external to DKGE"
    ),
    R = R.version.string,
    platform = R.version$platform
  ),
  file.path(run_dir, "source-test-manifest.json"),
  pretty = TRUE, auto_unbox = TRUE
)
if (failed > 0L) quit(status = 1L)
cat(run_dir, "\n")
