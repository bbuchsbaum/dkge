#!/usr/bin/env Rscript

# Assemble the final functional-alignment certification record. This script is
# intentionally strict: it only certifies an exact source-tree fingerprint and
# never turns descriptive power into an inferential promotion gate.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
} else {
  normalizePath(
    "inst/validation/functional-alignment-certification/collect-certification.R",
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
power_helper <- file.path(dirname(script), "known-warp-power.R")
script_env <- environment()
source(utils_script, local = TRUE)
collector_loader_paths <- c(
  script, utils_script, court_script, power_helper, provenance_script
)
input_snapshot <- function() {
  list(
    source_tree_sha256 = dkfa_candidate_source_tree_hash(root),
    test_validation = dkfa_test_validation_manifest(root),
    court_harness = dkfa_bundle_manifest(
      dkfa_court_harness_paths(root), root
    ),
    power_harness = dkfa_bundle_manifest(
      dkfa_power_harness_paths(root), root
    ),
    collector_loader = dkfa_bundle_manifest(collector_loader_paths, root)
  )
}
loaded_binding <- dkfa_bind_loaded_inputs(
  input_snapshot,
  function() {
    source(utils_script, local = script_env)
    source(court_script, local = script_env)
    source(power_helper, local = script_env)
    source(provenance_script, local = script_env)
    invisible(NULL)
  },
  "Certification collector loader"
)

parse_args <- function(x) {
  x <- x[grepl("^--[^=]+=", x)]
  out <- sub("^--", "", x)
  keys <- sub("=.*$", "", out)
  values <- sub("^[^=]+=", "", out)
  stats::setNames(as.list(values), keys)
}

cli <- parse_args(commandArgs(trailingOnly = TRUE))
required_args <- c("tarball", "check-dir", "pkgdown-dir", "review")
missing_args <- setdiff(required_args, names(cli))
if (length(missing_args)) {
  stop("Missing required argument(s): ", paste(missing_args, collapse = ", "),
       call. = FALSE)
}

must_file <- function(path, label) {
  path <- normalizePath(path, mustWork = FALSE)
  if (!file.exists(path) || dir.exists(path)) {
    stop(label, " is not a file: ", path, call. = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

must_dir <- function(path, label) {
  path <- normalizePath(path, mustWork = FALSE)
  if (!dir.exists(path)) {
    stop(label, " is not a directory: ", path, call. = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

tarball <- must_file(cli[["tarball"]], "Tarball")
check_dir <- must_dir(cli[["check-dir"]], "R CMD check directory")
pkgdown_dir <- must_dir(cli[["pkgdown-dir"]], "pkgdown directory")
review_path <- must_file(cli[["review"]], "Independent review")

source_hash <- loaded_binding$snapshot$source_tree_sha256
protocol_path <- file.path(dirname(script), "protocol.json")
protocol_v2_path <- file.path(dirname(script), "protocol-v2.json")
protocol_v3_path <- file.path(dirname(script), "protocol-v3.json")
v3_supersession_path <- file.path(
  dirname(script), "protocol-v3-supersession.json"
)
v2_invalidation_path <- file.path(
  dirname(script), "court-v2-invalidation.json"
)
sinkhorn_fixture_path <- file.path(
  dirname(script), "court-v2-failing-sinkhorn-fixture.json"
)
v2_historical_paths <- dkfa_v2_historical_paths(root)
v2_historical_bundle <- dkfa_review_evidence_bundle(v2_historical_paths)
v3_history_paths <- dkfa_v3_history_paths(root)
v3_history_bundle <- dkfa_review_evidence_bundle(v3_history_paths)
v4_history_paths <- dkfa_v4_history_paths(root)
v4_history_bundle <- dkfa_review_evidence_bundle(v4_history_paths)
v5_history_paths <- dkfa_v5_history_paths(root)
v5_history_bundle <- dkfa_review_evidence_bundle(v5_history_paths)
v6_history_paths <- dkfa_v6_history_paths(root)
v6_history_bundle <- dkfa_review_evidence_bundle(v6_history_paths)
v7_history_paths <- dkfa_v7_history_paths(root)
v7_history_bundle <- dkfa_review_evidence_bundle(v7_history_paths)
v8_history_paths <- dkfa_v8_history_paths(root)
v8_history_bundle <- dkfa_review_evidence_bundle(v8_history_paths)
protocol <- jsonlite::read_json(protocol_path, simplifyVector = TRUE)
protocol_hash <- fa_court_hash_file(protocol_path)
protocol_amendment <- dkfa_validate_current_protocol(
  root, candidate_source_tree_sha256 = source_hash
)
v2_invalidation_status <- dkfa_validate_court_v2_invalidation(
  v2_invalidation_path, sinkhorn_fixture_path, protocol_v2_path,
  historical_paths = v2_historical_paths,
  verify_solver = FALSE, verify_origin = TRUE
)
evidence_root <- file.path(root, "data-raw",
                           "functional-alignment-certification")

test_validation_hash <-
  loaded_binding$snapshot$test_validation$tree_sha256

court_runner <- file.path(dirname(script), "run-court.R")
court_analysis <- file.path(dirname(script), "analyze-court.R")
court_bundle <- loaded_binding$snapshot$court_harness
power_runner <- file.path(dirname(script), "run-known-warp-power.R")
template_helper <- file.path(
  root, "inst", "validation", "functional-alignment-template",
  "template-benefit.R"
)
determinism_script <- file.path(
  dirname(script), "verify-known-warp-determinism.R"
)
power_bundle <- loaded_binding$snapshot$power_harness

find_run <- function(prefix, required_file, label) {
  candidates <- list.dirs(evidence_root, recursive = FALSE, full.names = TRUE)
  candidates <- candidates[basename(candidates) == prefix]
  candidates <- candidates[file.exists(file.path(candidates, required_file))]
  if (length(candidates) != 1L) {
    stop("Expected exactly one ", label, " for prefix '", prefix,
         "'; found ", length(candidates), ".", call. = FALSE)
  }
  candidates[[1L]]
}

court_dir <- find_run(
  paste0("court-", substr(source_hash, 1L, 16L), "-",
         substr(court_bundle$bundle_sha256, 1L, 12L)),
  "court-verdict.json", "formal court run"
)
power_dir <- find_run(
  paste0("known-warp-", substr(source_hash, 1L, 12L), "-",
         substr(protocol_hash, 1L, 12L), "-",
         substr(power_bundle$bundle_sha256, 1L, 12L)),
  "known-warp-power-manifest.json", "known-warp run"
)
source_test_dir <- find_run(
  paste0("package-", substr(source_hash, 1L, 16L), "-",
         substr(test_validation_hash, 1L, 12L)),
  "source-test-manifest.json", "source-test run"
)

court <- jsonlite::read_json(file.path(court_dir, "court-verdict.json"),
                             simplifyVector = TRUE)
court_run_state <- jsonlite::read_json(
  file.path(court_dir, "court-run-state.json"), simplifyVector = TRUE
)
court_audit <- utils::read.csv(file.path(court_dir, "court-audit.csv"),
                               stringsAsFactors = FALSE)
court_raw <- utils::read.csv(
  file.path(court_dir, "formal-raw.csv"), stringsAsFactors = FALSE
)
court_recorded_summary <- utils::read.csv(
  file.path(court_dir, "formal-summary.csv"), stringsAsFactors = FALSE
)
court_evidence <- dkfa_validate_court_evidence(
  court_raw, court_recorded_summary, protocol
)
court_focus <- utils::read.csv(
  file.path(court_dir, "court-focus-summary.csv"),
  stringsAsFactors = FALSE
)
power_manifest <- jsonlite::read_json(
  file.path(power_dir, "known-warp-power-manifest.json"),
  simplifyVector = TRUE
)
power_run_state <- jsonlite::read_json(
  file.path(power_dir, "known-warp-power-run-state.json"),
  simplifyVector = TRUE
)
power_harness_record <- jsonlite::read_json(
  file.path(power_dir, "known-warp-power-harness-manifest.json"),
  simplifyVector = TRUE
)
power_raw <- utils::read.csv(
  file.path(power_dir, "known-warp-power-raw.csv"),
  stringsAsFactors = FALSE
)
power_summary <- utils::read.csv(
  file.path(power_dir, "known-warp-power-summary.csv"),
  stringsAsFactors = FALSE
)
power_evidence <- dkfa_validate_known_warp_evidence(
  power_raw, power_summary, protocol
)
if (isTRUE(power_evidence$passed)) {
  power_summary <- power_evidence$recomputed_summary
}
determinism <- jsonlite::read_json(
  file.path(power_dir, "serial-parallel-determinism.json"),
  simplifyVector = TRUE
)
source_test_evidence <- dkfa_validate_source_test_evidence(
  file.path(source_test_dir, "source-test-manifest.json"),
  file.path(source_test_dir, "source-test-results.csv"),
  file.path(source_test_dir, "session-info.txt")
)
source_tests <- source_test_evidence$manifest
review <- jsonlite::read_json(review_path, simplifyVector = TRUE)

court_harness_record <- jsonlite::read_json(
  file.path(court_dir, "court-harness-manifest.json"),
  simplifyVector = TRUE
)
dkfa_verify_bundle_manifest(
  power_harness_record,
  dkfa_power_harness_paths(root),
  root, "Known-warp run-state"
)
dkfa_verify_bundle_manifest(
  court_harness_record,
  dkfa_court_harness_paths(root),
  root, "Court"
)
dkfa_verify_bundle_manifest(
  list(
    bundle_sha256 = power_manifest$harness_bundle_sha256,
    files = power_manifest$harness_files
  ),
  dkfa_power_harness_paths(root),
  root, "Known-warp"
)
court_manifest_lines <- readLines(
  file.path(court_dir, "formal-manifest.txt"), warn = FALSE
)
court_manifest_value <- function(key) {
  lines <- court_manifest_lines[
    startsWith(court_manifest_lines, paste0(key, "="))
  ]
  if (!length(lines)) return(NA_character_)
  sub(paste0("^", key, "="), "", lines[[length(lines)]])
}

audit_value <- function(name) {
  hit <- court_audit$value[court_audit$check == name]
  if (length(hit) != 1L) stop("Missing court audit entry: ", name,
                              call. = FALSE)
  as.numeric(hit)
}
run_state_output <- function(stage, name) {
  value <- court_run_state$outputs[[stage]][[name]]
  if (is.null(value)) return(NA_character_)
  as.character(value)
}
court_formal_state_paths <- c(
  formal_budget = file.path(court_dir, "formal-budget.csv"),
  formal_grid = file.path(court_dir, "formal-grid.csv"),
  formal_raw = file.path(court_dir, "formal-raw.csv"),
  formal_summary = file.path(court_dir, "formal-summary.csv"),
  formal_inflation_flags = file.path(court_dir, "formal-inflation-flags.csv"),
  latent_span_model = file.path(court_dir, "latent-span-model.csv"),
  formal_manifest = file.path(court_dir, "formal-manifest.txt")
)
court_analysis_state_paths <- c(
  court_audit = file.path(court_dir, "court-audit.csv"),
  court_focus_summary = file.path(court_dir, "court-focus-summary.csv"),
  factor_uncertainty = file.path(court_dir, "factor-uncertainty.csv"),
  latent_span_excess = file.path(court_dir, "latent-span-excess-model.csv"),
  latent_span_logistic = file.path(
    court_dir, "latent-span-logistic-model.csv"
  ),
  court_verdict = file.path(court_dir, "court-verdict.json"),
  court_verdict_markdown = file.path(court_dir, "court-verdict.md")
)
run_state_hashes_match <- function(stage, paths) {
  recorded <- unlist(court_run_state$outputs[[stage]], use.names = TRUE)
  all(file.exists(paths)) && identical(
    recorded, vapply(paths, dkfa_hash_file, character(1))
  )
}
power_run_state_hashes_match <- function(stage, paths) {
  recorded <- unlist(power_run_state$outputs[[stage]], use.names = TRUE)
  all(file.exists(paths)) && identical(
    recorded, vapply(paths, dkfa_hash_file, character(1))
  )
}
power_state_paths <- c(
  power_manifest = file.path(power_dir, "known-warp-power-manifest.json"),
  power_raw = file.path(power_dir, "known-warp-power-raw.csv"),
  power_summary = file.path(power_dir, "known-warp-power-summary.csv")
)
determinism_state_path <- c(
  determinism = file.path(power_dir, "serial-parallel-determinism.json")
)

numeric_power <- power_raw[vapply(power_raw, is.numeric, logical(1))]
power_finite <- all(vapply(numeric_power, function(x) all(is.finite(x)),
                           logical(1)))
expected_power_rows <- protocol$known_warp_power$seeds$n_cohorts *
  length(protocol$known_warp_power$arms)
expected_power_arms <- sort(protocol$known_warp_power$arms)
template_convergence_by_seed <- vapply(
  split(power_raw$template_converged, power_raw$seed),
  function(x) {
    observed <- unique(x)
    if (length(observed) == 1L && !is.na(observed)) observed else NA
  },
  logical(1)
)

required_site_files <- c(
  "index.html",
  "articles/dkge-functional-alignment.html",
  "reference/index.html",
  "reference/dkge_prepare_alignment.html",
  "reference/dkge_fit_functional_template.html",
  "reference/dkge_infer_aligned.html",
  "reference/dkge_render_aligned.html"
)

extract_dir <- tempfile("dkge-certification-tarball-")
dir.create(extract_dir)
on.exit(unlink(extract_dir, recursive = TRUE), add = TRUE)
utils::untar(tarball, exdir = extract_dir)
top_dirs <- list.dirs(extract_dir, recursive = FALSE, full.names = TRUE)
if (length(top_dirs) != 1L) {
  stop("Tarball must contain exactly one package root.", call. = FALSE)
}
tarball_root <- top_dirs[[1L]]
tarball_sha256 <- dkfa_hash_file(tarball)
source_provenance <- dkfa_verify_built_source(
  root, tarball_root, certified_source_tree_sha256 = source_hash
)
docs_provenance <- dkfa_docs_projection(root, tarball_root)
tarball_source_hash <- source_provenance$certified_source_tree_sha256
tarball_raw_source_hash <- source_provenance$built_raw_source_tree_sha256
tarball_docs <- docs_provenance$certified_docs_input_sha256
current_docs <- dkfa_pkgdown_input_manifest(root)$tree_sha256

check_log <- must_file(file.path(check_dir, "00check.log"),
                       "R CMD check log")
check_lines <- readLines(check_log, warn = FALSE)
check_status <- dkfa_validate_r_cmd_check_log(check_lines)
check_ok <- isTRUE(check_status$checks[["status_ok"]])
check_policy_ok <- all(check_status$checks[names(check_status$checks) != "status_ok"])
required_check_options <- check_status$required_options
allowed_check_options <- check_status$allowed_options
package_name <- read.dcf(file.path(tarball_root, "DESCRIPTION"),
                         fields = "Package")[[1L]]
check_source_root <- file.path(check_dir, "00_pkg_src", package_name)
if (!dir.exists(check_source_root)) {
  stop("R CMD check did not retain its checked package source tree: ",
       check_source_root, call. = FALSE)
}
compiled_excludes <- c(
  "\\.(o|so|dll|dylib)$", "(^|/)src/symbols\\.rds$", "(^|/)\\.DS_Store$"
)
tarball_tree <- dkfa_tree_manifest(
  tarball_root, roots = ".", exclude = compiled_excludes
)
check_tree <- dkfa_tree_manifest(
  check_source_root, roots = ".", exclude = compiled_excludes
)

pkgdown_status <- dkfa_validate_certified_pkgdown_receipt(
  pkgdown_dir,
  tarball_sha256 = tarball_sha256,
  source_provenance = source_provenance,
  docs_provenance = docs_provenance,
  required_files = required_site_files
)
pkgdown_builder <- dkfa_bundle_manifest(
  c(file.path(dirname(script), "run-pkgdown.R"), utils_script, court_script,
    provenance_script),
  root
)
pkgdown_loaded_snapshot_sha256 <- dkfa_input_snapshot_sha256(list(
  source_tree_sha256 = source_hash,
  docs = dkfa_pkgdown_input_manifest(root),
  builder = pkgdown_builder
))
if (!is.null(pkgdown_status$receipt)) {
  dkfa_verify_bundle_manifest(
    list(
      bundle_sha256 = pkgdown_status$receipt$builder_bundle_sha256,
      files = pkgdown_status$receipt$builder_files
    ),
    c(file.path(dirname(script), "run-pkgdown.R"), utils_script, court_script,
      provenance_script),
    root, "Pkgdown builder"
  )
}

review_evidence_paths <- c(
  protocol = protocol_path,
  protocol_v2 = protocol_v2_path,
  protocol_v3 = protocol_v3_path,
  protocol_v3_supersession = v3_supersession_path,
  court_v2_invalidation = v2_invalidation_path,
  court_v2_sinkhorn_fixture = sinkhorn_fixture_path,
  v2_source_test_manifest =
    v2_historical_paths[["source_test_manifest"]],
  v2_source_test_results =
    v2_historical_paths[["source_test_results"]],
  v2_source_test_session_info =
    v2_historical_paths[["source_test_session_info"]],
  v2_failed_court_formal_budget =
    v2_historical_paths[["failed_court_formal_budget"]],
  court_run_state = file.path(court_dir, "court-run-state.json"),
  court_verdict = file.path(court_dir, "court-verdict.json"),
  court_verdict_markdown = file.path(court_dir, "court-verdict.md"),
  court_raw = file.path(court_dir, "formal-raw.csv"),
  court_summary = file.path(court_dir, "formal-summary.csv"),
  court_audit = file.path(court_dir, "court-audit.csv"),
  court_focus = file.path(court_dir, "court-focus-summary.csv"),
  court_factor_uncertainty = file.path(court_dir, "factor-uncertainty.csv"),
  court_formal_budget = file.path(court_dir, "formal-budget.csv"),
  court_formal_grid = file.path(court_dir, "formal-grid.csv"),
  court_inflation_flags = file.path(
    court_dir, "formal-inflation-flags.csv"
  ),
  court_latent_span_excess = file.path(
    court_dir, "latent-span-excess-model.csv"
  ),
  court_latent_span_logistic = file.path(
    court_dir, "latent-span-logistic-model.csv"
  ),
  court_latent_span_model = file.path(court_dir, "latent-span-model.csv"),
  court_formal_manifest = file.path(court_dir, "formal-manifest.txt"),
  court_harness_manifest = file.path(court_dir, "court-harness-manifest.json"),
  known_warp_manifest = file.path(
    power_dir, "known-warp-power-manifest.json"
  ),
  known_warp_run_state = file.path(
    power_dir, "known-warp-power-run-state.json"
  ),
  known_warp_harness_manifest = file.path(
    power_dir, "known-warp-power-harness-manifest.json"
  ),
  known_warp_raw = file.path(power_dir, "known-warp-power-raw.csv"),
  known_warp_summary = file.path(
    power_dir, "known-warp-power-summary.csv"
  ),
  known_warp_determinism = file.path(
    power_dir, "serial-parallel-determinism.json"
  ),
  source_test_manifest = file.path(
    source_test_dir, "source-test-manifest.json"
  ),
  source_test_results = file.path(
    source_test_dir, "source-test-results.csv"
  ),
  source_test_session_info = file.path(source_test_dir, "session-info.txt"),
  check_log = check_log,
  pkgdown_receipt = file.path(pkgdown_dir, ".dkge-pkgdown-receipt.json")
)
review_evidence_paths <- c(
  review_evidence_paths,
  v4_history_paths,
  v5_history_paths,
  v6_history_paths,
  v7_history_paths,
  v8_history_paths
)
review_evidence <- dkfa_review_evidence_bundle(review_evidence_paths)
review_evidence_hashes <- review_evidence$hashes
review_evidence_bundle_sha256 <- review_evidence$bundle_sha256

blockers <- review$blockers
if (is.null(blockers)) blockers <- character()
blockers <- as.character(blockers)
review_ok <- review$verdict %in% c("pass", "pass_with_limitations") &&
  length(blockers) == 0L &&
  identical(review$reviewed_source_tree_sha256, source_hash) &&
  identical(review$reviewed_test_validation_tree_sha256,
            test_validation_hash) &&
  identical(review$reviewed_tarball_sha256, tarball_sha256) &&
  identical(review$reviewed_check_source_tree_sha256,
            check_tree$tree_sha256) &&
  identical(review$reviewed_pkgdown_site_sha256,
            pkgdown_status$site_content_sha256) &&
  identical(review$reviewed_court_bundle_sha256,
            court_bundle$bundle_sha256) &&
  identical(review$reviewed_known_warp_bundle_sha256,
            power_bundle$bundle_sha256) &&
  identical(review$reviewed_evidence_bundle_sha256,
            review_evidence_bundle_sha256) &&
  identical(review$reviewed_collector_sha256, dkfa_hash_file(script))

protocol_history_bindings <- c(
  protocol_v2_sha256 = protocol_amendment$protocol_v2_sha256,
  protocol_v3_sha256 = protocol_amendment$protocol_v3_sha256,
  protocol_v3_supersession_sha256 =
    protocol_amendment$protocol_v3_supersession_sha256,
  protocol_v4_sha256 = protocol_amendment$protocol_v4_sha256,
  protocol_v4_invalidation_sha256 =
    protocol_amendment$protocol_v4_invalidation_sha256,
  v4_history_bundle_sha256 =
    protocol_amendment$v4_history_bundle$bundle_sha256,
  protocol_v5_sha256 = protocol_amendment$protocol_v5_sha256,
  protocol_v5_supersession_sha256 =
    protocol_amendment$protocol_v5_supersession_sha256,
  v5_history_bundle_sha256 =
    protocol_amendment$v5_history_bundle$bundle_sha256,
  protocol_v6_sha256 = protocol_amendment$protocol_v6_sha256,
  protocol_v6_supersession_sha256 =
    protocol_amendment$protocol_v6_supersession_sha256,
  v6_history_bundle_sha256 =
    protocol_amendment$v6_history_bundle$bundle_sha256,
  protocol_v7_sha256 = protocol_amendment$protocol_v7_sha256,
  protocol_v7_supersession_sha256 =
    protocol_amendment$protocol_v7_supersession_sha256,
  v7_history_bundle_sha256 =
    protocol_amendment$v7_history_bundle$bundle_sha256,
  protocol_v8_sha256 = protocol_amendment$protocol_v8_sha256,
  protocol_v8_supersession_sha256 =
    protocol_amendment$protocol_v8_supersession_sha256,
  v8_history_bundle_sha256 =
    protocol_amendment$v8_history_bundle$bundle_sha256,
  protocol_v9_sha256 = protocol_amendment$protocol_v9_sha256
)
court_loaded_snapshot_sha256 <- dkfa_input_snapshot_sha256(list(
  source_tree_sha256 = source_hash,
  harness = court_bundle
))
power_loaded_snapshot_sha256 <- dkfa_input_snapshot_sha256(list(
  source_tree_sha256 = source_hash,
  harness = power_bundle
))
source_test_loaded_snapshot_sha256 <- dkfa_input_snapshot_sha256(list(
  source_tree_sha256 = source_hash,
  test_validation = loaded_binding$snapshot$test_validation,
  protocol_history = loaded_binding$snapshot$court_harness
))
object_history_bound <- function(object) {
  all(vapply(names(protocol_history_bindings), function(name) {
    identical(object[[name]], unname(protocol_history_bindings[[name]]))
  }, logical(1)))
}
manifest_history_bound <- function() {
  all(vapply(names(protocol_history_bindings), function(name) {
    identical(court_manifest_value(name),
              unname(protocol_history_bindings[[name]]))
  }, logical(1)))
}

gates <- list(
  v9_protocol_history_valid = isTRUE(protocol_amendment$passed),
  retained_history_evidence_bound = identical(
    unname(v3_history_bundle$hashes[["protocol_v3"]]),
    protocol_amendment$protocol_v3_sha256
  ) && identical(
    unname(v3_history_bundle$hashes[["protocol_v3_supersession"]]),
    protocol_amendment$protocol_v3_supersession_sha256
  ) && identical(
    power_manifest$protocol_v3_sha256,
    protocol_amendment$protocol_v3_sha256
  ) && identical(
    power_manifest$protocol_v3_supersession_sha256,
    protocol_amendment$protocol_v3_supersession_sha256
  ) && identical(
    v4_history_bundle$bundle_sha256,
    protocol_amendment$v4_history_bundle$bundle_sha256
  ) && identical(
    v5_history_bundle$bundle_sha256,
    protocol_amendment$v5_history_bundle$bundle_sha256
  ) && identical(
    v6_history_bundle$bundle_sha256,
    protocol_amendment$v6_history_bundle$bundle_sha256
  ) && identical(
    v7_history_bundle$bundle_sha256,
    protocol_amendment$v7_history_bundle$bundle_sha256
  ) && identical(
    v8_history_bundle$bundle_sha256,
    protocol_amendment$v8_history_bundle$bundle_sha256
  ),
  v2_invalidation_chain_valid = isTRUE(v2_invalidation_status$passed),
  v2_historical_evidence_bound = identical(
    power_manifest$v2_historical_evidence_bundle_sha256,
    v2_historical_bundle$bundle_sha256
  ) && identical(
    unlist(power_manifest$v2_historical_evidence_hashes, use.names = TRUE),
    v2_historical_bundle$hashes
  ) && identical(
    determinism$v2_historical_evidence_bundle_sha256,
    v2_historical_bundle$bundle_sha256
  ) && identical(
    court$v2_historical_evidence_bundle_sha256,
    v2_historical_bundle$bundle_sha256
  ),
  source_tree_matches_court = identical(court$source_tree_sha256, source_hash),
  protocol_matches_court = identical(
    court$protocol_id, protocol$protocol_id
  ) && object_history_bound(court) && identical(
    court$court_v2_invalidation_sha256,
    dkfa_hash_file(v2_invalidation_path)
  ) && identical(
    court$court_v2_sinkhorn_fixture_sha256,
    dkfa_hash_file(sinkhorn_fixture_path)
  ),
  source_tree_matches_power = identical(
    power_manifest$source_tree_sha256, source_hash
  ),
  source_tree_matches_source_tests = identical(
    source_tests$source_tree_sha256, source_hash
  ),
  test_validation_tree_matches_source_tests = identical(
    source_tests$test_validation_tree_sha256, test_validation_hash
  ),
  tarball_contains_certified_source = isTRUE(source_provenance$passed) &&
    identical(tarball_source_hash, source_hash),
  tarball_docs_match_live_docs = isTRUE(docs_provenance$passed) &&
    identical(tarball_docs, current_docs),
  protocol_matches_power = identical(
    power_manifest$protocol_sha256, protocol_hash
  ) && identical(
    power_manifest$protocol_id, protocol$protocol_id
  ) && object_history_bound(power_manifest) && identical(
    power_manifest$court_v2_invalidation_sha256,
    dkfa_hash_file(v2_invalidation_path)
  ) && identical(
    power_manifest$court_v2_sinkhorn_fixture_sha256,
    dkfa_hash_file(sinkhorn_fixture_path)
  ) && identical(
    power_manifest$v2_historical_evidence_bundle_sha256,
    v2_historical_bundle$bundle_sha256
  ),
  court_harness_matches_record = identical(
    court_harness_record$bundle_sha256, court_bundle$bundle_sha256
  ) && identical(
    court$court_harness_bundle_sha256, court_bundle$bundle_sha256
  ) && identical(
    court_manifest_value("court_harness_bundle_sha256"),
    court_bundle$bundle_sha256
  ),
  court_recorded_helpers_match = identical(
    court_manifest_value("script_sha256"), dkfa_hash_file(court_script)
  ) && identical(
    court_manifest_value("historical_court_script_sha256"),
    dkfa_hash_file(court_script)
  ) && identical(
    court_manifest_value("certification_wrapper_sha256"),
    dkfa_hash_file(court_runner)
  ) && identical(
    court_manifest_value("protocol_sha256"), protocol_hash
  ) && identical(
    court_manifest_value("certification_protocol_id"), protocol$protocol_id
  ) && manifest_history_bound() && identical(
    court_manifest_value("court_v2_invalidation_sha256"),
    dkfa_hash_file(v2_invalidation_path)
  ) && identical(
    court_manifest_value("court_v2_sinkhorn_fixture_sha256"),
    dkfa_hash_file(sinkhorn_fixture_path)
  ) && identical(
    court_manifest_value("v2_historical_evidence_bundle_sha256"),
    v2_historical_bundle$bundle_sha256
  ) && identical(
    court_manifest_value("v2_source_test_manifest_sha256"),
    v2_historical_bundle$hashes[["source_test_manifest"]]
  ) && identical(
    court_manifest_value("v2_source_test_results_sha256"),
    v2_historical_bundle$hashes[["source_test_results"]]
  ) && identical(
    court_manifest_value("v2_source_test_session_info_sha256"),
    v2_historical_bundle$hashes[["source_test_session_info"]]
  ) && identical(
    court_manifest_value("v2_failed_court_formal_budget_sha256"),
    v2_historical_bundle$hashes[["failed_court_formal_budget"]]
  ) && identical(
    court_manifest_value("sinkhorn_max_iter"), "50000"
  ) && identical(
    court_manifest_value("v3_only_numerical_amendment_verified"), "TRUE"
  ) && identical(
    court_manifest_value(
      "v9_no_numerical_or_statistical_amendment_verified"
    ),
    "TRUE"
  ) && identical(
    court_manifest_value("loaded_input_snapshot_sha256"),
    court_loaded_snapshot_sha256
  ) && identical(
    court_manifest_value("n_sim"), as.character(protocol$court$n_sim)
  ) && identical(
    court_manifest_value("n_perm"), as.character(protocol$court$n_perm)
  ) && identical(
    court_manifest_value("requested_formal_cores"),
    as.character(protocol$court$formal_cores)
  ),
  certification_source_loaders_bound = identical(
    court_manifest_value("source_loader"), "pkgload::load_all"
  ) && identical(
    court_manifest_value("source_loader_compile"), "TRUE"
  ) && identical(
    power_manifest$source_loader$method, "pkgload::load_all"
  ) && isTRUE(
    power_manifest$source_loader$compile
  ) && identical(
    determinism$source_loader$method, "pkgload::load_all"
  ) && isTRUE(
    determinism$source_loader$compile
  ) && identical(
    source_tests$source_loader$method, "pkgload::load_all"
  ) && isTRUE(
    source_tests$source_loader$compile
  ) && identical(
    power_manifest$loaded_input_snapshot_sha256,
    power_loaded_snapshot_sha256
  ) && identical(
    power_manifest$source_loader$bound_input_snapshot_sha256,
    power_loaded_snapshot_sha256
  ) && identical(
    determinism$loaded_input_snapshot_sha256,
    power_loaded_snapshot_sha256
  ) && identical(
    determinism$source_loader$bound_input_snapshot_sha256,
    power_loaded_snapshot_sha256
  ) && identical(
    source_tests$loaded_input_snapshot_sha256,
    source_test_loaded_snapshot_sha256
  ) && identical(
    source_tests$source_loader$bound_input_snapshot_sha256,
    source_test_loaded_snapshot_sha256
  ),
  known_warp_harness_matches_record = identical(
    power_manifest$harness_bundle_sha256, power_bundle$bundle_sha256
  ) && identical(
    power_harness_record$bundle_sha256, power_bundle$bundle_sha256
  ),
  known_warp_run_state_complete = identical(
    power_run_state$status, "determinism_complete"
  ) && identical(
    power_run_state$run_id, basename(power_dir)
  ) && identical(
    power_run_state$source_tree_sha256, source_hash
  ) && identical(
    power_run_state$protocol_sha256, protocol_hash
  ) && object_history_bound(power_run_state) && identical(
    power_run_state$loaded_input_snapshot_sha256,
    power_loaded_snapshot_sha256
  ) && identical(
    power_run_state$court_v2_invalidation_sha256,
    dkfa_hash_file(v2_invalidation_path)
  ) && identical(
    power_run_state$court_v2_sinkhorn_fixture_sha256,
    dkfa_hash_file(sinkhorn_fixture_path)
  ) && identical(
    power_run_state$power_harness_bundle_sha256,
    power_bundle$bundle_sha256
  ) && power_run_state_hashes_match(
    "power_complete", power_state_paths
  ) && power_run_state_hashes_match(
    "determinism_complete", determinism_state_path
  ),
  court_run_state_analyzed_complete = identical(
    court_run_state$status, "analyzed_complete"
  ) && identical(
    court_run_state$run_id, basename(court_dir)
  ) && identical(
    court_run_state$source_tree_sha256, source_hash
  ) && identical(
    court_run_state$protocol_sha256, protocol_hash
  ) && object_history_bound(court_run_state) && identical(
    court_run_state$loaded_input_snapshot_sha256,
    court_loaded_snapshot_sha256
  ) && identical(
    court_run_state$court_v2_invalidation_sha256,
    dkfa_hash_file(v2_invalidation_path)
  ) && identical(
    court_run_state$court_v2_sinkhorn_fixture_sha256,
    dkfa_hash_file(sinkhorn_fixture_path)
  ) && identical(
    court_run_state$court_harness_bundle_sha256,
    court_bundle$bundle_sha256
  ) && run_state_hashes_match(
    "formal_complete", court_formal_state_paths
  ) && run_state_hashes_match(
    "analyzed_complete", court_analysis_state_paths
  ),
  court_schedule_and_summary_recomputed = isTRUE(court_evidence$passed),
  court_verdict_recomputed_from_raw = isTRUE(court_evidence$passed) &&
    identical(
      court$exact_oracle_gate,
      court_evidence$recomputed_verdict$exact_oracle_gate
    ) && identical(
      court$negative_control_promotion_gate,
      court_evidence$recomputed_verdict$negative_control_gate
    ),
  court_complete = identical(
    court$status,
    "complete_approximate_only_inferential_promotion_blocked"
  ),
  exact_oracle_passed = isTRUE(court$exact_oracle_gate) &&
    audit_value("exact_oracle_gate") == 1,
  inferential_promotion_blocked = !isTRUE(
    court$negative_control_promotion_gate
  ) && audit_value("negative_control_promotion_gate") == 0,
  court_source_unchanged = isTRUE(court$source_tree_unchanged) &&
    audit_value("source_tree_unchanged") == 1,
  no_posthoc_threshold_changes = !isTRUE(
    court$thresholds_changed_after_results
  ) && audit_value("posthoc_threshold_changes") == 0,
  court_rows_complete = audit_value("raw_rows") ==
    audit_value("expected_raw_rows") &&
    audit_value("unique_cohorts") == audit_value("expected_unique_cohorts"),
  court_required_values_finite = audit_value("nonfinite_required_raw") == 0,
  known_warp_budget_complete = nrow(power_raw) == expected_power_rows &&
    length(unique(power_raw$seed)) ==
      protocol$known_warp_power$seeds$n_cohorts &&
    identical(sort(unique(power_raw$arm)), expected_power_arms),
  known_warp_schedule_and_summary_recomputed = isTRUE(power_evidence$passed),
  known_warp_manifest_matches_frozen_schedule = identical(
    power_manifest$seeds$first,
    protocol$known_warp_power$seeds$first
  ) && identical(
    power_manifest$seeds$last,
    protocol$known_warp_power$seeds$last
  ) && identical(
    power_manifest$seeds$n,
    protocol$known_warp_power$seeds$n_cohorts
  ) && identical(
    power_manifest$n_perm, protocol$known_warp_power$n_perm
  ) && identical(
    power_manifest$heldout_signal_scale,
    protocol$known_warp_power$heldout_signal_scale
  ) && identical(
    power_manifest$heldout_value_noise,
    protocol$known_warp_power$heldout_value_noise
  ),
  known_warp_values_finite = power_finite,
  known_warp_numerically_valid = all(power_raw$solver_converged),
  known_warp_template_status_recorded =
    length(template_convergence_by_seed) ==
      protocol$known_warp_power$seeds$n_cohorts &&
    !anyNA(template_convergence_by_seed),
  known_warp_serial_parallel_deterministic = isTRUE(determinism$passed),
  known_warp_determinism_matches_run = identical(
    determinism$parallel_run_id, basename(power_dir)
  ) && object_history_bound(determinism) && identical(
    determinism$loaded_input_snapshot_sha256,
    power_loaded_snapshot_sha256
  ) && identical(
    determinism$court_v2_invalidation_sha256,
    dkfa_hash_file(v2_invalidation_path)
  ) && identical(
    determinism$court_v2_sinkhorn_fixture_sha256,
    dkfa_hash_file(sinkhorn_fixture_path)
  ) && identical(
    determinism$power_harness_bundle_sha256,
    power_bundle$bundle_sha256
  ) && identical(
    unlist(determinism$power_output_hashes, use.names = TRUE),
    vapply(power_state_paths, dkfa_hash_file, character(1))
  ),
  source_test_evidence_bound = isTRUE(source_test_evidence$passed),
  source_test_protocol_history_bound = identical(
    source_tests$protocol_id, protocol$protocol_id
  ) && object_history_bound(source_tests) && identical(
    source_tests$protocol_history_bundle_sha256,
    court_bundle$bundle_sha256
  ),
  full_source_suite_passed =
    identical(source_tests$gate, "full_source_test_suite") &&
    source_tests$test_blocks > 0 && source_tests$expectations > 0 &&
    isTRUE(source_tests$passed) &&
    source_tests$failed == 0 && source_tests$errors == 0 &&
    source_tests$expectation_failures == 0,
  built_tarball_check_passed = check_ok,
  built_tarball_check_policy_bound = check_policy_ok,
  checked_source_matches_tarball = identical(
    check_tree$tree_sha256, tarball_tree$tree_sha256
  ),
  pkgdown_receipt_matches_tarball = isTRUE(pkgdown_status$passed),
  pkgdown_utf8_locale_recorded = !is.null(pkgdown_status$receipt) &&
    isTRUE(pkgdown_status$receipt$utf8_locale) &&
    is.character(pkgdown_status$receipt$locale) &&
    length(pkgdown_status$receipt$locale) == 1L &&
    !is.na(pkgdown_status$receipt$locale) &&
    nzchar(pkgdown_status$receipt$locale),
  pkgdown_builder_matches_record = !is.null(pkgdown_status$receipt) &&
    identical(pkgdown_status$receipt$builder_bundle_sha256,
              pkgdown_builder$bundle_sha256) && identical(
      pkgdown_status$receipt$loaded_input_snapshot_sha256,
      pkgdown_loaded_snapshot_sha256
    ),
  independent_review_passed = review_ok
)

failed_gates <- names(gates)[!unlist(gates, use.names = FALSE)]
if (length(failed_gates)) {
  stop("Certification gate(s) failed: ", paste(failed_gates, collapse = ", "),
       call. = FALSE)
}

find_test <- function(path, name) {
  absolute <- file.path(root, path)
  lines <- readLines(absolute, warn = FALSE)
  needle <- paste0("test_that(\"", name, "\"")
  line <- which(grepl(needle, lines, fixed = TRUE))
  if (length(line) != 1L) {
    stop("Could not resolve unique test anchor '", name, "' in ", path,
         call. = FALSE)
  }
  list(file = path, line = unname(line), test = name)
}

evidence_map <- list(
  algebraic_oracles = list(
    fold_weight_renormalization = find_test(
      "tests/testthat/test-fit.R",
      "held-out beta energy cannot change training-fold MFA weights or Chat"
    ),
    legacy_weight_inversion = find_test(
      "tests/testthat/test-fit.R",
      "legacy fits recover fold weights from final shrunken weights"
    ),
    value_reconstruction_and_covariance = find_test(
      "tests/testthat/test-alignment-features.R",
      "conditional residualization matches a hand-matrix oracle"
    ),
    whole_field_covariance = find_test(
      "tests/testthat/test-alignment-features.R",
      "the separable oracle covers every cross-parcel covariance"
    ),
    constant_and_mass_preservation = find_test(
      "tests/testthat/test-transport-contracts.R",
      "joint plans and application operators obey distinct conservation laws"
    ),
    reference_self_map = find_test(
      "tests/testthat/test-alignment-operator-symmetry.R",
      "reference and non-reference subjects use one mapper policy"
    ),
    template_mass_update = find_test(
      "tests/testthat/test-functional-template-fit.R",
      "mass-aware template update matches a direct plan oracle"
    )
  ),
  regression_contracts = list(
    parametric_transport = find_test(
      "tests/testthat/test-inference.R",
      "parametric helper consumes transported rows with unequal cluster counts"
    ),
    center_mode = find_test(
      "tests/testthat/test-inference.R",
      "mean max-T rejects retired center modes at the public boundary"
    ),
    mixed_cache_provenance = find_test(
      "tests/testthat/test-transport-prepare.R",
      "fitted alignment rejects every mutated structural input"
    ),
    full_fit_loading_fallback = find_test(
      "tests/testthat/test-inference-transport.R",
      "inferential transport never falls back to full-fit loadings"
    ),
    identity_reference_asymmetry = find_test(
      "tests/testthat/test-alignment-operator-symmetry.R",
      "reference and non-reference subjects use one mapper policy"
    ),
    feature_collapse_and_singular_span = find_test(
      "tests/testthat/test-alignment-features.R",
      "rank, condition, energy, and covariance gates fail closed"
    ),
    bare_grid_semantics = find_test(
      "tests/testthat/test-functional-template-fit.R",
      "a bare MNI grid is support, not functional correspondence"
    )
  ),
  statistical_evidence = list(
    frozen_protocol = "inst/validation/functional-alignment-certification/protocol.json",
    superseded_v3_protocol = paste0(
      "inst/validation/functional-alignment-certification/protocol-v3.json"
    ),
    v3_supersession = paste0(
      "inst/validation/functional-alignment-certification/",
      "protocol-v3-supersession.json"
    ),
    invalidated_v4_protocol = paste0(
      "inst/validation/functional-alignment-certification/protocol-v4.json"
    ),
    v4_invalidation = paste0(
      "inst/validation/functional-alignment-certification/",
      "protocol-v4-invalidation.json"
    ),
    superseded_v5_protocol = paste0(
      "inst/validation/functional-alignment-certification/protocol-v5.json"
    ),
    v5_supersession = paste0(
      "inst/validation/functional-alignment-certification/",
      "protocol-v5-supersession.json"
    ),
    superseded_v6_protocol = paste0(
      "inst/validation/functional-alignment-certification/protocol-v6.json"
    ),
    v6_supersession = paste0(
      "inst/validation/functional-alignment-certification/",
      "protocol-v6-supersession.json"
    ),
    superseded_v2_protocol = paste0(
      "inst/validation/functional-alignment-certification/protocol-v2.json"
    ),
    v2_invalidation = paste0(
      "inst/validation/functional-alignment-certification/",
      "court-v2-invalidation.json"
    ),
    v2_sinkhorn_fixture = paste0(
      "inst/validation/functional-alignment-certification/",
      "court-v2-failing-sinkhorn-fixture.json"
    ),
    v2_source_test_manifest = dkfa_relative_path(
      v2_historical_paths[["source_test_manifest"]], root
    ),
    v2_source_test_results = dkfa_relative_path(
      v2_historical_paths[["source_test_results"]], root
    ),
    v2_source_test_session_info = dkfa_relative_path(
      v2_historical_paths[["source_test_session_info"]], root
    ),
    v2_failed_court_formal_budget = dkfa_relative_path(
      v2_historical_paths[["failed_court_formal_budget"]], root
    ),
    formal_court = file.path("data-raw", "functional-alignment-certification",
                            basename(court_dir), "court-verdict.json"),
    formal_court_run_state = file.path(
      "data-raw", "functional-alignment-certification", basename(court_dir),
      "court-run-state.json"
    ),
    formal_court_harness = file.path(
      "data-raw", "functional-alignment-certification", basename(court_dir),
      "court-harness-manifest.json"
    ),
    known_warp = file.path("data-raw", "functional-alignment-certification",
                          basename(power_dir), "known-warp-power-summary.csv"),
    known_warp_run_state = file.path(
      "data-raw", "functional-alignment-certification", basename(power_dir),
      "known-warp-power-run-state.json"
    ),
    known_warp_harness = file.path(
      "data-raw", "functional-alignment-certification", basename(power_dir),
      "known-warp-power-harness-manifest.json"
    ),
    source_test_manifest = file.path(
      "data-raw", "functional-alignment-certification",
      basename(source_test_dir), "source-test-manifest.json"
    ),
    source_test_results = file.path(
      "data-raw", "functional-alignment-certification",
      basename(source_test_dir), "source-test-results.csv"
    ),
    source_test_session_info = file.path(
      "data-raw", "functional-alignment-certification",
      basename(source_test_dir), "session-info.txt"
    ),
    pkgdown_receipt = file.path(pkgdown_dir, ".dkge-pkgdown-receipt.json"),
    independent_review = file.path("data-raw",
                                   "functional-alignment-certification",
                                   basename(review_path))
  )
)

hash_record <- function(path, role) {
  data.frame(
    role = role,
    path = path,
    bytes = unname(file.info(path)$size),
    sha256 = fa_court_hash_file(path),
    stringsAsFactors = FALSE
  )
}

artifact_group <- function(paths, role) {
  data.frame(path = paths, role = rep(role, length(paths)),
             stringsAsFactors = FALSE)
}
artifact_spec <- do.call(rbind, list(
  artifact_group(tarball, "built_tarball"),
  artifact_group(check_log, "r_cmd_check_log"),
  artifact_group(protocol_path, "frozen_protocol"),
  artifact_group(protocol_v2_path, "superseded_v2_protocol"),
  artifact_group(protocol_v3_path, "superseded_v3_protocol"),
  artifact_group(v3_supersession_path, "v3_supersession"),
  artifact_group(unname(v4_history_paths), "v4_invalidated_history"),
  artifact_group(unname(v5_history_paths), "v5_superseded_history"),
  artifact_group(unname(v6_history_paths), "v6_superseded_history"),
  artifact_group(unname(v7_history_paths), "v7_superseded_history"),
  artifact_group(unname(v8_history_paths), "v8_superseded_history"),
  artifact_group(v2_invalidation_path, "v2_court_invalidation"),
  artifact_group(sinkhorn_fixture_path, "v2_sinkhorn_fixture"),
  artifact_group(unname(v2_historical_paths), "v2_historical_evidence"),
  artifact_group(review_path, "independent_review"),
  artifact_group(file.path(court_dir, c(
    "court-verdict.json", "court-verdict.md", "court-audit.csv",
    "court-focus-summary.csv", "court-harness-manifest.json",
    "court-run-state.json",
    "factor-uncertainty.csv", "formal-budget.csv", "formal-grid.csv",
    "formal-inflation-flags.csv", "formal-manifest.txt", "formal-raw.csv",
    "formal-summary.csv", "latent-span-excess-model.csv",
    "latent-span-logistic-model.csv", "latent-span-model.csv"
  )), "formal_court"),
  artifact_group(file.path(power_dir, c(
    "known-warp-power-harness-manifest.json",
    "known-warp-power-manifest.json", "known-warp-power-raw.csv",
    "known-warp-power-run-state.json", "known-warp-power-summary.csv",
    "serial-parallel-determinism.json"
  )), "known_warp"),
  artifact_group(file.path(source_test_dir, c(
    "source-test-manifest.json", "source-test-results.csv", "session-info.txt"
  )), "source_test_suite"),
  artifact_group(file.path(pkgdown_dir, required_site_files), "pkgdown_page"),
  artifact_group(file.path(pkgdown_dir, ".dkge-pkgdown-receipt.json"),
                 "pkgdown_receipt")
))
artifact_paths <- artifact_spec$path
if (any(!file.exists(artifact_paths))) {
  stop("An expected certification artifact disappeared during collection.",
       call. = FALSE)
}
inventory <- do.call(
  rbind,
  Map(hash_record, artifact_spec$path, artifact_spec$role)
)

git_output <- function(...) {
  out <- suppressWarnings(system2("git", c("-C", root, ...),
                                  stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0L) {
    return(NA_character_)
  }
  paste(out, collapse = "\n")
}
git_status <- git_output("status", "--porcelain")

final_dir <- file.path(
  evidence_root, paste0(
    "final-", substr(source_hash, 1L, 16L), "-",
    substr(tarball_sha256, 1L, 12L)
  )
)
if (dir.exists(final_dir)) {
  stop("Final certification directory already exists: ", final_dir,
       call. = FALSE)
}
staging_dir <- tempfile(
  paste0(".", basename(final_dir), "-staging-"), tmpdir = evidence_root
)
dir.create(staging_dir, recursive = TRUE)
on.exit(if (dir.exists(staging_dir)) unlink(staging_dir, recursive = TRUE),
        add = TRUE)
utils::write.csv(inventory, file.path(staging_dir, "artifact-inventory.csv"),
                 row.names = FALSE)
jsonlite::write_json(evidence_map, file.path(staging_dir, "evidence-map.json"),
                     pretty = TRUE, auto_unbox = TRUE)

manifest <- list(
  protocol_id = protocol$protocol_id,
  protocol_amendment = list(
    passed = protocol_amendment$passed,
    protocol_v2_sha256 = protocol_amendment$protocol_v2_sha256,
    protocol_v3_sha256 = protocol_amendment$protocol_v3_sha256,
    protocol_v3_supersession_sha256 =
      protocol_amendment$protocol_v3_supersession_sha256,
    protocol_v4_sha256 = protocol_amendment$protocol_v4_sha256,
    protocol_v4_invalidation_sha256 =
      protocol_amendment$protocol_v4_invalidation_sha256,
    v4_history_bundle_sha256 = v4_history_bundle$bundle_sha256,
    v4_history_hashes = as.list(v4_history_bundle$hashes),
    protocol_v5_sha256 = protocol_amendment$protocol_v5_sha256,
    protocol_v5_supersession_sha256 =
      protocol_amendment$protocol_v5_supersession_sha256,
    v5_history_bundle_sha256 = v5_history_bundle$bundle_sha256,
    v5_history_hashes = as.list(v5_history_bundle$hashes),
    protocol_v6_sha256 = protocol_amendment$protocol_v6_sha256,
    protocol_v6_supersession_sha256 =
      protocol_amendment$protocol_v6_supersession_sha256,
    v6_history_bundle_sha256 = v6_history_bundle$bundle_sha256,
    v6_history_hashes = as.list(v6_history_bundle$hashes),
    protocol_v7_sha256 = protocol_amendment$protocol_v7_sha256,
    protocol_v7_supersession_sha256 =
      protocol_amendment$protocol_v7_supersession_sha256,
    v7_history_bundle_sha256 = v7_history_bundle$bundle_sha256,
    v7_history_hashes = as.list(v7_history_bundle$hashes),
    protocol_v8_sha256 = protocol_amendment$protocol_v8_sha256,
    protocol_v8_supersession_sha256 =
      protocol_amendment$protocol_v8_supersession_sha256,
    v8_history_bundle_sha256 = v8_history_bundle$bundle_sha256,
    v8_history_hashes = as.list(v8_history_bundle$hashes),
    protocol_v9_sha256 = protocol_amendment$protocol_v9_sha256,
    v2_invalidation_sha256 = dkfa_hash_file(v2_invalidation_path),
    v2_sinkhorn_fixture_sha256 = dkfa_hash_file(sinkhorn_fixture_path),
    v2_historical_evidence_bundle_sha256 =
      v2_historical_bundle$bundle_sha256,
    v2_historical_evidence_hashes =
      as.list(v2_historical_bundle$hashes),
    v3_history_bundle_sha256 = v3_history_bundle$bundle_sha256,
    v3_history_hashes = as.list(v3_history_bundle$hashes),
    collector_loaded_input_snapshot_sha256 =
      loaded_binding$snapshot_sha256,
    inherited_numerical_amendment =
      protocol$amendment$inherited_numerical_amendment,
    protocol_change = protocol$amendment$protocol_change,
    no_second_numerical_amendment =
      protocol$amendment$no_second_numerical_amendment
  ),
  status = "complete_approximate_only",
  inferential_promotion = FALSE,
  claim = paste(
    "The exact comparator and package-integrity gates passed, but the frozen",
    "negative-control promotion gate failed. Functional alignment is available",
    "only through the fail-closed approximate workflow with explicit override."
  ),
  source_tree_sha256 = source_hash,
  test_validation_tree_sha256 = test_validation_hash,
  tarball_source_tree_sha256 = tarball_source_hash,
  tarball_raw_source_tree_sha256 = tarball_raw_source_hash,
  package_source_projection_sha256 =
    source_provenance$package_projection_sha256,
  tarball_shipped_payload_tree_sha256 =
    source_provenance$shipped_payload_tree_sha256,
  tarball_sha256 = tarball_sha256,
  live_docs_input_sha256 = current_docs,
  tarball_docs_input_sha256 = tarball_docs,
  docs_overlay_bundle_sha256 = docs_provenance$overlay_bundle_sha256,
  tarball_package_tree_sha256 = tarball_tree$tree_sha256,
  checked_package_tree_sha256 = check_tree$tree_sha256,
  pkgdown_site_content_sha256 = pkgdown_status$site_content_sha256,
  pkgdown_builder_bundle_sha256 = pkgdown_builder$bundle_sha256,
  court_harness_bundle_sha256 = court_bundle$bundle_sha256,
  known_warp_harness_bundle_sha256 = power_bundle$bundle_sha256,
  review_evidence_bundle_sha256 = review_evidence_bundle_sha256,
  review_evidence_hashes = as.list(review_evidence_hashes),
  collector_sha256 = dkfa_hash_file(script),
  protocol_sha256 = protocol_hash,
  git_head = git_output("rev-parse", "HEAD"),
  git_branch = git_output("branch", "--show-current"),
  git_worktree_dirty = !is.na(git_status) && nzchar(git_status),
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  gates = gates,
  court = list(
    run_id = court$run_id,
    run_state = court_run_state$status,
    run_state_sha256 = dkfa_hash_file(
      file.path(court_dir, "court-run-state.json")
    ),
    status = court$status,
    exact_oracle_gate = court$exact_oracle_gate,
    negative_control_promotion_gate = court$negative_control_promotion_gate,
    rows = audit_value("raw_rows"),
    cohorts = audit_value("unique_cohorts"),
    minimum_convergence = audit_value("minimum_cell_arm_convergence")
  ),
  known_warp = list(
    run_id = power_manifest$run_id,
    run_state = power_run_state$status,
    run_state_sha256 = dkfa_hash_file(
      file.path(power_dir, "known-warp-power-run-state.json")
    ),
    cohorts = length(unique(power_raw$seed)),
    rows = nrow(power_raw),
    workers = power_manifest$workers,
    n_perm = power_manifest$n_perm,
    runtime_seconds = power_manifest$runtime_seconds,
    solver_converged = all(power_raw$solver_converged),
    template_convergence = mean(template_convergence_by_seed),
    summary = unname(split(power_summary, seq_len(nrow(power_summary))))
  ),
  source_tests = source_tests,
  source_test_evidence = list(
    checks = as.list(source_test_evidence$checks),
    counts = as.list(source_test_evidence$counts),
    hashes = as.list(source_test_evidence$hashes)
  ),
  package = list(
    tarball = basename(tarball),
    check_log_sha256 = fa_court_hash_file(check_log),
    check_required_options = as.list(required_check_options),
    check_allowed_options = as.list(allowed_check_options),
    check_recorded_options = as.list(check_status$recorded_options),
    check_policy_passed = check_policy_ok,
    checked_package_tree_sha256 = check_tree$tree_sha256,
    pkgdown_required_pages = required_site_files,
    pkgdown_receipt = ".dkge-pkgdown-receipt.json",
    pkgdown_site_content_sha256 = pkgdown_status$site_content_sha256,
    pkgdown_utf8_locale = pkgdown_status$receipt$utf8_locale,
    pkgdown_locale = pkgdown_status$receipt$locale
  ),
  independent_review = review,
  environment = list(R = R.version.string, platform = R.version$platform),
  scope = list(
    certified = "the exact source-tree and tarball fingerprints recorded here",
    excluded = c("merge", "tag", "hosted publication", "public release")
  )
)
jsonlite::write_json(
  manifest, file.path(staging_dir, "certification-manifest.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = 16
)

baseline <- court_focus[
  court_focus$cell == "baseline" & court_focus$arm %in% c(
    "full_permutation_reestimation", "geometry_only",
    "independent_alignment", "kernel_image_residual_prototype"
  ), c("arm", "fwer_rate", "fwer_ci_lower", "fwer_ci_upper")
]
baseline_lines <- apply(baseline, 1L, function(x) {
  sprintf("- `%s`: FWER %.3f (Wilson 95%% CI %.3f-%.3f)",
          x[["arm"]], as.numeric(x[["fwer_rate"]]),
          as.numeric(x[["fwer_ci_lower"]]),
          as.numeric(x[["fwer_ci_upper"]]))
})
power_lines <- apply(power_summary, 1L, function(x) {
  sprintf(paste0("- `%s`: power %.3f (Wilson 95%% CI %.3f-%.3f), ",
                 "correlation %.3f, RMSE %.3f, amplitude ratio %.3f, ",
                 "point spread %.3f"),
          x[["arm"]], as.numeric(x[["power"]]),
          as.numeric(x[["power_ci_lower"]]),
          as.numeric(x[["power_ci_upper"]]),
          as.numeric(x[["mean_correlation"]]), as.numeric(x[["mean_rmse"]]),
          as.numeric(x[["mean_amplitude_ratio"]]),
          as.numeric(x[["mean_point_spread"]]))
})

report <- c(
  "# Functional-alignment certification",
  "",
  paste0("- Status: `", manifest$status, "`"),
  paste0("- Source-tree SHA-256: `", source_hash, "`"),
  paste0("- Tarball SHA-256: `", manifest$tarball_sha256, "`"),
  paste0("- Protocol SHA-256: `", protocol_hash, "`"),
  "- Inferential promotion: **blocked by the frozen court**",
  "- Scope: exact candidate only; no merge, tag, hosted publication, or release claim",
  "",
  "## Gate result",
  "",
  paste(
    "All provenance, deterministic-execution, algebraic/regression, source-test,",
    "built-tarball, documentation, and independent-review gates passed. The exact",
    "full-pipeline re-estimation exact null-action comparator passed. The",
    "negative-control promotion gate did",
    "not pass, so aligned inference remains explicitly approximate and fail-closed."
  ),
  "",
  "## Formal court: baseline",
  "",
  paste0(
    "The superseded v2 protocol, failed-run budget, source-test manifest, ",
    "results, and session receipt are retained under evidence bundle `",
    v2_historical_bundle$bundle_sha256, "`. V3 changed only the Sinkhorn ",
    "iteration ceiling and its directory label. V3 is preserved under history ",
    "bundle `", v3_history_bundle$bundle_sha256, "`. V4 completed computation ",
    "after its candidate drifted, but its analyzer failed before reading ",
    "statistical outputs; the complete inadmissible history is bound under `",
    v4_history_bundle$bundle_sha256, "`. V5 and V6 were superseded before ",
    "execution under `", v5_history_bundle$bundle_sha256, "` and `",
    v6_history_bundle$bundle_sha256, "`. V7 completed source, court, and power ",
    "execution but was superseded before collection after the exact package ",
    "check exposed nonstatistical documentation and test-harness defects; its ",
    "protocol and supersession are bound under `",
    v7_history_bundle$bundle_sha256, "`. V8 completed source, court, power, ",
    "and exact package-check execution but failed the separately required ",
    "pkgdown gate before collection; every V8 output and the failure ",
    "supersession are bound under `", v8_history_bundle$bundle_sha256,
    "`. V9 changes only the 19 vignette comparator namespaces, validation ",
    "and supersession provenance, and the output-directory label."
  ),
  "",
  baseline_lines,
  "",
  "The complete factor grid, uncertainty intervals, adverse cells, and latent-span",
  "models are retained in the court run directory; no observed result was removed",
  "and no threshold was changed after inspection.",
  "",
  "## Held-out known-warp efficacy",
  "",
  power_lines,
  "",
  "These power and recovery metrics are descriptive. They cannot override the",
  "failed Type-I promotion gate.",
  paste0(
    "Template outer-loop convergence was ",
    sum(template_convergence_by_seed), "/",
    length(template_convergence_by_seed),
    "; every scheduled cohort, including non-convergence, is retained."
  ),
  "",
  "## Package evidence",
  "",
  paste0("- Full source suite: ", source_tests$expectations,
         " expectations, 0 failures; ", source_tests$warnings,
         " warning events recorded with generic classifications in the ",
         "source-test manifest; ", source_tests$skipped,
         " manifest-counted skips"),
  "- Known-warp serial/parallel comparison: identical within 1e-12",
  "- Built-tarball R CMD check: `Status: OK`",
  "- R CMD check policy: vignettes excluded only through `--ignore-vignettes`; installation, tests, examples, and codoc retained",
  "- Checked package source is byte-bound to the tarball contents",
  paste0("- pkgdown required pages: ", length(required_site_files), "/",
         length(required_site_files),
         " present and byte-bound to the exact tarball"),
  paste0("- Independent review verdict: `", review$verdict, "` with 0 blockers"),
  "",
  "See `certification-manifest.json`, `artifact-inventory.csv`, and",
  "`evidence-map.json` for the machine-readable record."
)
writeLines(report, file.path(staging_dir, "certification-report.md"),
           useBytes = TRUE)

# Close the collection-time provenance window immediately before the atomic
# publish. Every mutable live input used above must still have exactly the bytes
# that were reviewed and recorded.
final_source_hash <- dkfa_certified_source_tree_hash(root)
final_test_validation_hash <- dkfa_test_validation_manifest(root)$tree_sha256
final_docs_hash <- dkfa_pkgdown_input_manifest(root)$tree_sha256
final_check_tree <- dkfa_tree_manifest(
  check_source_root, roots = ".", exclude = compiled_excludes
)$tree_sha256
final_site_hash <- dkfa_site_manifest(pkgdown_dir)$tree_sha256
final_inventory <- do.call(
  rbind,
  Map(hash_record, artifact_spec$path, artifact_spec$role)
)
final_review_evidence_hash <-
  dkfa_review_evidence_bundle(review_evidence_paths)$bundle_sha256
final_input_snapshot <- input_snapshot()
dkfa_assert_input_snapshot(
  loaded_binding$snapshot, final_input_snapshot,
  "Certification collection"
)
final_checks <- c(
  loaded_inputs = identical(final_input_snapshot, loaded_binding$snapshot),
  source_tree = identical(final_source_hash, source_hash),
  test_validation_tree = identical(
    final_test_validation_hash, test_validation_hash
  ),
  docs_input = identical(final_docs_hash, current_docs),
  checked_package_tree = identical(final_check_tree, check_tree$tree_sha256),
  pkgdown_site = identical(final_site_hash, pkgdown_status$site_content_sha256),
  artifact_inventory = identical(final_inventory, inventory),
  review_evidence = identical(
    final_review_evidence_hash, review_evidence_bundle_sha256
  )
)
if (!all(final_checks)) {
  stop(
    "Certification inputs changed during collection: ",
    paste(names(final_checks)[!final_checks], collapse = ", "),
    call. = FALSE
  )
}
if (!file.rename(staging_dir, final_dir)) {
  stop("Could not atomically publish the final certification directory.",
       call. = FALSE)
}
cat(final_dir, "\n")
