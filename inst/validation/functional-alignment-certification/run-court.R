#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
} else {
  normalizePath(
    "inst/validation/functional-alignment-certification/run-court.R",
    mustWork = TRUE
  )
}
root <- normalizePath(file.path(dirname(script), "..", "..", ".."),
                      mustWork = TRUE)
utils_script <- file.path(dirname(script), "evidence-utils.R")
court_script <- file.path(
  root, "inst", "validation", "functional-alignment-court", "court.R"
)
protocol <- file.path(dirname(script), "protocol.json")
protocol_v2 <- file.path(dirname(script), "protocol-v2.json")
v2_invalidation <- file.path(dirname(script), "court-v2-invalidation.json")
sinkhorn_fixture <- file.path(
  dirname(script), "court-v2-failing-sinkhorn-fixture.json"
)
analysis_script <- file.path(dirname(script), "analyze-court.R")
script_env <- environment()
source(utils_script, local = TRUE)
input_snapshot <- function() {
  list(
    source_tree_sha256 = dkfa_candidate_source_tree_hash(root),
    harness = dkfa_bundle_manifest(dkfa_court_harness_paths(root), root)
  )
}
loaded_binding <- dkfa_bind_loaded_inputs(
  input_snapshot,
  function() {
    source(utils_script, local = script_env)
    source(court_script, local = script_env)
    dkfa_load_source_package(root, quiet = TRUE)
  },
  "Court loader"
)
source_loader <- loaded_binding$value
source_loader$bound_input_snapshot_sha256 <-
  loaded_binding$snapshot_sha256
source_tree_hash <- loaded_binding$snapshot$source_tree_sha256
harness <- loaded_binding$snapshot$harness
v2_historical_paths <- dkfa_v2_historical_paths(root)
v2_historical_bundle <- dkfa_review_evidence_bundle(v2_historical_paths)

# V3's sole numerical amendment is inherited unchanged by v9 and remains
# scoped to this certification wrapper. The historical v1/v2 helper stays
# byte-semantically fixed at 5000.
fa_court_mapper_v2 <- fa_court_mapper
fa_court_mapper <- function(cfg, geometry_only = FALSE) {
  mapper <- fa_court_mapper_v2(cfg, geometry_only = geometry_only)
  mapper$params$max_iter <- 50000L
  mapper
}

protocol_amendment <- dkfa_validate_current_protocol(
  root, candidate_source_tree_sha256 = source_tree_hash
)
v2_invalidation_status <- dkfa_validate_court_v2_invalidation(
  v2_invalidation, sinkhorn_fixture, protocol_v2,
  historical_paths = v2_historical_paths,
  verify_solver = TRUE, verify_origin = TRUE
)
if (!isTRUE(protocol_amendment$passed) ||
    !isTRUE(v2_invalidation_status$passed)) {
  stop(
    "The v9 protocol history or retained v2 invalidation chain is invalid.",
    call. = FALSE
  )
}

run_id <- paste0(
  "court-", substr(source_tree_hash, 1L, 16L), "-",
  substr(harness$bundle_sha256, 1L, 12L)
)
run_dir <- file.path(
  root, "data-raw", "functional-alignment-certification", run_id
)
run_claim <- dkfa_claim_run_directory(
  run_dir, harness,
  metadata = list(
    protocol_id = "dkge-functional-alignment-certification-v9",
    run_id = run_id,
    source_tree_sha256 = source_tree_hash,
    protocol_sha256 = protocol_amendment$protocol_v9_sha256,
    protocol_v2_sha256 = protocol_amendment$protocol_v2_sha256,
    court_v2_invalidation_sha256 = dkfa_hash_file(v2_invalidation),
    court_v2_sinkhorn_fixture_sha256 = dkfa_hash_file(sinkhorn_fixture),
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
    protocol_v9_sha256 = protocol_amendment$protocol_v9_sha256,
    protocol_v3_sha256 = protocol_amendment$protocol_v3_sha256,
    protocol_v3_supersession_sha256 =
      protocol_amendment$protocol_v3_supersession_sha256,
    court_harness_bundle_sha256 = harness$bundle_sha256,
    v2_historical_evidence_bundle_sha256 =
      v2_historical_bundle$bundle_sha256,
    loaded_input_snapshot_sha256 = loaded_binding$snapshot_sha256
  )
)

fa_court_output_dir <- function() run_dir
fa_court_protocol_path <- function() protocol
Sys.setenv(DKGE_COURT_CORES = "4")

manifest_path <- file.path(run_dir, "formal-manifest.txt")
budget <- data.frame(
  n_sim = 40L,
  n_perm = 31L,
  predicted_seconds = NA_real_,
  seconds_per_exact_perm = NA_real_,
  resource_limited = TRUE,
  exact_cells = 6L,
  formal_cores = 4L
)
extra <- c(
  paste0("certification_protocol_id=dkge-functional-alignment-certification-v9"),
  paste0("certification_wrapper_sha256=", fa_court_hash_file(script)),
  paste0("historical_court_script_sha256=", fa_court_hash_file(court_script)),
  paste0("court_harness_bundle_sha256=", harness$bundle_sha256),
  paste0("protocol_v2_sha256=", protocol_amendment$protocol_v2_sha256),
  paste0("protocol_v3_sha256=", protocol_amendment$protocol_v3_sha256),
  paste0(
    "protocol_v3_supersession_sha256=",
    protocol_amendment$protocol_v3_supersession_sha256
  ),
  paste0("protocol_v4_sha256=", protocol_amendment$protocol_v4_sha256),
  paste0(
    "protocol_v4_invalidation_sha256=",
    protocol_amendment$protocol_v4_invalidation_sha256
  ),
  paste0(
    "v4_history_bundle_sha256=",
    protocol_amendment$v4_history_bundle$bundle_sha256
  ),
  paste0("protocol_v5_sha256=", protocol_amendment$protocol_v5_sha256),
  paste0(
    "protocol_v5_supersession_sha256=",
    protocol_amendment$protocol_v5_supersession_sha256
  ),
  paste0(
    "v5_history_bundle_sha256=",
    protocol_amendment$v5_history_bundle$bundle_sha256
  ),
  paste0("protocol_v6_sha256=", protocol_amendment$protocol_v6_sha256),
  paste0(
    "protocol_v6_supersession_sha256=",
    protocol_amendment$protocol_v6_supersession_sha256
  ),
  paste0(
    "v6_history_bundle_sha256=",
    protocol_amendment$v6_history_bundle$bundle_sha256
  ),
  paste0("protocol_v7_sha256=", protocol_amendment$protocol_v7_sha256),
  paste0(
    "protocol_v7_supersession_sha256=",
    protocol_amendment$protocol_v7_supersession_sha256
  ),
  paste0(
    "v7_history_bundle_sha256=",
    protocol_amendment$v7_history_bundle$bundle_sha256
  ),
  paste0("protocol_v8_sha256=", protocol_amendment$protocol_v8_sha256),
  paste0(
    "protocol_v8_supersession_sha256=",
    protocol_amendment$protocol_v8_supersession_sha256
  ),
  paste0(
    "v8_history_bundle_sha256=",
    protocol_amendment$v8_history_bundle$bundle_sha256
  ),
  paste0("protocol_v9_sha256=", protocol_amendment$protocol_v9_sha256),
  paste0("court_v2_invalidation_sha256=", dkfa_hash_file(v2_invalidation)),
  paste0("court_v2_sinkhorn_fixture_sha256=", dkfa_hash_file(sinkhorn_fixture)),
  paste0(
    "v2_historical_evidence_bundle_sha256=",
    v2_historical_bundle$bundle_sha256
  ),
  paste0(
    "v2_source_test_manifest_sha256=",
    v2_historical_bundle$hashes[["source_test_manifest"]]
  ),
  paste0(
    "v2_source_test_results_sha256=",
    v2_historical_bundle$hashes[["source_test_results"]]
  ),
  paste0(
    "v2_source_test_session_info_sha256=",
    v2_historical_bundle$hashes[["source_test_session_info"]]
  ),
  paste0(
    "v2_failed_court_formal_budget_sha256=",
    v2_historical_bundle$hashes[["failed_court_formal_budget"]]
  ),
  paste0("sinkhorn_max_iter=50000"),
  paste0("v3_only_numerical_amendment_verified=TRUE"),
  paste0("v9_no_numerical_or_statistical_amendment_verified=TRUE"),
  paste0("source_loader=", source_loader$method),
  paste0("source_loader_version=", source_loader$loader_version),
  paste0("source_loader_compile=", source_loader$compile),
  paste0(
    "loaded_input_snapshot_sha256=",
    loaded_binding$snapshot_sha256
  ),
  paste0("run_id=", run_id),
  paste0("source_tree_hash_at_start=", source_tree_hash)
)
result <- tryCatch({
  utils::write.csv(
    budget, file.path(run_dir, "formal-budget.csv"), row.names = FALSE
  )
  value <- fa_court_run_formal()
  writeLines(
    c(readLines(manifest_path, warn = FALSE), extra),
    manifest_path,
    useBytes = TRUE
  )
  formal_outputs <- c(
    formal_budget = file.path(run_dir, "formal-budget.csv"),
    formal_grid = file.path(run_dir, "formal-grid.csv"),
    formal_raw = file.path(run_dir, "formal-raw.csv"),
    formal_summary = file.path(run_dir, "formal-summary.csv"),
    formal_inflation_flags = file.path(
      run_dir, "formal-inflation-flags.csv"
    ),
    latent_span_model = file.path(run_dir, "latent-span-model.csv"),
    formal_manifest = manifest_path
  )
  if (!all(file.exists(formal_outputs))) {
    stop("The formal court returned without all required outputs.",
         call. = FALSE)
  }
  dkfa_assert_input_snapshot(
    loaded_binding$snapshot, input_snapshot(), "Court run"
  )
  dkfa_update_run_state(
    run_claim$state_path, "formal_complete",
    outputs = as.list(vapply(formal_outputs, dkfa_hash_file, character(1)))
  )
  value
}, error = function(e) {
  current <- tryCatch(
    jsonlite::read_json(run_claim$state_path, simplifyVector = TRUE),
    error = function(...) NULL
  )
  if (!is.null(current) && identical(current$status, "running")) {
    try(dkfa_update_run_state(
      run_claim$state_path, "failed",
      condition = list(
        classes = as.list(class(e)),
        message = conditionMessage(e)
      )
    ), silent = TRUE)
  }
  stop(e)
})
cat(run_dir, "\n")
invisible(result)
