#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
} else {
  normalizePath(
    paste0(
      "inst/validation/functional-alignment-certification/",
      "verify-known-warp-determinism.R"
    ),
    mustWork = TRUE
  )
}
root <- normalizePath(file.path(dirname(script), "..", "..", ".."),
                      mustWork = TRUE)
utils_script <- file.path(dirname(script), "evidence-utils.R")
template_script <- file.path(
  root, "inst", "validation", "functional-alignment-template",
  "template-benefit.R"
)
power_script <- file.path(dirname(script), "known-warp-power.R")
court_script <- file.path(
  root, "inst", "validation", "functional-alignment-court", "court.R"
)
protocol_path <- file.path(dirname(script), "protocol.json")
protocol_v2_path <- file.path(dirname(script), "protocol-v2.json")
v2_invalidation_path <- file.path(
  dirname(script), "court-v2-invalidation.json"
)
sinkhorn_fixture_path <- file.path(
  dirname(script), "court-v2-failing-sinkhorn-fixture.json"
)
runner_script <- file.path(dirname(script), "run-known-warp-power.R")
script_env <- environment()
source(utils_script, local = TRUE)
input_snapshot <- function() {
  list(
    source_tree_sha256 = dkfa_candidate_source_tree_hash(root),
    harness = dkfa_bundle_manifest(dkfa_power_harness_paths(root), root)
  )
}
loaded_binding <- dkfa_bind_loaded_inputs(
  input_snapshot,
  function() {
    source(utils_script, local = script_env)
    source(template_script, local = script_env)
    source(power_script, local = script_env)
    source(court_script, local = script_env)
    dkfa_load_source_package(root, quiet = TRUE)
  },
  "Known-warp determinism loader"
)
source_loader <- loaded_binding$value
source_loader$bound_input_snapshot_sha256 <-
  loaded_binding$snapshot_sha256
source_hash <- loaded_binding$snapshot$source_tree_sha256
harness <- loaded_binding$snapshot$harness
v2_historical_paths <- dkfa_v2_historical_paths(root)
v2_historical_bundle <- dkfa_review_evidence_bundle(v2_historical_paths)
protocol_amendment <- dkfa_validate_current_protocol(
  root, candidate_source_tree_sha256 = source_hash
)
protocol <- jsonlite::read_json(protocol_path, simplifyVector = TRUE)
v2_invalidation_status <- dkfa_validate_court_v2_invalidation(
  v2_invalidation_path, sinkhorn_fixture_path, protocol_v2_path,
  historical_paths = v2_historical_paths,
  verify_solver = FALSE, verify_origin = TRUE
)
if (!isTRUE(protocol_amendment$passed) ||
    !isTRUE(v2_invalidation_status$passed)) {
  stop(
    "The v9 protocol history or retained v2 invalidation chain is invalid.",
    call. = FALSE
  )
}

evidence_root <- file.path(
  root, "data-raw", "functional-alignment-certification"
)
candidates <- list.dirs(evidence_root, recursive = FALSE, full.names = TRUE)
protocol_hash <- dkfa_hash_file(protocol_path)
prefix <- paste0(
  "known-warp-", substr(source_hash, 1L, 12L), "-",
  substr(protocol_hash, 1L, 12L), "-",
  substr(harness$bundle_sha256, 1L, 12L)
)
candidates <- candidates[basename(candidates) == prefix]
candidates <- candidates[file.exists(file.path(
  candidates, "known-warp-power-raw.csv"
))]
if (!length(candidates)) stop("No complete known-warp power run.", call. = FALSE)
run_dir <- candidates[[which.max(file.info(candidates)$mtime)]]
run_state_path <- file.path(run_dir, "known-warp-power-run-state.json")
eligible_state <- jsonlite::read_json(
  run_state_path, simplifyVector = TRUE
)
if (!identical(eligible_state$status, "power_complete")) {
  stop(
    "Determinism verification requires a power-complete run and cannot rerun ",
    "a failed or terminal verification.",
    call. = FALSE
  )
}
verification_passed <- dkfa_with_run_failure_receipt(
  run_state_path, "power_complete", {
power_state_paths <- c(
  power_manifest = file.path(run_dir, "known-warp-power-manifest.json"),
  power_raw = file.path(run_dir, "known-warp-power-raw.csv"),
  power_summary = file.path(run_dir, "known-warp-power-summary.csv")
)
power_state_validation <- dkfa_validate_run_state_outputs(
  run_state_path, "power_complete",
  list(power_complete = power_state_paths)
)
if (!isTRUE(power_state_validation$passed)) {
  stop(
    "Determinism verification requires an intact power-complete receipt: ",
    paste(
      names(power_state_validation$checks)[!power_state_validation$checks],
      collapse = ", "
    ),
    call. = FALSE
  )
}
run_manifest <- jsonlite::read_json(
  file.path(run_dir, "known-warp-power-manifest.json"),
  simplifyVector = TRUE
)
run_state <- power_state_validation$state
binding_checks <- c(
  state_run_id = identical(run_state$run_id, basename(run_dir)),
  state_source = identical(run_state$source_tree_sha256, source_hash),
  state_protocol = identical(run_state$protocol_sha256, protocol_hash),
  state_loaded_inputs = identical(
    run_state$loaded_input_snapshot_sha256,
    run_manifest$loaded_input_snapshot_sha256
  ),
  state_v2 = identical(
    run_state$protocol_v2_sha256,
    protocol_amendment$protocol_v2_sha256
  ),
  state_v2_invalidation = identical(
    run_state$court_v2_invalidation_sha256,
    dkfa_hash_file(v2_invalidation_path)
  ),
  state_v2_fixture = identical(
    run_state$court_v2_sinkhorn_fixture_sha256,
    dkfa_hash_file(sinkhorn_fixture_path)
  ),
  state_v2_history = identical(
    run_state$v2_historical_evidence_bundle_sha256,
    v2_historical_bundle$bundle_sha256
  ),
  state_v4 = identical(
    run_state$protocol_v4_sha256,
    protocol_amendment$protocol_v4_sha256
  ),
  state_v4_invalidation = identical(
    run_state$protocol_v4_invalidation_sha256,
    protocol_amendment$protocol_v4_invalidation_sha256
  ),
  state_v4_history = identical(
    run_state$v4_history_bundle_sha256,
    protocol_amendment$v4_history_bundle$bundle_sha256
  ),
  state_v5 = identical(
    run_state$protocol_v5_sha256,
    protocol_amendment$protocol_v5_sha256
  ),
  state_v5_supersession = identical(
    run_state$protocol_v5_supersession_sha256,
    protocol_amendment$protocol_v5_supersession_sha256
  ),
  state_v5_history = identical(
    run_state$v5_history_bundle_sha256,
    protocol_amendment$v5_history_bundle$bundle_sha256
  ),
  state_v6 = identical(
    run_state$protocol_v6_sha256,
    protocol_amendment$protocol_v6_sha256
  ),
  state_v6_supersession = identical(
    run_state$protocol_v6_supersession_sha256,
    protocol_amendment$protocol_v6_supersession_sha256
  ),
  state_v6_history = identical(
    run_state$v6_history_bundle_sha256,
    protocol_amendment$v6_history_bundle$bundle_sha256
  ),
  state_v7 = identical(
    run_state$protocol_v7_sha256,
    protocol_amendment$protocol_v7_sha256
  ),
  state_v7_supersession = identical(
    run_state$protocol_v7_supersession_sha256,
    protocol_amendment$protocol_v7_supersession_sha256
  ),
  state_v7_history = identical(
    run_state$v7_history_bundle_sha256,
    protocol_amendment$v7_history_bundle$bundle_sha256
  ),
  state_v8 = identical(
    run_state$protocol_v8_sha256,
    protocol_amendment$protocol_v8_sha256
  ),
  state_v8_supersession = identical(
    run_state$protocol_v8_supersession_sha256,
    protocol_amendment$protocol_v8_supersession_sha256
  ),
  state_v8_history = identical(
    run_state$v8_history_bundle_sha256,
    protocol_amendment$v8_history_bundle$bundle_sha256
  ),
  state_v9 = identical(
    run_state$protocol_v9_sha256,
    protocol_amendment$protocol_v9_sha256
  ),
  state_v3 = identical(
    run_state$protocol_v3_sha256,
    protocol_amendment$protocol_v3_sha256
  ),
  state_v3_supersession = identical(
    run_state$protocol_v3_supersession_sha256,
    protocol_amendment$protocol_v3_supersession_sha256
  ),
  manifest_run_id = identical(run_manifest$run_id, basename(run_dir)),
  manifest_source = identical(run_manifest$source_tree_sha256, source_hash),
  manifest_protocol = identical(run_manifest$protocol_sha256, protocol_hash),
  manifest_loaded_inputs = identical(
    run_manifest$loaded_input_snapshot_sha256,
    run_manifest$source_loader$bound_input_snapshot_sha256
  ),
  manifest_v2 = identical(
    run_manifest$protocol_v2_sha256,
    protocol_amendment$protocol_v2_sha256
  ),
  manifest_v2_invalidation = identical(
    run_manifest$court_v2_invalidation_sha256,
    dkfa_hash_file(v2_invalidation_path)
  ),
  manifest_v2_fixture = identical(
    run_manifest$court_v2_sinkhorn_fixture_sha256,
    dkfa_hash_file(sinkhorn_fixture_path)
  ),
  manifest_v2_history = identical(
    run_manifest$v2_historical_evidence_bundle_sha256,
    v2_historical_bundle$bundle_sha256
  ),
  manifest_v3 = identical(
    run_manifest$protocol_v3_sha256,
    protocol_amendment$protocol_v3_sha256
  ),
  manifest_v3_supersession = identical(
    run_manifest$protocol_v3_supersession_sha256,
    protocol_amendment$protocol_v3_supersession_sha256
  ),
  manifest_v4 = identical(
    run_manifest$protocol_v4_sha256,
    protocol_amendment$protocol_v4_sha256
  ),
  manifest_v4_invalidation = identical(
    run_manifest$protocol_v4_invalidation_sha256,
    protocol_amendment$protocol_v4_invalidation_sha256
  ),
  manifest_v4_history = identical(
    run_manifest$v4_history_bundle_sha256,
    protocol_amendment$v4_history_bundle$bundle_sha256
  ),
  manifest_v5 = identical(
    run_manifest$protocol_v5_sha256,
    protocol_amendment$protocol_v5_sha256
  ),
  manifest_v5_supersession = identical(
    run_manifest$protocol_v5_supersession_sha256,
    protocol_amendment$protocol_v5_supersession_sha256
  ),
  manifest_v5_history = identical(
    run_manifest$v5_history_bundle_sha256,
    protocol_amendment$v5_history_bundle$bundle_sha256
  ),
  manifest_v6 = identical(
    run_manifest$protocol_v6_sha256,
    protocol_amendment$protocol_v6_sha256
  ),
  manifest_v6_supersession = identical(
    run_manifest$protocol_v6_supersession_sha256,
    protocol_amendment$protocol_v6_supersession_sha256
  ),
  manifest_v6_history = identical(
    run_manifest$v6_history_bundle_sha256,
    protocol_amendment$v6_history_bundle$bundle_sha256
  ),
  manifest_v7 = identical(
    run_manifest$protocol_v7_sha256,
    protocol_amendment$protocol_v7_sha256
  ),
  manifest_v7_supersession = identical(
    run_manifest$protocol_v7_supersession_sha256,
    protocol_amendment$protocol_v7_supersession_sha256
  ),
  manifest_v7_history = identical(
    run_manifest$v7_history_bundle_sha256,
    protocol_amendment$v7_history_bundle$bundle_sha256
  ),
  manifest_v8 = identical(
    run_manifest$protocol_v8_sha256,
    protocol_amendment$protocol_v8_sha256
  ),
  manifest_v8_supersession = identical(
    run_manifest$protocol_v8_supersession_sha256,
    protocol_amendment$protocol_v8_supersession_sha256
  ),
  manifest_v8_history = identical(
    run_manifest$v8_history_bundle_sha256,
    protocol_amendment$v8_history_bundle$bundle_sha256
  ),
  manifest_v9 = identical(
    run_manifest$protocol_v9_sha256,
    protocol_amendment$protocol_v9_sha256
  )
)
if (!all(binding_checks)) {
  stop(
    "Known-warp run provenance is invalid: ",
    paste(names(binding_checks)[!binding_checks], collapse = ", "),
    call. = FALSE
  )
}
power_harness_record <- jsonlite::read_json(
  file.path(run_dir, "known-warp-power-harness-manifest.json"),
  simplifyVector = TRUE
)
dkfa_verify_bundle_manifest(
  power_harness_record,
  dkfa_power_harness_paths(root),
  root,
  "Known-warp run-state"
)
dkfa_verify_bundle_manifest(
  list(
    bundle_sha256 = run_manifest$harness_bundle_sha256,
    files = run_manifest$harness_files
  ),
  dkfa_power_harness_paths(root),
  root,
  "Known-warp"
)
observed <- utils::read.csv(
  file.path(run_dir, "known-warp-power-raw.csv"),
  stringsAsFactors = FALSE
)
recorded_summary <- utils::read.csv(
  file.path(run_dir, "known-warp-power-summary.csv"),
  stringsAsFactors = FALSE
)
power_evidence <- dkfa_validate_known_warp_evidence(
  observed, recorded_summary, protocol
)
if (!isTRUE(power_evidence$passed)) {
  stop(
    "Known-warp raw schedule or recomputed summary is invalid: ",
    paste(names(power_evidence$checks)[!power_evidence$checks],
          collapse = ", "),
    call. = FALSE
  )
}
seeds <- 91201:91202
serial <- do.call(rbind, lapply(seeds, dktf_power_run_one))
expected <- observed[observed$seed %in% seeds, names(serial), drop = FALSE]
serial <- serial[order(serial$seed, serial$arm), , drop = FALSE]
expected <- expected[order(expected$seed, expected$arm), , drop = FALSE]
rownames(serial) <- rownames(expected) <- NULL
compared <- setdiff(names(serial), "runtime_seconds")
comparison <- all.equal(
  serial[compared], expected[compared],
  tolerance = 1e-12, check.attributes = TRUE
)
passed <- isTRUE(comparison)
determinism_path <- file.path(
  run_dir, "serial-parallel-determinism.json"
)
dkfa_assert_input_snapshot(
  loaded_binding$snapshot, input_snapshot(),
  "Known-warp determinism verification"
)
dkfa_write_json_atomic(
  list(
    check = "known_warp_serial_parallel_determinism",
    passed = passed,
    seeds = as.list(seeds),
    tolerance = 1e-12,
    excluded_nondeterministic_fields = "runtime_seconds",
    compared_fields = as.list(compared),
    parallel_run_id = basename(run_dir),
    source_loader = source_loader,
    loaded_input_snapshot_sha256 = loaded_binding$snapshot_sha256,
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
    protocol_v9_sha256 = protocol_amendment$protocol_v9_sha256,
    court_v2_invalidation_sha256 = dkfa_hash_file(v2_invalidation_path),
    court_v2_sinkhorn_fixture_sha256 = dkfa_hash_file(sinkhorn_fixture_path),
    power_harness_bundle_sha256 = harness$bundle_sha256,
    power_output_hashes = as.list(
      power_state_validation$live_hashes$power_complete
    ),
    v2_historical_evidence_bundle_sha256 =
      v2_historical_bundle$bundle_sha256,
    detail = if (passed) "identical within tolerance" else comparison
  ),
  determinism_path
)
dkfa_update_run_state(
  run_state_path, "determinism_complete",
  outputs = list(determinism = dkfa_hash_file(determinism_path))
)
passed
})
if (!isTRUE(verification_passed)) quit(status = 1L)
cat(run_dir, "\n")
