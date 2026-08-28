#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
} else {
  normalizePath(
    "inst/validation/functional-alignment-certification/run-known-warp-power.R",
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
protocol <- file.path(dirname(script), "protocol.json")
protocol_v2 <- file.path(dirname(script), "protocol-v2.json")
v2_invalidation <- file.path(dirname(script), "court-v2-invalidation.json")
sinkhorn_fixture <- file.path(
  dirname(script), "court-v2-failing-sinkhorn-fixture.json"
)
determinism_script <- file.path(
  dirname(script), "verify-known-warp-determinism.R"
)
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
  "Known-warp loader"
)
source_loader <- loaded_binding$value
source_loader$bound_input_snapshot_sha256 <-
  loaded_binding$snapshot_sha256
source_tree_hash <- loaded_binding$snapshot$source_tree_sha256
harness <- loaded_binding$snapshot$harness
v2_historical_paths <- dkfa_v2_historical_paths(root)
v2_historical_bundle <- dkfa_review_evidence_bundle(v2_historical_paths)
protocol_amendment <- dkfa_validate_current_protocol(
  root, candidate_source_tree_sha256 = source_tree_hash
)
v2_invalidation_status <- dkfa_validate_court_v2_invalidation(
  v2_invalidation, sinkhorn_fixture, protocol_v2,
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

protocol_hash <- digest::digest(file = protocol, algo = "sha256",
                                serialize = FALSE)
run_id <- paste0(
  "known-warp-", substr(source_tree_hash, 1L, 12L), "-",
  substr(protocol_hash, 1L, 12L), "-",
  substr(harness$bundle_sha256, 1L, 12L)
)
run_dir <- file.path(
  root, "data-raw", "functional-alignment-certification", run_id
)
output_path <- file.path(run_dir, "known-warp-power-raw.csv")
summary_path <- file.path(run_dir, "known-warp-power-summary.csv")
power_manifest_path <- file.path(
  run_dir, "known-warp-power-manifest.json"
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
    power_harness_bundle_sha256 = harness$bundle_sha256,
    v2_historical_evidence_bundle_sha256 =
      v2_historical_bundle$bundle_sha256,
    loaded_input_snapshot_sha256 = loaded_binding$snapshot_sha256
  ),
  harness_filename = "known-warp-power-harness-manifest.json",
  state_filename = "known-warp-power-run-state.json"
)

seeds <- 91201:91240
workers <- if (.Platform$OS.type == "windows") 1L else 4L
started <- proc.time()[["elapsed"]]
run_seed <- function(seed) {
  tryCatch(
    list(ok = TRUE, value = dktf_power_run_one(seed)),
    error = function(e) list(ok = FALSE, seed = seed,
                             message = conditionMessage(e))
  )
}
result <- tryCatch({
  results <- if (.Platform$OS.type != "windows" && workers > 1L) {
    parallel::mclapply(
      seeds, run_seed, mc.cores = workers,
      mc.preschedule = TRUE, mc.set.seed = FALSE
    )
  } else {
    lapply(seeds, run_seed)
  }
  failed <- which(!vapply(results, `[[`, logical(1), "ok"))
  if (length(failed)) {
    details <- vapply(results[failed], function(x) {
      sprintf("seed=%d: %s", x$seed, x$message)
    }, character(1))
    stop("Known-warp task failure(s):\n", paste(details, collapse = "\n"),
         call. = FALSE)
  }
  raw <- do.call(rbind, lapply(results, `[[`, "value"))
  summary <- dktf_power_summarize(raw)
  elapsed <- proc.time()[["elapsed"]] - started
  dkfa_assert_input_snapshot(
    loaded_binding$snapshot, input_snapshot(), "Known-warp run"
  )
  utils::write.csv(raw, output_path, row.names = FALSE)
  utils::write.csv(summary, summary_path, row.names = FALSE)
  dkfa_write_json_atomic(
    list(
      protocol_id = "dkge-functional-alignment-certification-v9",
      run_id = run_id,
      source_tree_sha256 = source_tree_hash,
      protocol_sha256 = protocol_hash,
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
      court_v2_invalidation_sha256 = dkfa_hash_file(v2_invalidation),
      court_v2_sinkhorn_fixture_sha256 = dkfa_hash_file(sinkhorn_fixture),
      v2_historical_evidence_bundle_sha256 =
        v2_historical_bundle$bundle_sha256,
      v2_historical_evidence_hashes =
        as.list(v2_historical_bundle$hashes),
      harness_bundle_sha256 = harness$bundle_sha256,
      harness_files = as.list(harness$files),
      runner_sha256 = digest::digest(file = script, algo = "sha256",
                                     serialize = FALSE),
      source_loader = source_loader,
      loaded_input_snapshot_sha256 = loaded_binding$snapshot_sha256,
      seeds = list(first = min(seeds), last = max(seeds), n = length(seeds)),
      workers = workers,
      n_perm = 127L,
      heldout_signal_scale = 0.35,
      heldout_value_noise = 0.8,
      runtime_seconds = elapsed,
      role = paste(
        "Descriptive held-out efficacy only; power cannot override court",
        "or provenance gates."
      )
    ),
    power_manifest_path
  )
  power_outputs <- c(
    power_manifest = power_manifest_path,
    power_raw = output_path,
    power_summary = summary_path
  )
  dkfa_update_run_state(
    run_claim$state_path, "power_complete",
    outputs = as.list(vapply(power_outputs, dkfa_hash_file, character(1)))
  )
  raw
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
