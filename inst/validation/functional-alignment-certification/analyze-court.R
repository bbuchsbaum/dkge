#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
} else {
  normalizePath(
    "inst/validation/functional-alignment-certification/analyze-court.R",
    mustWork = TRUE
  )
}
root <- normalizePath(file.path(dirname(script), "..", "..", ".."),
                      mustWork = TRUE)
utils_script <- file.path(dirname(script), "evidence-utils.R")
court_script <- file.path(
  root, "inst", "validation", "functional-alignment-court", "court.R"
)
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
    invisible(NULL)
  },
  "Court analyzer loader"
)
candidate_source_hash <- loaded_binding$snapshot$source_tree_sha256
v2_historical_paths <- dkfa_v2_historical_paths(root)
v2_historical_bundle <- dkfa_review_evidence_bundle(v2_historical_paths)
args <- commandArgs(trailingOnly = TRUE)
evidence_root <- file.path(
  root, "data-raw", "functional-alignment-certification"
)
run_dir <- if (length(args)) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  candidates <- list.dirs(evidence_root, recursive = FALSE, full.names = TRUE)
  candidates <- candidates[grepl("^court-", basename(candidates))]
  completed <- file.exists(file.path(candidates, "formal-manifest.txt"))
  candidates <- candidates[completed]
  if (!length(candidates)) stop("No completed v9 court run found.", call. = FALSE)
  candidates[[which.max(file.info(candidates)$mtime)]]
}

runner_script <- file.path(dirname(script), "run-court.R")
protocol_path <- file.path(dirname(script), "protocol.json")
protocol_v2_path <- file.path(dirname(script), "protocol-v2.json")
v2_invalidation_path <- file.path(
  dirname(script), "court-v2-invalidation.json"
)
sinkhorn_fixture_path <- file.path(
  dirname(script), "court-v2-failing-sinkhorn-fixture.json"
)
protocol_amendment <- dkfa_validate_current_protocol(
  root, candidate_source_tree_sha256 = candidate_source_hash
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
run_state_path <- file.path(run_dir, "court-run-state.json")
run_state <- jsonlite::read_json(run_state_path, simplifyVector = TRUE)
if (!identical(run_state$status, "formal_complete")) {
  stop(
    "Court analysis requires a formally completed, non-failed run state.",
    call. = FALSE
  )
}
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
state_binding_checks <- c(
  protocol_id = identical(
    run_state$protocol_id,
    "dkge-functional-alignment-certification-v9"
  ),
  source = identical(run_state$source_tree_sha256, candidate_source_hash),
  protocol = identical(
    run_state$protocol_sha256, protocol_amendment$protocol_v9_sha256
  ) && identical(
    run_state$protocol_v9_sha256, protocol_amendment$protocol_v9_sha256
  ),
  loaded_inputs = identical(
    run_state$loaded_input_snapshot_sha256,
    loaded_binding$snapshot_sha256
  ),
  protocol_history = all(vapply(
    names(protocol_history_bindings),
    function(name) identical(
      run_state[[name]], unname(protocol_history_bindings[[name]])
    ),
    logical(1)
  )),
  v2_invalidation = identical(
    run_state$court_v2_invalidation_sha256,
    dkfa_hash_file(v2_invalidation_path)
  ) && identical(
    run_state$court_v2_sinkhorn_fixture_sha256,
    dkfa_hash_file(sinkhorn_fixture_path)
  ) && identical(
    run_state$v2_historical_evidence_bundle_sha256,
    v2_historical_bundle$bundle_sha256
  )
)
if (!all(state_binding_checks)) {
  stop(
    "Court run state is not bound to the live v9 candidate: ",
    paste(names(state_binding_checks)[!state_binding_checks], collapse = ", "),
    call. = FALSE
  )
}
analysis_passed <- dkfa_with_run_failure_receipt(
  run_state_path, "formal_complete", {
formal_outputs <- c(
  formal_budget = file.path(run_dir, "formal-budget.csv"),
  formal_grid = file.path(run_dir, "formal-grid.csv"),
  formal_raw = file.path(run_dir, "formal-raw.csv"),
  formal_summary = file.path(run_dir, "formal-summary.csv"),
  formal_inflation_flags = file.path(run_dir, "formal-inflation-flags.csv"),
  latent_span_model = file.path(run_dir, "latent-span-model.csv"),
  formal_manifest = file.path(run_dir, "formal-manifest.txt")
)
recorded_formal_hashes <- unlist(
  run_state$outputs$formal_complete, use.names = TRUE
)
current_formal_hashes <- if (all(file.exists(formal_outputs))) {
  vapply(formal_outputs, dkfa_hash_file, character(1))
} else {
  setNames(rep(NA_character_, length(formal_outputs)), names(formal_outputs))
}
if (!identical(recorded_formal_hashes, current_formal_hashes)) {
  stop("Formal court outputs do not match their pre-analysis run-state receipt.",
       call. = FALSE)
}
raw <- utils::read.csv(file.path(run_dir, "formal-raw.csv"),
                       stringsAsFactors = FALSE)
recorded_summary <- utils::read.csv(
  file.path(run_dir, "formal-summary.csv"), stringsAsFactors = FALSE
)
court_evidence <- dkfa_validate_court_evidence(
  raw, recorded_summary, protocol
)
if (!isTRUE(court_evidence$passed)) {
  stop(
    "Court raw schedule or recomputed summary is invalid: ",
    paste(names(court_evidence$checks)[!court_evidence$checks],
          collapse = ", "),
    call. = FALSE
  )
}
summary_rows <- court_evidence$recomputed_summary
manifest <- readLines(file.path(run_dir, "formal-manifest.txt"), warn = FALSE)
harness_record <- jsonlite::read_json(
  file.path(run_dir, "court-harness-manifest.json"), simplifyVector = TRUE
)
court_harness <- dkfa_verify_bundle_manifest(
  harness_record,
  dkfa_court_harness_paths(root),
  root,
  "Court"
)
manifest_value <- function(key) {
  line <- manifest[startsWith(manifest, paste0(key, "="))]
  if (!length(line)) return(NA_character_)
  sub(paste0("^", key, "="), "", line[[length(line)]])
}

focus_arms <- c(
  "geometry_only", "independent_alignment",
  "kernel_image_residual_prototype", "full_permutation_reestimation"
)
focus_summary <- summary_rows[
  summary_rows$arm %in% focus_arms, , drop = FALSE
]
utils::write.csv(
  focus_summary, file.path(run_dir, "court-focus-summary.csv"),
  row.names = FALSE
)

factor_cells <- rbind(
  data.frame(factor = "S", cell = c("small_S", "medium_S", "large_S"),
             setting = c("8", "16", "32")),
  data.frame(factor = "estimation_rank", cell = c("rank_2", "rank_4"),
             setting = c("2", "4")),
  data.frame(factor = "kernel_rank", cell = c("baseline", "kernel_rank_4"),
             setting = c("6", "4")),
  data.frame(factor = "contrast_family_dimension",
             cell = c("baseline", "family_2"), setting = c("1", "2")),
  data.frame(factor = "eigengap", cell = c("baseline", "small_gap"),
             setting = c("large", "small")),
  stringsAsFactors = FALSE
)
factor_uncertainty <- merge(
  factor_cells, focus_summary,
  by = "cell", all.x = FALSE, all.y = FALSE, sort = FALSE
)
factor_uncertainty <- factor_uncertainty[
  order(match(factor_uncertainty$factor, unique(factor_cells$factor)),
        factor_uncertainty$setting,
        match(factor_uncertainty$arm, focus_arms)), , drop = FALSE
]
utils::write.csv(
  factor_uncertainty, file.path(run_dir, "factor-uncertainty.csv"),
  row.names = FALSE
)

same_data <- raw[raw$arm %in% c(
  "l2_fold_loading", "kernel_image_residual_prototype"
), , drop = FALSE]
same_data$excess_rejection <- as.numeric(same_data$reject_fwer) - 0.05
excess_model <- stats::lm(
  excess_rejection ~ latent_span_predictor + arm,
  data = same_data
)
logistic_model <- stats::glm(
  reject_fwer ~ latent_span_predictor + arm,
  family = stats::binomial(), data = same_data
)
write_coefficients <- function(model, path) {
  coefficients <- as.data.frame(summary(model)$coefficients)
  coefficients$term <- rownames(coefficients)
  rownames(coefficients) <- NULL
  coefficients <- coefficients[, c(
    "term", setdiff(names(coefficients), "term")
  ), drop = FALSE]
  utils::write.csv(coefficients, path, row.names = FALSE)
}
write_coefficients(
  excess_model, file.path(run_dir, "latent-span-excess-model.csv")
)
write_coefficients(
  logistic_model, file.path(run_dir, "latent-span-logistic-model.csv")
)

verdict <- court_evidence$recomputed_verdict
required_numeric <- setdiff(
  names(raw)[vapply(raw, is.numeric, logical(1))],
  "exact_elapsed_seconds"
)
expected_rows <- 21L * 40L * 6L + 6L * 40L
source_at_start <- manifest_value("source_tree_hash_at_start")
source_at_end <- fa_court_source_tree_hash(root)
audit <- data.frame(
  check = c(
    "raw_rows", "expected_raw_rows", "unique_cohorts",
    "expected_unique_cohorts", "unique_seeds", "nonfinite_required_raw",
    "minimum_cell_arm_convergence", "exact_oracle_gate",
    "negative_control_promotion_gate", "source_tree_unchanged",
    "posthoc_threshold_changes"
  ),
  value = c(
    nrow(raw), expected_rows,
    length(unique(paste(raw$cell, raw$seed))), 21L * 40L,
    length(unique(raw$seed)),
    sum(!is.finite(as.matrix(raw[, required_numeric, drop = FALSE]))),
    min(summary_rows$convergence),
    as.integer(verdict$exact_oracle_gate),
    as.integer(verdict$negative_control_gate),
    as.integer(identical(source_at_start, source_at_end)),
    0L
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(audit, file.path(run_dir, "court-audit.csv"),
                 row.names = FALSE)

complete <- nrow(raw) == expected_rows &&
  length(unique(paste(raw$cell, raw$seed))) == 21L * 40L &&
  audit$value[audit$check == "nonfinite_required_raw"] == 0 &&
  min(summary_rows$convergence) >= 0.95 &&
  isTRUE(verdict$exact_oracle_gate) &&
  identical(source_at_start, source_at_end)
status <- if (!complete) {
  "invalid_or_incomplete_court"
} else if (isTRUE(verdict$negative_control_gate)) {
  "complete_controls_passed_public_modes_remain_contractually_approximate"
} else {
  "complete_approximate_only_inferential_promotion_blocked"
}

jsonlite::write_json(
  list(
    protocol_id = "dkge-functional-alignment-certification-v9",
    run_id = basename(run_dir),
    status = status,
    exact_oracle_gate = verdict$exact_oracle_gate,
    negative_control_promotion_gate = verdict$negative_control_gate,
    materially_inflated = verdict$materially_inflated,
    source_tree_sha256 = source_at_end,
    court_harness_bundle_sha256 = court_harness$bundle_sha256,
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
    v2_historical_evidence_bundle_sha256 =
      v2_historical_bundle$bundle_sha256,
    source_tree_unchanged = identical(source_at_start, source_at_end),
    thresholds_changed_after_results = FALSE,
    claim = if (complete && !verdict$negative_control_gate) {
      paste(
        "The exact comparator and evidence-integrity gates passed, but the",
        "frozen-plan promotion gate failed. Public alignment remains",
        "approximate and requires an explicit inference override."
      )
    } else if (complete) {
      paste(
        "The court controls passed for this simulated scope. Public",
        "independent and same-data modes nevertheless retain their",
        "predeclared approximate labels because the estimator contract is",
        "broader than this court."
      )
    } else {
      "The court is invalid or incomplete; no calibration claim is licensed."
    }
  ),
  file.path(run_dir, "court-verdict.json"), pretty = TRUE,
  auto_unbox = TRUE
)

format_factor_row <- function(i) {
  x <- factor_uncertainty[i, , drop = FALSE]
  sprintf(
    "| %s | %s | %s | %.3f | [%.3f, %.3f] | %.3f |",
    x$factor[[1]], x$setting[[1]], x$arm[[1]], x$fwer_rate[[1]],
    x$fwer_ci_lower[[1]], x$fwer_ci_upper[[1]], x$convergence[[1]]
  )
}
factor_lines <- vapply(
  seq_len(nrow(factor_uncertainty)), format_factor_row, character(1)
)
report <- c(
  "# Functional-alignment v9 court verdict",
  "",
  paste0("Status: `", status, "`."),
  "",
  paste0("Exact-oracle gate: **", verdict$exact_oracle_gate, "**. ",
         "Frozen-plan promotion gate: **", verdict$negative_control_gate,
         "**."),
  "",
  paste0(
    "The run contains ", nrow(raw), " arm-level results from ",
    length(unique(paste(raw$cell, raw$seed))),
    paste0(
      " fixed cohorts. No seed, threshold, epsilon, or adverse cell was ",
      "changed. V9 inherits v3's sole predeclared Sinkhorn iteration-ceiling ",
      "amendment after v2 stopped without statistical outputs; v4-v9 change ",
      "only candidate/documentation binding, package-gate and supersession ",
      "provenance, and directory labels."
    )
  ),
  "",
  "The table below reports the predeclared S, estimation-rank, kernel-rank, contrast-span, and eigengap comparisons with Wilson 95% intervals. It is an OAT/fractional grid, so differences are descriptive rather than causal factor effects.",
  "",
  "| Factor | Setting | Arm | FWER | Wilson 95% CI | Convergence |",
  "|---|---:|---|---:|---:|---:|",
  factor_lines,
  "",
  "Same-data residualized alignment remains approximate regardless of favorable cells. Independent features remove direct reuse in correspondence but do not remove the rank-truncated estimator's latent-span dependence. Geometry-only is an ancillary bound, not evidence of functional alignment.",
  "",
  "Held-out known-warp power is reported separately and cannot override this verdict."
)
writeLines(report, file.path(run_dir, "court-verdict.md"), useBytes = TRUE)
analysis_outputs <- c(
  court_audit = file.path(run_dir, "court-audit.csv"),
  court_focus_summary = file.path(run_dir, "court-focus-summary.csv"),
  factor_uncertainty = file.path(run_dir, "factor-uncertainty.csv"),
  latent_span_excess = file.path(run_dir, "latent-span-excess-model.csv"),
  latent_span_logistic = file.path(run_dir, "latent-span-logistic-model.csv"),
  court_verdict = file.path(run_dir, "court-verdict.json"),
  court_verdict_markdown = file.path(run_dir, "court-verdict.md")
)
if (!all(file.exists(analysis_outputs))) {
  stop("Court analysis did not materialize every required derivative.",
       call. = FALSE)
}
dkfa_assert_input_snapshot(
  loaded_binding$snapshot, input_snapshot(), "Court analysis"
)
dkfa_update_run_state(
  run_state_path, "analyzed_complete",
  outputs = as.list(vapply(analysis_outputs, dkfa_hash_file, character(1)))
)
TRUE
})
if (!isTRUE(analysis_passed)) quit(status = 1L)
cat(run_dir, "\n")
