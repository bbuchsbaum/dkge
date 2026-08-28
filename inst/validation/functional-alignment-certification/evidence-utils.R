#!/usr/bin/env Rscript

# Hashing, receipt, and exact-source loading helpers shared by the certification
# runners and collector. Hash/receipt functions remain side-effect free; the
# loader is explicit and fail-closed so it can never substitute an installed
# namespace for the live candidate.

dkfa_load_source_package <- function(package_root, quiet = TRUE,
                                     source_loader = "pkgload",
                                     load_all = NULL) {
  package_root <- normalizePath(package_root, mustWork = TRUE)
  if (!identical(source_loader, "pkgload")) {
    stop(
      paste0(
        "Exact-source validation requires `pkgload::load_all()`; an installed ",
        "dkge namespace is not an admissible fallback."
      ),
      call. = FALSE
    )
  }
  if (is.null(load_all)) {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop(
        paste0(
          "Exact-source validation requires `pkgload::load_all()`; an installed ",
          "dkge namespace is not an admissible fallback."
        ),
        call. = FALSE
      )
    }
    load_all <- pkgload::load_all
    loader_version <- as.character(utils::packageVersion("pkgload"))
  } else {
    if (!is.function(load_all)) {
      stop("`load_all` must be a function when supplied.", call. = FALSE)
    }
    loader_version <- "test-double"
  }
  load_all(package_root, quiet = quiet, compile = TRUE)
  list(
    method = "pkgload::load_all",
    package_root = package_root,
    loader_version = loader_version,
    compile = TRUE
  )
}

dkfa_hash_file <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  if (dir.exists(path)) stop("Expected a file, not a directory: ", path,
                             call. = FALSE)
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

dkfa_relative_path <- function(path, root) {
  path <- normalizePath(path, mustWork = TRUE)
  root <- normalizePath(root, mustWork = TRUE)
  prefix <- paste0(root, .Platform$file.sep)
  if (!startsWith(path, prefix)) {
    stop("Bundle file is outside its declared root: ", path, call. = FALSE)
  }
  substring(path, nchar(prefix) + 1L)
}

dkfa_bundle_manifest <- function(paths, root) {
  paths <- unique(vapply(paths, normalizePath, character(1), mustWork = TRUE))
  if (!length(paths) || any(file.info(paths)$isdir %in% TRUE)) {
    stop("A dependency bundle requires one or more files.", call. = FALSE)
  }
  relative <- vapply(paths, dkfa_relative_path, character(1), root = root)
  order_index <- order(relative, method = "radix")
  relative <- relative[order_index]
  paths <- paths[order_index]
  hashes <- vapply(paths, dkfa_hash_file, character(1))
  names(hashes) <- relative
  records <- paste(relative, hashes, sep = "=")
  list(
    files = hashes,
    bundle_sha256 = digest::digest(
      paste(records, collapse = "\n"), algo = "sha256", serialize = FALSE
    )
  )
}

dkfa_verify_bundle_manifest <- function(recorded, paths, root, label) {
  current <- dkfa_bundle_manifest(paths, root)
  recorded_files <- unlist(recorded$files, use.names = TRUE)
  if (!identical(as.character(recorded$bundle_sha256),
                 current$bundle_sha256) ||
      !identical(recorded_files, current$files)) {
    stop(label, " dependency bundle does not match the recorded manifest.",
         call. = FALSE)
  }
  invisible(current)
}

dkfa_review_evidence_bundle <- function(paths) {
  if (is.null(names(paths)) || any(!nzchar(names(paths))) ||
      anyDuplicated(names(paths))) {
    stop("Reviewed evidence paths require unique non-empty role names.",
         call. = FALSE)
  }
  missing <- names(paths)[!file.exists(paths)]
  if (length(missing)) {
    stop("Missing reviewed evidence: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  roles <- sort(names(paths), method = "radix")
  hashes <- vapply(paths[roles], dkfa_hash_file, character(1))
  records <- paste(names(hashes), hashes, sep = "=")
  list(
    hashes = hashes,
    bundle_sha256 = digest::digest(
      paste(records, collapse = "\n"), algo = "sha256", serialize = FALSE
    )
  )
}

dkfa_v2_historical_paths <- function(root) {
  relative <- c(
    source_test_manifest = paste0(
      "data-raw/functional-alignment-certification/",
      "package-0be19103609aa074-b143f64a3f6e/source-test-manifest.json"
    ),
    source_test_results = paste0(
      "data-raw/functional-alignment-certification/",
      "package-0be19103609aa074-b143f64a3f6e/source-test-results.csv"
    ),
    source_test_session_info = paste0(
      "data-raw/functional-alignment-certification/",
      "package-0be19103609aa074-b143f64a3f6e/session-info.txt"
    ),
    failed_court_formal_budget = paste0(
      "data-raw/functional-alignment-certification/",
      "court-0be19103609aa074-30c9732ef157/formal-budget.csv"
    )
  )
  stats::setNames(file.path(root, unname(relative)), names(relative))
}

dkfa_v3_history_paths <- function(root) {
  certification <- file.path(
    root, "inst", "validation", "functional-alignment-certification"
  )
  c(
    protocol_v3 = file.path(certification, "protocol-v3.json"),
    protocol_v3_supersession = file.path(
      certification, "protocol-v3-supersession.json"
    )
  )
}

dkfa_v4_history_paths <- function(root) {
  certification <- file.path(
    root, "inst", "validation", "functional-alignment-certification"
  )
  source_suite <- file.path(
    root, "data-raw", "functional-alignment-certification",
    "package-17e10ec404ab7758-555e000abf33"
  )
  court <- file.path(
    root, "data-raw", "functional-alignment-certification",
    "court-17e10ec404ab7758-20e29b85364e"
  )
  c(
    protocol_v4 = file.path(certification, "protocol-v4.json"),
    protocol_v4_invalidation = file.path(
      certification, "protocol-v4-invalidation.json"
    ),
    v4_source_test_manifest = file.path(
      source_suite, "source-test-manifest.json"
    ),
    v4_source_test_results = file.path(source_suite, "source-test-results.csv"),
    v4_source_test_session_info = file.path(source_suite, "session-info.txt"),
    v4_court_run_state = file.path(court, "court-run-state.json"),
    v4_court_harness_manifest = file.path(court, "court-harness-manifest.json"),
    v4_formal_budget = file.path(court, "formal-budget.csv"),
    v4_formal_grid = file.path(court, "formal-grid.csv"),
    v4_formal_inflation_flags = file.path(court, "formal-inflation-flags.csv"),
    v4_formal_manifest = file.path(court, "formal-manifest.txt"),
    v4_formal_raw = file.path(court, "formal-raw.csv"),
    v4_formal_summary = file.path(court, "formal-summary.csv"),
    v4_latent_span_model = file.path(court, "latent-span-model.csv")
  )
}

dkfa_v5_history_paths <- function(root) {
  certification <- file.path(
    root, "inst", "validation", "functional-alignment-certification"
  )
  c(
    protocol_v5 = file.path(certification, "protocol-v5.json"),
    protocol_v5_supersession = file.path(
      certification, "protocol-v5-supersession.json"
    )
  )
}

dkfa_v6_history_paths <- function(root) {
  certification <- file.path(
    root, "inst", "validation", "functional-alignment-certification"
  )
  c(
    protocol_v6 = file.path(certification, "protocol-v6.json"),
    protocol_v6_supersession = file.path(
      certification, "protocol-v6-supersession.json"
    )
  )
}

dkfa_named_directory_files <- function(root, relative_directory) {
  directory <- file.path(root, relative_directory)
  paths <- list.files(
    directory, recursive = TRUE, full.names = TRUE,
    all.files = TRUE, no.. = TRUE
  )
  paths <- paths[file.info(paths)$isdir %in% FALSE]
  if (!length(paths)) return(stats::setNames(character(), character()))
  roles <- substring(paths, nchar(directory) + 2L)
  stats::setNames(paths, roles)
}

dkfa_v7_source_test_paths <- function(root) {
  dkfa_named_directory_files(
    root,
    file.path(
      "data-raw", "functional-alignment-certification",
      "package-7bee7084bc4c66d2-5363594a0a95"
    )
  )
}

dkfa_v7_court_paths <- function(root) {
  dkfa_named_directory_files(
    root,
    file.path(
      "data-raw", "functional-alignment-certification",
      "court-7bee7084bc4c66d2-56f8cf4b5eb3"
    )
  )
}

dkfa_v7_power_paths <- function(root) {
  dkfa_named_directory_files(
    root,
    file.path(
      "data-raw", "functional-alignment-certification",
      "known-warp-7bee7084bc4c-c6d3a282b967-71c6137c8826"
    )
  )
}

dkfa_v7_history_paths <- function(root) {
  certification <- file.path(
    root, "inst", "validation", "functional-alignment-certification"
  )
  source_test <- dkfa_v7_source_test_paths(root)
  court <- dkfa_v7_court_paths(root)
  power <- dkfa_v7_power_paths(root)
  c(
    protocol_v7 = file.path(certification, "protocol-v7.json"),
    protocol_v7_supersession = file.path(
      certification, "protocol-v7-supersession.json"
    ),
    stats::setNames(unname(source_test), paste0("v7_source_test_", names(source_test))),
    stats::setNames(unname(court), paste0("v7_court_", names(court))),
    stats::setNames(unname(power), paste0("v7_power_", names(power)))
  )
}

dkfa_v8_source_test_paths <- function(root) {
  dkfa_named_directory_files(
    root,
    file.path(
      "data-raw", "functional-alignment-certification",
      "package-e38585975999d2bc-d881948e5df7"
    )
  )
}

dkfa_v8_court_paths <- function(root) {
  dkfa_named_directory_files(
    root,
    file.path(
      "data-raw", "functional-alignment-certification",
      "court-e38585975999d2bc-ce04e4dbf8cc"
    )
  )
}

dkfa_v8_power_paths <- function(root) {
  dkfa_named_directory_files(
    root,
    file.path(
      "data-raw", "functional-alignment-certification",
      "known-warp-e38585975999-0e68dd7558d3-1b550d962536"
    )
  )
}

dkfa_v8_history_paths <- function(root) {
  certification <- file.path(
    root, "inst", "validation", "functional-alignment-certification"
  )
  source_test <- dkfa_v8_source_test_paths(root)
  court <- dkfa_v8_court_paths(root)
  power <- dkfa_v8_power_paths(root)
  c(
    protocol_v8 = file.path(certification, "protocol-v8.json"),
    protocol_v8_supersession = file.path(
      certification, "protocol-v8-supersession.json"
    ),
    stats::setNames(
      unname(source_test), paste0("v8_source_test_", names(source_test))
    ),
    stats::setNames(unname(court), paste0("v8_court_", names(court))),
    stats::setNames(unname(power), paste0("v8_power_", names(power)))
  )
}

dkfa_court_harness_paths <- function(root) {
  certification <- file.path(
    root, "inst", "validation", "functional-alignment-certification"
  )
  c(
    file.path(certification, "run-court.R"),
    file.path(certification, "evidence-utils.R"),
    file.path(
      root, "inst", "validation", "functional-alignment-court", "court.R"
    ),
    file.path(certification, "protocol.json"),
    file.path(certification, "protocol-v2.json"),
    unname(dkfa_v3_history_paths(root)),
    unname(dkfa_v4_history_paths(root)),
    unname(dkfa_v5_history_paths(root)),
    unname(dkfa_v6_history_paths(root)),
    unname(dkfa_v7_history_paths(root)),
    unname(dkfa_v8_history_paths(root)),
    file.path(certification, "court-v2-invalidation.json"),
    file.path(certification, "court-v2-failing-sinkhorn-fixture.json"),
    file.path(certification, "analyze-court.R"),
    unname(dkfa_v2_historical_paths(root))
  )
}

dkfa_power_harness_paths <- function(root) {
  certification <- file.path(
    root, "inst", "validation", "functional-alignment-certification"
  )
  c(
    file.path(certification, "run-known-warp-power.R"),
    file.path(certification, "evidence-utils.R"),
    file.path(
      root, "inst", "validation", "functional-alignment-template",
      "template-benefit.R"
    ),
    file.path(certification, "known-warp-power.R"),
    file.path(
      root, "inst", "validation", "functional-alignment-court", "court.R"
    ),
    file.path(certification, "protocol.json"),
    file.path(certification, "protocol-v2.json"),
    unname(dkfa_v3_history_paths(root)),
    unname(dkfa_v4_history_paths(root)),
    unname(dkfa_v5_history_paths(root)),
    unname(dkfa_v6_history_paths(root)),
    unname(dkfa_v7_history_paths(root)),
    unname(dkfa_v8_history_paths(root)),
    file.path(certification, "court-v2-invalidation.json"),
    file.path(certification, "court-v2-failing-sinkhorn-fixture.json"),
    file.path(certification, "verify-known-warp-determinism.R"),
    unname(dkfa_v2_historical_paths(root))
  )
}

dkfa_write_json_atomic <- function(value, path) {
  parent <- dirname(path)
  if (!dir.exists(parent)) {
    stop("Atomic JSON destination directory does not exist: ", parent,
         call. = FALSE)
  }
  temporary <- tempfile(paste0(".", basename(path), "-"), tmpdir = parent)
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  jsonlite::write_json(
    value, temporary, pretty = TRUE, auto_unbox = TRUE, null = "null"
  )
  if (!file.rename(temporary, path)) {
    stop("Could not atomically publish JSON receipt: ", path, call. = FALSE)
  }
  invisible(path)
}

dkfa_claim_run_directory <- function(
    run_dir, harness, metadata = list(),
    harness_filename = "court-harness-manifest.json",
    state_filename = "court-run-state.json") {
  if (dir.exists(run_dir)) {
    stop(
      "Certification run directory is already claimed and cannot be rerun: ",
      run_dir,
      call. = FALSE
    )
  }
  if (!dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not atomically claim certification run directory: ", run_dir,
         call. = FALSE)
  }
  if (length(harness_filename) != 1L || !nzchar(harness_filename) ||
      basename(harness_filename) != harness_filename ||
      length(state_filename) != 1L || !nzchar(state_filename) ||
      basename(state_filename) != state_filename ||
      identical(harness_filename, state_filename)) {
    stop("Run receipt filenames must be distinct plain filenames.",
         call. = FALSE)
  }
  harness_path <- file.path(run_dir, harness_filename)
  state_path <- file.path(run_dir, state_filename)
  dkfa_write_json_atomic(
    list(
      schema_version = "1.0.0",
      bundle_sha256 = harness$bundle_sha256,
      files = as.list(harness$files)
    ),
    harness_path
  )
  dkfa_write_json_atomic(
    c(
      list(
        schema_version = "1.0.0",
        status = "running",
        started_at_utc = format(
          Sys.time(), tz = "UTC", usetz = TRUE
        )
      ),
      metadata
    ),
    state_path
  )
  list(harness_path = harness_path, state_path = state_path)
}

dkfa_update_run_state <- function(state_path, status,
                                  condition = NULL, outputs = NULL) {
  allowed <- c(
    "failed", "formal_complete", "analyzed_complete",
    "power_complete", "determinism_complete", "invalidated"
  )
  if (!status %in% allowed) {
    stop("Invalid certification run status: ", status, call. = FALSE)
  }
  state <- jsonlite::read_json(state_path, simplifyVector = FALSE)
  terminal <- c(
    "failed", "analyzed_complete", "determinism_complete", "invalidated"
  )
  if (state$status %in% terminal) {
    stop("A terminal certification run state cannot be changed.", call. = FALSE)
  }
  valid_transition <- (
    identical(state$status, "running") &&
      status %in% c("failed", "formal_complete", "power_complete")
  ) || (
    identical(state$status, "formal_complete") &&
      identical(status, "analyzed_complete")
  ) || (
    identical(state$status, "power_complete") &&
      status %in% c("failed", "determinism_complete")
  ) || (
    identical(status, "invalidated") &&
      state$status %in% c("running", "formal_complete", "power_complete")
  )
  if (!isTRUE(valid_transition)) {
    stop(
      "Invalid certification run-state transition: ", state$status,
      " -> ", status,
      call. = FALSE
    )
  }
  state$status <- status
  timestamp_field <- switch(
    status,
    failed = "failed_at_utc",
    formal_complete = "formal_completed_at_utc",
    analyzed_complete = "analyzed_completed_at_utc",
    power_complete = "power_completed_at_utc",
    determinism_complete = "determinism_completed_at_utc",
    invalidated = "invalidated_at_utc"
  )
  state[[timestamp_field]] <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  if (!is.null(condition)) state$condition <- condition
  if (!is.null(outputs)) {
    if (is.null(state$outputs)) state$outputs <- list()
    state$outputs[[status]] <- outputs
  }
  dkfa_write_json_atomic(state, state_path)
  invisible(state)
}

dkfa_with_run_failure_receipt <- function(
    state_path, active_status, code) {
  tryCatch(
    force(code),
    error = function(e) {
      current <- tryCatch(
        jsonlite::read_json(state_path, simplifyVector = TRUE),
        error = function(...) NULL
      )
      if (!is.null(current) && identical(current$status, active_status)) {
        try(dkfa_update_run_state(
          state_path, "failed",
          condition = list(
            classes = as.list(class(e)),
            message = conditionMessage(e)
          )
        ), silent = TRUE)
      }
      stop(e)
    }
  )
}

dkfa_validate_run_state_outputs <- function(
    state_path, expected_status, stage_paths) {
  state <- jsonlite::read_json(state_path, simplifyVector = TRUE)
  if (is.null(names(stage_paths)) || any(!nzchar(names(stage_paths))) ||
      anyDuplicated(names(stage_paths)) ||
      !all(vapply(stage_paths, function(x) {
        is.character(x) && length(x) > 0L && !is.null(names(x)) &&
          all(nzchar(names(x))) && !anyDuplicated(names(x))
      }, logical(1)))) {
    stop("Run-state validation requires named stages of named file paths.",
         call. = FALSE)
  }
  checks <- c(status = identical(state$status, expected_status))
  live_hashes <- list()
  for (stage in names(stage_paths)) {
    paths <- stage_paths[[stage]]
    present <- all(file.exists(paths) & !dir.exists(paths))
    current <- if (present) {
      vapply(paths, dkfa_hash_file, character(1))
    } else {
      stats::setNames(rep(NA_character_, length(paths)), names(paths))
    }
    recorded <- unlist(state$outputs[[stage]], use.names = TRUE)
    checks <- c(
      checks,
      stats::setNames(present, paste0(stage, "_files_present")),
      stats::setNames(
        identical(recorded, current), paste0(stage, "_hashes_match")
      )
    )
    live_hashes[[stage]] <- current
  }
  list(
    passed = all(checks), checks = checks,
    state = state, live_hashes = live_hashes
  )
}

dkfa_data_frames_equal <- function(observed, expected, keys,
                                   tolerance = 1e-12) {
  if (!identical(names(observed), names(expected)) ||
      !all(keys %in% names(observed))) {
    return(FALSE)
  }
  order_rows <- function(x) {
    arguments <- c(unname(x[keys]), list(method = "radix"))
    do.call(order, arguments)
  }
  observed <- observed[order_rows(observed), , drop = FALSE]
  expected <- expected[order_rows(expected), , drop = FALSE]
  rownames(observed) <- rownames(expected) <- NULL
  isTRUE(all.equal(
    observed, expected,
    tolerance = tolerance, check.attributes = TRUE
  ))
}

dkfa_court_expected_schedule <- function(protocol, court_functions = NULL) {
  if (is.null(court_functions)) {
    caller <- parent.frame()
    court_functions <- list(
      grid = get0("fa_court_grid", envir = caller, mode = "function",
                  inherits = TRUE),
      as_config = get0(
        "fa_court_as_config", envir = caller, mode = "function",
        inherits = TRUE
      )
    )
  }
  if (!all(vapply(court_functions[c("grid", "as_config")], is.function,
                  logical(1)))) {
    stop("Frozen court schedule functions are unavailable.", call. = FALSE)
  }
  grid <- court_functions$grid()
  base_arms <- c(
    protocol$court$additional_historical_arms,
    setdiff(protocol$court$arms_of_record, "full_permutation_reestimation")
  )
  base_arms <- unique(as.character(base_arms))
  n_sim <- as.integer(protocol$court$n_sim)
  n_perm <- as.integer(protocol$court$n_perm)
  expected_rows <- lapply(seq_len(nrow(grid)), function(cell_index) {
    cfg <- court_functions$as_config(grid[cell_index, , drop = FALSE])
    arms <- base_arms
    if (isTRUE(cfg$exact_oracle)) {
      arms <- c(arms, "full_permutation_reestimation")
    }
    do.call(rbind, lapply(seq_len(n_sim), function(sim) {
      data.frame(
        cell = rep(cfg$cell, length(arms)),
        seed = rep(7300000L + 10000L * cell_index + sim, length(arms)),
        arm = arms,
        cell_index = rep(as.integer(cell_index), length(arms)),
        S = rep(cfg$S, length(arms)),
        estimation_rank = rep(cfg$estimation_rank, length(arms)),
        kernel_rank = rep(cfg$kernel_rank, length(arms)),
        contrast_family_dimension = rep(
          cfg$contrast_family_dimension, length(arms)
        ),
        eigengap_setting = rep(cfg$eigengap, length(arms)),
        nuisance_magnitude = rep(cfg$nuisance_magnitude, length(arms)),
        spatial_correlation = rep(cfg$spatial_correlation, length(arms)),
        spatial_heteroscedasticity = rep(
          cfg$spatial_heteroscedasticity, length(arms)
        ),
        parcel_count_heterogeneity = rep(
          cfg$parcel_count_heterogeneity, length(arms)
        ),
        epsilon = rep(cfg$epsilon, length(arms)),
        functional_cost_weight = rep(
          cfg$functional_cost_weight, length(arms)
        ),
        dgp = rep(cfg$dgp, length(arms)),
        n_perm = rep(n_perm, length(arms)),
        stringsAsFactors = FALSE
      )
    }))
  })
  expected <- do.call(rbind, expected_rows)
  rownames(expected) <- NULL
  expected
}

dkfa_validate_court_evidence <- function(raw, recorded_summary, protocol) {
  caller <- parent.frame()
  court_functions <- list(
    grid = get0("fa_court_grid", envir = caller, mode = "function",
                inherits = TRUE),
    as_config = get0(
      "fa_court_as_config", envir = caller, mode = "function",
      inherits = TRUE
    ),
    summarize = get0(
      "fa_court_summarize", envir = caller, mode = "function",
      inherits = TRUE
    ),
    verdict = get0(
      "fa_court_verdict", envir = caller, mode = "function",
      inherits = TRUE
    )
  )
  available <- vapply(court_functions, is.function, logical(1))
  if (!all(available)) {
    return(list(
      passed = FALSE,
      checks = c(court_functions_available = FALSE)
    ))
  }
  expected <- dkfa_court_expected_schedule(protocol, court_functions)
  schedule_fields <- names(expected)
  keys <- c("cell", "seed", "arm")
  fields_present <- all(schedule_fields %in% names(raw))
  observed_schedule <- if (fields_present) {
    raw[, schedule_fields, drop = FALSE]
  } else {
    data.frame()
  }
  exact_schedule <- fields_present && dkfa_data_frames_equal(
    observed_schedule, expected, keys, tolerance = 0
  )
  recomputed_summary <- court_functions$summarize(raw)
  summary_fields_match <- setequal(
    names(recorded_summary), names(recomputed_summary)
  )
  summary_equal <- FALSE
  if (summary_fields_match) {
    recomputed_summary <- recomputed_summary[
      names(recorded_summary)
    ]
    summary_equal <- dkfa_data_frames_equal(
      recorded_summary, recomputed_summary,
      c("cell", "arm"), tolerance = 1e-12
    )
  }
  checks <- c(
    court_functions_available = all(available),
    exact_schedule = exact_schedule,
    no_duplicate_schedule_keys = fields_present &&
      !anyDuplicated(raw[keys]),
    summary_recomputed_from_raw = summary_equal
  )
  list(
    passed = all(checks), checks = checks,
    expected_schedule = expected,
    recomputed_summary = recomputed_summary,
    recomputed_verdict = court_functions$verdict(recomputed_summary)
  )
}

dkfa_validate_known_warp_evidence <- function(raw, recorded_summary,
                                              protocol) {
  summarize <- get0(
    "dktf_power_summarize", envir = parent.frame(), mode = "function",
    inherits = TRUE
  )
  if (!is.function(summarize)) {
    return(list(
      passed = FALSE,
      checks = c(power_summary_function_available = FALSE)
    ))
  }
  expected <- dkfa_known_warp_expected_schedule(protocol)
  keys <- c("seed", "arm")
  keys_present <- all(keys %in% names(raw))
  observed_keys <- if (keys_present) raw[keys] else data.frame()
  exact_schedule <- keys_present && dkfa_data_frames_equal(
    observed_keys, expected, keys, tolerance = 0
  )
  recomputed_summary <- summarize(raw)
  summary_fields_match <- setequal(
    names(recorded_summary), names(recomputed_summary)
  )
  summary_equal <- FALSE
  if (summary_fields_match) {
    recomputed_summary <- recomputed_summary[names(recorded_summary)]
    summary_equal <- dkfa_data_frames_equal(
      recorded_summary, recomputed_summary, "arm", tolerance = 1e-12
    )
  }
  checks <- c(
    exact_schedule = exact_schedule,
    no_duplicate_schedule_keys = keys_present && !anyDuplicated(raw[keys]),
    summary_recomputed_from_raw = summary_equal
  )
  list(
    passed = all(checks), checks = checks,
    expected_schedule = expected,
    recomputed_summary = recomputed_summary
  )
}

dkfa_known_warp_expected_schedule <- function(protocol) {
  seeds <- seq.int(
    as.integer(protocol$known_warp_power$seeds$first),
    as.integer(protocol$known_warp_power$seeds$last)
  )
  arms <- as.character(protocol$known_warp_power$arms)
  expand.grid(
    seed = seeds, arm = arms,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
}

dkfa_tree_manifest <- function(root, roots = ".", exclude = character()) {
  root <- normalizePath(root, mustWork = TRUE)
  paths <- unlist(lapply(roots, function(item) {
    absolute <- file.path(root, item)
    if (dir.exists(absolute)) {
      list.files(absolute, recursive = TRUE, full.names = TRUE,
                 all.files = TRUE, no.. = TRUE)
    } else if (file.exists(absolute)) {
      absolute
    } else {
      character()
    }
  }), use.names = FALSE)
  paths <- unique(paths[file.info(paths)$isdir %in% FALSE])
  relative <- substring(paths, nchar(root) + 2L)
  if (length(exclude) && length(relative)) {
    drop <- Reduce(`|`, lapply(exclude, grepl, x = relative,
                                perl = TRUE, ignore.case = TRUE))
    paths <- paths[!drop]
  }
  if (!length(paths)) {
    return(list(files = setNames(character(), character()),
                tree_sha256 = digest::digest(
                  "", algo = "sha256", serialize = FALSE
                )))
  }
  manifest <- dkfa_bundle_manifest(paths, root)
  list(files = manifest$files, tree_sha256 = manifest$bundle_sha256)
}

dkfa_candidate_source_manifest <- function(root) {
  dkfa_tree_manifest(
    root,
    roots = c("DESCRIPTION", "NAMESPACE", "R", "src"),
    exclude = c("\\.(o|so|dll|dylib)$")
  )
}

dkfa_candidate_source_tree_hash <- function(root) {
  dkfa_candidate_source_manifest(root)$tree_sha256
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

dkfa_input_snapshot_sha256 <- function(snapshot) {
  digest::digest(snapshot, algo = "sha256", serialize = TRUE)
}

dkfa_assert_input_snapshot <- function(expected, observed, label) {
  if (!identical(expected, observed)) {
    stop(
      label,
      " inputs changed while executable code was being loaded or run; ",
      "the evidence is inadmissible and must not be reused.",
      call. = FALSE
    )
  }
  invisible(observed)
}

dkfa_bind_loaded_inputs <- function(snapshot, reload, label) {
  if (!is.function(snapshot) || !is.function(reload)) {
    stop("`snapshot` and `reload` must be functions.", call. = FALSE)
  }
  before <- snapshot()
  value <- reload()
  after <- snapshot()
  dkfa_assert_input_snapshot(before, after, label)
  list(
    value = value,
    snapshot = after,
    snapshot_sha256 = dkfa_input_snapshot_sha256(after)
  )
}

dkfa_pkgdown_input_manifest <- function(package_root) {
  roots <- c(
    "DESCRIPTION", "NAMESPACE", "R", "man", "vignettes", "README.md",
    "_pkgdown.yml", "pkgdown"
  )
  dkfa_tree_manifest(package_root, roots = roots,
                     exclude = c("(^|/)\\.DS_Store$"))
}

dkfa_site_manifest <- function(site_dir,
                               receipt = ".dkge-pkgdown-receipt.json") {
  dkfa_tree_manifest(
    site_dir, roots = ".",
    exclude = c(
      paste0("(^|/)\\Q", receipt, "\\E$"),
      "(^|/)\\.DS_Store$"
    )
  )
}

dkfa_validate_pkgdown_receipt <- function(
    site_dir, tarball_sha256, source_tree_sha256, docs_input_sha256,
    required_files) {
  receipt_path <- file.path(site_dir, ".dkge-pkgdown-receipt.json")
  if (!file.exists(receipt_path)) {
    return(list(passed = FALSE, reason = "missing pkgdown receipt"))
  }
  receipt <- jsonlite::read_json(receipt_path, simplifyVector = TRUE)
  missing <- required_files[!file.exists(file.path(site_dir, required_files))]
  current_site <- dkfa_site_manifest(site_dir)
  checks <- c(
    schema = identical(receipt$schema_version, "1.0.0"),
    tarball = identical(receipt$tarball_sha256, tarball_sha256),
    source = identical(receipt$tarball_source_tree_sha256,
                       source_tree_sha256),
    docs = identical(receipt$docs_input_sha256, docs_input_sha256),
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

dkfa_canonicalize_json <- function(x) {
  if (!is.list(x)) return(x)
  if (!is.null(names(x))) x <- x[order(names(x), method = "radix")]
  lapply(x, dkfa_canonicalize_json)
}

dkfa_validate_v3_protocol_amendment <- function(
    protocol_v2_path, protocol_v3_path) {
  v2 <- jsonlite::read_json(protocol_v2_path, simplifyVector = FALSE)
  v3 <- jsonlite::read_json(protocol_v3_path, simplifyVector = FALSE)
  expected <- v2
  expected$protocol_id <- "dkge-functional-alignment-certification-v3"
  expected$status <- "frozen_before_v3_results"
  expected$amendment <- v3$amendment
  expected$court$sinkhorn_max_iter <- 50000L
  expected$court$implementation <- v3$court$implementation
  checks <- c(
    v2_protocol_hash = identical(
      dkfa_hash_file(protocol_v2_path),
      "9ce7d46a4e2efd98f9319f157dd477a78125972e7a4e961aa772a8377541b25e"
    ),
    v3_protocol_hash = identical(
      dkfa_hash_file(protocol_v3_path),
      "f2b843c630ec2166a152b306938b105e42a9e06203fc0cd50e3855d343d78dc7"
    ),
    v2_id = identical(
      v2$protocol_id, "dkge-functional-alignment-certification-v2"
    ),
    v3_id = identical(
      v3$protocol_id, "dkge-functional-alignment-certification-v3"
    ),
    v3_status = identical(v3$status, "frozen_before_v3_results"),
    no_candidate_source_binding = is.null(
      v3$candidate_source_tree_sha256
    ),
    v3_iteration_cap = identical(v3$court$sinkhorn_max_iter, 50000L),
    supersedes_v2 = identical(
      v3$amendment$supersedes,
      "dkge-functional-alignment-certification-v2"
    ),
    one_numerical_change = grepl(
      "5000 to 50000", v3$amendment$only_numerical_change,
      fixed = TRUE
    ),
    metadata_only_change = identical(
      v3$amendment$metadata_only_change,
      paste(
        "Correct the court output-directory version label from v2 to v3.",
        "This is a provenance-label correction and changes no data generation,",
        "estimator, solver, statistic, gate, or stopping rule."
      )
    ) && identical(
      v3$court$implementation,
      paste(
        "The v1 court implementation and OAT grid are reused without changing",
        "estimands, signs, seeds, arms, or gates; outputs are written to a",
        "source-tree-addressed v3 directory."
      )
    ),
    no_second_amendment = grepl(
      "will not be raised again", v3$amendment$no_second_amendment,
      fixed = TRUE
    ),
    exact_v2_to_v3_delta = identical(
      dkfa_canonicalize_json(v3), dkfa_canonicalize_json(expected)
    )
  )
  list(
    passed = all(checks),
    checks = checks,
    protocol_v2_sha256 = dkfa_hash_file(protocol_v2_path),
    protocol_v3_sha256 = dkfa_hash_file(protocol_v3_path)
  )
}

dkfa_validate_v3_supersession <- function(
    supersession_path, protocol_v3_path,
    candidate_source_tree_sha256 = NULL, spatial_path = NULL) {
  receipt <- jsonlite::read_json(supersession_path, simplifyVector = FALSE)
  expected_source <-
    "17e10ec404ab7758ada4591522ec415a98cc6c93e7a203d0ba6b42ac9aee066f"
  expected_spatial <-
    "f3e7c84558c34ebe2b88d8427faf2deb1316b9071c214bf419435d2765d8dca2"
  checks <- c(
    receipt_hash = identical(
      dkfa_hash_file(supersession_path),
      "9e29eeb8d58da95e30eafe4818361bad7b1e44b51dbac5ee074b46dce4be0f16"
    ),
    status = identical(
      receipt$status,
      "superseded_before_execution_due_candidate_source_drift"
    ),
    v3_protocol_id = identical(
      receipt$superseded_protocol$protocol_id,
      "dkge-functional-alignment-certification-v3"
    ),
    v3_protocol_path = identical(
      receipt$superseded_protocol$path, "protocol-v3.json"
    ),
    v3_protocol_hash = identical(
      receipt$superseded_protocol$sha256,
      "f2b843c630ec2166a152b306938b105e42a9e06203fc0cd50e3855d343d78dc7"
    ) && identical(
      dkfa_hash_file(protocol_v3_path),
      receipt$superseded_protocol$sha256
    ),
    pre_drift_source = identical(
      receipt$superseded_protocol$candidate_source_tree_sha256_before_drift,
      "0be19103609aa074d2ff86c15c27f6f701442623800a10cbe599b01c141ae87f"
    ),
    no_v3_execution = !isTRUE(receipt$materialization$court_started) &&
      !isTRUE(receipt$materialization$known_warp_power_started) &&
      !isTRUE(receipt$materialization$statistical_outputs_materialized) &&
      !isTRUE(
        receipt$materialization$statistical_results_used_for_supersession
      ),
    changed_path = identical(receipt$source_change$path, "R/dkge-spatial.R"),
    changed_file_hash = identical(
      receipt$source_change$sha256_after_drift, expected_spatial
    ) && (
      is.null(spatial_path) || identical(
        dkfa_hash_file(spatial_path), expected_spatial
      )
    ),
    successor_id = identical(
      receipt$successor$protocol_id,
      "dkge-functional-alignment-certification-v4"
    ),
    successor_path = identical(receipt$successor$path, "protocol.json"),
    candidate_source_binding = identical(
      receipt$source_change$candidate_source_tree_sha256_after_drift,
      expected_source
    ) && identical(
      receipt$successor$candidate_source_tree_sha256, expected_source
    ) && (
      is.null(candidate_source_tree_sha256) || identical(
        candidate_source_tree_sha256, expected_source
      )
    ),
    no_numerical_change = identical(
      receipt$successor$numerical_protocol_change, FALSE
    )
  )
  list(
    passed = all(checks), checks = checks,
    supersession_sha256 = dkfa_hash_file(supersession_path),
    receipt = receipt
  )
}

dkfa_validate_v4_protocol_amendment <- function(
    protocol_v2_path, protocol_v3_path, supersession_path,
    protocol_v4_path, candidate_source_tree_sha256 = NULL,
    spatial_path = NULL) {
  v3_validation <- dkfa_validate_v3_protocol_amendment(
    protocol_v2_path, protocol_v3_path
  )
  supersession <- dkfa_validate_v3_supersession(
    supersession_path, protocol_v3_path,
    candidate_source_tree_sha256 = candidate_source_tree_sha256,
    spatial_path = spatial_path
  )
  v3 <- jsonlite::read_json(protocol_v3_path, simplifyVector = FALSE)
  v4 <- jsonlite::read_json(protocol_v4_path, simplifyVector = FALSE)
  expected_source <-
    "17e10ec404ab7758ada4591522ec415a98cc6c93e7a203d0ba6b42ac9aee066f"
  expected <- v3
  expected$protocol_id <- "dkge-functional-alignment-certification-v4"
  expected$status <- "frozen_before_v4_results"
  expected$candidate_source_tree_sha256 <- expected_source
  expected$amendment <- v4$amendment
  expected$court$implementation <- v4$court$implementation
  checks <- c(
    v3_chain = isTRUE(v3_validation$passed),
    v3_supersession = isTRUE(supersession$passed),
    v4_protocol_hash = identical(
      dkfa_hash_file(protocol_v4_path),
      "728a40762eeeca517760452d854cf2a870ba299cad05d77366f7afb43acd4cba"
    ),
    v4_id = identical(
      v4$protocol_id, "dkge-functional-alignment-certification-v4"
    ),
    v4_status = identical(v4$status, "frozen_before_v4_results"),
    candidate_source_binding = identical(
      v4$candidate_source_tree_sha256, expected_source
    ) && (
      is.null(candidate_source_tree_sha256) || identical(
        candidate_source_tree_sha256, expected_source
      )
    ),
    supersedes_v3 = identical(
      v4$amendment$supersedes,
      "dkge-functional-alignment-certification-v3"
    ),
    v3_receipts = identical(
      v4$amendment$v3_protocol, "protocol-v3.json"
    ) && identical(
      v4$amendment$v3_protocol_sha256,
      v3_validation$protocol_v3_sha256
    ) && identical(
      v4$amendment$v3_supersession,
      "protocol-v3-supersession.json"
    ) && identical(
      v4$amendment$v3_supersession_sha256,
      supersession$supersession_sha256
    ),
    inherited_iteration_cap = identical(
      v4$court$sinkhorn_max_iter, 50000L
    ) && identical(
      v4$court$sinkhorn_max_iter, v3$court$sinkhorn_max_iter
    ),
    no_second_numerical_amendment = grepl(
      "V4 makes no numerical protocol amendment",
      v4$amendment$no_second_numerical_amendment,
      fixed = TRUE
    ) && grepl(
      "will not be raised again",
      v4$amendment$no_second_numerical_amendment,
      fixed = TRUE
    ),
    v4_directory_label = identical(
      v4$court$implementation,
      paste(
        "The v1 court implementation and OAT grid are reused without changing",
        "estimands, signs, seeds, arms, or gates; outputs are written to a",
        "source-tree-addressed v4 directory."
      )
    ),
    exact_v3_to_v4_delta = identical(
      dkfa_canonicalize_json(v4), dkfa_canonicalize_json(expected)
    )
  )
  list(
    passed = all(checks), checks = checks,
    v3_checks = v3_validation$checks,
    supersession_checks = supersession$checks,
    protocol_v2_sha256 = v3_validation$protocol_v2_sha256,
    protocol_v3_sha256 = v3_validation$protocol_v3_sha256,
    protocol_v3_supersession_sha256 = supersession$supersession_sha256,
    protocol_v4_sha256 = dkfa_hash_file(protocol_v4_path)
  )
}

dkfa_validate_v4_invalidation <- function(
    invalidation_path, protocol_v4_path, root) {
  receipt <- jsonlite::read_json(invalidation_path, simplifyVector = FALSE)
  state_path <- file.path(
    root, receipt$court_run$files$run_state$path
  )
  state <- jsonlite::read_json(state_path, simplifyVector = FALSE)
  source_manifest_path <- file.path(
    root,
    receipt$validated_source_suite$files$source_test_manifest$path
  )
  source_manifest <- jsonlite::read_json(
    source_manifest_path, simplifyVector = FALSE
  )
  file_entries <- c(
    receipt$validated_source_suite$files,
    receipt$court_run$files
  )
  file_checks <- vapply(file_entries, function(entry) {
    path <- file.path(root, entry$path)
    file.exists(path) &&
      identical(dkfa_hash_file(path), entry$sha256) &&
      identical(as.numeric(file.info(path)$size), as.numeric(entry$bytes))
  }, logical(1))
  names(file_checks) <- paste0("retained_", names(file_entries))

  formal_state_hashes <- state$outputs$formal_complete
  formal_receipt_names <- c(
    formal_budget = "formal_budget",
    formal_grid = "formal_grid",
    formal_raw = "formal_raw",
    formal_summary = "formal_summary",
    formal_inflation_flags = "formal_inflation_flags",
    latent_span_model = "latent_span_model",
    formal_manifest = "formal_manifest"
  )
  state_output_checks <- vapply(names(formal_receipt_names), function(name) {
    receipt_name <- formal_receipt_names[[name]]
    identical(
      formal_state_hashes[[name]],
      receipt$court_run$files[[receipt_name]]$sha256
    )
  }, logical(1))
  names(state_output_checks) <- paste0("state_output_", names(state_output_checks))

  checks <- c(
    receipt_hash = identical(
      dkfa_hash_file(invalidation_path),
      "470ff4158791bb5bc1e1c598f676201913cc0534532154bbd5a52f1de321a369"
    ),
    status = identical(
      receipt$status,
      "invalidated_before_analysis_due_candidate_source_drift"
    ),
    protocol_path = identical(
      receipt$protocol$path,
      paste0(
        "inst/validation/functional-alignment-certification/",
        "protocol-v4.json"
      )
    ),
    protocol_hash = identical(
      receipt$protocol$sha256,
      "728a40762eeeca517760452d854cf2a870ba299cad05d77366f7afb43acd4cba"
    ) && identical(
      dkfa_hash_file(protocol_v4_path), receipt$protocol$sha256
    ),
    v4_source_binding = identical(
      receipt$protocol$candidate_source_tree_sha256,
      "17e10ec404ab7758ada4591522ec415a98cc6c93e7a203d0ba6b42ac9aee066f"
    ),
    v4_test_binding = identical(
      receipt$protocol$test_validation_tree_sha256,
      "555e000abf338344763bfdb4cc718d530c6b86c4a5fd1f403a245f7ca1d0a7d7"
    ),
    v4_harness_binding = identical(
      receipt$protocol$court_harness_bundle_sha256,
      "20e29b85364ec73df7a0ba218f984908d333b07e9f7905399ef269b6c3c18e2b"
    ),
    source_suite_passed = isTRUE(receipt$validated_source_suite$passed) &&
      isFALSE(receipt$validated_source_suite$statistical_court_output) &&
      isTRUE(source_manifest$passed) &&
      identical(
        source_manifest$source_tree_sha256,
        receipt$protocol$candidate_source_tree_sha256
      ) && identical(
        source_manifest$test_validation_tree_sha256,
        receipt$protocol$test_validation_tree_sha256
      ),
    state_terminal = identical(state$status, "invalidated") &&
      identical(
        state$condition$classes,
        list("dkfa_candidate_source_drift")
      ),
    state_protocol = identical(
      state$protocol_v4_sha256, receipt$protocol$sha256
    ) && identical(
      state$source_tree_sha256,
      receipt$protocol$candidate_source_tree_sha256
    ),
    drift_binding = identical(
      receipt$source_drift$candidate_source_tree_sha256_at_start,
      receipt$protocol$candidate_source_tree_sha256
    ) && identical(
      receipt$source_drift$candidate_source_tree_sha256_at_detection,
      "1cff0f30bf445c21078cbf64adc316ed50b7aa995059c043f2f10f9fdb017136"
    ),
    no_statistical_inspection = isTRUE(
      receipt$materialization$formal_outputs_materialized
    ) && isTRUE(receipt$materialization$analyzer_started) &&
      isFALSE(receipt$materialization$analyzer_read_statistical_outputs) &&
      isFALSE(receipt$materialization$analysis_outputs_materialized) &&
      isFALSE(receipt$materialization$statistical_results_inspected) &&
      isFALSE(receipt$materialization$known_warp_power_started) &&
      isFALSE(receipt$materialization$court_rerun_under_v4) &&
      isFALSE(state$condition$analysis_started) &&
      isFALSE(state$condition$statistical_outputs_inspected),
    successor = identical(
      receipt$successor$protocol_id,
      "dkge-functional-alignment-certification-v5"
    ) && identical(
      receipt$successor$path,
      paste0(
        "inst/validation/functional-alignment-certification/",
        "protocol.json"
      )
    ) && isFALSE(receipt$successor$numerical_protocol_change) &&
      isFALSE(receipt$successor$statistical_protocol_change) &&
      isTRUE(receipt$successor$fresh_source_suite_required) &&
      isTRUE(receipt$successor$fresh_court_required) &&
      isTRUE(receipt$successor$fresh_known_warp_power_required),
    file_checks,
    state_output_checks
  )
  list(
    passed = all(checks), checks = checks,
    protocol_v4_invalidation_sha256 = dkfa_hash_file(invalidation_path),
    v4_history_bundle = dkfa_review_evidence_bundle(
      dkfa_v4_history_paths(root)
    ),
    receipt = receipt
  )
}

dkfa_validate_v5_protocol_amendment <- function(
    protocol_v2_path, protocol_v3_path, supersession_path,
    protocol_v4_path, v4_invalidation_path, protocol_v5_path, root,
    candidate_source_tree_sha256 = NULL) {
  v4_validation <- dkfa_validate_v4_protocol_amendment(
    protocol_v2_path, protocol_v3_path, supersession_path,
    protocol_v4_path,
    candidate_source_tree_sha256 = NULL,
    spatial_path = NULL
  )
  invalidation <- dkfa_validate_v4_invalidation(
    v4_invalidation_path, protocol_v4_path, root
  )
  v4 <- jsonlite::read_json(protocol_v4_path, simplifyVector = FALSE)
  v5 <- jsonlite::read_json(protocol_v5_path, simplifyVector = FALSE)
  expected_source <-
    "1904edc4832aa028636fba9a57ef74497db62b6c9cee5cf38c917513499832e4"
  expected <- v4
  expected$protocol_id <- "dkge-functional-alignment-certification-v5"
  expected$status <- "frozen_before_v5_results"
  expected$candidate_source_tree_sha256 <- expected_source
  expected$amendment <- v5$amendment
  expected$court$implementation <- v5$court$implementation
  checks <- c(
    v4_chain = isTRUE(v4_validation$passed),
    v4_invalidation = isTRUE(invalidation$passed),
    v5_protocol_hash = identical(
      dkfa_hash_file(protocol_v5_path),
      "18b19d49ee853b62a0ae1a59e95c9d82f1b56207a5eb1c3e9cf097a81f576612"
    ),
    v5_id = identical(
      v5$protocol_id, "dkge-functional-alignment-certification-v5"
    ),
    v5_status = identical(v5$status, "frozen_before_v5_results"),
    candidate_source_binding = identical(
      v5$candidate_source_tree_sha256, expected_source
    ) && (
      is.null(candidate_source_tree_sha256) || identical(
        candidate_source_tree_sha256, expected_source
      )
    ),
    supersedes_v4 = identical(
      v5$amendment$supersedes,
      "dkge-functional-alignment-certification-v4"
    ),
    v4_receipts = identical(
      v5$amendment$v4_protocol, "protocol-v4.json"
    ) && identical(
      v5$amendment$v4_protocol_sha256,
      v4_validation$protocol_v4_sha256
    ) && identical(
      v5$amendment$v4_invalidation, "protocol-v4-invalidation.json"
    ) && identical(
      v5$amendment$v4_invalidation_sha256,
      invalidation$protocol_v4_invalidation_sha256
    ),
    inherited_iteration_cap = identical(
      v5$court$sinkhorn_max_iter, 50000L
    ) && identical(
      v5$court$sinkhorn_max_iter, v4$court$sinkhorn_max_iter
    ),
    no_new_protocol_amendment = grepl(
      "V5 makes no numerical or statistical protocol amendment",
      v5$amendment$no_second_numerical_amendment,
      fixed = TRUE
    ) && grepl(
      "will not be raised again",
      v5$amendment$no_second_numerical_amendment,
      fixed = TRUE
    ),
    v5_directory_label = identical(
      v5$court$implementation,
      paste(
        "The v1 court implementation and OAT grid are reused without changing",
        "estimands, signs, seeds, arms, or gates; outputs are written to a",
        "source-tree-addressed v5 directory."
      )
    ),
    exact_v4_to_v5_delta = identical(
      dkfa_canonicalize_json(v5), dkfa_canonicalize_json(expected)
    )
  )
  list(
    passed = all(checks), checks = checks,
    v4_checks = v4_validation$checks,
    v4_invalidation_checks = invalidation$checks,
    protocol_v2_sha256 = v4_validation$protocol_v2_sha256,
    protocol_v3_sha256 = v4_validation$protocol_v3_sha256,
    protocol_v3_supersession_sha256 =
      v4_validation$protocol_v3_supersession_sha256,
    protocol_v4_sha256 = v4_validation$protocol_v4_sha256,
    protocol_v4_invalidation_sha256 =
      invalidation$protocol_v4_invalidation_sha256,
    v4_history_bundle = invalidation$v4_history_bundle,
    protocol_v5_sha256 = dkfa_hash_file(protocol_v5_path)
  )
}

dkfa_validate_v5_supersession <- function(
    supersession_path, protocol_v5_path,
    candidate_source_tree_sha256 = NULL) {
  receipt <- jsonlite::read_json(supersession_path, simplifyVector = FALSE)
  expected_source <-
    "9054ca2e51ea11ea88cca447c15536194090d21732dd7c50f95462d7b98f3887"
  checks <- c(
    receipt_hash = identical(
      dkfa_hash_file(supersession_path),
      "cf3d487d81dfb163589d8f65ae597344c08301b046d5862f89c8623ef6665168"
    ),
    status = identical(
      receipt$status,
      "superseded_before_execution_due_candidate_source_drift"
    ),
    v5_protocol_id = identical(
      receipt$superseded_protocol$protocol_id,
      "dkge-functional-alignment-certification-v5"
    ),
    v5_protocol_path = identical(
      receipt$superseded_protocol$path, "protocol-v5.json"
    ),
    v5_protocol_hash = identical(
      receipt$superseded_protocol$sha256,
      "18b19d49ee853b62a0ae1a59e95c9d82f1b56207a5eb1c3e9cf097a81f576612"
    ) && identical(
      dkfa_hash_file(protocol_v5_path),
      receipt$superseded_protocol$sha256
    ),
    pre_drift_source = identical(
      receipt$superseded_protocol$candidate_source_tree_sha256_before_drift,
      "1904edc4832aa028636fba9a57ef74497db62b6c9cee5cf38c917513499832e4"
    ),
    no_v5_execution = isFALSE(
      receipt$materialization$independent_freeze_approved
    ) && isFALSE(receipt$materialization$source_test_suite_started) &&
      isFALSE(receipt$materialization$court_started) &&
      isFALSE(receipt$materialization$known_warp_power_started) &&
      isFALSE(receipt$materialization$statistical_outputs_materialized) &&
      isFALSE(
        receipt$materialization$statistical_results_used_for_supersession
      ),
    reviewed_source = identical(
      receipt$source_change$candidate_source_tree_sha256_after_drift,
      expected_source
    ) && identical(
      receipt$source_change$review_verdict, "FREEZE"
    ) && identical(
      receipt$source_change$reviewed_source_tree_sha256, expected_source
    ) && identical(
      receipt$source_change$reviewed_test_validation_tree_sha256_before_protocol_tests,
      "dc5373ed9038bc89308124033e0c556abb9e0b1404f822530e8ae3308121ece8"
    ) && (
      is.null(candidate_source_tree_sha256) || identical(
        candidate_source_tree_sha256, expected_source
      )
    ),
    successor = identical(
      receipt$successor$protocol_id,
      "dkge-functional-alignment-certification-v6"
    ) && identical(receipt$successor$path, "protocol.json") &&
      identical(
        receipt$successor$candidate_source_tree_sha256, expected_source
      ) && isFALSE(receipt$successor$numerical_protocol_change) &&
      isFALSE(receipt$successor$statistical_protocol_change)
  )
  list(
    passed = all(checks), checks = checks,
    protocol_v5_supersession_sha256 = dkfa_hash_file(supersession_path),
    receipt = receipt
  )
}

dkfa_validate_v6_protocol_amendment <- function(
    protocol_v2_path, protocol_v3_path, v3_supersession_path,
    protocol_v4_path, v4_invalidation_path, protocol_v5_path,
    v5_supersession_path, protocol_v6_path, root,
    candidate_source_tree_sha256 = NULL) {
  v5_validation <- dkfa_validate_v5_protocol_amendment(
    protocol_v2_path, protocol_v3_path, v3_supersession_path,
    protocol_v4_path, v4_invalidation_path, protocol_v5_path, root,
    candidate_source_tree_sha256 = NULL
  )
  supersession <- dkfa_validate_v5_supersession(
    v5_supersession_path, protocol_v5_path,
    candidate_source_tree_sha256 = candidate_source_tree_sha256
  )
  v5 <- jsonlite::read_json(protocol_v5_path, simplifyVector = FALSE)
  v6 <- jsonlite::read_json(protocol_v6_path, simplifyVector = FALSE)
  expected_source <-
    "9054ca2e51ea11ea88cca447c15536194090d21732dd7c50f95462d7b98f3887"
  expected <- v5
  expected$protocol_id <- "dkge-functional-alignment-certification-v6"
  expected$status <- "frozen_before_v6_results"
  expected$candidate_source_tree_sha256 <- expected_source
  expected$amendment <- v6$amendment
  expected$court$implementation <- v6$court$implementation
  checks <- c(
    v5_chain = isTRUE(v5_validation$passed),
    v5_supersession = isTRUE(supersession$passed),
    v6_protocol_hash = identical(
      dkfa_hash_file(protocol_v6_path),
      "6c40c99387fd157e875871a242256d1ab6bea9d2fcc8f2ba21309f26b150f5ef"
    ),
    v6_id = identical(
      v6$protocol_id, "dkge-functional-alignment-certification-v6"
    ),
    v6_status = identical(v6$status, "frozen_before_v6_results"),
    candidate_source_binding = identical(
      v6$candidate_source_tree_sha256, expected_source
    ) && (
      is.null(candidate_source_tree_sha256) || identical(
        candidate_source_tree_sha256, expected_source
      )
    ),
    supersedes_v5 = identical(
      v6$amendment$supersedes,
      "dkge-functional-alignment-certification-v5"
    ),
    v5_receipts = identical(
      v6$amendment$v5_protocol, "protocol-v5.json"
    ) && identical(
      v6$amendment$v5_protocol_sha256,
      v5_validation$protocol_v5_sha256
    ) && identical(
      v6$amendment$v5_supersession,
      "protocol-v5-supersession.json"
    ) && identical(
      v6$amendment$v5_supersession_sha256,
      supersession$protocol_v5_supersession_sha256
    ),
    inherited_iteration_cap = identical(
      v6$court$sinkhorn_max_iter, 50000L
    ) && identical(
      v6$court$sinkhorn_max_iter, v5$court$sinkhorn_max_iter
    ),
    no_new_protocol_amendment = grepl(
      "V6 makes no numerical or statistical protocol amendment",
      v6$amendment$no_second_numerical_amendment,
      fixed = TRUE
    ) && grepl(
      "will not be raised again",
      v6$amendment$no_second_numerical_amendment,
      fixed = TRUE
    ),
    v6_directory_label = identical(
      v6$court$implementation,
      paste(
        "The v1 court implementation and OAT grid are reused without changing",
        "estimands, signs, seeds, arms, or gates; outputs are written to a",
        "source-tree-addressed v6 directory."
      )
    ),
    exact_v5_to_v6_delta = identical(
      dkfa_canonicalize_json(v6), dkfa_canonicalize_json(expected)
    )
  )
  list(
    passed = all(checks), checks = checks,
    v5_checks = v5_validation$checks,
    v5_supersession_checks = supersession$checks,
    protocol_v2_sha256 = v5_validation$protocol_v2_sha256,
    protocol_v3_sha256 = v5_validation$protocol_v3_sha256,
    protocol_v3_supersession_sha256 =
      v5_validation$protocol_v3_supersession_sha256,
    protocol_v4_sha256 = v5_validation$protocol_v4_sha256,
    protocol_v4_invalidation_sha256 =
      v5_validation$protocol_v4_invalidation_sha256,
    v4_history_bundle = v5_validation$v4_history_bundle,
    protocol_v5_sha256 = v5_validation$protocol_v5_sha256,
    protocol_v5_supersession_sha256 =
      supersession$protocol_v5_supersession_sha256,
    v5_history_bundle = dkfa_review_evidence_bundle(
      dkfa_v5_history_paths(root)
    ),
    protocol_v6_sha256 = dkfa_hash_file(protocol_v6_path)
  )
}

dkfa_validate_v6_supersession <- function(
    supersession_path, protocol_v6_path,
    candidate_source_tree_sha256 = NULL) {
  receipt <- jsonlite::read_json(supersession_path, simplifyVector = FALSE)
  expected_source <-
    "7bee7084bc4c66d29880efd3006c16d4118e38b08c3242127ab9d1b4e9f4a2b6"
  checks <- c(
    receipt_hash = identical(
      dkfa_hash_file(supersession_path),
      "1a6858b01b8b7a62d7155bd3b3cf9ef331717c9fba63f156e7a0817ea5db4463"
    ),
    status = identical(
      receipt$status,
      "superseded_before_execution_due_candidate_source_drift"
    ),
    v6_protocol_id = identical(
      receipt$superseded_protocol$protocol_id,
      "dkge-functional-alignment-certification-v6"
    ),
    v6_protocol_path = identical(
      receipt$superseded_protocol$path, "protocol-v6.json"
    ),
    v6_protocol_hash = identical(
      receipt$superseded_protocol$sha256,
      "6c40c99387fd157e875871a242256d1ab6bea9d2fcc8f2ba21309f26b150f5ef"
    ) && identical(
      dkfa_hash_file(protocol_v6_path),
      receipt$superseded_protocol$sha256
    ),
    pre_drift_source = identical(
      receipt$superseded_protocol$candidate_source_tree_sha256_before_drift,
      "9054ca2e51ea11ea88cca447c15536194090d21732dd7c50f95462d7b98f3887"
    ),
    freeze_revoked_before_execution = isTRUE(
      receipt$materialization$independent_freeze_issued
    ) && isTRUE(
      receipt$materialization$independent_freeze_subsequently_revoked
    ) && isFALSE(receipt$materialization$source_test_suite_started) &&
      isFALSE(receipt$materialization$court_started) &&
      isFALSE(receipt$materialization$known_warp_power_started) &&
      isFALSE(receipt$materialization$statistical_outputs_materialized) &&
      isFALSE(
        receipt$materialization$statistical_results_used_for_supersession
      ),
    reviewed_source = identical(
      receipt$source_change$candidate_source_tree_sha256_after_drift,
      expected_source
    ) && identical(
      receipt$source_change$review_verdict, "FREEZE"
    ) && identical(
      receipt$source_change$reviewed_source_tree_sha256, expected_source
    ) && identical(
      receipt$source_change$reviewed_test_validation_tree_sha256_before_protocol_tests,
      "cc001f434e191516a673bcda063d79f2adb5923fa7107c2b38734b28c6a24860"
    ) && (
      is.null(candidate_source_tree_sha256) || identical(
        candidate_source_tree_sha256, expected_source
      )
    ),
    successor = identical(
      receipt$successor$protocol_id,
      "dkge-functional-alignment-certification-v7"
    ) && identical(receipt$successor$path, "protocol.json") &&
      identical(
        receipt$successor$candidate_source_tree_sha256, expected_source
      ) && isFALSE(receipt$successor$numerical_protocol_change) &&
      isFALSE(receipt$successor$statistical_protocol_change)
  )
  list(
    passed = all(checks), checks = checks,
    protocol_v6_supersession_sha256 = dkfa_hash_file(supersession_path),
    receipt = receipt
  )
}

dkfa_validate_v7_protocol_amendment <- function(
    protocol_v2_path, protocol_v3_path, v3_supersession_path,
    protocol_v4_path, v4_invalidation_path, protocol_v5_path,
    v5_supersession_path, protocol_v6_path, v6_supersession_path,
    protocol_v7_path, root, candidate_source_tree_sha256 = NULL) {
  v6_validation <- dkfa_validate_v6_protocol_amendment(
    protocol_v2_path, protocol_v3_path, v3_supersession_path,
    protocol_v4_path, v4_invalidation_path, protocol_v5_path,
    v5_supersession_path, protocol_v6_path, root,
    candidate_source_tree_sha256 = NULL
  )
  supersession <- dkfa_validate_v6_supersession(
    v6_supersession_path, protocol_v6_path,
    candidate_source_tree_sha256 = candidate_source_tree_sha256
  )
  v6 <- jsonlite::read_json(protocol_v6_path, simplifyVector = FALSE)
  v7 <- jsonlite::read_json(protocol_v7_path, simplifyVector = FALSE)
  expected_source <-
    "7bee7084bc4c66d29880efd3006c16d4118e38b08c3242127ab9d1b4e9f4a2b6"
  expected <- v6
  expected$protocol_id <- "dkge-functional-alignment-certification-v7"
  expected$status <- "frozen_before_v7_results"
  expected$candidate_source_tree_sha256 <- expected_source
  expected$amendment <- v7$amendment
  expected$court$implementation <- v7$court$implementation
  checks <- c(
    v6_chain = isTRUE(v6_validation$passed),
    v6_supersession = isTRUE(supersession$passed),
    v7_protocol_hash = identical(
      dkfa_hash_file(protocol_v7_path),
      "c6d3a282b967bacef0f7da76a9192389e624e91075eccd3de4de648386f28487"
    ),
    v7_id = identical(
      v7$protocol_id, "dkge-functional-alignment-certification-v7"
    ),
    v7_status = identical(v7$status, "frozen_before_v7_results"),
    candidate_source_binding = identical(
      v7$candidate_source_tree_sha256, expected_source
    ) && (
      is.null(candidate_source_tree_sha256) || identical(
        candidate_source_tree_sha256, expected_source
      )
    ),
    supersedes_v6 = identical(
      v7$amendment$supersedes,
      "dkge-functional-alignment-certification-v6"
    ),
    v6_receipts = identical(
      v7$amendment$v6_protocol, "protocol-v6.json"
    ) && identical(
      v7$amendment$v6_protocol_sha256,
      v6_validation$protocol_v6_sha256
    ) && identical(
      v7$amendment$v6_supersession,
      "protocol-v6-supersession.json"
    ) && identical(
      v7$amendment$v6_supersession_sha256,
      supersession$protocol_v6_supersession_sha256
    ),
    inherited_iteration_cap = identical(
      v7$court$sinkhorn_max_iter, 50000L
    ) && identical(
      v7$court$sinkhorn_max_iter, v6$court$sinkhorn_max_iter
    ),
    no_new_protocol_amendment = grepl(
      "V7 makes no numerical or statistical protocol amendment",
      v7$amendment$no_second_numerical_amendment,
      fixed = TRUE
    ) && grepl(
      "will not be raised again",
      v7$amendment$no_second_numerical_amendment,
      fixed = TRUE
    ),
    v7_directory_label = identical(
      v7$court$implementation,
      paste(
        "The v1 court implementation and OAT grid are reused without changing",
        "estimands, signs, seeds, arms, or gates; outputs are written to a",
        "source-tree-addressed v7 directory."
      )
    ),
    exact_v6_to_v7_delta = identical(
      dkfa_canonicalize_json(v7), dkfa_canonicalize_json(expected)
    )
  )
  list(
    passed = all(checks), checks = checks,
    v6_checks = v6_validation$checks,
    v6_supersession_checks = supersession$checks,
    protocol_v2_sha256 = v6_validation$protocol_v2_sha256,
    protocol_v3_sha256 = v6_validation$protocol_v3_sha256,
    protocol_v3_supersession_sha256 =
      v6_validation$protocol_v3_supersession_sha256,
    protocol_v4_sha256 = v6_validation$protocol_v4_sha256,
    protocol_v4_invalidation_sha256 =
      v6_validation$protocol_v4_invalidation_sha256,
    v4_history_bundle = v6_validation$v4_history_bundle,
    protocol_v5_sha256 = v6_validation$protocol_v5_sha256,
    protocol_v5_supersession_sha256 =
      v6_validation$protocol_v5_supersession_sha256,
    v5_history_bundle = v6_validation$v5_history_bundle,
    protocol_v6_sha256 = v6_validation$protocol_v6_sha256,
    protocol_v6_supersession_sha256 =
      supersession$protocol_v6_supersession_sha256,
    v6_history_bundle = dkfa_review_evidence_bundle(
      dkfa_v6_history_paths(root)
    ),
    protocol_v7_sha256 = dkfa_hash_file(protocol_v7_path)
  )
}

dkfa_validate_v7_supersession <- function(
    supersession_path, protocol_v7_path,
    root = NULL, candidate_source_tree_sha256 = NULL,
    verify_retained_evidence = TRUE) {
  receipt <- jsonlite::read_json(supersession_path, simplifyVector = FALSE)
  expected_source <-
    "e38585975999d2bc6bdc1b76747e90a3ab89328940c92d676592494255831400"
  if (isTRUE(verify_retained_evidence) && is.null(root)) {
    stop(
      "V7 supersession validation requires the repository root so retained ",
      "execution evidence can be rehashed. Set `verify_retained_evidence = ",
      "FALSE` only for installed-layout structural checks.",
      call. = FALSE
    )
  }
  retained <- if (!isTRUE(verify_retained_evidence)) {
    NULL
  } else {
    list(
      source_test = dkfa_review_evidence_bundle(
        dkfa_v7_source_test_paths(root)
      ),
      court = dkfa_review_evidence_bundle(dkfa_v7_court_paths(root)),
      power = dkfa_review_evidence_bundle(dkfa_v7_power_paths(root))
    )
  }
  checks <- c(
    receipt_hash = identical(
      dkfa_hash_file(supersession_path),
      "ceea7aed678eabcd7683531dc295685d489e3223b1b5cb4f9d13aef393a024f6"
    ),
    status = identical(
      receipt$status,
      "superseded_after_execution_due_nonstatistical_package_check_repair"
    ),
    v7_protocol = identical(
      receipt$superseded_protocol$protocol_id,
      "dkge-functional-alignment-certification-v7"
    ) && identical(receipt$superseded_protocol$path, "protocol-v7.json") &&
      identical(
        receipt$superseded_protocol$sha256,
        "c6d3a282b967bacef0f7da76a9192389e624e91075eccd3de4de648386f28487"
      ) && identical(
        dkfa_hash_file(protocol_v7_path),
        receipt$superseded_protocol$sha256
      ) && identical(
        receipt$superseded_protocol$candidate_source_tree_sha256,
        "7bee7084bc4c66d29880efd3006c16d4118e38b08c3242127ab9d1b4e9f4a2b6"
      ),
    executed_before_supersession = isTRUE(
      receipt$materialization$independent_freeze_issued
    ) && isTRUE(receipt$materialization$source_test_suite_completed) &&
      isTRUE(receipt$materialization$court_analyzed_complete) &&
      isTRUE(receipt$materialization$known_warp_determinism_complete) &&
      isTRUE(receipt$materialization$statistical_outputs_materialized) &&
      isTRUE(receipt$materialization$statistical_results_observed) &&
      isFALSE(receipt$materialization$certification_collector_started) &&
      isFALSE(receipt$materialization$certification_record_materialized) &&
      isFALSE(
        receipt$materialization$statistical_results_used_to_change_successor_design
      ),
    retained_outcomes = identical(
      receipt$retained_evidence$source_test_bundle_sha256,
      "59cee5c123670d9db6bc51c9ae5cb18ecd00221d7d45007d7bb2862130f81984"
    ) && identical(
      receipt$retained_evidence$court_bundle_sha256,
      "3f19285f4f707abb56147106824e2fc50c79b03e3f8228b0eb54d713d73abd63"
    ) && identical(
      receipt$retained_evidence$known_warp_bundle_sha256,
      "326598b219a13ae4cadc7663bcd1d2e357e05a997e4dd7def9e52688ef1def1b"
    ) && isTRUE(receipt$retained_evidence$court_exact_oracle_gate) &&
      isFALSE(receipt$retained_evidence$court_negative_control_promotion_gate) &&
      identical(
        receipt$retained_evidence$court_status,
        "complete_approximate_only_inferential_promotion_blocked"
      ) && (
        is.null(retained) || (
          identical(
            retained$source_test$bundle_sha256,
            receipt$retained_evidence$source_test_bundle_sha256
          ) && identical(
            retained$court$bundle_sha256,
            receipt$retained_evidence$court_bundle_sha256
          ) && identical(
            retained$power$bundle_sha256,
            receipt$retained_evidence$known_warp_bundle_sha256
          )
        )
      ),
    package_check_failure_retained = identical(
      receipt$package_check_discovery$default_build_tarball_sha256_not_certified,
      "56d33163d384b4e7df1ee8b1a1898467e40b70493fb5c62ec43bc83af99c672a"
    ) && identical(
      receipt$package_check_discovery$closed_world_tarball_sha256_not_certified,
      "517ca26a38b5b3a24af7ac325f345db2dd378263130340bc4e309cce40083ec4"
    ),
    repaired_source = identical(
      receipt$source_change$candidate_source_tree_sha256_after_repair,
      expected_source
    ) && isFALSE(receipt$source_change$package_runtime_change) && (
      is.null(candidate_source_tree_sha256) || identical(
        candidate_source_tree_sha256, expected_source
      )
    ),
    successor = identical(
      receipt$successor$protocol_id,
      "dkge-functional-alignment-certification-v8"
    ) && identical(receipt$successor$path, "protocol.json") &&
      identical(
        receipt$successor$candidate_source_tree_sha256, expected_source
      ) && isFALSE(receipt$successor$numerical_protocol_change) &&
      isFALSE(receipt$successor$statistical_protocol_change)
  )
  list(
    passed = all(checks), checks = checks,
    protocol_v7_supersession_sha256 = dkfa_hash_file(supersession_path),
    receipt = receipt,
    retained_evidence = retained
  )
}

dkfa_validate_v8_protocol_amendment <- function(
    protocol_v2_path, protocol_v3_path, v3_supersession_path,
    protocol_v4_path, v4_invalidation_path, protocol_v5_path,
    v5_supersession_path, protocol_v6_path, v6_supersession_path,
    protocol_v7_path, v7_supersession_path, protocol_v8_path, root,
    candidate_source_tree_sha256 = NULL) {
  v7_validation <- dkfa_validate_v7_protocol_amendment(
    protocol_v2_path, protocol_v3_path, v3_supersession_path,
    protocol_v4_path, v4_invalidation_path, protocol_v5_path,
    v5_supersession_path, protocol_v6_path, v6_supersession_path,
    protocol_v7_path, root, candidate_source_tree_sha256 = NULL
  )
  supersession <- dkfa_validate_v7_supersession(
    v7_supersession_path, protocol_v7_path,
    root = root,
    candidate_source_tree_sha256 = candidate_source_tree_sha256
  )
  v7 <- jsonlite::read_json(protocol_v7_path, simplifyVector = FALSE)
  v8 <- jsonlite::read_json(protocol_v8_path, simplifyVector = FALSE)
  expected_source <-
    "e38585975999d2bc6bdc1b76747e90a3ab89328940c92d676592494255831400"
  expected <- v7
  expected$protocol_id <- "dkge-functional-alignment-certification-v8"
  expected$status <- "frozen_before_v8_results"
  expected$candidate_source_tree_sha256 <- expected_source
  expected$amendment <- v8$amendment
  expected$artifact_validation <- v8$artifact_validation
  expected$court$implementation <- v8$court$implementation
  checks <- c(
    v7_chain = isTRUE(v7_validation$passed),
    v7_supersession = isTRUE(supersession$passed),
    v8_protocol_hash = identical(
      dkfa_hash_file(protocol_v8_path),
      "0e68dd7558d3df28878b8ac9f0edb49147620d9dd5c518898ed6243f305876c3"
    ),
    v8_id = identical(
      v8$protocol_id, "dkge-functional-alignment-certification-v8"
    ),
    v8_status = identical(v8$status, "frozen_before_v8_results"),
    candidate_source_binding = identical(
      v8$candidate_source_tree_sha256, expected_source
    ) && (
      is.null(candidate_source_tree_sha256) || identical(
        candidate_source_tree_sha256, expected_source
      )
    ),
    supersedes_v7 = identical(
      v8$amendment$supersedes,
      "dkge-functional-alignment-certification-v7"
    ),
    v7_receipts = identical(
      v8$amendment$v7_protocol, "protocol-v7.json"
    ) && identical(
      v8$amendment$v7_protocol_sha256,
      v7_validation$protocol_v7_sha256
    ) && identical(
      v8$amendment$v7_supersession,
      "protocol-v7-supersession.json"
    ) && identical(
      v8$amendment$v7_supersession_sha256,
      supersession$protocol_v7_supersession_sha256
    ),
    inherited_iteration_cap = identical(
      v8$court$sinkhorn_max_iter, 50000L
    ) && identical(
      v8$court$sinkhorn_max_iter, v7$court$sinkhorn_max_iter
    ),
    no_new_protocol_amendment = grepl(
      "V8 makes no numerical or statistical protocol amendment",
      v8$amendment$no_second_numerical_amendment,
      fixed = TRUE
    ) && grepl(
      "V7 outcomes were observed but did not alter V8's design",
      v8$amendment$no_second_numerical_amendment,
      fixed = TRUE
    ) && grepl(
      "will not be raised again",
      v8$amendment$no_second_numerical_amendment,
      fixed = TRUE
    ),
    closed_world_artifact_policy = identical(
      v8$artifact_validation$certified_tarball_build,
      "R CMD build --no-build-vignettes"
    ) && identical(
      v8$artifact_validation$built_tarball_check,
      "R CMD check --no-manual --no-clean --ignore-vignettes"
    ),
    v8_directory_label = identical(
      v8$court$implementation,
      paste(
        "The v1 court implementation and OAT grid are reused without changing",
        "estimands, signs, seeds, arms, or gates; outputs are written to a",
        "source-tree-addressed v8 directory."
      )
    ),
    exact_v7_to_v8_delta = identical(
      dkfa_canonicalize_json(v8), dkfa_canonicalize_json(expected)
    )
  )
  list(
    passed = all(checks), checks = checks,
    v7_checks = v7_validation$checks,
    v7_supersession_checks = supersession$checks,
    protocol_v2_sha256 = v7_validation$protocol_v2_sha256,
    protocol_v3_sha256 = v7_validation$protocol_v3_sha256,
    protocol_v3_supersession_sha256 =
      v7_validation$protocol_v3_supersession_sha256,
    protocol_v4_sha256 = v7_validation$protocol_v4_sha256,
    protocol_v4_invalidation_sha256 =
      v7_validation$protocol_v4_invalidation_sha256,
    v4_history_bundle = v7_validation$v4_history_bundle,
    protocol_v5_sha256 = v7_validation$protocol_v5_sha256,
    protocol_v5_supersession_sha256 =
      v7_validation$protocol_v5_supersession_sha256,
    v5_history_bundle = v7_validation$v5_history_bundle,
    protocol_v6_sha256 = v7_validation$protocol_v6_sha256,
    protocol_v6_supersession_sha256 =
      v7_validation$protocol_v6_supersession_sha256,
    v6_history_bundle = v7_validation$v6_history_bundle,
    protocol_v7_sha256 = v7_validation$protocol_v7_sha256,
    protocol_v7_supersession_sha256 =
      supersession$protocol_v7_supersession_sha256,
    v7_history_bundle = dkfa_review_evidence_bundle(
      dkfa_v7_history_paths(root)
    ),
    protocol_v8_sha256 = dkfa_hash_file(protocol_v8_path)
  )
}

dkfa_validate_v8_supersession <- function(
    supersession_path, protocol_v8_path, root = NULL,
    candidate_source_tree_sha256 = NULL,
    verify_retained_evidence = TRUE) {
  receipt <- jsonlite::read_json(supersession_path, simplifyVector = FALSE)
  expected_source <-
    "e38585975999d2bc6bdc1b76747e90a3ab89328940c92d676592494255831400"
  if (isTRUE(verify_retained_evidence) && is.null(root)) {
    stop(
      "V8 supersession validation requires the repository root so retained ",
      "execution evidence can be rehashed. Set `verify_retained_evidence = ",
      "FALSE` only for installed-layout structural checks.",
      call. = FALSE
    )
  }
  retained <- if (!isTRUE(verify_retained_evidence)) {
    NULL
  } else {
    list(
      source_test = dkfa_review_evidence_bundle(
        dkfa_v8_source_test_paths(root)
      ),
      court = dkfa_review_evidence_bundle(dkfa_v8_court_paths(root)),
      power = dkfa_review_evidence_bundle(dkfa_v8_power_paths(root))
    )
  }
  candidate <- NULL
  if (!is.null(root) && dir.exists(file.path(root, "vignettes"))) {
    vignette_paths <- list.files(
      file.path(root, "vignettes"), pattern = "[.]Rmd$", full.names = TRUE
    )
    vignette_paths <- Filter(function(path) {
      any(grepl("^has_albers_v2 <-", readLines(path, warn = FALSE)))
    }, vignette_paths)
    candidate <- list(
      guard_count = length(vignette_paths),
      guard_bundle = dkfa_bundle_manifest(vignette_paths, root),
      docs = dkfa_tree_manifest(
        root, c("man", "vignettes", "README.md"),
        exclude = "(^|/)[.]DS_Store$"
      ),
      forbidden_comparator = any(vapply(
        vignette_paths,
        function(path) any(grepl(
          "utils::package_version", readLines(path, warn = FALSE),
          fixed = TRUE
        )),
        logical(1)
      ))
    )
  }
  checks <- c(
    receipt_hash = identical(
      dkfa_hash_file(supersession_path),
      "e55cba258474256e9d8d320a649a399ba1ef7e55855b91ad21494b3a6eec1e73"
    ),
    status = identical(
      receipt$status,
      "superseded_after_execution_due_documentation_gate_failure"
    ),
    v8_protocol = identical(
      receipt$superseded_protocol$protocol_id,
      "dkge-functional-alignment-certification-v8"
    ) && identical(
      receipt$superseded_protocol$path, "protocol-v8.json"
    ) && identical(
      receipt$superseded_protocol$sha256,
      "0e68dd7558d3df28878b8ac9f0edb49147620d9dd5c518898ed6243f305876c3"
    ) && identical(
      dkfa_hash_file(protocol_v8_path),
      receipt$superseded_protocol$sha256
    ) && identical(
      receipt$superseded_protocol$candidate_source_tree_sha256,
      expected_source
    ) && identical(
      receipt$superseded_protocol$test_validation_tree_sha256,
      "d881948e5df70f5c03bcb4880e503f64c9ec35e93c1bdbee19a3e88a3e5b41db"
    ),
    executed_before_supersession = isTRUE(
      receipt$materialization$independent_freeze_issued
    ) && isTRUE(receipt$materialization$source_test_suite_completed) &&
      isTRUE(receipt$materialization$court_analyzed_complete) &&
      isTRUE(receipt$materialization$known_warp_determinism_complete) &&
      isTRUE(receipt$materialization$statistical_outputs_materialized) &&
      isTRUE(receipt$materialization$statistical_results_observed) &&
      isTRUE(receipt$materialization$closed_world_tarball_built) &&
      isTRUE(receipt$materialization$exact_tarball_check_status_ok) &&
      isTRUE(receipt$materialization$pkgdown_gate_started) &&
      isFALSE(receipt$materialization$pkgdown_gate_completed) &&
      isFALSE(receipt$materialization$certification_collector_started) &&
      isFALSE(receipt$materialization$certification_record_materialized) &&
      isFALSE(receipt$materialization$inferential_promotion_granted) &&
      isFALSE(
        receipt$materialization$statistical_results_used_to_change_successor_design
      ),
    retained_outcomes = identical(
      receipt$retained_evidence$source_test_bundle_sha256,
      "f71926dc8fee41e1d9c7ae65c52e31161ecd8091ca5eaa83b102b4bfbf24cc06"
    ) && identical(
      receipt$retained_evidence$court_bundle_sha256,
      "63a5d7830cbccb11e210386bf5eab943a0961a643c690faadc1ba6f9d10fae3e"
    ) && identical(
      receipt$retained_evidence$known_warp_bundle_sha256,
      "007f7ecd9e534d38abe70ad87fbe0a0de30e3a8df10a1e727868da514a82c9b1"
    ) && isTRUE(receipt$retained_evidence$court_exact_oracle_gate) &&
      isFALSE(receipt$retained_evidence$court_negative_control_promotion_gate) &&
      identical(
        receipt$retained_evidence$court_status,
        "complete_approximate_only_inferential_promotion_blocked"
      ) && isTRUE(receipt$retained_evidence$v8_outputs_are_historical_only) &&
      isFALSE(receipt$retained_evidence$v8_outputs_reused_as_successor_evidence) &&
      (
        is.null(retained) || (
          identical(
            retained$source_test$bundle_sha256,
            receipt$retained_evidence$source_test_bundle_sha256
          ) && identical(
            retained$court$bundle_sha256,
            receipt$retained_evidence$court_bundle_sha256
          ) && identical(
            retained$power$bundle_sha256,
            receipt$retained_evidence$known_warp_bundle_sha256
          )
        )
      ),
    package_documentation_failure = identical(
      receipt$package_documentation_gate_failure$closed_world_tarball_sha256_not_certified,
      "ca062be141e79fad6bb667c954d6998e3dd7e27a379590c81c29944cbacddb57"
    ) && identical(
      receipt$package_documentation_gate_failure$exact_check_status, "OK"
    ) && isFALSE(
      receipt$package_documentation_gate_failure$pkgdown_receipt_materialized
    ) && identical(
      receipt$package_documentation_gate_failure$failure_signature,
      "'package_version' is not an exported object from 'namespace:utils'"
    ) && isFALSE(
      receipt$package_documentation_gate_failure$durable_wrapper_log_or_receipt_retained
    ),
    documentation_only_repair = identical(
      receipt$candidate_change$candidate_source_tree_sha256_before_repair,
      expected_source
    ) && identical(
      receipt$candidate_change$candidate_source_tree_sha256_after_repair,
      expected_source
    ) && identical(
      receipt$candidate_change$packaged_docs_tree_sha256_before_repair,
      "4ab1c477cdb878d7ada00db9e51de0fcc4243792e014f715c2d5b402f113f28d"
    ) && identical(
      receipt$candidate_change$packaged_docs_tree_sha256_after_repair,
      "eec80acfc2102b6f97d2c6b1f388bb173c90616d32681c6487a91383fb3ada5a"
    ) && identical(receipt$candidate_change$vignette_guard_count, 19L) &&
      identical(
        receipt$candidate_change$vignette_guard_bundle_sha256_before_repair,
        "39ae42b8e2727bf0fb98257a631ea6ca0a46abfa6f9353e8c99021f83bbb81f3"
      ) && identical(
        receipt$candidate_change$vignette_guard_bundle_sha256_after_repair,
        "cedbc8ed92a8cfce3d4ac004feddba2b66e9e814b628a7a26f05b54b1ee42099"
      ) && isFALSE(receipt$candidate_change$package_runtime_change) &&
      isFALSE(receipt$candidate_change$executable_r_or_cpp_change) &&
      isFALSE(receipt$candidate_change$court_or_power_statistical_helper_change) &&
      isFALSE(receipt$candidate_change$numerical_protocol_change) &&
      isFALSE(receipt$candidate_change$statistical_protocol_change) &&
      (
        is.null(candidate) || (
          identical(candidate$guard_count, 19L) &&
          identical(
            candidate$guard_bundle$bundle_sha256,
            receipt$candidate_change$vignette_guard_bundle_sha256_after_repair
          ) && identical(
            candidate$docs$tree_sha256,
            receipt$candidate_change$packaged_docs_tree_sha256_after_repair
          ) && !isTRUE(candidate$forbidden_comparator)
        )
      ) && (
        is.null(candidate_source_tree_sha256) || identical(
          candidate_source_tree_sha256, expected_source
        )
      ),
    successor = identical(
      receipt$successor$protocol_id,
      "dkge-functional-alignment-certification-v9"
    ) && identical(receipt$successor$path, "protocol.json") &&
      identical(receipt$successor$candidate_source_tree_sha256, expected_source) &&
      isTRUE(receipt$successor$fresh_source_court_power_required) &&
      isFALSE(receipt$successor$numerical_protocol_change) &&
      isFALSE(receipt$successor$statistical_protocol_change)
  )
  list(
    passed = all(checks), checks = checks,
    protocol_v8_supersession_sha256 = dkfa_hash_file(supersession_path),
    receipt = receipt, retained_evidence = retained,
    candidate = candidate
  )
}

dkfa_validate_v9_protocol_amendment <- function(
    protocol_v2_path, protocol_v3_path, v3_supersession_path,
    protocol_v4_path, v4_invalidation_path, protocol_v5_path,
    v5_supersession_path, protocol_v6_path, v6_supersession_path,
    protocol_v7_path, v7_supersession_path, protocol_v8_path,
    v8_supersession_path, protocol_v9_path, root,
    candidate_source_tree_sha256 = NULL) {
  v8_validation <- dkfa_validate_v8_protocol_amendment(
    protocol_v2_path, protocol_v3_path, v3_supersession_path,
    protocol_v4_path, v4_invalidation_path, protocol_v5_path,
    v5_supersession_path, protocol_v6_path, v6_supersession_path,
    protocol_v7_path, v7_supersession_path, protocol_v8_path, root,
    candidate_source_tree_sha256 = NULL
  )
  supersession <- dkfa_validate_v8_supersession(
    v8_supersession_path, protocol_v8_path, root = root,
    candidate_source_tree_sha256 = candidate_source_tree_sha256
  )
  v8 <- jsonlite::read_json(protocol_v8_path, simplifyVector = FALSE)
  v9 <- jsonlite::read_json(protocol_v9_path, simplifyVector = FALSE)
  expected_source <-
    "e38585975999d2bc6bdc1b76747e90a3ab89328940c92d676592494255831400"
  expected <- v8
  expected$protocol_id <- "dkge-functional-alignment-certification-v9"
  expected$status <- "frozen_before_v9_results"
  expected$candidate_source_tree_sha256 <- expected_source
  expected$amendment <- v9$amendment
  expected$court$implementation <- v9$court$implementation
  checks <- c(
    v8_chain = isTRUE(v8_validation$passed),
    v8_supersession = isTRUE(supersession$passed),
    v9_protocol_hash = identical(
      dkfa_hash_file(protocol_v9_path),
      "4e6b828f463d9ae57de5a100f440dc7e96a376da3fdcae34509a42362223ac0a"
    ),
    v9_id = identical(
      v9$protocol_id, "dkge-functional-alignment-certification-v9"
    ),
    v9_status = identical(v9$status, "frozen_before_v9_results"),
    candidate_source_binding = identical(
      v9$candidate_source_tree_sha256, expected_source
    ) && (
      is.null(candidate_source_tree_sha256) || identical(
        candidate_source_tree_sha256, expected_source
      )
    ),
    supersedes_v8 = identical(
      v9$amendment$supersedes,
      "dkge-functional-alignment-certification-v8"
    ),
    v8_receipts = identical(
      v9$amendment$v8_protocol, "protocol-v8.json"
    ) && identical(
      v9$amendment$v8_protocol_sha256,
      v8_validation$protocol_v8_sha256
    ) && identical(
      v9$amendment$v8_supersession,
      "protocol-v8-supersession.json"
    ) && identical(
      v9$amendment$v8_supersession_sha256,
      supersession$protocol_v8_supersession_sha256
    ),
    inherited_iteration_cap = identical(
      v9$court$sinkhorn_max_iter, 50000L
    ) && identical(
      v9$court$sinkhorn_max_iter, v8$court$sinkhorn_max_iter
    ),
    no_new_protocol_amendment = grepl(
      "V9 makes no numerical or statistical protocol amendment",
      v9$amendment$no_second_numerical_amendment,
      fixed = TRUE
    ) && grepl(
      "V8 outcomes were observed but did not alter V9's design",
      v9$amendment$no_second_numerical_amendment,
      fixed = TRUE
    ) && grepl(
      "will not be reused",
      v9$amendment$no_second_numerical_amendment,
      fixed = TRUE
    ) && grepl(
      "will not be raised again",
      v9$amendment$no_second_numerical_amendment,
      fixed = TRUE
    ),
    exact_documentation_repair = grepl(
      "utils::package_version",
      v9$amendment$candidate_change_before_results,
      fixed = TRUE
    ) && grepl(
      "base::package_version",
      v9$amendment$candidate_change_before_results,
      fixed = TRUE
    ) && grepl(
      "exactly 19 vignette",
      v9$amendment$candidate_change_before_results,
      fixed = TRUE
    ),
    v9_directory_label = identical(
      v9$court$implementation,
      paste(
        "The v1 court implementation and OAT grid are reused without changing",
        "estimands, signs, seeds, arms, or gates; outputs are written to a",
        "source-tree-addressed v9 directory."
      )
    ),
    exact_v8_to_v9_delta = identical(
      dkfa_canonicalize_json(v9), dkfa_canonicalize_json(expected)
    )
  )
  list(
    passed = all(checks), checks = checks,
    v8_checks = v8_validation$checks,
    v8_supersession_checks = supersession$checks,
    protocol_v2_sha256 = v8_validation$protocol_v2_sha256,
    protocol_v3_sha256 = v8_validation$protocol_v3_sha256,
    protocol_v3_supersession_sha256 =
      v8_validation$protocol_v3_supersession_sha256,
    protocol_v4_sha256 = v8_validation$protocol_v4_sha256,
    protocol_v4_invalidation_sha256 =
      v8_validation$protocol_v4_invalidation_sha256,
    v4_history_bundle = v8_validation$v4_history_bundle,
    protocol_v5_sha256 = v8_validation$protocol_v5_sha256,
    protocol_v5_supersession_sha256 =
      v8_validation$protocol_v5_supersession_sha256,
    v5_history_bundle = v8_validation$v5_history_bundle,
    protocol_v6_sha256 = v8_validation$protocol_v6_sha256,
    protocol_v6_supersession_sha256 =
      v8_validation$protocol_v6_supersession_sha256,
    v6_history_bundle = v8_validation$v6_history_bundle,
    protocol_v7_sha256 = v8_validation$protocol_v7_sha256,
    protocol_v7_supersession_sha256 =
      v8_validation$protocol_v7_supersession_sha256,
    v7_history_bundle = v8_validation$v7_history_bundle,
    protocol_v8_sha256 = v8_validation$protocol_v8_sha256,
    protocol_v8_supersession_sha256 =
      supersession$protocol_v8_supersession_sha256,
    v8_history_bundle = dkfa_review_evidence_bundle(
      dkfa_v8_history_paths(root)
    ),
    protocol_v9_sha256 = dkfa_hash_file(protocol_v9_path)
  )
}

dkfa_validate_current_protocol <- function(
    root, candidate_source_tree_sha256 = NULL) {
  certification <- file.path(
    root, "inst", "validation", "functional-alignment-certification"
  )
  dkfa_validate_v9_protocol_amendment(
    file.path(certification, "protocol-v2.json"),
    file.path(certification, "protocol-v3.json"),
    file.path(certification, "protocol-v3-supersession.json"),
    file.path(certification, "protocol-v4.json"),
    file.path(certification, "protocol-v4-invalidation.json"),
    file.path(certification, "protocol-v5.json"),
    file.path(certification, "protocol-v5-supersession.json"),
    file.path(certification, "protocol-v6.json"),
    file.path(certification, "protocol-v6-supersession.json"),
    file.path(certification, "protocol-v7.json"),
    file.path(certification, "protocol-v7-supersession.json"),
    file.path(certification, "protocol-v8.json"),
    file.path(certification, "protocol-v8-supersession.json"),
    file.path(certification, "protocol.json"),
    root,
    candidate_source_tree_sha256 = candidate_source_tree_sha256
  )
}

dkfa_validate_court_v2_invalidation <- function(
    invalidation_path, fixture_path, protocol_v2_path,
    historical_paths = NULL, verify_solver = FALSE,
    verify_origin = FALSE) {
  invalidation <- jsonlite::read_json(
    invalidation_path, simplifyVector = TRUE
  )
  fixture <- jsonlite::read_json(fixture_path, simplifyVector = TRUE)
  official_invalidation_hash <-
    "08e4ff679ebf34cb2d1a8e9408691920a4ae7cfc6ecb93b02e2b10ea8472866c"
  official_fixture_hash <-
    "628031d4361e734a5dbbf1860efe55c8081913882f47f914ee81fa0ec8aabb13"
  official_protocol_hash <-
    "9ce7d46a4e2efd98f9319f157dd477a78125972e7a4e961aa772a8377541b25e"
  expected_historical <- c(
    source_test_manifest =
      "868f9a4bedb92bb2e9d4e1be1fd9827d8dc01137bfeef127effba6aed8413cc6",
    source_test_results =
      "1eff5db2e88af04d3b59dc99824b756f310d78f1ce3513841fc877b53a73d141",
    source_test_session_info =
      "7cda8fe33c3c20fb2eaf70959496c69c9bca21a7ea3e865dc88a54429f433b1a",
    failed_court_formal_budget =
      "c553e88d2e8a4876063b92fab74d1ebc4aa88dc5feea3ac54fb2bee523034052"
  )
  expected_historical_paths <- c(
    source_test_manifest = paste0(
      "data-raw/functional-alignment-certification/",
      "package-0be19103609aa074-b143f64a3f6e/source-test-manifest.json"
    ),
    source_test_results = paste0(
      "data-raw/functional-alignment-certification/",
      "package-0be19103609aa074-b143f64a3f6e/source-test-results.csv"
    ),
    source_test_session_info = paste0(
      "data-raw/functional-alignment-certification/",
      "package-0be19103609aa074-b143f64a3f6e/session-info.txt"
    ),
    failed_court_formal_budget = paste0(
      "data-raw/functional-alignment-certification/",
      "court-0be19103609aa074-30c9732ef157/formal-budget.csv"
    )
  )
  historical_receipts <- invalidation$historical_evidence
  receipt_value <- function(role, field) {
    record <- historical_receipts[[role]]
    if (is.null(record) || is.null(record[[field]])) return(NA_character_)
    as.character(record[[field]])
  }
  receipt_hashes <- vapply(
    names(expected_historical),
    receipt_value, character(1), field = "sha256"
  )
  receipt_paths <- vapply(
    names(expected_historical),
    receipt_value, character(1), field = "path"
  )
  checks <- c(
    invalidation_file = identical(
      dkfa_hash_file(invalidation_path), official_invalidation_hash
    ),
    fixture_file = identical(
      dkfa_hash_file(fixture_path), official_fixture_hash
    ),
    protocol_v2_file = identical(
      dkfa_hash_file(protocol_v2_path), official_protocol_hash
    ),
    status = identical(
      invalidation$status, "invalidated_before_statistical_outputs"
    ),
    protocol_id = identical(
      invalidation$v2$protocol_id,
      "dkge-functional-alignment-certification-v2"
    ),
    protocol_receipt = identical(
      invalidation$v2$protocol_sha256, official_protocol_hash
    ),
    fixture_receipt = identical(
      invalidation$fixture$sha256, official_fixture_hash
    ),
    source_receipt = identical(
      invalidation$v2$source_tree_sha256,
      "0be19103609aa074d2ff86c15c27f6f701442623800a10cbe599b01c141ae87f"
    ),
    test_validation_receipt = identical(
      invalidation$v2$test_validation_tree_sha256,
      "b143f64a3f6e1d17af0d1d30153eae4f885e800c7a1b63da3a3b8e5a466f539d"
    ),
    v2_harness_receipt = identical(
      invalidation$v2$court_harness_bundle_sha256,
      "30c9732ef157fc6ba7c543865c3b6e0d3eb7c38ede5ed931a5026314704ce4b7"
    ),
    historical_receipt_hashes = identical(
      unname(receipt_hashes), unname(expected_historical)
    ),
    historical_receipt_paths = identical(
      unname(receipt_paths), unname(expected_historical_paths)
    ),
    source_manifest_receipt = identical(
      invalidation$v2$source_test_manifest_sha256,
      unname(expected_historical[["source_test_manifest"]])
    ),
    failure_identity = identical(
      unname(unlist(invalidation$representative_failure[
        c("cell", "sim", "seed", "arm", "subject_index")
      ])),
      unname(unlist(list(
        cell = "baseline", sim = 9L, seed = 7310009L,
        arm = "geometry_only", subject_index = 6L
      )))
    ),
    numerical_failure = identical(
      invalidation$representative_failure$max_iter, 5000L
    ) && identical(
      invalidation$representative_failure$tolerance, 1e-4
    ) && invalidation$representative_failure$marginal_error >
      invalidation$representative_failure$tolerance &&
      !isTRUE(invalidation$representative_failure$converged),
    bitwise_cap_equivalence = isTRUE(
      invalidation$numerical_amendment_evidence$plans_bitwise_identical
    ) && isTRUE(
      invalidation$numerical_amendment_evidence$operators_bitwise_identical
    ) && identical(
      invalidation$numerical_amendment_evidence$max_iter_20000$plan_sha256,
      invalidation$numerical_amendment_evidence$max_iter_50000$plan_sha256
    ) && identical(
      invalidation$numerical_amendment_evidence$max_iter_20000$operator_sha256,
      invalidation$numerical_amendment_evidence$max_iter_50000$operator_sha256
    ),
    no_statistical_outputs = !isTRUE(
      invalidation$materialization$statistical_outputs_materialized
    ) && !isTRUE(invalidation$materialization$formal_raw_present) &&
      !isTRUE(invalidation$materialization$formal_summary_present) &&
      !isTRUE(invalidation$materialization$formal_manifest_present) &&
      !isTRUE(
        invalidation$materialization$statistical_results_used_for_amendment
      ),
    fixture_identity = identical(fixture$source$cell, "baseline") &&
      identical(fixture$source$sim, 9L) &&
      identical(fixture$source$seed, 7310009L) &&
      identical(fixture$source$arm, "geometry_only") &&
      identical(fixture$source$subject_index, 6L) &&
      identical(fixture$solver$failed_max_iter, 5000L) &&
      identical(fixture$solver$amended_max_iter, 50000L)
  )

  historical_checks <- logical()
  if (!is.null(historical_paths)) {
    if (is.null(names(historical_paths)) ||
        !setequal(names(historical_paths), names(expected_historical))) {
      historical_checks <- c(historical_path_roles = FALSE)
    } else {
      historical_paths <- historical_paths[names(expected_historical)]
      present <- file.exists(historical_paths) & !dir.exists(historical_paths)
      live_hashes <- rep(NA_character_, length(historical_paths))
      names(live_hashes) <- names(historical_paths)
      if (all(present)) {
        live_hashes <- vapply(historical_paths, dkfa_hash_file, character(1))
      }
      content_ok <- FALSE
      if (all(present) && identical(
          unname(live_hashes), unname(expected_historical))) {
        source_manifest <- jsonlite::read_json(
          historical_paths[["source_test_manifest"]], simplifyVector = TRUE
        )
        failed_budget <- utils::read.csv(
          historical_paths[["failed_court_formal_budget"]],
          stringsAsFactors = FALSE
        )
        content_ok <- isTRUE(source_manifest$passed) &&
          identical(
            source_manifest$source_tree_sha256,
            invalidation$v2$source_tree_sha256
          ) && identical(
            source_manifest$test_validation_tree_sha256,
            invalidation$v2$test_validation_tree_sha256
          ) && identical(
            source_manifest$results_sha256,
            expected_historical[["source_test_results"]]
          ) && identical(
            source_manifest$session_info_sha256,
            expected_historical[["source_test_session_info"]]
          ) && nrow(failed_budget) == 1L &&
          identical(failed_budget$n_sim[[1L]], 40L) &&
          identical(failed_budget$n_perm[[1L]], 31L) &&
          isTRUE(failed_budget$resource_limited[[1L]]) &&
          identical(failed_budget$exact_cells[[1L]], 6L) &&
          identical(failed_budget$formal_cores[[1L]], 4L)
      }
      historical_checks <- c(
        historical_files_present = all(present),
        historical_file_hashes = identical(
          unname(live_hashes), unname(expected_historical)
        ),
        historical_file_contents = content_ok
      )
    }
  }

  solver_checks <- logical()
  if (isTRUE(verify_solver)) {
    sinkhorn <- getFromNamespace(".dkge_sinkhorn_plan", "dkge")
    make_operator <- getFromNamespace(".dkge_transport_operator", "dkge")
    C <- as.matrix(fixture$cost)
    mu <- as.numeric(fixture$mu)
    nu <- as.numeric(fixture$nu)
    solve_at <- function(max_iter) suppressWarnings(sinkhorn(
      C, mu, nu,
      epsilon = fixture$solver$epsilon,
      max_iter = max_iter,
      tol = fixture$solver$tolerance,
      warm_start = FALSE,
      return_diagnostics = TRUE
    ))
    low <- solve_at(5000L)
    medium <- solve_at(20000L)
    high <- solve_at(50000L)
    plan_hash <- function(x) digest::digest(x$plan, algo = "sha256")
    operator_hash <- function(x) digest::digest(
      make_operator(x$plan, mu, nu, "intensive"), algo = "sha256"
    )
    amendment <- invalidation$numerical_amendment_evidence
    solver_checks <- c(
      fails_at_v2_cap = !isTRUE(low$diagnostics$converged) &&
        low$diagnostics$marginal_error > fixture$solver$tolerance,
      converges_at_20000 = isTRUE(medium$diagnostics$converged),
      converges_at_v3_cap = isTRUE(high$diagnostics$converged),
      plans_bitwise_identical = identical(medium$plan, high$plan),
      # The official invalidation file hash above binds the historical
      # platform's bitwise plan and operator hashes.  A fresh solve must be
      # bitwise stable across the two iteration caps on *this* platform, but
      # its floating-point bytes are not required to match another libm/compiler
      # build.  Requiring both recorded cap hashes to agree retains the receipt
      # mutation check without misrepresenting cross-platform reproducibility.
      plan_hashes = identical(plan_hash(medium), plan_hash(high)) &&
        identical(
          amendment$max_iter_20000$plan_sha256,
          amendment$max_iter_50000$plan_sha256
        ),
      operator_hashes = identical(operator_hash(medium), operator_hash(high)) &&
        identical(
          amendment$max_iter_20000$operator_sha256,
          amendment$max_iter_50000$operator_sha256
        )
    )
  }
  origin_checks <- logical()
  if (isTRUE(verify_origin)) {
    caller <- parent.frame()
    base_config <- get0(
      "fa_court_base_config", envir = caller, mode = "function",
      inherits = TRUE
    )
    generate <- get0(
      "fa_court_generate", envir = caller, mode = "function",
      inherits = TRUE
    )
    available <- c(
      fa_court_base_config = is.function(base_config),
      fa_court_generate = is.function(generate)
    )
    reconstructed <- NULL
    if (all(available)) {
      cfg <- base_config()
      dat <- generate(cfg, fixture$source$seed)
      source_index <- fixture$source$subject_index
      reference_index <- cfg$reference_subject
      normalize_mass <- function(x) as.numeric(x) / sum(x)
      reconstructed_mu <- normalize_mass(dat$sizes[[source_index]])
      reconstructed_nu <- normalize_mass(dat$sizes[[reference_index]])
      source_xyz <- dat$centroids[[source_index]] / cfg$sigma_spatial
      target_xyz <- dat$centroids[[reference_index]] / cfg$sigma_spatial
      reconstructed_cost <- matrix(
        rowSums(source_xyz * source_xyz),
        nrow(source_xyz), nrow(target_xyz)
      ) + matrix(
        rowSums(target_xyz * target_xyz),
        nrow(source_xyz), nrow(target_xyz), byrow = TRUE
      ) - 2 * tcrossprod(source_xyz, target_xyz)
      reconstructed_cost[reconstructed_cost < 0] <- 0
      reconstructed <- list(
        cost = reconstructed_cost,
        mu = reconstructed_mu,
        nu = reconstructed_nu
      )
    }
    numerically_equal <- function(x, y) isTRUE(all.equal(
      as.numeric(x), as.numeric(y),
      tolerance = 64 * .Machine$double.eps,
      check.attributes = FALSE
    ))
    origin_checks <- c(
      fixture_generators_available = all(available),
      fixture_cost_reconstructed = !is.null(reconstructed) &&
        numerically_equal(reconstructed$cost, as.matrix(fixture$cost)),
      fixture_source_mass_reconstructed = !is.null(reconstructed) &&
        numerically_equal(reconstructed$mu, fixture$mu),
      fixture_target_mass_reconstructed = !is.null(reconstructed) &&
        numerically_equal(reconstructed$nu, fixture$nu)
    )
  }
  checks <- c(checks, historical_checks, solver_checks, origin_checks)
  list(
    passed = all(checks), checks = checks,
    invalidation = invalidation, fixture = fixture
  )
}
