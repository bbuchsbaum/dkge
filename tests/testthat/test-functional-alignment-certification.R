library(testthat)

certification_root <- dkge_validation_path(
  "functional-alignment-certification"
)
certification_repo_root <- normalizePath(
  file.path(certification_root, "..", "..", ".."),
  mustWork = FALSE
)
court_source <- file.path(
  dirname(certification_root), "functional-alignment-court", "court.R"
)
if (file.exists(court_source)) source(court_source, local = TRUE)
power_source <- file.path(certification_root, "known-warp-power.R")
if (file.exists(power_source)) source(power_source, local = TRUE)
evidence_utils_source <- file.path(certification_root, "evidence-utils.R")
if (file.exists(evidence_utils_source)) {
  source(evidence_utils_source, local = TRUE)
}
package_provenance_source <- file.path(
  certification_root, "package-provenance.R"
)
if (file.exists(package_provenance_source)) {
  source(package_provenance_source, local = TRUE)
}

test_that("certification executors fail closed without a source loader", {
  skip_if_not(file.exists(evidence_utils_source))
  expect_error(
    dkfa_load_source_package(tempdir(), source_loader = "unavailable"),
    "installed dkge namespace is not an admissible fallback"
  )
  captured <- NULL
  source_loader <- dkfa_load_source_package(
    tempdir(),
    load_all = function(path, quiet, compile) {
      captured <<- list(path = path, quiet = quiet, compile = compile)
      invisible(NULL)
    }
  )
  expect_true(captured$compile)
  expect_true(source_loader$compile)
  expect_identical(source_loader$loader_version, "test-double")
  executor_paths <- c(
    file.path(certification_root, "run-court.R"),
    file.path(certification_root, "run-known-warp-power.R"),
    file.path(certification_root, "verify-known-warp-determinism.R"),
    file.path(certification_root, "run-source-tests.R"),
    court_source,
    file.path(
      dirname(certification_root), "functional-alignment-template",
      "run-template-benefit.R"
    )
  )
  expect_true(all(file.exists(executor_paths)))
  for (path in executor_paths) {
    code <- paste(readLines(path, warn = FALSE), collapse = "\n")
    expect_match(code, "dkfa_load_source_package", fixed = TRUE)
    expect_false(grepl("library(dkge)", code, fixed = TRUE))
  }
  bound_paths <- c(
    file.path(certification_root, "run-court.R"),
    file.path(certification_root, "analyze-court.R"),
    file.path(certification_root, "run-known-warp-power.R"),
    file.path(certification_root, "verify-known-warp-determinism.R"),
    file.path(certification_root, "run-source-tests.R"),
    file.path(certification_root, "run-pkgdown.R"),
    file.path(certification_root, "collect-certification.R")
  )
  for (path in bound_paths) {
    code <- paste(readLines(path, warn = FALSE), collapse = "\n")
    expect_match(code, "dkfa_bind_loaded_inputs", fixed = TRUE)
    expect_match(code, "dkfa_assert_input_snapshot", fixed = TRUE)
  }
})

test_that("two-sided loading rejects bytes changed by the loader", {
  root <- tempfile("dkfa-loader-race-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  covered <- file.path(root, "covered.R")
  writeLines("value <- 1", covered)
  snapshot <- function() dkfa_bundle_manifest(covered, root)
  expect_error(
    dkfa_bind_loaded_inputs(
      snapshot,
      function() {
        dkfa_load_source_package(
          root,
          load_all = function(path, quiet, compile) {
            expect_true(compile)
            writeLines("value <- 2", covered)
          }
        )
      },
      "Adversarial loader"
    ),
    "inputs changed while executable code was being loaded"
  )
})

test_that("v4 preserves the byte-frozen v3 amendment and binds source drift", {
  protocol <- file.path(certification_root, "protocol-v4.json")
  protocol_v2 <- file.path(certification_root, "protocol-v2.json")
  protocol_v3 <- file.path(certification_root, "protocol-v3.json")
  supersession <- file.path(
    certification_root, "protocol-v3-supersession.json"
  )
  expect_true(all(file.exists(c(
    protocol, protocol_v2, protocol_v3, supersession
  ))))
  expect_identical(
    dkfa_hash_file(protocol_v3),
    "f2b843c630ec2166a152b306938b105e42a9e06203fc0cd50e3855d343d78dc7"
  )
  parsed <- jsonlite::read_json(protocol, simplifyVector = TRUE)
  expect_identical(
    parsed$protocol_id, "dkge-functional-alignment-certification-v4"
  )
  amendment <- dkfa_validate_v4_protocol_amendment(
    protocol_v2, protocol_v3, supersession, protocol,
    candidate_source_tree_sha256 =
      "17e10ec404ab7758ada4591522ec415a98cc6c93e7a203d0ba6b42ac9aee066f",
    spatial_path = NULL
  )
  expect_true(amendment$passed, info = paste(
    names(amendment$checks)[!amendment$checks], collapse = ", "
  ))
  expect_true(all(amendment$v3_checks))
  expect_true(all(amendment$supersession_checks))
  expect_identical(parsed$court$n_sim, 40L)
  expect_identical(parsed$court$n_perm, 31L)
  expect_identical(parsed$court$sinkhorn_max_iter, 50000L)
  expect_identical(
    parsed$candidate_source_tree_sha256,
    "17e10ec404ab7758ada4591522ec415a98cc6c93e7a203d0ba6b42ac9aee066f"
  )
  expect_match(parsed$court$implementation,
               "source-tree-addressed v4 directory")
  expect_match(parsed$amendment$protocol_change,
               "version label from v3 to v4")
  expect_match(parsed$amendment$no_second_numerical_amendment,
               "V4 makes no numerical protocol amendment")
  expect_identical(
    fa_court_mapper(fa_court_base_config())$params$max_iter, 5000L
  )
  court_runner <- paste(readLines(
    file.path(certification_root, "run-court.R"), warn = FALSE
  ), collapse = "\n")
  expect_match(court_runner, "mapper$params$max_iter <- 50000L", fixed = TRUE)
  expect_identical(parsed$known_warp_power$seeds$n_cohorts, 40L)
  expect_setequal(parsed$court$arms_of_record, c(
    "geometry_only", "independent_alignment",
    "kernel_image_residual_prototype", "full_permutation_reestimation"
  ))
  expect_match(parsed$known_warp_power$power_role,
               "No power superiority threshold")
  expect_match(parsed$candidate_contract$release_claim,
               "Approximate inference requires an explicit override")

  for (mutation in c("analytical", "metadata", "source_binding")) {
    payload <- jsonlite::read_json(protocol, simplifyVector = FALSE)
    if (identical(mutation, "analytical")) {
      payload$court$n_perm <- 63L
    } else if (identical(mutation, "metadata")) {
      payload$amendment$protocol_change <- paste(
        payload$amendment$protocol_change, "mutated"
      )
    } else {
      payload$candidate_source_tree_sha256 <- paste0(
        "0", payload$candidate_source_tree_sha256
      )
    }
    candidate <- tempfile("dkfa-v4-protocol-", fileext = ".json")
    jsonlite::write_json(
      payload, candidate, pretty = TRUE, auto_unbox = TRUE
    )
    expect_false(dkfa_validate_v4_protocol_amendment(
      protocol_v2, protocol_v3, supersession, candidate
    )$passed, info = paste("protocol mutation must fail:", mutation))
  }

  bad_v3 <- tempfile("dkfa-v3-protocol-", fileext = ".json")
  v3_payload <- jsonlite::read_json(protocol_v3, simplifyVector = FALSE)
  v3_payload$court$n_perm <- 63L
  jsonlite::write_json(v3_payload, bad_v3, pretty = TRUE, auto_unbox = TRUE)
  expect_false(dkfa_validate_v3_protocol_amendment(
    protocol_v2, bad_v3
  )$passed)

  bad_supersession <- tempfile("dkfa-v3-supersession-", fileext = ".json")
  supersession_payload <- jsonlite::read_json(
    supersession, simplifyVector = FALSE
  )
  supersession_payload$materialization$court_started <- TRUE
  jsonlite::write_json(
    supersession_payload, bad_supersession,
    pretty = TRUE, auto_unbox = TRUE
  )
  expect_false(dkfa_validate_v4_protocol_amendment(
    protocol_v2, protocol_v3, bad_supersession, protocol
  )$passed)
  expect_false(dkfa_validate_v4_protocol_amendment(
    protocol_v2, protocol_v3, supersession, protocol,
    candidate_source_tree_sha256 = "wrong-source"
  )$passed)
})

test_that("v4 invalidation and v5-v9 succession are immutable", {
  paths <- list(
    v2 = file.path(certification_root, "protocol-v2.json"),
    v3 = file.path(certification_root, "protocol-v3.json"),
    v3_sup = file.path(certification_root, "protocol-v3-supersession.json"),
    v4 = file.path(certification_root, "protocol-v4.json"),
    v4_inv = file.path(certification_root, "protocol-v4-invalidation.json"),
    v5 = file.path(certification_root, "protocol-v5.json"),
    v5_sup = file.path(certification_root, "protocol-v5-supersession.json"),
    v6 = file.path(certification_root, "protocol-v6.json"),
    v6_sup = file.path(certification_root, "protocol-v6-supersession.json"),
    v7 = file.path(certification_root, "protocol-v7.json"),
    v7_sup = file.path(certification_root, "protocol-v7-supersession.json"),
    v8 = file.path(certification_root, "protocol-v8.json"),
    v8_sup = file.path(certification_root, "protocol-v8-supersession.json"),
    v9 = file.path(certification_root, "protocol.json")
  )
  expect_true(all(file.exists(unlist(paths, use.names = FALSE))))

  v5_sup <- dkfa_validate_v5_supersession(paths$v5_sup, paths$v5)
  v6_sup <- dkfa_validate_v6_supersession(paths$v6_sup, paths$v6)
  v7_sup <- dkfa_validate_v7_supersession(
    paths$v7_sup, paths$v7,
    candidate_source_tree_sha256 =
      "e38585975999d2bc6bdc1b76747e90a3ab89328940c92d676592494255831400",
    verify_retained_evidence = FALSE
  )
  v8_sup <- dkfa_validate_v8_supersession(
    paths$v8_sup, paths$v8,
    candidate_source_tree_sha256 =
      "e38585975999d2bc6bdc1b76747e90a3ab89328940c92d676592494255831400",
    verify_retained_evidence = FALSE
  )
  expect_true(v5_sup$passed)
  expect_true(v6_sup$passed)
  expect_true(v7_sup$passed)
  expect_true(v8_sup$passed)

  skip_if_not(
    dir.exists(file.path(
      certification_repo_root, "data-raw",
      "functional-alignment-certification"
    )),
    "repository-only historical evidence is not shipped in the package tarball"
  )

  invalidation <- dkfa_validate_v4_invalidation(
    paths$v4_inv, paths$v4, certification_repo_root
  )
  expect_true(invalidation$passed, info = paste(
    names(invalidation$checks)[!invalidation$checks], collapse = ", "
  ))
  expect_identical(
    invalidation$protocol_v4_invalidation_sha256,
    "470ff4158791bb5bc1e1c598f676201913cc0534532154bbd5a52f1de321a369"
  )

  current <- dkfa_validate_current_protocol(
    certification_repo_root,
    candidate_source_tree_sha256 =
      "e38585975999d2bc6bdc1b76747e90a3ab89328940c92d676592494255831400"
  )
  # The historical V9 protocol is source-bound. This successor changes
  # executable source and packaged documentation, so the old certification
  # must fail closed rather than silently certify the new candidate.
  expect_false(current$passed)
  expect_false(current$checks[["v8_supersession"]])
  expect_identical(
    current$protocol_v7_sha256,
    "c6d3a282b967bacef0f7da76a9192389e624e91075eccd3de4de648386f28487"
  )
  expect_identical(
    current$protocol_v8_sha256,
    "0e68dd7558d3df28878b8ac9f0edb49147620d9dd5c518898ed6243f305876c3"
  )
  expect_identical(
    current$protocol_v9_sha256,
    "4e6b828f463d9ae57de5a100f440dc7e96a376da3fdcae34509a42362223ac0a"
  )
  expect_false(dkfa_validate_current_protocol(
    certification_repo_root,
    candidate_source_tree_sha256 = "wrong-source"
  )$passed)

  bad_invalidation <- tempfile("dkfa-v4-invalidation-", fileext = ".json")
  payload <- jsonlite::read_json(paths$v4_inv, simplifyVector = FALSE)
  payload$status <- "mutated"
  jsonlite::write_json(payload, bad_invalidation, pretty = TRUE,
                       auto_unbox = TRUE)
  expect_false(dkfa_validate_v4_invalidation(
    bad_invalidation, paths$v4, certification_repo_root
  )$passed)

  for (version in c("v5", "v6")) {
    supersession <- paths[[paste0(version, "_sup")]]
    protocol <- paths[[version]]
    payload <- jsonlite::read_json(supersession, simplifyVector = FALSE)
    payload$materialization$source_test_suite_started <- TRUE
    candidate <- tempfile(paste0("dkfa-", version, "-supersession-"),
                          fileext = ".json")
    jsonlite::write_json(payload, candidate, pretty = TRUE, auto_unbox = TRUE)
    validator <- if (identical(version, "v5")) {
      dkfa_validate_v5_supersession
    } else {
      dkfa_validate_v6_supersession
    }
    expect_false(validator(candidate, protocol)$passed)
  }

  bad_v8 <- tempfile("dkfa-v8-protocol-", fileext = ".json")
  payload <- jsonlite::read_json(paths$v8, simplifyVector = FALSE)
  payload$court$n_perm <- 63L
  jsonlite::write_json(payload, bad_v8, pretty = TRUE, auto_unbox = TRUE)
  expect_false(dkfa_validate_v8_protocol_amendment(
    paths$v2, paths$v3, paths$v3_sup, paths$v4, paths$v4_inv,
    paths$v5, paths$v5_sup, paths$v6, paths$v6_sup,
    paths$v7, paths$v7_sup, bad_v8,
    certification_repo_root
  )$passed)

  bad_v9 <- tempfile("dkfa-v9-protocol-", fileext = ".json")
  payload <- jsonlite::read_json(paths$v9, simplifyVector = FALSE)
  payload$court$n_perm <- 63L
  jsonlite::write_json(payload, bad_v9, pretty = TRUE, auto_unbox = TRUE)
  expect_false(dkfa_validate_v9_protocol_amendment(
    paths$v2, paths$v3, paths$v3_sup, paths$v4, paths$v4_inv,
    paths$v5, paths$v5_sup, paths$v6, paths$v6_sup,
    paths$v7, paths$v7_sup, paths$v8, paths$v8_sup, bad_v9,
    certification_repo_root
  )$passed)
})

test_that("v7 supersession binds every retained execution byte", {
  skip_if_not(
    dir.exists(file.path(
      certification_repo_root, "data-raw",
      "functional-alignment-certification"
    )),
    "repository-only historical evidence is not shipped in the package tarball"
  )
  protocol_v7 <- file.path(certification_root, "protocol-v7.json")
  supersession <- file.path(
    certification_root, "protocol-v7-supersession.json"
  )
  expected_source <-
    "e38585975999d2bc6bdc1b76747e90a3ab89328940c92d676592494255831400"
  expect_true(dkfa_validate_v7_supersession(
    supersession, protocol_v7, root = certification_repo_root,
    candidate_source_tree_sha256 = expected_source
  )$passed)

  path_functions <- list(
    source_test = dkfa_v7_source_test_paths,
    court = dkfa_v7_court_paths,
    power = dkfa_v7_power_paths
  )
  for (role in names(path_functions)) {
    sandbox <- tempfile(paste0("dkfa-v7-", role, "-"))
    dir.create(sandbox)
    on.exit(unlink(sandbox, recursive = TRUE), add = TRUE)
    all_paths <- c(
      dkfa_v7_source_test_paths(certification_repo_root),
      dkfa_v7_court_paths(certification_repo_root),
      dkfa_v7_power_paths(certification_repo_root)
    )
    for (path in unname(all_paths)) {
      relative <- substring(path, nchar(certification_repo_root) + 2L)
      destination <- file.path(sandbox, relative)
      dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
      expect_true(file.copy(path, destination, overwrite = FALSE))
    }
    victim <- unname(path_functions[[role]](sandbox))[[1L]]
    writeLines(c(readLines(victim, warn = FALSE), "v7 mutation"), victim)
    expect_false(dkfa_validate_v7_supersession(
      supersession, protocol_v7, root = sandbox,
      candidate_source_tree_sha256 = expected_source
    )$passed, info = role)
  }
})

test_that("v8 supersession binds every retained execution and failure role", {
  skip_if_not(
    dir.exists(file.path(
      certification_repo_root, "data-raw",
      "functional-alignment-certification"
    )),
    "repository-only historical evidence is not shipped in the package tarball"
  )
  protocol_v8 <- file.path(certification_root, "protocol-v8.json")
  supersession <- file.path(
    certification_root, "protocol-v8-supersession.json"
  )
  expected_source <-
    "e38585975999d2bc6bdc1b76747e90a3ab89328940c92d676592494255831400"
  current <- dkfa_validate_v8_supersession(
    supersession, protocol_v8, root = certification_repo_root,
    candidate_source_tree_sha256 = expected_source
  )
  expect_false(current$passed)
  expect_false(current$checks[["documentation_only_repair"]])
  expect_true(all(current$checks[
    names(current$checks) != "documentation_only_repair"
  ]))

  path_functions <- list(
    source_test = dkfa_v8_source_test_paths,
    court = dkfa_v8_court_paths,
    power = dkfa_v8_power_paths
  )
  for (role in names(path_functions)) {
    sandbox <- tempfile(paste0("dkfa-v8-", role, "-"))
    dir.create(sandbox)
    on.exit(unlink(sandbox, recursive = TRUE), add = TRUE)
    all_paths <- c(
      dkfa_v8_source_test_paths(certification_repo_root),
      dkfa_v8_court_paths(certification_repo_root),
      dkfa_v8_power_paths(certification_repo_root)
    )
    for (path in unname(all_paths)) {
      relative <- substring(path, nchar(certification_repo_root) + 2L)
      destination <- file.path(sandbox, relative)
      dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
      expect_true(file.copy(path, destination, overwrite = FALSE))
    }
    expect_true(dkfa_validate_v8_supersession(
      supersession, protocol_v8, root = sandbox,
      candidate_source_tree_sha256 = expected_source
    )$passed, info = paste("clean retained evidence must pass:", role))
    victim <- unname(path_functions[[role]](sandbox))[[1L]]
    writeLines(c(readLines(victim, warn = FALSE), "v8 mutation"), victim)
    expect_false(dkfa_validate_v8_supersession(
      supersession, protocol_v8, root = sandbox,
      candidate_source_tree_sha256 = expected_source
    )$passed, info = role)
  }

  mutated_supersession <- tempfile(
    "dkfa-v8-package-failure-", fileext = ".json"
  )
  payload <- jsonlite::read_json(supersession, simplifyVector = FALSE)
  payload$package_documentation_gate_failure$failure_signature <- "mutated"
  jsonlite::write_json(
    payload, mutated_supersession, pretty = TRUE, auto_unbox = TRUE
  )
  expect_false(dkfa_validate_v8_supersession(
    mutated_supersession, protocol_v8,
    candidate_source_tree_sha256 = expected_source,
    verify_retained_evidence = FALSE
  )$passed)
})

test_that("v2 invalidation and the one-change Sinkhorn fixture are immutable", {
  protocol <- file.path(certification_root, "protocol.json")
  protocol_v2 <- file.path(certification_root, "protocol-v2.json")
  protocol_v3 <- file.path(certification_root, "protocol-v3.json")
  invalidation <- file.path(
    certification_root, "court-v2-invalidation.json"
  )
  fixture <- file.path(
    certification_root, "court-v2-failing-sinkhorn-fixture.json"
  )
  expect_true(all(file.exists(c(
    protocol, protocol_v2, invalidation, fixture
  ))))
  historical_paths <- dkfa_v2_historical_paths(certification_repo_root)
  live_historical <- all(file.exists(historical_paths))
  validated <- dkfa_validate_court_v2_invalidation(
    invalidation, fixture, protocol_v2,
    historical_paths = if (live_historical) historical_paths else NULL,
    verify_solver = TRUE, verify_origin = TRUE
  )
  expect_true(validated$passed, info = paste(
    names(validated$checks)[!validated$checks], collapse = ", "
  ))
  expect_true(validated$checks[["fixture_cost_reconstructed"]])
  expect_true(validated$checks[["fixture_source_mass_reconstructed"]])
  expect_true(validated$checks[["fixture_target_mass_reconstructed"]])

  parsed <- jsonlite::read_json(invalidation, simplifyVector = FALSE)
  mutations <- list(
    protocol_hash = function(x) {
      x$v2$protocol_sha256 <- paste0("0", x$v2$protocol_sha256)
      x
    },
    harness_hash = function(x) {
      x$v2$court_harness_bundle_sha256 <- "mutated"
      x
    },
    source_hash = function(x) {
      x$v2$source_tree_sha256 <- "mutated"
      x
    },
    test_hash = function(x) {
      x$v2$test_validation_tree_sha256 <- "mutated"
      x
    },
    historical_manifest_hash = function(x) {
      x$historical_evidence$source_test_manifest$sha256 <- "mutated"
      x
    },
    historical_results_path = function(x) {
      x$historical_evidence$source_test_results$path <- "mutated"
      x
    },
    failed_budget_hash = function(x) {
      x$historical_evidence$failed_court_formal_budget$sha256 <- "mutated"
      x
    },
    failure_identity = function(x) {
      x$representative_failure$seed <- 1L
      x
    },
    condition_class = function(x) {
      x$representative_failure$condition_class <- list("error")
      x
    },
    condition_message = function(x) {
      x$representative_failure$condition_message <- "mutated"
      x
    },
    iteration_cap = function(x) {
      x$representative_failure$max_iter <- 5001L
      x
    },
    marginal_tolerance = function(x) {
      x$representative_failure$tolerance <- 2e-4
      x
    },
    marginal_error = function(x) {
      x$representative_failure$marginal_error <- 0
      x
    },
    fixture_hash = function(x) {
      x$fixture$sha256 <- "mutated"
      x
    },
    converged_plan_hash = function(x) {
      x$numerical_amendment_evidence$max_iter_50000$plan_sha256 <- "mutated"
      x
    },
    converged_operator_hash = function(x) {
      x$numerical_amendment_evidence$max_iter_50000$operator_sha256 <-
        "mutated"
      x
    },
    materialization = function(x) {
      x$materialization$statistical_outputs_materialized <- TRUE
      x
    }
  )
  for (name in names(mutations)) {
    bad_path <- tempfile(paste0("dkfa-v2-invalidation-", name),
                         fileext = ".json")
    jsonlite::write_json(
      mutations[[name]](parsed), bad_path, pretty = TRUE, auto_unbox = TRUE
    )
    expect_false(dkfa_validate_court_v2_invalidation(
      bad_path, fixture, protocol_v2
    )$passed, info = paste("receipt mutation must fail:", name))
  }

  bad_fixture <- tempfile("dkfa-v2-fixture-", fileext = ".json")
  fixture_payload <- jsonlite::read_json(fixture, simplifyVector = FALSE)
  fixture_payload$cost[[1L]][[1L]] <-
    fixture_payload$cost[[1L]][[1L]] + 0.01
  jsonlite::write_json(
    fixture_payload, bad_fixture, pretty = TRUE, auto_unbox = TRUE
  )
  expect_false(dkfa_validate_court_v2_invalidation(
    invalidation, bad_fixture, protocol_v2
  )$passed)

  bad_protocol <- tempfile("dkfa-v2-protocol-", fileext = ".json")
  protocol_payload <- jsonlite::read_json(protocol_v2, simplifyVector = FALSE)
  protocol_payload$court$n_perm <- 63L
  jsonlite::write_json(
    protocol_payload, bad_protocol, pretty = TRUE, auto_unbox = TRUE
  )
  expect_false(dkfa_validate_v3_protocol_amendment(
    bad_protocol, protocol_v3
  )$passed)

  if (live_historical) {
    copied_root <- tempfile("dkfa-v2-historical-")
    dir.create(copied_root)
    copied <- stats::setNames(
      file.path(copied_root, basename(historical_paths)),
      names(historical_paths)
    )
    expect_true(all(file.copy(historical_paths, copied)))
    expect_true(dkfa_validate_court_v2_invalidation(
      invalidation, fixture, protocol_v2, historical_paths = copied
    )$passed)
    writeLines(
      c(readLines(copied[["failed_court_formal_budget"]]), "mutated"),
      copied[["failed_court_formal_budget"]]
    )
    expect_false(dkfa_validate_court_v2_invalidation(
      invalidation, fixture, protocol_v2, historical_paths = copied
    )$passed)
  }
})

test_that("court and power run claims fail closed and bind output bytes", {
  root <- tempfile("dkfa-court-run-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dependency <- file.path(root, "dependency.txt")
  writeLines("frozen", dependency)
  harness <- dkfa_bundle_manifest(dependency, root)
  empty_dir <- file.path(root, "preexisting-empty")
  dir.create(empty_dir)
  expect_error(
    dkfa_claim_run_directory(empty_dir, harness),
    "already claimed and cannot be rerun"
  )
  run_dir <- file.path(root, "court-fixed-id")
  claim <- dkfa_claim_run_directory(
    run_dir, harness, metadata = list(run_id = "court-fixed-id")
  )
  expect_true(all(file.exists(unlist(claim))))
  dkfa_update_run_state(
    claim$state_path, "failed",
    condition = list(classes = "error", message = "deliberate failure")
  )
  state <- jsonlite::read_json(claim$state_path, simplifyVector = TRUE)
  expect_identical(state$status, "failed")
  expect_error(
    dkfa_claim_run_directory(run_dir, harness),
    "already claimed and cannot be rerun"
  )
  expect_error(
    dkfa_update_run_state(claim$state_path, "formal_complete"),
    "terminal certification run state"
  )

  power_dir <- file.path(root, "power-fixed-id")
  power_claim <- dkfa_claim_run_directory(
    power_dir, harness, metadata = list(run_id = "power-fixed-id"),
    harness_filename = "known-warp-power-harness-manifest.json",
    state_filename = "known-warp-power-run-state.json"
  )
  dkfa_update_run_state(
    power_claim$state_path, "failed",
    condition = list(classes = "error", message = "power failure")
  )
  expect_error(
    dkfa_claim_run_directory(
      power_dir, harness,
      harness_filename = "known-warp-power-harness-manifest.json",
      state_filename = "known-warp-power-run-state.json"
    ),
    "already claimed and cannot be rerun"
  )

  complete_dir <- file.path(root, "power-complete-id")
  complete_claim <- dkfa_claim_run_directory(
    complete_dir, harness, metadata = list(run_id = "power-complete-id"),
    harness_filename = "known-warp-power-harness-manifest.json",
    state_filename = "known-warp-power-run-state.json"
  )
  raw_path <- file.path(complete_dir, "raw.csv")
  writeLines("value\n1", raw_path)
  dkfa_update_run_state(
    complete_claim$state_path, "power_complete",
    outputs = list(raw = dkfa_hash_file(raw_path))
  )
  intact <- dkfa_validate_run_state_outputs(
    complete_claim$state_path, "power_complete",
    list(power_complete = c(raw = raw_path))
  )
  expect_true(intact$passed)
  writeLines(c(readLines(raw_path), "2"), raw_path)
  mutated <- dkfa_validate_run_state_outputs(
    complete_claim$state_path, "power_complete",
    list(power_complete = c(raw = raw_path))
  )
  expect_false(mutated$passed)
  expect_false(unname(mutated$checks[["power_complete_hashes_match"]]))
  expect_error(
    dkfa_with_run_failure_receipt(
      complete_claim$state_path, "power_complete", {
        verified <- dkfa_validate_run_state_outputs(
          complete_claim$state_path, "power_complete",
          list(power_complete = c(raw = raw_path))
        )
        if (!isTRUE(verified$passed)) {
          stop("deliberate verifier failure", call. = FALSE)
        }
      }
    ),
    "deliberate verifier failure"
  )
  failed_verifier <- jsonlite::read_json(
    complete_claim$state_path, simplifyVector = TRUE
  )
  expect_identical(failed_verifier$status, "failed")
  expect_error(
    dkfa_claim_run_directory(
      complete_dir, harness,
      harness_filename = "known-warp-power-harness-manifest.json",
      state_filename = "known-warp-power-run-state.json"
    ),
    "already claimed and cannot be rerun"
  )
})

test_that("court and known-warp harnesses bind the complete history", {
  skip_if_not(
    dir.exists(file.path(
      certification_repo_root, "data-raw",
      "functional-alignment-certification"
    )),
    "repository-only historical evidence is not shipped in the package tarball"
  )
  court_paths <- dkfa_court_harness_paths(certification_repo_root)
  power_paths <- dkfa_power_harness_paths(certification_repo_root)
  historical <- unname(dkfa_v2_historical_paths(certification_repo_root))
  v3_history <- unname(dkfa_v3_history_paths(certification_repo_root))
  v4_history <- unname(dkfa_v4_history_paths(certification_repo_root))
  v5_history <- unname(dkfa_v5_history_paths(certification_repo_root))
  v6_history <- unname(dkfa_v6_history_paths(certification_repo_root))
  v7_history <- unname(dkfa_v7_history_paths(certification_repo_root))
  v8_history <- unname(dkfa_v8_history_paths(certification_repo_root))
  expect_true(all(historical %in% court_paths))
  expect_true(all(historical %in% power_paths))
  expect_true(all(v3_history %in% court_paths))
  expect_true(all(v3_history %in% power_paths))
  expect_true(all(v4_history %in% court_paths))
  expect_true(all(v4_history %in% power_paths))
  expect_true(all(v5_history %in% court_paths))
  expect_true(all(v5_history %in% power_paths))
  expect_true(all(v6_history %in% court_paths))
  expect_true(all(v6_history %in% power_paths))
  expect_true(all(v7_history %in% court_paths))
  expect_true(all(v7_history %in% power_paths))
  expect_true(all(v8_history %in% court_paths))
  expect_true(all(v8_history %in% power_paths))
  expect_identical(
    sum(basename(court_paths) == "evidence-utils.R"), 1L
  )
  expect_identical(
    sum(basename(power_paths) == "evidence-utils.R"), 1L
  )
})

test_that("court evidence requires the exact frozen schedule and raw summary", {
  protocol <- jsonlite::read_json(
    file.path(certification_root, "protocol.json"), simplifyVector = TRUE
  )
  raw <- dkfa_court_expected_schedule(protocol)
  raw$reject_fwer <- FALSE
  raw$reject_unadjusted <- FALSE
  raw$effect_bias <- 0
  raw$effect_abs_bias <- 0
  raw$operator_sensitivity <- 0
  raw$convergence <- 1
  raw$relative_eigengap <- 1
  raw$latent_span_predictor <- 0
  summary <- fa_court_summarize(raw)
  valid <- dkfa_validate_court_evidence(raw, summary, protocol)
  expect_true(valid$passed, info = paste(
    names(valid$checks)[!valid$checks], collapse = ", "
  ))

  duplicate_key <- raw
  duplicate_key[nrow(duplicate_key), c("cell", "seed", "arm")] <-
    duplicate_key[1L, c("cell", "seed", "arm")]
  expect_false(dkfa_validate_court_evidence(
    duplicate_key, fa_court_summarize(duplicate_key), protocol
  )$passed)

  relabeled <- raw
  relabeled$arm[[1L]] <- "geometry_only"
  expect_false(dkfa_validate_court_evidence(
    relabeled, fa_court_summarize(relabeled), protocol
  )$passed)

  favorable_summary <- summary
  favorable_summary$fwer_rate[[1L]] <- 0.5
  expect_false(dkfa_validate_court_evidence(
    raw, favorable_summary, protocol
  )$passed)
})

test_that("known-warp evidence requires every frozen seed-arm pair", {
  protocol <- jsonlite::read_json(
    file.path(certification_root, "protocol.json"), simplifyVector = TRUE
  )
  raw <- dkfa_known_warp_expected_schedule(protocol)
  raw$reject <- FALSE
  raw$correlation <- 0
  raw$rmse <- 1
  raw$amplitude_ratio <- 1
  raw$positive_peak_error <- 1
  raw$negative_peak_error <- 1
  raw$latent_error <- 1
  raw$point_spread <- 1
  raw$template_converged <- TRUE
  summary <- dktf_power_summarize(raw)
  expect_true(dkfa_validate_known_warp_evidence(
    raw, summary, protocol
  )$passed)

  duplicate_key <- raw
  duplicate_key[nrow(duplicate_key), c("seed", "arm")] <-
    duplicate_key[1L, c("seed", "arm")]
  expect_false(dkfa_validate_known_warp_evidence(
    duplicate_key, dktf_power_summarize(duplicate_key), protocol
  )$passed)

  favorable_summary <- summary
  favorable_summary$power[[1L]] <- 1
  expect_false(dkfa_validate_known_warp_evidence(
    raw, favorable_summary, protocol
  )$passed)
})

test_that("known-warp maxT signs are exhaustive modulo global sign", {
  skip_if_not(file.exists(power_source))
  signs <- dktf_power_signs(8L)
  expect_equal(dim(signs), c(8L, 127L))
  expect_true(all(signs[1L, ] == 1))
  expect_equal(ncol(unique(t(signs))), 8L)
  # The 127 sign columns themselves must be unique.
  expect_equal(nrow(unique(t(signs))), 127L)
  expect_false(any(colSums(signs == 1) == nrow(signs)))
  expect_true(any(apply(signs, 2L, identical, c(1, rep(-1, 7L)))))
  complete <- cbind(identity = rep(1, 8L), signs)
  expect_equal(nrow(unique(t(complete))), 2^(8L - 1L))

  set.seed(918)
  Y <- matrix(rnorm(8L * 6L), 8L, 6L)
  first <- dktf_power_test(Y, signs)
  second <- dktf_power_test(Y, signs)
  expect_identical(first, second)
  expect_true(first[["p_value"]] >= 0 && first[["p_value"]] <= 1)
  expect_identical(
    as.logical(first[["reject"]]), first[["p_value"]] <= 0.05
  )
})

test_that("certification runners are source-addressed and refuse overwrite", {
  court_runner <- readLines(
    file.path(certification_root, "run-court.R"), warn = FALSE
  )
  power_runner <- readLines(
    file.path(certification_root, "run-known-warp-power.R"), warn = FALSE
  )
  pkgdown_runner <- readLines(
    file.path(certification_root, "run-pkgdown.R"), warn = FALSE
  )
  source_test_runner <- readLines(
    file.path(certification_root, "run-source-tests.R"), warn = FALSE
  )
  collector <- readLines(
    file.path(certification_root, "collect-certification.R"), warn = FALSE
  )
  provenance <- readLines(
    file.path(certification_root, "package-provenance.R"), warn = FALSE
  )
  expect_true(any(grepl("source_tree_hash", court_runner, fixed = TRUE)))
  expect_true(any(grepl("harness", court_runner, fixed = TRUE)))
  expect_true(any(grepl("dkfa_claim_run_directory", court_runner,
                        fixed = TRUE)))
  expect_true(any(grepl("dkfa_update_run_state", court_runner,
                        fixed = TRUE)))
  expect_true(any(grepl("source_tree_hash", power_runner, fixed = TRUE)))
  expect_true(any(grepl("harness", power_runner, fixed = TRUE)))
  expect_true(any(grepl("dkfa_claim_run_directory", power_runner,
                        fixed = TRUE)))
  expect_true(any(grepl("known-warp-power-run-state.json", power_runner,
                        fixed = TRUE)))
  expect_true(any(grepl("dkfa_update_run_state", power_runner,
                        fixed = TRUE)))
  expect_true(any(grepl("tarball_source_tree_sha256", pkgdown_runner,
                        fixed = TRUE)))
  expect_true(any(grepl("l10n_info", pkgdown_runner, fixed = TRUE)))
  expect_true(any(grepl("utf8_locale", pkgdown_runner, fixed = TRUE)))
  expect_true(any(grepl("dkfa_test_validation_manifest", source_test_runner,
                        fixed = TRUE)))
  expect_true(any(grepl("load_package = \"none\"", source_test_runner,
                        fixed = TRUE)))
  expect_true(any(grepl("protocol_history_bundle_sha256", source_test_runner,
                        fixed = TRUE)))
  expect_true(any(grepl("pkgdown_site_content_sha256", collector,
                        fixed = TRUE)))
  expect_true(any(grepl("checked_package_tree_sha256", collector,
                        fixed = TRUE)))
  expect_true(any(grepl("reviewed_evidence_bundle_sha256", collector,
                        fixed = TRUE)))
  expect_true(any(grepl("--ignore-vignettes", collector, fixed = TRUE)))
  expect_true(any(grepl("no_unexpected_options", provenance, fixed = TRUE)))
  expect_true(any(grepl("built_tarball_check_policy_bound", collector,
                        fixed = TRUE)))
  expect_true(any(grepl("pkgdown_utf8_locale_recorded", collector,
                        fixed = TRUE)))
  expect_true(any(grepl("v7_superseded_history", collector, fixed = TRUE)))
  expect_true(any(grepl("v8_superseded_history", collector, fixed = TRUE)))
  expect_true(any(grepl("staging_dir", collector, fixed = TRUE)))
})

test_that("certified R CMD check policy cannot skip substantive checks", {
  good <- c(
    "* using options '--no-manual --ignore-vignettes'",
    "Status: OK"
  )
  expect_true(dkfa_validate_r_cmd_check_log(good)$passed)
  expect_true(dkfa_validate_r_cmd_check_log(c(
    "* using options '--no-manual --ignore-vignettes --no-clean'",
    "Status: OK"
  ))$passed)
  expect_true(dkfa_validate_r_cmd_check_log(c(
    "* using options ‘--no-manual --ignore-vignettes’",
    "Status: OK"
  ))$passed)
  expect_false(dkfa_validate_r_cmd_check_log(
    sub("--ignore-vignettes", "--no-vignettes", good, fixed = TRUE)
  )$passed)
  for (option in c(
    "--no-install", "--install=skip", "--no-tests",
    "--test-dir=/private/tmp/empty", "--no-examples", "--no-codoc",
    "--check-subdirs=no"
  )) {
    weakened <- good
    weakened[[1L]] <- sub(
      "--ignore-vignettes",
      paste("--ignore-vignettes", option),
      weakened[[1L]], fixed = TRUE
    )
    expect_false(
      dkfa_validate_r_cmd_check_log(weakened)$passed,
      info = option
    )
  }
  expect_false(dkfa_validate_r_cmd_check_log(
    sub("Status: OK", "Status: 1 WARNING", good, fixed = TRUE)
  )$passed)
})

test_that("all albersdown vignette version guards use an exported comparator", {
  vignette_dir <- file.path(certification_repo_root, "vignettes")
  skip_if_not(
    dir.exists(vignette_dir),
    paste0(
      "source vignette guards are verified by the source suite and exact ",
      "pkgdown gate"
    )
  )
  paths <- list.files(
    vignette_dir, pattern = "\\.Rmd$", full.names = TRUE
  )
  guarded <- Filter(function(path) {
    any(grepl(
      "^has_albers_v2 <-", readLines(path, warn = FALSE)
    ))
  }, paths)
  expect_length(guarded, 19L)

  for (path in guarded) {
    lines <- readLines(path, warn = FALSE)
    expect_false(
      any(grepl("utils::package_version", lines, fixed = TRUE)),
      info = basename(path)
    )
    start <- grep("^has_albers_v2 <-", lines)
    expression <- paste(lines[start:(start + 1L)], collapse = "\n")
    expect_silent(
      eval(parse(text = expression), envir = new.env(parent = baseenv()))
    )
  }
})

test_that("dependency bundle receipts fail after any harness mutation", {
  skip_if_not(file.exists(evidence_utils_source))
  root <- tempfile("dkfa-bundle-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  first <- file.path(root, "runner.R")
  second <- file.path(root, "helper.R")
  writeLines("runner <- TRUE", first)
  writeLines("helper <- TRUE", second)

  receipt <- dkfa_bundle_manifest(c(first, second), root)
  expect_silent(dkfa_verify_bundle_manifest(
    receipt, c(second, first), root, "Adversarial"
  ))
  writeLines("helper <- FALSE", second)
  expect_error(
    dkfa_verify_bundle_manifest(receipt, c(first, second), root,
                                "Adversarial"),
    "does not match"
  )
})

test_that("review evidence digest binds every consumed derivative", {
  skip_if_not(file.exists(evidence_utils_source))
  root <- tempfile("dkfa-reviewed-evidence-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  roles <- c(
    "protocol", "protocol_v2", "court_v2_invalidation",
    "court_v2_sinkhorn_fixture", "v2_source_test_manifest",
    "v2_source_test_results", "v2_source_test_session_info",
    "v2_failed_court_formal_budget", "court_run_state", "court_verdict",
    "court_verdict_markdown", "court_raw",
    "court_summary", "court_audit", "court_focus",
    "court_factor_uncertainty", "court_formal_budget", "court_formal_grid",
    "court_inflation_flags", "court_latent_span_excess",
    "court_latent_span_logistic", "court_latent_span_model",
    "court_formal_manifest", "court_harness_manifest",
    "known_warp_manifest", "known_warp_raw", "known_warp_summary",
    "known_warp_determinism", "source_test_manifest",
    "source_test_results", "source_test_session_info", "check_log",
    "pkgdown_receipt"
  )
  paths <- stats::setNames(file.path(root, paste0(roles, ".txt")), roles)
  for (role in roles) writeLines(paste("original", role), paths[[role]])
  baseline <- dkfa_review_evidence_bundle(paths)
  expect_setequal(names(baseline$hashes), roles)
  for (role in roles) {
    writeLines(paste("mutated", role), paths[[role]])
    changed <- dkfa_review_evidence_bundle(paths)
    expect_false(
      identical(changed$bundle_sha256, baseline$bundle_sha256),
      info = paste("mutation must invalidate review evidence role", role)
    )
    writeLines(paste("original", role), paths[[role]])
  }
  unlink(paths[["court_audit"]])
  expect_error(
    dkfa_review_evidence_bundle(paths),
    "Missing reviewed evidence"
  )
})

test_that("checked-source and pkgdown receipts are mutation-sensitive", {
  skip_if_not(file.exists(evidence_utils_source))
  root <- tempfile("dkfa-artifacts-")
  tar_tree <- file.path(root, "tarball")
  check_tree <- file.path(root, "checked")
  site <- file.path(root, "site")
  dir.create(tar_tree, recursive = TRUE)
  dir.create(check_tree, recursive = TRUE)
  dir.create(file.path(site, "reference"), recursive = TRUE)
  writeLines("same source", file.path(tar_tree, "R-code.R"))
  writeLines("same source", file.path(check_tree, "R-code.R"))
  expect_identical(
    dkfa_tree_manifest(tar_tree)$tree_sha256,
    dkfa_tree_manifest(check_tree)$tree_sha256
  )
  writeLines("stale checked source", file.path(check_tree, "R-code.R"))
  expect_false(identical(
    dkfa_tree_manifest(tar_tree)$tree_sha256,
    dkfa_tree_manifest(check_tree)$tree_sha256
  ))

  required <- c("index.html", "reference/index.html")
  writeLines("home", file.path(site, required[[1L]]))
  writeLines("reference", file.path(site, required[[2L]]))
  site_hash <- dkfa_site_manifest(site)$tree_sha256
  jsonlite::write_json(
    list(
      schema_version = "1.0.0",
      tarball_sha256 = "tarball-hash",
      tarball_source_tree_sha256 = "source-hash",
      docs_input_sha256 = "docs-hash",
      site_content_sha256 = site_hash,
      required_files = as.list(required)
    ),
    file.path(site, ".dkge-pkgdown-receipt.json"),
    auto_unbox = TRUE
  )
  valid <- dkfa_validate_pkgdown_receipt(
    site, "tarball-hash", "source-hash", "docs-hash", required
  )
  expect_true(valid$passed)
  writeLines("mutated after build", file.path(site, required[[1L]]))
  invalid <- dkfa_validate_pkgdown_receipt(
    site, "tarball-hash", "source-hash", "docs-hash", required
  )
  expect_false(invalid$passed)
  expect_false(unname(invalid$checks[["site"]]))
})

test_that("certified pkgdown receipts require UTF-8 locale provenance", {
  skip_if_not(file.exists(package_provenance_source))
  site <- tempfile("dkfa-certified-pkgdown-")
  dir.create(file.path(site, "reference"), recursive = TRUE)
  on.exit(unlink(site, recursive = TRUE), add = TRUE)
  required <- c("index.html", "reference/index.html")
  writeLines("home", file.path(site, required[[1L]]))
  writeLines("reference", file.path(site, required[[2L]]))

  source_provenance <- list(
    certified_source_tree_sha256 = "source-hash",
    built_raw_source_tree_sha256 = "raw-source-hash",
    package_projection_sha256 = "projection-hash",
    shipped_payload_tree_sha256 = "payload-hash",
    description = list(semantic_sha256 = "description-hash")
  )
  docs_provenance <- list(
    certified_docs_input_sha256 = "docs-hash",
    overlay_bundle_sha256 = "overlay-hash"
  )
  receipt <- list(
    schema_version = "1.4.0",
    utf8_locale = TRUE,
    locale = "LC_CTYPE=en_US.UTF-8;LC_COLLATE=en_US.UTF-8",
    tarball_sha256 = "tarball-hash",
    tarball_source_tree_sha256 = "source-hash",
    tarball_raw_source_tree_sha256 = "raw-source-hash",
    package_source_projection_sha256 = "projection-hash",
    tarball_shipped_payload_tree_sha256 = "payload-hash",
    description_semantic_sha256 = "description-hash",
    docs_input_sha256 = "docs-hash",
    docs_overlay_bundle_sha256 = "overlay-hash",
    site_content_sha256 = dkfa_site_manifest(site)$tree_sha256,
    required_files = as.list(required)
  )
  receipt_path <- file.path(site, ".dkge-pkgdown-receipt.json")
  write_receipt <- function(value) {
    jsonlite::write_json(value, receipt_path, auto_unbox = TRUE)
  }
  validate <- function() {
    dkfa_validate_certified_pkgdown_receipt(
      site,
      tarball_sha256 = "tarball-hash",
      source_provenance = source_provenance,
      docs_provenance = docs_provenance,
      required_files = required
    )
  }

  write_receipt(receipt)
  expect_true(validate()$passed)

  false_locale <- receipt
  false_locale$utf8_locale <- FALSE
  write_receipt(false_locale)
  invalid_false <- validate()
  expect_false(invalid_false$passed)
  expect_false(unname(invalid_false$checks[["utf8_locale"]]))

  missing_locale_gate <- receipt
  missing_locale_gate$utf8_locale <- NULL
  write_receipt(missing_locale_gate)
  invalid_missing <- validate()
  expect_false(invalid_missing$passed)
  expect_false(unname(invalid_missing$checks[["utf8_locale"]]))

  missing_locale_string <- receipt
  missing_locale_string$locale <- NULL
  write_receipt(missing_locale_string)
  invalid_string <- validate()
  expect_false(invalid_string$passed)
  expect_false(unname(invalid_string$checks[["locale_recorded"]]))
})

test_that("R build projection and docs overlay are explicit and fail closed", {
  skip_if_not(file.exists(package_provenance_source))
  root <- tempfile("dkfa-package-projection-")
  live <- file.path(root, "live")
  built <- file.path(root, "built")
  for (path in c(
    file.path(live, c(
      "R", "src", "man", "vignettes", "pkgdown", "inst/validation",
      "inst/extdata", "tests/testthat"
    )),
    file.path(built, c(
      "R", "src", "man", "vignettes", "inst/validation",
      "inst/extdata", "tests/testthat"
    ))
  )) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  live_description <- c(
    "Package: projectionfixture",
    "Version: 1.0.0",
    "Title: Projection Fixture",
    "Description: A fixture for package projection tests.",
    "Authors@R: person(\"Test\", \"Author\", role = c(\"aut\", \"cre\"),",
    "    email = \"test@example.org\")",
    "License: MIT"
  )
  built_description <- c(
    "Package: projectionfixture",
    "Version: 1.0.0",
    "Title: Projection Fixture",
    "Description: A fixture for package projection tests.",
    paste0("Authors@R: person(\"Test\", \"Author\", role = ",
           "c(\"aut\", \"cre\"), email = \"test@example.org\")"),
    "License: MIT",
    "NeedsCompilation: yes",
    "Packaged: 2026-08-27 00:00:00 UTC; test",
    "Author: Test Author [aut, cre]",
    "Maintainer: Test Author <test@example.org>"
  )
  writeLines(live_description, file.path(live, "DESCRIPTION"))
  writeLines(built_description, file.path(built, "DESCRIPTION"))
  for (relative in c(
    "NAMESPACE", "R/code.R", "src/code.cpp", "man/topic.Rd",
    "vignettes/article.Rmd", "README.md", "inst/validation/helper.R",
    "inst/extdata/calibration.csv", "tests/testthat/test-projection.R",
    "LICENSE", "NEWS.md"
  )) {
    writeLines(paste("content", relative), file.path(live, relative))
    writeLines(paste("content", relative), file.path(built, relative))
  }
  # Repository guidance is deliberately excluded by .Rbuildignore and is not
  # part of the closed-world built-package projection.
  writeLines("repository guidance", file.path(live, "CONTRIBUTING.md"))
  writeLines("url: https://example.org", file.path(live, "_pkgdown.yml"))
  writeLines("body {}", file.path(live, "pkgdown", "extra.css"))

  live_source_hash <- dkfa_certified_source_tree_hash(live)
  source_projection <- dkfa_verify_built_source(
    live, built, certified_source_tree_sha256 = live_source_hash
  )
  expect_true(
    source_projection$passed,
    info = paste0(
      paste(
        names(source_projection$checks)[!source_projection$checks],
        collapse = ", "
      ),
      "; live-only: ",
      paste(source_projection$live_only_package_payload_files,
            collapse = ", "),
      "; built-only: ",
      paste(source_projection$built_only_package_payload_files,
            collapse = ", ")
    )
  )
  expect_identical(
    source_projection$certified_source_tree_sha256, live_source_hash
  )
  wrong_source <- dkfa_verify_built_source(
    live, built, certified_source_tree_sha256 = "wrong-certified-source"
  )
  expect_false(wrong_source$passed)
  expect_false(unname(wrong_source$checks[["certified_source_tree"]]))
  expect_identical(
    wrong_source$certified_source_tree_sha256, live_source_hash
  )
  validation_path <- file.path(built, "inst", "validation", "helper.R")
  writeLines("mutated shipped validation", validation_path)
  invalid_validation <- dkfa_verify_built_source(
    live, built, certified_source_tree_sha256 = live_source_hash
  )
  expect_false(invalid_validation$passed)
  expect_false(unname(
    invalid_validation$checks[["shipped_payload_file_bytes"]]
  ))
  writeLines("content inst/validation/helper.R", validation_path)

  test_path <- file.path(built, "tests", "testthat", "test-projection.R")
  writeLines("mutated shipped tests", test_path)
  invalid_tests <- dkfa_verify_built_source(
    live, built, certified_source_tree_sha256 = live_source_hash
  )
  expect_false(invalid_tests$passed)
  expect_false(unname(
    invalid_tests$checks[["shipped_payload_file_bytes"]]
  ))
  writeLines("content tests/testthat/test-projection.R", test_path)

  extdata_path <- file.path(built, "inst", "extdata", "calibration.csv")
  writeLines("mutated shipped extdata", extdata_path)
  invalid_extdata <- dkfa_verify_built_source(
    live, built, certified_source_tree_sha256 = live_source_hash
  )
  expect_false(invalid_extdata$passed)
  expect_false(unname(
    invalid_extdata$checks[["shipped_payload_file_bytes"]]
  ))
  writeLines("content inst/extdata/calibration.csv", extdata_path)

  configure_path <- file.path(built, "configure")
  writeLines("#!/bin/sh", configure_path)
  invalid_hook <- dkfa_verify_built_source(
    live, built, certified_source_tree_sha256 = live_source_hash
  )
  expect_false(invalid_hook$passed)
  expect_false(unname(
    invalid_hook$checks[["closed_world_package_file_names"]]
  ))
  unlink(configure_path)

  live_extdata_path <- file.path(
    live, "inst", "extdata", "calibration.csv"
  )
  validation_before <- dkfa_test_validation_manifest(live)$tree_sha256
  writeLines("mutated live extdata", live_extdata_path)
  validation_after <- dkfa_test_validation_manifest(live)$tree_sha256
  expect_false(identical(validation_before, validation_after))
  writeLines("content inst/extdata/calibration.csv", live_extdata_path)
  old_collate <- Sys.getlocale("LC_COLLATE")
  on.exit(suppressWarnings(Sys.setlocale("LC_COLLATE", old_collate)),
          add = TRUE)
  locales <- c("C", "POSIX", "en_US.UTF-8", "nl_NL.UTF-8", "fr_FR.UTF-8")
  locale_hashes <- vapply(locales, function(locale) {
    selected <- suppressWarnings(Sys.setlocale("LC_COLLATE", locale))
    if (is.na(selected)) return(NA_character_)
    dkfa_certified_source_tree_hash(live)
  }, character(1))
  expect_length(unique(stats::na.omit(locale_hashes)), 1L)
  court_locale_hashes <- vapply(locales, function(locale) {
    selected <- suppressWarnings(Sys.setlocale("LC_COLLATE", locale))
    if (is.na(selected)) return(NA_character_)
    fa_court_source_tree_hash(live)
  }, character(1))
  expect_length(unique(stats::na.omit(court_locale_hashes)), 1L)
  expect_false(identical(
    digest::digest(file = file.path(live, "DESCRIPTION"), algo = "sha256",
                   serialize = FALSE),
    digest::digest(file = file.path(built, "DESCRIPTION"), algo = "sha256",
                   serialize = FALSE)
  ))
  docs_projection <- dkfa_docs_projection(live, built)
  expect_true(docs_projection$passed)
  overlaid <- dkfa_apply_docs_overlay(live, built)
  expect_identical(
    overlaid$overlaid_docs_input_sha256,
    dkfa_pkgdown_input_manifest(live)$tree_sha256
  )

  writeLines("mutated source", file.path(built, "R", "code.R"))
  expect_false(dkfa_verify_built_source(
    live, built, certified_source_tree_sha256 = live_source_hash
  )$passed)
})

test_that("source-test evidence binds results, counts, and session provenance", {
  skip_if_not(file.exists(package_provenance_source))
  root <- tempfile("dkfa-source-tests-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  results_path <- file.path(root, "source-test-results.csv")
  manifest_path <- file.path(root, "source-test-manifest.json")
  session_path <- file.path(root, "session-info.txt")
  results <- data.frame(
    file = c("a.R", "b.R"),
    nb = c(3L, 2L),
    failed = c(0L, 0L),
    skipped = c(FALSE, TRUE),
    error = c(FALSE, FALSE),
    warning = c(0L, 1L)
  )
  utils::write.csv(results, results_path, row.names = FALSE)
  writeLines("R session fixture", session_path)
  valid_results <- results
  write_manifest <- function(expectations = 5L, test_blocks = 2L,
                             warnings = 1L, skipped = 1L,
                             gate = "full_source_test_suite",
                             source_compile = TRUE) {
    jsonlite::write_json(
      list(
        schema_version = "1.2.0",
        gate = gate,
        passed = TRUE,
        source_loader = list(
          method = "pkgload::load_all", compile = source_compile
        ),
        test_blocks = test_blocks,
        expectations = expectations,
        expectation_failures = 0L,
        errors = 0L,
        failed = 0L,
        warnings = warnings,
        skipped = skipped,
        results_sha256 = dkfa_hash_file(results_path),
        session_info_sha256 = dkfa_hash_file(session_path)
      ),
      manifest_path,
      auto_unbox = TRUE
    )
  }
  write_manifest()
  evidence <- dkfa_validate_source_test_evidence(
    manifest_path, results_path, session_path
  )
  expect_true(evidence$passed)
  expect_equal(unname(evidence$counts[["expectations"]]), 5)

  write_manifest(source_compile = FALSE)
  expect_false(dkfa_validate_source_test_evidence(
    manifest_path, results_path, session_path
  )$passed)

  write_manifest(expectations = 6L)
  expect_false(dkfa_validate_source_test_evidence(
    manifest_path, results_path, session_path
  )$passed)

  write_manifest()
  results$nb[[1L]] <- 4L
  utils::write.csv(results, results_path, row.names = FALSE)
  expect_false(dkfa_validate_source_test_evidence(
    manifest_path, results_path, session_path
  )$passed)

  utils::write.csv(valid_results, results_path, row.names = FALSE)
  write_manifest(gate = "not_the_full_source_suite")
  expect_false(dkfa_validate_source_test_evidence(
    manifest_path, results_path, session_path
  )$passed)

  empty_results <- valid_results[0L, , drop = FALSE]
  utils::write.csv(empty_results, results_path, row.names = FALSE)
  write_manifest(
    expectations = 0L, test_blocks = 0L, warnings = 0L, skipped = 0L
  )
  empty_evidence <- dkfa_validate_source_test_evidence(
    manifest_path, results_path, session_path
  )
  expect_false(empty_evidence$passed)
  expect_false(unname(empty_evidence$checks[["nonempty_test_blocks"]]))
  expect_false(unname(empty_evidence$checks[["nonempty_expectations"]]))

  unlink(session_path)
  expect_error(
    dkfa_validate_source_test_evidence(
      manifest_path, results_path, session_path
    ),
    "Missing source-test evidence"
  )
})
