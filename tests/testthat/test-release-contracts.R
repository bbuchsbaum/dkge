test_that("core dependency metadata is immutable and beta-resolvable", {
  desc_path <- testthat::test_path("..", "..", "DESCRIPTION")
  desc <- if (file.exists(desc_path)) {
    as.list(read.dcf(desc_path)[1, ])
  } else {
    as.list(utils::packageDescription("dkge"))
  }
  fields <- names(desc)
  suggests <- trimws(strsplit(desc[["Suggests"]], ",", fixed = TRUE)[[1]])
  enhances <- trimws(strsplit(desc[["Enhances"]], ",", fixed = TRUE)[[1]])

  expect_false("Remotes" %in% fields,
               info = "core checks must not resolve mutable GitHub remotes")
  expect_false("T4transport" %in% suggests,
               info = "the experimental oracle belongs in its opt-in workflow")
  expect_true(all(c("future", "withr") %in% suggests))
  expect_true("neuralign" %in% enhances)
  expect_match(desc[["Additional_repositories"]], "bbuchsbaum[.]r-universe[.]dev")

  lock_path <- testthat::test_path("..", "..", "tools", "release",
                                   "noncran-lock.csv")
  if (file.exists(lock_path)) {
    lock <- utils::read.csv(lock_path, stringsAsFactors = FALSE)
    expect_setequal(lock$package, c("fmridesign", "fmrireg", "neuralign"))
    expect_true(all(grepl("^[0-9a-f]{40}$", lock$sha)))
    expect_identical(lock$role, c("Import", "Import", "Enhances"))
  }
})

test_that("release workflows cover checks sanitizers and the independent oracle", {
  root <- testthat::test_path("..", "..")
  workflows <- file.path(root, ".github", "workflows")
  skip_if_not(dir.exists(workflows), "source-only workflow contract")
  expect_true(file.exists(file.path(workflows, "R-CMD-check.yaml")))
  expect_true(file.exists(file.path(workflows, "sanitizers.yaml")))
  expect_true(file.exists(file.path(workflows, "sinkhorn-oracle.yaml")))

  check_yaml <- paste(readLines(file.path(workflows, "R-CMD-check.yaml"),
                               warn = FALSE), collapse = "\n")
  expect_match(check_yaml, "macos-latest")
  expect_match(check_yaml, "windows-latest")
  expect_match(check_yaml, "ubuntu-latest")
  expect_match(check_yaml, "oldrel")
  expect_match(check_yaml, "devel")
  expect_match(check_yaml, "build-source")
  expect_match(check_yaml, "actions/download-artifact@v4")
  expect_match(check_yaml, "Verify exact candidate identity")
  expect_match(check_yaml, "rcmdcheck::rcmdcheck")
  expect_match(check_yaml, "source-tarball[.]sha256")
  expect_match(check_yaml, "e627d4d1861e5cc4aea5a5668d6052e394429b92")
  expect_match(check_yaml, "20c9f07a225f21c4eea1732ea31418ccedd1c056")

  sanitizer_yaml <- paste(readLines(file.path(workflows, "sanitizers.yaml"),
                                    warn = FALSE), collapse = "\n")
  expect_match(sanitizer_yaml, "r-devel-san")
  expect_match(sanitizer_yaml, "ASAN_OPTIONS")
  expect_match(sanitizer_yaml, "UBSAN_OPTIONS")

  oracle_yaml <- paste(readLines(file.path(workflows, "sinkhorn-oracle.yaml"),
                                 warn = FALSE), collapse = "\n")
  expect_false(grepl("T4transport", oracle_yaml, fixed = TRUE))
  oracle_test <- paste(readLines(file.path(root, "tests", "testthat",
                                          "test-transport-sinkhorn.R"),
                                warn = FALSE), collapse = "\n")
  expect_match(oracle_test, "dense_sinkhorn_reference")

  integration_yaml <- paste(readLines(
    file.path(workflows, "neuralign-integration.yaml"), warn = FALSE
  ), collapse = "\n")
  expect_match(integration_yaml, "ae1468340497ef89a32b148946cebdf6dfa42c47")
})
