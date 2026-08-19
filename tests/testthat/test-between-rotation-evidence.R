library(testthat)

.audit0817_file_sha256 <- function(path) {
  expect_true(file.exists(path), info = path)
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(file = path, algo = "sha256", serialize = FALSE))
  }
  hashed <- suppressWarnings(system2(
    "shasum", c("-a", "256", path), stdout = TRUE, stderr = FALSE
  ))
  if (is.null(attr(hashed, "status")) || identical(attr(hashed, "status"), 0L)) {
    return(sub("\\s.*$", "", hashed[[1]]))
  }
  hashed <- system2("sha256sum", path, stdout = TRUE, stderr = FALSE)
  sub("\\s.*$", "", hashed[[1]])
}

test_that("frozen rotation calibration evidence is complete and decision-stable", {
  pkg_root <- normalizePath(testthat::test_path("../.."))

  artifact <- function(name) {
    path <- system.file("extdata", name, package = "dkge")
    if (!nzchar(path)) {
      path <- system.file("extdata", paste0(name, ".gz"), package = "dkge")
    }
    expect_true(nzchar(path), info = name)
    path
  }
  replicates <- read.csv(
    artifact("dkge-between-rotation-replicates.csv"),
    stringsAsFactors = FALSE
  )
  summary <- read.csv(
    artifact("dkge-between-rotation-summary.csv"),
    stringsAsFactors = FALSE
  )
  metadata <- read.csv(
    artifact("dkge-between-rotation-metadata.csv"),
    stringsAsFactors = FALSE
  )
  metadata <- stats::setNames(metadata$value, metadata$key)

  expect_equal(nrow(replicates), 8100L)
  expect_equal(nrow(summary), 14L)
  expect_true(all(is.finite(replicates$p)))
  expect_true(all(replicates$p > 0 & replicates$p <= 1))
  expect_equal(sort(unique(replicates$B)), c(199L, 399L))
  plan_path <- file.path(pkg_root, "data-raw/dkge-between-rotation-plan.md")
  runner_path <- file.path(pkg_root, "dev/calibrate-dkge-between-rotation.R")
  if (file.exists(plan_path) && file.exists(runner_path)) {
    expect_equal(
      unname(metadata[["plan_sha256"]]),
      .audit0817_file_sha256(plan_path)
    )
    expect_equal(
      unname(metadata[["runner_sha256"]]),
      .audit0817_file_sha256(runner_path)
    )
  }

  null <- replicates[replicates$arm != "strong_interaction_power", ]
  key <- interaction(null$arm, null$error, null$term, null$method,
                     drop = TRUE, lex.order = TRUE)
  rejection <- vapply(split(null$reject_005, key), mean, numeric(1))
  summary_null <- summary[summary$gate_class != "power", ]
  summary_key <- interaction(
    summary_null$arm, summary_null$error, summary_null$term,
    summary_null$method, drop = TRUE, lex.order = TRUE
  )
  expect_equal(unname(rejection[as.character(summary_key)]),
               summary_null$rejection_005, tolerance = 1e-15)

  global <- summary[summary$arm == "gaussian_global_null", ]
  expect_equal(
    global$rejection_005[match(c("group", "trait", "group:trait"), global$term)],
    c(0.049, 0.055, 0.058)
  )
  expect_true(all(global$gate_pass))

  partial <- summary[summary$arm == "gaussian_partial_null", ]
  expect_equal(
    partial$rejection_005[match(c("group", "trait", "group:trait"),
                                partial$term)],
    c(0.040, 0.042, 0.066)
  )
  expect_equal(
    partial$gate_pass[match(c("group", "trait", "group:trait"), partial$term)],
    c(TRUE, TRUE, FALSE)
  )

  robustness <- summary[summary$gate_class == "non_gaussian_robustness", ]
  expect_equal(nrow(robustness), 6L)
  expect_true(all(robustness$gate_pass))

  power <- summary[summary$gate_class == "power", ]
  expect_equal(power$rejection_005[match("rotation", power$method)], 0.43)
  expect_equal(power$rejection_005[match("freedman_lane", power$method)], 0.58)
  expect_false(power$gate_pass[match("rotation", power$method)])
  expect_equal(
    unname(metadata[c("gaussian_primary_all_pass",
                      "non_gaussian_all_pass", "power_pass")]),
    c("FALSE", "TRUE", "FALSE")
  )

  report_path <- file.path(pkg_root, "data-raw/dkge-between-rotation-report.md")
  if (file.exists(report_path)) {
    report <- paste(readLines(report_path), collapse = "\n")
    expect_true(grepl("(?i)size-adjusted comparison", report, perl = TRUE))
    expect_true(grepl("0\\.086", report))
    expect_true(grepl("0\\.49", report))
  }

  vignette_path <- file.path(pkg_root, "vignettes/dkge-between-subjects.Rmd")
  if (file.exists(vignette_path)) {
    vignette <- paste(readLines(vignette_path), collapse = "\n")
    expect_true(grepl("B = 199", vignette, fixed = TRUE))
    expect_true(grepl("(?i)size inflation", vignette, perl = TRUE))
  }
})
