library(testthat)

test_that("frozen null-calibration evidence is complete and internally consistent", {
  artifact <- function(name) {
    path <- system.file("extdata", name, package = "dkge")
    expect_true(nzchar(path), info = name)
    path
  }
  replicates <- read.csv(artifact("dkge-null-calibration-replicates.csv"),
                         stringsAsFactors = FALSE)
  summary <- read.csv(artifact("dkge-null-calibration-summary.csv"),
                      stringsAsFactors = FALSE)
  metadata <- read.csv(artifact("dkge-null-calibration-metadata.csv"),
                       stringsAsFactors = FALSE)

  expect_equal(nrow(replicates), 9000L)
  expect_equal(nrow(summary), 12L)
  expect_equal(metadata$calibration_version, "null-calibration-v1")
  expect_equal(
    metadata$plan_sha256,
    "d671401031967ce8260ee4399c2cfe592af1ca4b0dd17289fbdd013fc31c7a74"
  )
  expect_true(all(is.finite(replicates$p_value)))
  expect_true(all(replicates$p_value > 0 & replicates$p_value <= 1))
  expect_equal(sort(unique(replicates$legacy_ks_threshold)), 0.05)
  expect_equal(sort(unique(replicates$legacy_quantile_threshold)),
               c(0.08, 0.09))

  key <- interaction(replicates$family, replicates$setting,
                     replicates$estimand, drop = TRUE, lex.order = TRUE)
  recomputed <- vapply(split(replicates$p_value, key),
                       function(p) mean(p <= 0.05), numeric(1))
  summary_key <- interaction(summary$family, summary$setting,
                             summary$estimand, drop = TRUE, lex.order = TRUE)
  expect_equal(unname(recomputed[as.character(summary_key)]),
               summary$empirical_size_0.05, tolerance = 1e-15)

  adaptive <- summary[summary$family == "adaptive_loso_signflip", ]
  expect_equal(sort(adaptive$estimand), c("kenergy", "kenergy_prec", "none"))
  expect_true(all(adaptive$replicates == 1000L))
  expect_true(all(adaptive$legacy_gate_passed))
  expect_true(all(adaptive$classification == "calibrated"))

  primary <- summary[summary$setting == "primary_n20", ]
  expect_equal(sort(primary$empirical_size_0.05), c(0.074, 0.076, 0.077))
  expect_true(all(primary$size_wilson_lower > 0.05))
  expect_false(any(primary$legacy_gate_passed))
  expect_true(all(primary$classification == "anti_conservative"))

  diagnostic80 <- summary[summary$setting == "diagnostic_n80", ]
  expect_equal(sort(diagnostic80$empirical_size_0.05), c(0.030, 0.048, 0.052))
  expect_true(all(diagnostic80$legacy_gate_passed))
  expect_equal(sort(diagnostic80$classification),
               c("calibrated", "calibrated", "conservative"))
})
