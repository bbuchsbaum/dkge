library(testthat)

test_that("tracked q = 60 performance evidence matches the frozen contract", {
  path <- system.file("extdata", "dkge-q60-performance.csv", package = "dkge")
  expect_true(nzchar(path))
  evidence <- utils::read.csv(path)
  runs <- evidence[evidence$row_type == "run", ]
  summaries <- evidence[evidence$row_type == "summary", ]

  workloads <- c(
    "trial_chunks", "group_legacy", "group_advanced", "fold_refit",
    "poisson_bootstrap"
  )
  expect_equal(nrow(runs), 15L)
  expect_equal(nrow(summaries), 5L)
  expect_setequal(unique(runs$workload), workloads)
  expect_equal(sort(unique(runs$seed)), 19421:19423)
  expect_true(all(runs$q == 60L))
  expect_true(all(is.finite(runs$peak_rss_bytes)))
  expect_false(any(runs$p_by_p_detected))

  expected_P <- c(
    trial_chunks = 100000L,
    group_legacy = 10000L,
    group_advanced = 10000L,
    fold_refit = 10000L,
    poisson_bootstrap = 2000L
  )
  expect_equal(
    stats::setNames(summaries$P, summaries$workload)[names(expected_P)],
    expected_P,
    ignore_attr = TRUE
  )
  forbidden_bytes <- as.numeric(runs$P)^2 * 8
  expect_true(all(runs$largest_allocation_bytes < forbidden_bytes))

  imaging <- summaries[summaries$workload == "trial_chunks", ]
  expect_lte(imaging$object_bytes, 100e6)
  expect_lte(imaging$design_bytes, 0.25e6)
  expect_lte(imaging$largest_allocation_bytes, 100e6)
  expect_lte(imaging$largest_allocation_bytes,
             as.numeric(imaging$P)^2 * 8 / 1000)

  limits <- data.frame(
    workload = workloads,
    elapsed = c(20, 20, 30, 30, 60),
    rss = c(1.25e9, 2e9, 2e9, 2e9, 1.5e9)
  )
  recomputed <- vapply(seq_len(nrow(summaries)), function(i) {
    row <- summaries[i, ]
    limit <- limits[limits$workload == row$workload, ]
    extra <- if (row$workload == "trial_chunks") {
      row$object_bytes <= 100e6 && row$design_bytes <= 0.25e6 &&
        row$largest_allocation_bytes <= 100e6
    } else {
      TRUE
    }
    row$elapsed_seconds <= limit$elapsed && row$peak_rss_bytes <= limit$rss &&
      row$elapsed_iqr <= row$elapsed_seconds && !row$p_by_p_detected && extra
  }, logical(1))
  expect_identical(as.logical(summaries$passed), recomputed)

  # The tracked host run preserves one scheduling-dispersion miss rather than
  # rerunning until green. Its time, RSS, and allocation-safety gates pass.
  expect_equal(sum(summaries$passed), 4L)
  expect_setequal(summaries$workload[!summaries$passed], "poisson_bootstrap")
  poisson <- summaries[summaries$workload == "poisson_bootstrap", ]
  expect_lte(poisson$elapsed_seconds, 60)
  expect_lte(poisson$peak_rss_bytes, 1.5e9)
  expect_gt(poisson$elapsed_iqr, poisson$elapsed_seconds)
})
