testthat::local_edition(3)

test_that("Between-subject permutation p-values are approximately Uniform(0,1) under the null (quick)", {
  pvals <- dkge_test_between_null_uniformity(
    nrep = 24,
    n = 18,
    p = 20,
    rank = 2,
    B = 149,
    seed = 400
  )

  for (term in names(pvals)) {
    p <- pvals[[term]]
    expect_true(abs(mean(p) - 0.5) < 0.14, info = paste("mean p for", term))
    expect_true(abs(stats::median(p) - 0.5) < 0.16, info = paste("median p for", term))
    ks <- suppressWarnings(stats::ks.test(p, "punif"))
    expect_true(ks$p.value > 0.01, info = paste("KS test for", term))
  }
})

test_that("Between-subject permutation p-values remain calibrated in a longer null run (opt-in)", {
  testthat::skip_on_cran()
  if (!identical(Sys.getenv("DKGE_LONG_TESTS", "false"), "true")) {
    testthat::skip("Set DKGE_LONG_TESTS=true to run the long between-subject null test.")
  }

  pvals <- dkge_test_between_null_uniformity(
    nrep = 120,
    n = 20,
    p = 24,
    rank = 2,
    B = 399,
    seed = 800
  )

  for (term in names(pvals)) {
    p <- pvals[[term]]
    ks <- suppressWarnings(stats::ks.test(p, "punif"))
    expect_gt(ks$p.value, 0.05, info = paste("KS test for", term))
    q_emp <- as.numeric(stats::quantile(p, probs = c(0.1, 0.25, 0.5, 0.75, 0.9)))
    expect_lt(max(abs(q_emp - c(0.1, 0.25, 0.5, 0.75, 0.9))),
              0.09,
              info = paste("quantiles for", term))
  }
})
