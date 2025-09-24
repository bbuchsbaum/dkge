testthat::local_edition(3)

quick_configs <- c("none", "kenergy", "kenergy_prec")

test_that("Null p-values are approximately Uniform(0,1) with LOSO and adaptive weights (quick)", {
  for (adapt in quick_configs) {
    p <- dkge_test_null_uniformity(nrep = 24,
                                   nsub = 12,
                                   V = 96,
                                   adapt = adapt,
                                   n_perm = 149,
                                   seed = switch(adapt,
                                     none = 11L,
                                     kenergy = 21L,
                                     kenergy_prec = 31L))

    expect_true(abs(mean(p) - 0.5) < 0.12, info = paste("mean p for", adapt))
    expect_true(abs(stats::median(p) - 0.5) < 0.12, info = paste("median p for", adapt))
    ks <- suppressWarnings(stats::ks.test(p, "punif"))
    expect_true(ks$p.value > 0.01, info = paste("KS test for", adapt))
  }
})

test_that("Null p-values ~ Uniform(0,1) in a longer run (opt-in)", {
  testthat::skip_on_cran()
  if (!identical(Sys.getenv("DKGE_LONG_TESTS", "false"), "true")) {
    testthat::skip("Set DKGE_LONG_TESTS=true to run the long null simulation test.")
  }

  p_none  <- dkge_test_null_uniformity(nrep = 120, nsub = 16, V = 128,
                                       adapt = "none",          n_perm = 499, seed = 101)
  p_ke    <- dkge_test_null_uniformity(nrep = 120, nsub = 16, V = 128,
                                       adapt = "kenergy",       n_perm = 499, seed = 201)
  p_kprec <- dkge_test_null_uniformity(nrep = 120, nsub = 16, V = 128,
                                       adapt = "kenergy_prec",  n_perm = 499, seed = 301)

  for (p in list(p_none, p_ke, p_kprec)) {
    ks <- suppressWarnings(stats::ks.test(p, "punif"))
    expect_gt(ks$p.value, 0.05)
    q_emp <- as.numeric(stats::quantile(p, probs = c(0.1, 0.25, 0.5, 0.75, 0.9)))
    expect_lt(max(abs(q_emp - c(0.1, 0.25, 0.5, 0.75, 0.9))), 0.08)
  }
})
