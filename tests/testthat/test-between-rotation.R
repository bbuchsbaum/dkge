test_that("Haar draws are orthogonal and isolated seeds are reproducible", {
  for (d in c(1L, 2L, 7L, 15L)) {
    G <- dkge:::.dkge_seeded_haar_rotation(1000L + d, d)
    expect_equal(crossprod(G), diag(d), tolerance = 1e-12)
    expect_equal(G, dkge:::.dkge_seeded_haar_rotation(1000L + d, d))
  }
  expect_error(dkge:::.dkge_haar_orthogonal(0), "positive integer")
})

test_that("compressed rotation matches dense reconstruction invariants", {
  fixture <- dkge_test_rotation_simulate(seed = 6101, n = 14, p = 9)
  term <- "trait"
  Xred <- dkge:::.dkge_between_reduced_X(
    fixture$design, fixture$X, term
  )
  full_setup <- dkge:::.dkge_between_rrr_fast_setup(fixture$X)
  red_setup <- dkge:::.dkge_between_rrr_fast_setup(Xred)
  setup <- dkge:::.dkge_between_rotation_setup(
    full_setup, red_setup, fixture$Y
  )
  G <- dkge:::.dkge_seeded_haar_rotation(6102, setup$residual_dimension)
  compressed <- dkge:::.dkge_between_rotation_eval_input(setup, G)
  Yb <- setup$Q0 %*% setup$fitted_coordinates +
    setup$Qe %*% G %*% setup$residual_coordinates

  expect_equal(compressed$C_full, crossprod(full_setup$Q, Yb),
               tolerance = 1e-11)
  expect_equal(compressed$C_red, crossprod(red_setup$Q, Yb),
               tolerance = 1e-11)
  expect_equal(compressed$yss, sum(Yb^2), tolerance = 1e-11)
  expect_equal(crossprod(Yb), crossprod(fixture$Y), tolerance = 1e-10)
  expect_equal(crossprod(setup$Q0, Yb), setup$fitted_coordinates,
               tolerance = 1e-11)
})

test_that("compressed rotation null matches an explicit dense refit loop", {
  fixture <- dkge_test_rotation_simulate(seed = 6111, n = 16, p = 8)
  B <- 17L
  seed <- 6112L
  term <- "group:trait"
  result <- dkge_between_permute(
    fixture$fit, terms = term, method = "rotation", B = B, seed = seed,
    scope = "both", feature_adjust = "maxT"
  )
  oracle <- dkge_test_rotation_dense_result(fixture, term, B, seed)

  expect_equal(result$tests[[term]]$null, oracle$null, tolerance = 1e-10)
  observed_red <- dkge:::.dkge_rrr_fit_core(
    dkge:::.dkge_between_reduced_X(fixture$design, fixture$X, term),
    fixture$Y, rank = fixture$fit$rank, tol = fixture$fit$tol,
    warn_rank = FALSE
  )
  expect_equal(result$tests[[term]]$statistic,
               dkge:::.dkge_between_stat(fixture$fit, observed_red),
               tolerance = 1e-11)
  expect_equal(result$summary$p, oracle$p, tolerance = 1e-12)
  expect_equal(result$feature_tests[[term]]$statistic,
               oracle$feature_statistic, tolerance = 1e-10)
  expect_equal(result$feature_tests[[term]]$p,
               oracle$feature_p, tolerance = 1e-12)
  expect_equal(result$feature_tests[[term]]$null_max,
               oracle$null_max, tolerance = 1e-10)
  expect_equal(result$feature_tests[[term]]$p_adjusted,
               oracle$maxT_p, tolerance = 1e-12)
})

test_that("rotation supports feature inference and records assumptions", {
  fixture <- dkge_test_rotation_simulate(
    seed = 6121, n = 18, p = 7,
    power_term = "group:trait", effect_size = 1
  )
  result <- dkge_between_permute(
    fixture$fit, terms = "group:trait", method = "rotation",
    B = 39, seed = 6122, scope = "both", feature_adjust = "maxT"
  )

  expect_equal(result$method, "rotation")
  expect_match(result$resampling_assumptions$exact_under, "matrix-normal")
  expect_equal(length(result$feature_tests[["group:trait"]]$p),
               ncol(fixture$Y))
  expect_identical(result$feature_tests[["group:trait"]]$feature_ids,
                   colnames(fixture$Y))
  expect_true(all(result$feature_tests[["group:trait"]]$p_adjusted >= 0))
  expect_true(all(result$feature_tests[["group:trait"]]$p_adjusted <= 1))

  fdr <- dkge_between_permute(
    fixture$fit, terms = "group:trait", method = "rotation",
    B = 19, seed = 6123, scope = "features", feature_adjust = "fdr"
  )
  expect_equal(
    fdr$feature_tests[["group:trait"]]$p_adjusted,
    stats::p.adjust(fdr$feature_tests[["group:trait"]]$p, method = "fdr")
  )
})

test_that("rotation setup avoids a complete subject-square QR basis", {
  implementation <- paste(deparse(body(dkge:::.dkge_between_rotation_setup)),
                          collapse = "\n")
  expect_false(grepl("complete = TRUE", implementation, fixed = TRUE))
  fixture <- dkge_test_rotation_simulate(seed = 6124, n = 17, p = 4)
  Xred <- dkge:::.dkge_between_reduced_X(
    fixture$design, fixture$X, "trait"
  )
  setup <- dkge:::.dkge_between_rotation_setup(
    dkge:::.dkge_between_rrr_fast_setup(fixture$X),
    dkge:::.dkge_between_rrr_fast_setup(Xred),
    fixture$Y
  )
  expect_equal(dim(setup$Qe),
               c(nrow(fixture$X), setup$residual_dimension))
})

test_that("seeded rotation restores caller RNG and is reproducible", {
  fixture <- dkge_test_rotation_simulate(seed = 6131, n = 14, p = 5)
  set.seed(6132)
  expected <- runif(5)
  set.seed(6132)
  first <- dkge_between_permute(
    fixture$fit, terms = "trait", method = "rotation", B = 11, seed = 6133
  )
  observed <- runif(5)
  second <- dkge_between_permute(
    fixture$fit, terms = "trait", method = "rotation", B = 11, seed = 6133
  )

  expect_equal(observed, expected)
  expect_equal(first$tests$trait$null, second$tests$trait$null)
  expect_equal(first$summary$p, second$summary$p)
})

test_that("rotation serial and parallel paths are bit-identical", {
  fixture <- dkge_test_rotation_simulate(seed = 6141, n = 14, p = 6)
  serial <- dkge_between_permute(
    fixture$fit, terms = "trait", method = "rotation", B = 9, seed = 6142,
    scope = "both", feature_adjust = "maxT"
  )
  parallel <- suppressWarnings(dkge_between_permute(
    fixture$fit, terms = "trait", method = "rotation", B = 9, seed = 6142,
    scope = "both", feature_adjust = "maxT", parallel = TRUE
  ))

  expect_identical(serial$tests$trait$null, parallel$tests$trait$null)
  expect_identical(serial$summary$p, parallel$summary$p)
  expect_identical(serial$feature_tests$trait$p,
                   parallel$feature_tests$trait$p)
  expect_identical(serial$feature_tests$trait$null_max,
                   parallel$feature_tests$trait$null_max)
})

test_that("rotation fails closed outside its validated scope", {
  fixture <- dkge_test_rotation_simulate(seed = 6151, n = 12, p = 5)
  weighted <- dkge_between_rrr(
    dkge_make_target(Y = fixture$Y, subject_ids = fixture$data$subject_id),
    fixture$design, rank = 2,
    weights = list(subject = seq(0.5, 1.5, length.out = 12))
  )
  expect_error(
    dkge_between_permute(weighted, terms = "trait", method = "rotation", B = 5),
    "requires an unweighted"
  )
  expect_error(
    dkge_between_permute(
      fixture$fit, terms = "trait", method = "rotation", B = 5,
      blocks = rep(1:2, each = 6)
    ),
    "one exchangeability block"
  )
  expect_error(
    dkge_between_permute(
      fixture$fit, terms = "trait", method = "rotation", B = 5,
      blocks = c(rep(1, 11), NA)
    ),
    "cannot contain missing"
  )

  X <- cbind("(Intercept)" = 1, x = c(-1, 1))
  rownames(X) <- c("s1", "s2")
  Y <- matrix(c(-1, 1, 1, -1), 2, 2,
              dimnames = list(rownames(X), c("f1", "f2")))
  saturated <- dkge_between_rrr(dkge_make_target(Y = Y), X, rank = 1)
  expect_error(
    dkge_between_permute(
      saturated, terms = "x", method = "rotation", B = 5
    ),
    "at least two residual dimensions"
  )
})

test_that("rotation detects a strong rank-one interaction", {
  fixture <- dkge_test_rotation_simulate(
    seed = 6161, n = 20, p = 12,
    power_term = "group:trait", effect_size = 1.5
  )
  result <- dkge_between_permute(
    fixture$fit, terms = "group:trait", method = "rotation",
    B = 99, seed = 6162
  )
  expect_lte(result$summary$p, 0.05)
})

test_that("two-term serial and parallel nulls match and term 2 is order-independent", {
  skip_if_not_installed("future")
  fixture <- dkge_test_rotation_simulate(seed = 5150, n = 14, p = 6)
  terms <- c("trait", "group:trait")
  run <- function(method, which_terms, parallel = FALSE) {
    dkge_between_permute(
      fixture$fit, terms = which_terms, method = method,
      B = 25, seed = 5150, parallel = parallel
    )
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::sequential)

  for (method in c("rotation", "freedman_lane")) {
    serial <- run(method, terms, parallel = FALSE)
    parallel_seq <- suppressWarnings(run(method, terms, parallel = TRUE))
    for (term in terms) {
      expect_identical(serial$tests[[term]]$null,
                       parallel_seq$tests[[term]]$null,
                       info = paste(method, "sequential", term))
    }

    only_term2 <- run(method, "group:trait", parallel = FALSE)
    expect_identical(serial$tests[["group:trait"]]$null,
                     only_term2$tests[["group:trait"]]$null,
                     info = paste(method, "term-2 alone"))
  }

  n_cores <- tryCatch(future::availableCores(), error = function(e) 1L)
  if (!is.finite(n_cores) || n_cores < 2L) {
    skip("fewer than 2 cores for future::multisession")
  }
  workers <- tryCatch({
    future::plan(future::multisession, workers = 2)
    future::value(future::future({
      requireNamespace("dkge", quietly = TRUE)
    }))
  }, error = function(e) FALSE)
  if (!isTRUE(workers)) {
    skip("future::multisession workers cannot load dkge")
  }
  for (method in c("rotation", "freedman_lane")) {
    serial <- run(method, terms, parallel = FALSE)
    parallel_ms <- suppressWarnings(run(method, terms, parallel = TRUE))
    for (term in terms) {
      expect_identical(serial$tests[[term]]$null,
                       parallel_ms$tests[[term]]$null,
                       info = paste(method, "multisession", term))
    }
  }
})
