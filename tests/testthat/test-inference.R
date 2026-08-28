# test-inference.R
# Simple diagnostics for dkge_infer helpers

library(testthat)

test_that("aligned inference help states its exchangeability and calibration limits", {
  source_rd <- testthat::test_path(
    "..", "..", "man", "dkge_infer_aligned.Rd"
  )
  rd <- if (file.exists(source_rd)) {
    tools::parse_Rd(source_rd)
  } else {
    topic <- utils::help("dkge_infer_aligned", package = "dkge")
    expect_gt(length(topic), 0L)
    if (!length(topic)) return(invisible(NULL))
    utils:::.getHelpFile(topic)
  }
  rendered <- paste(
    capture.output(tools::Rd2txt(rd)),
    collapse = "\n"
  )
  expect_match(rendered, "subject is the sampling and exchangeability unit",
               fixed = TRUE)
  expect_match(rendered, "contrast-by-location", fixed = TRUE)
  expect_match(rendered, "flattening that complete family", fixed = TRUE)
  expect_match(rendered, "negative-control promotion gate did not",
               fixed = TRUE)
  expect_match(rendered, "equal[[:space:]]+subject weighting only")
})
library(dkge)

make_inference_fixture <- function(S = 5, q = 3, P = 4, T = 60, seed = 5151) {
  set.seed(seed)
  effects <- paste0("eff", seq_len(q))
  betas <- replicate(S, {
    mat <- matrix(rnorm(q * P), q, P)
    rownames(mat) <- effects
    mat
  }, simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(rnorm(T * q), T, q)
    X <- qr.Q(qr(X))
    colnames(X) <- effects
    X
  }, simplify = FALSE)
  fit <- dkge_fit(dkge_data(betas, designs = designs), K = diag(q), rank = 2)
  list(fit = fit, betas = betas, effects = effects)
}

test_that("dkge_infer returns expected structure", {
  fixture <- make_inference_fixture()
  expect_error(
    dkge_infer(fixture$fit, c(1, -1, 0)),
    "approximate|rank-truncated",
    class = "dkge_alignment_ineligible_error"
  )
  expect_error(
    dkge_infer(fixture$fit, c(1, -1, 0), transport = list()),
    "Transport inside.*retired",
    class = "dkge_alignment_ineligible_error"
  )
  res <- dkge_infer(
    fixture$fit, c(1, -1, 0), allow_approximate_alignment = TRUE
  )

  expect_s3_class(res, "dkge_inference")
  expect_equal(res$method, "loso")
  expect_equal(res$inference, "signflip")
  expect_equal(res$correction, "maxT")
  expect_equal(length(res$statistics), 1)
  expect_equal(length(res$p_values), 1)
  expect_false(anyNA(res$p_values[[1]]))
  expect_length(res$significant[[1]], ncol(fixture$betas[[1]]))
  expect_identical(res$metadata$estimator$status, "approximate")
  expect_true(res$metadata$estimator$approximate_override)

  expect_error(
    suppressWarnings(dkge_infer(
      fixture$fit, c(1, -1, 0), method = "analytic", n_perm = 100
    )),
    "approximate|analytic|rank-truncated",
    class = "dkge_alignment_ineligible_error"
  )
  analytic <- suppressWarnings(dkge_infer(
    fixture$fit,
    c(1, -1, 0),
    method = "analytic",
    n_perm = 100,
    allow_approximate_alignment = TRUE
  ))
  expect_s3_class(analytic, "dkge_inference")
  expect_identical(analytic$method, "analytic")
  expect_identical(analytic$metadata$estimator$status, "approximate")
  expect_match(analytic$metadata$estimator$reason, "analytic LOSO")
  expect_true(analytic$metadata$estimator$approximate_override)
})

test_that("unsupported parametric max-T combinations fail closed", {
  fixture <- make_inference_fixture()
  expect_error(
    dkge_infer(
      fixture$fit,
      c(1, -1, 0),
      inference = "parametric",
      correction = "maxT",
      allow_approximate_alignment = TRUE
    ),
    "does not define a max-T null distribution",
    class = "dkge_inference_correction_error"
  )
})

test_that("multi-contrast max-T uses one family and shared subject signs", {
  set.seed(260854L)
  subjects <- paste0("s", 1:8)
  Y1 <- matrix(rnorm(8 * 3), 8, 3,
               dimnames = list(subjects, paste0("a", 1:3)))
  Y2 <- matrix(rnorm(8 * 2), 8, 2,
               dimnames = list(rev(subjects), paste0("b", 1:2)))
  stub <- list(contrasts = c("first", "second"), method = "fixture")

  set.seed(17)
  public <- dkge:::.infer_signflip(
    stub, 200L, "maxT", mapped_values = list(Y1, Y2)
  )
  family <- cbind(Y1, Y2[subjects, , drop = FALSE])
  colnames(family) <- c(paste0("first::", colnames(Y1)),
                        paste0("second::", colnames(Y2)))
  set.seed(17)
  oracle <- dkge_signflip_maxT(family, B = 200L)

  expect_equal(unlist(public$statistics, use.names = FALSE),
               unname(oracle$stat), tolerance = 0)
  expect_equal(unlist(public$p_adjusted, use.names = FALSE),
               unname(oracle$p), tolerance = 0)
  expect_identical(public$metadata$family_scope,
                   "all_contrasts_and_support_locations")
  expect_true(public$metadata$shared_subject_signs)
  expect_equal(public$metadata$family_dimensions, c(first = 3L, second = 2L))

  set.seed(17)
  reordered <- dkge:::.infer_signflip(
    list(contrasts = c("second", "first"), method = "fixture"),
    200L, "maxT", mapped_values = list(Y2, Y1)
  )
  expect_equal(reordered$p_adjusted[[2L]], public$p_adjusted[[1L]], tolerance = 0)
  expect_equal(reordered$p_adjusted[[1L]], public$p_adjusted[[2L]], tolerance = 0)
})

test_that("marginal sign flips and corrections use one unequal-width family", {
  set.seed(260855L)
  subjects <- paste0("s", 1:8)
  Y1 <- matrix(rnorm(8 * 3), 8, 3,
               dimnames = list(subjects, paste0("a", 1:3)))
  Y2 <- matrix(rnorm(8 * 2), 8, 2,
               dimnames = list(rev(subjects), paste0("b", 1:2)))
  stub <- list(contrasts = c("first", "second"), method = "fixture")

  set.seed(23)
  raw <- dkge:::.infer_signflip(
    stub, 200L, "bonferroni",
    mapped_values = list(first = Y1, second = Y2)
  )
  flips <- raw$metadata$sign_flips
  oracle_p <- function(Y) {
    Y <- Y[rownames(flips), , drop = FALSE]
    n_subjects <- nrow(Y)
    statistic <- function(values) {
      colMeans(values) /
        (apply(values, 2, stats::sd) / sqrt(n_subjects) + 1e-12)
    }
    observed <- abs(statistic(Y))
    null <- vapply(seq_len(ncol(flips)), function(b) {
      abs(statistic(flips[, b] * Y))
    }, numeric(ncol(Y)))
    (1 + rowSums(null >= observed)) / (ncol(flips) + 1)
  }
  expect_equal(raw$p_values[[1L]], oracle_p(Y1), tolerance = 0,
               ignore_attr = TRUE)
  expect_equal(raw$p_values[[2L]], oracle_p(Y2), tolerance = 0,
               ignore_attr = TRUE)
  expect_true(raw$metadata$shared_subject_signs)
  expect_identical(raw$metadata$family_scope,
                   "all_contrasts_and_support_locations")
  expect_equal(raw$metadata$family_dimensions,
               c(first = 3L, second = 2L))

  family_p <- unlist(raw$p_values, use.names = FALSE)
  bonferroni <- dkge:::.apply_correction(raw, "bonferroni", 0.05)
  expect_equal(
    unlist(bonferroni$p_adjusted, use.names = FALSE),
    stats::p.adjust(family_p, method = "bonferroni"), tolerance = 0
  )
  expect_identical(
    bonferroni$metadata$adjustment_scope,
    "all_contrasts_and_support_locations"
  )
  expect_identical(bonferroni$metadata$n_family_tests, 5L)
  fdr <- dkge:::.apply_correction(raw, "fdr", 0.05)
  expect_equal(
    unlist(fdr$p_adjusted, use.names = FALSE),
    stats::p.adjust(family_p, method = "fdr"), tolerance = 0
  )

  set.seed(23)
  row_reordered <- dkge:::.infer_signflip(
    stub, 200L, "bonferroni",
    mapped_values = list(
      first = Y1[rev(subjects), , drop = FALSE],
      second = Y2[c("s3", "s8", "s1", "s6", "s2", "s7", "s4", "s5"),
                  , drop = FALSE]
    )
  )
  expect_equal(row_reordered$p_values, raw$p_values, tolerance = 0)

  set.seed(23)
  contrast_reordered <- dkge:::.infer_signflip(
    list(contrasts = c("second", "first"), method = "fixture"),
    200L, "bonferroni",
    mapped_values = list(second = Y2, first = Y1)
  )
  expect_equal(contrast_reordered$p_values[[2L]], raw$p_values[[1L]],
               tolerance = 0)
  expect_equal(contrast_reordered$p_values[[1L]], raw$p_values[[2L]],
               tolerance = 0)
})

test_that("dkge_infer errors when cluster counts differ without transport", {
  data <- create_mismatched_data()
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 2)

  expect_error(
    suppressWarnings(dkge_infer(
      fit, c(1, -1, 0), allow_approximate_alignment = TRUE
    )),
    "Subject cluster counts differ",
    fixed = FALSE
  )
})

test_that("dkge_infer refuses current static mapper-based transport", {
  data <- create_mismatched_data()
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 2)

  transport_cfg <- list(
    centroids = data$centroids,
    medoid = 1L,
    mapper = dkge_mapper_spec("ridge", lambda = 1e-2)
  )

  expect_error(
    suppressWarnings(dkge_infer(
      fit, c(1, -1, 0), transport = transport_cfg,
      allow_approximate_alignment = TRUE
    )),
    "followed by `dkge_infer_aligned",
    class = "dkge_alignment_ineligible_error"
  )
})

test_that("parametric helper consumes transported rows with unequal cluster counts", {
  data <- create_mismatched_data()
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 2)
  transport_cfg <- list(
    centroids = data$centroids,
    medoid = 1L,
    mapper = dkge_mapper_spec("ridge", lambda = 1e-2)
  )
  contrast <- suppressWarnings(
    dkge_contrast(fit, c(1, -1, 0), method = "loso")
  )
  transported <- suppressWarnings(dkge_transport_contrasts_to_medoid(
    fit, contrast, medoid = 1L, centroids = data$centroids,
    mapper = transport_cfg$mapper
  ))
  Y <- transported[[1]]$subj_values
  res <- dkge:::.infer_parametric(contrast, "none", list(Y))
  expected <- colMeans(Y) / (apply(Y, 2, stats::sd) / sqrt(nrow(Y)) + 1e-12)
  expect_equal(res$statistics[[1]], expected, tolerance = 1e-12)
  expect_length(res$p_adjusted[[1]], ncol(Y))
  expect_equal(res$p_adjusted[[1]], res$p_values[[1]])
})

test_that("mean max-T rejects retired center modes at the public boundary", {
  Y <- matrix(rnorm(6 * 3), 6, 3)
  expect_error(
    suppressWarnings(dkge_signflip_maxT(Y, B = 100, center = "median")),
    "only.*mean", class = "dkge_inference_center_error"
  )
  expect_error(
    suppressWarnings(dkge_inference_spec(B = 100, center = "none")),
    "only.*mean", class = "dkge_inference_center_error"
  )
  expect_no_error(dkge_inference_spec(B = 100, center = "mean"))
})

test_that("sign-flip inputs and permutation counts fail closed", {
  Y <- matrix(
    seq_len(18), 6, 3,
    dimnames = list(paste0("s", 1:6), paste0("f", 1:3))
  )

  for (nonfinite in c(NA_real_, NaN, Inf, -Inf)) {
    bad <- Y
    bad[1, 1] <- nonfinite
    expect_error(
      dkge_signflip_maxT(bad, B = 100),
      "finite numeric matrix",
      class = "dkge_inference_data_error"
    )
  }
  expect_error(
    dkge_signflip_maxT(matrix(character(), 6, 0), B = 100),
    class = "dkge_inference_data_error"
  )
  expect_error(
    dkge_signflip_maxT(matrix(letters[1:18], 6, 3), B = 100),
    class = "dkge_inference_data_error"
  )
  for (bad_count in list(99, 100.5, NA_real_, Inf, c(100, 101))) {
    expect_error(
      dkge_signflip_maxT(Y, B = bad_count),
      "finite whole number",
      class = "dkge_inference_permutation_error"
    )
  }

  duplicate_features <- Y
  colnames(duplicate_features) <- c("same", "same", "other")
  expect_error(
    dkge_signflip_maxT(duplicate_features, B = 100),
    class = "dkge_inference_feature_error"
  )
  expect_no_error(dkge_signflip_maxT(Y, B = 100L))
})

test_that("all sign-flip correction routes reject fractional n_perm", {
  Y <- matrix(seq_len(18), 6, 3)
  stub <- list(contrasts = "effect", method = "fixture")
  expect_error(
    dkge:::.infer_signflip(
      stub, n_perm = 100.5, correction = "bonferroni",
      mapped_values = list(effect = Y)
    ),
    class = "dkge_inference_permutation_error"
  )
})

test_that("one-sided sign-flip max-T uses a signed (not absolute) null", {
  set.seed(1)
  S <- 12L; Q <- 6L
  Y <- matrix(rnorm(S * Q, sd = 0.3), S, Q)
  Y[, 1] <- Y[, 1] + 2   # strong positive effect in cluster 1

  set.seed(42); res_two <- dkge_signflip_maxT(Y, B = 500, tail = "two.sided")
  set.seed(42); res_gt  <- dkge_signflip_maxT(Y, B = 500, tail = "greater")

  # Same seed => identical sign flips. The greater-tail null is max(t_b), which
  # must be <= the two-sided null max(|t_b|), and strictly smaller for some draws
  # (the buggy version made them identical).
  expect_true(all(res_gt$maxnull <= res_two$maxnull + 1e-12))
  expect_false(isTRUE(all.equal(res_gt$maxnull, res_two$maxnull)))

  # One-sided test is at least as powerful for the positive effect.
  expect_lte(res_gt$p[1], res_two$p[1] + 1e-12)
})

test_that("sign-flip max-T exposes uncorrected p-values bounded by the adjusted ones", {
  set.seed(3)
  Y <- matrix(rnorm(12 * 5, sd = 0.5), 12, 5)
  Y[, 1] <- Y[, 1] + 1.5
  res <- dkge_signflip_maxT(Y, B = 400, tail = "two.sided")
  expect_length(res$p_unadj, ncol(Y))
  expect_true(all(res$p_unadj >= 0 & res$p_unadj <= 1))
  # Uncorrected p is never larger than the max-T (FWER) adjusted p.
  expect_true(all(res$p_unadj <= res$p + 1e-9))
})

test_that("sign-flip adjusted and unadjusted p-values match the sampled null", {
  set.seed(2718)
  Y <- matrix(c(
    2,  0,  1,
    1, -1,  0,
    3,  1, -1,
    2,  0,  2,
    4, -2,  1
  ), nrow = 5, byrow = TRUE,
  dimnames = list(paste0("person ", 1:5), c("signal", "null-x", "zero-ish")))
  B <- 100L
  result <- dkge_signflip_maxT(Y, B = B, tail = "two.sided")

  observed <- colMeans(Y) /
    (apply(Y, 2, stats::sd) / sqrt(nrow(Y)) + 1e-12)
  null_stats <- vapply(seq_len(B), function(b) {
    Yb <- result$flips[, b] * Y
    abs(colMeans(Yb) /
          (apply(Yb, 2, stats::sd) / sqrt(nrow(Yb)) + 1e-12))
  }, numeric(ncol(Y)))
  max_null <- apply(null_stats, 2, max)
  expected_adjusted <- vapply(abs(observed), function(x) {
    (1 + sum(max_null >= x)) / (B + 1)
  }, numeric(1))
  expected_unadjusted <- vapply(seq_along(observed), function(j) {
    (1 + sum(null_stats[j, ] >= abs(observed[j]))) / (B + 1)
  }, numeric(1))
  names(expected_unadjusted) <- colnames(Y)
  names(max_null) <- paste0("perm", seq_len(B))

  expect_named(result, c("stat", "p", "p_unadj", "maxnull", "flips"))
  expect_equal(result$stat, observed, tolerance = 1e-12)
  expect_equal(result$p, expected_adjusted, tolerance = 1e-12)
  expect_equal(result$p_unadj, expected_unadjusted, tolerance = 1e-12)
  expect_equal(result$maxnull, max_null, tolerance = 1e-12)
  expect_equal(names(result$stat), colnames(Y))
  expect_equal(names(result$p), colnames(Y))
  expect_equal(names(result$p_unadj), colnames(Y))
  expect_equal(rownames(result$flips), rownames(Y))
  expect_equal(colnames(result$flips), names(result$maxnull))
})

test_that("sign-flip schema is seeded and stable for minimum and degenerate inputs", {
  Y <- matrix(0, 5, 1,
              dimnames = list(paste0("s", 1:5), "constant-zero"))
  set.seed(99)
  first <- dkge_signflip_maxT(Y, B = 100)
  set.seed(99)
  second <- dkge_signflip_maxT(Y, B = 100)

  expect_identical(first, second)
  expect_named(first, c("stat", "p", "p_unadj", "maxnull", "flips"))
  expect_equal(first$stat, c("constant-zero" = 0))
  expect_equal(first$p, c("constant-zero" = 1))
  expect_equal(first$p_unadj, c("constant-zero" = 1))
  expect_true(all(first$maxnull == 0))
  expect_equal(dim(first$flips), c(5, 100))
})
