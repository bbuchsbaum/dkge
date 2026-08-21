# test-inference.R
# Simple diagnostics for dkge_infer helpers

library(testthat)
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
  res <- dkge_infer(fixture$fit, c(1, -1, 0))

  expect_s3_class(res, "dkge_inference")
  expect_equal(res$method, "loso")
  expect_equal(res$inference, "signflip")
  expect_equal(res$correction, "maxT")
  expect_equal(length(res$statistics), 1)
  expect_equal(length(res$p_values), 1)
  expect_false(anyNA(res$p_values[[1]]))
  expect_length(res$significant[[1]], ncol(fixture$betas[[1]]))
})

test_that("dkge_infer errors when cluster counts differ without transport", {
  data <- create_mismatched_data()
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 2)

  expect_error(
    suppressWarnings(dkge_infer(fit, c(1, -1, 0))),
    "Subject cluster counts differ",
    fixed = FALSE
  )
})

test_that("dkge_infer applies mapper-based transport", {
  data <- create_mismatched_data()
  fit <- dkge_fit(data$betas, data$designs, K = data$K, rank = 2)

  transport_cfg <- list(
    centroids = data$centroids,
    medoid = 1L,
    mapper = dkge_mapper_spec("ridge", lambda = 1e-2),
    betas = data$betas
  )

  res <- suppressWarnings(dkge_infer(fit, c(1, -1, 0), transport = transport_cfg))

  expect_true(is.list(res$transport))
  expect_equal(ncol(res$transport[[1]]$subj_values), nrow(data$centroids[[1]]))
  expect_false(anyNA(res$p_values[[1]]))
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
