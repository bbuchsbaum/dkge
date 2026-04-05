# test-reliability-weights.R
# Tests that verify the mathematical correctness of adapt="reliability" weights.

# Helper: build B_list / B_list2 with known per-voxel correlation structure.
# Voxels in `signal_voxels` share a subject-level random effect across both runs;
# voxels in `noise_voxels` are independent noise in both runs.
.make_reliability_data <- function(S = 6, q = 4, V = 30,
                                   signal_voxels = 1:10,
                                   noise_snr = 0.1,
                                   seed = 42) {
  set.seed(seed)
  subject_signal <- rnorm(S)  # true subject-level latent

  make_run <- function() {
    lapply(seq_len(S), function(s) {
      B <- matrix(rnorm(q * V, sd = 1), q, V)           # pure noise everywhere
      B[, signal_voxels] <- B[, signal_voxels] +
        matrix(rep(subject_signal[s] / noise_snr, q * length(signal_voxels)),
               q, length(signal_voxels))
      B
    })
  }
  list(run1 = make_run(), run2 = make_run(),
       subject_signal = subject_signal,
       signal_voxels = signal_voxels,
       noise_voxels = setdiff(seq_len(V), signal_voxels))
}

test_that("reliability weights rank high-correlation voxels above noise voxels", {
  dat <- .make_reliability_data(S = 8, V = 40, signal_voxels = 1:15, noise_snr = 0.05)

  w <- .dkge_adapt_weights(dat$run1, adapt = "reliability",
                            B_list2 = dat$run2, winsor = 0.99)

  expect_length(w, 40)
  expect_true(all(is.finite(w)))
  expect_true(all(w > 0))

  mean_signal <- mean(w[dat$signal_voxels])
  mean_noise  <- mean(w[dat$noise_voxels])

  # Signal voxels (high cross-run correlation) must have higher average weight
  expect_gt(mean_signal, mean_noise)
})

test_that("reliability weight is near 1 for perfectly correlated voxels", {
  S <- 6
  q <- 3
  V <- 5
  subject_signal <- seq_len(S) * 10  # deterministic, large relative to noise

  # Identical runs (zero noise) → r = 1, weight = 1 + 1e-6
  B_perfect <- lapply(seq_len(S), function(s) {
    matrix(rep(subject_signal[s], q * V), q, V)
  })
  w <- .dkge_adapt_weights(B_perfect, adapt = "reliability",
                            B_list2 = B_perfect, winsor = 1.0)

  # All correlations should be 1 → weights ≈ 1
  expect_true(all(abs(w - 1) < 1e-4))
})

test_that("reliability weight falls back to uniform when B_list2 is NULL", {
  set.seed(1)
  B <- lapply(1:4, function(s) matrix(rnorm(20), 4, 5))
  w <- withCallingHandlers(
    .dkge_adapt_weights(B, adapt = "reliability", B_list2 = NULL),
    warning = function(w) invokeRestart("muffleWarning")
  )
  expect_equal(w, rep(1, 5))
})

test_that("reliability weights handle zero-variance voxels gracefully", {
  S <- 5
  q <- 3
  V <- 6
  # All subjects have identical values in voxel 1 → zero variance → r = NaN → weight = 1e-6
  B1 <- lapply(seq_len(S), function(s) {
    B <- matrix(rnorm(q * V), q, V)
    B[, 1] <- 0  # zero variance voxel in run 1
    B
  })
  B2 <- lapply(seq_len(S), function(s) {
    B <- matrix(rnorm(q * V), q, V)
    B[, 1] <- 0  # also zero in run 2
    B
  })
  w <- .dkge_adapt_weights(B1, adapt = "reliability", B_list2 = B2, winsor = 1.0)
  expect_true(all(is.finite(w)))
  expect_true(all(w > 0))
  # Zero-variance voxel gets the minimum weight
  expect_equal(w[1], 1e-6)
})

test_that("dkge_weights stores B_list2 and reports it in print", {
  B2 <- lapply(1:4, function(s) matrix(rnorm(20), 4, 5))
  w <- dkge_weights(adapt = "reliability", B_list2 = B2)
  expect_identical(w$B_list2, B2)
  expect_equal(w$adapt, "reliability")

  output <- capture.output(print(w))
  expect_true(any(grepl("B_list2=S4", output)))
})

test_that("dkge_weights warns when reliability requested without B_list2", {
  expect_warning(
    dkge_weights(adapt = "reliability"),
    regexp = "requires B_list2"
  )
})

test_that("reliability weights integrate into full fit pipeline", {
  toy  <- dkge_sim_toy(factors = list(A = list(L = 2)),
                       active_terms = "A", S = 4, P = 20, seed = 7)
  # Second run: add small noise to betas
  B2   <- lapply(toy$B_list, function(B) B + matrix(rnorm(length(B), sd = 0.05), nrow(B)))
  w    <- dkge_weights(adapt = "reliability", B_list2 = B2)
  fit2 <- dkge_update_weights(
    dkge(toy$B_list, toy$X_list, K = toy$K, rank = 1),
    weights = w
  )
  expect_s3_class(fit2, "dkge")
  expect_equal(fit2$rank, 1L)
})
