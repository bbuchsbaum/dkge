library(testthat)

test_that("chunked trial reduction matches the dense constructor", {
  set.seed(19401)
  cell <- rep(1:3, each = 8)
  X <- model.matrix(~ 0 + factor(cell, levels = 1:3))
  colnames(X) <- paste0("e", 1:3)
  truth <- matrix(rnorm(3 * 9), 3, 9)
  Y <- X %*% truth + matrix(rnorm(nrow(X) * 9, sd = 0.4), nrow(X), 9)
  colnames(Y) <- paste0("v", 1:9)
  covariance <- outer(seq_len(nrow(X)), seq_len(nrow(X)),
                      function(i, j) 0.35^abs(i - j))

  dense <- dkge_trial_subject(Y, X, trial_covariance = covariance)
  chunked <- dkge_trial_subject_chunks(
    list(Y[, 1:2], Y[, 3:6], Y[, 7:9]), X,
    trial_covariance = covariance
  )
  expect_equal(chunked$beta, dense$beta, tolerance = 1e-11,
               ignore_attr = TRUE)
  expect_equal(chunked$effect_information, dense$effect_information,
               tolerance = 1e-11, ignore_attr = TRUE)
  expect_null(chunked$effect_score)
  expect_equal(
    chunked$effect_information %*% chunked$beta,
    dense$effect_score,
    tolerance = 1e-11,
    ignore_attr = TRUE
  )
  expect_equal(chunked$residual_variance, dense$residual_variance,
               tolerance = 1e-11, ignore_attr = TRUE)
  expect_equal(chunked$residual_sum_squares, dense$residual_sum_squares,
               tolerance = 1e-11, ignore_attr = TRUE)
  expect_equal(chunked$noise_trace, dense$noise_trace, tolerance = 1e-11)
  expect_equal(chunked$input_type, "trialwise_chunks")
  expect_equal(chunked$error_model$n_chunks, 3L)
  expect_false(chunked$error_model$effect_score_retained)
})

test_that("function chunk sources terminate without retaining trialwise Y", {
  set.seed(19403)
  X <- model.matrix(~ 0 + factor(rep(1:2, each = 8)))
  colnames(X) <- c("e1", "e2")
  blocks <- lapply(c(3, 4, 2), function(width) {
    matrix(rnorm(nrow(X) * width), nrow(X), width)
  })
  names_seen <- 0L
  source <- function(i) {
    if (i > length(blocks)) return(NULL)
    block <- blocks[[i]]
    colnames(block) <- paste0("v", names_seen + seq_len(ncol(block)))
    names_seen <<- names_seen + ncol(block)
    block
  }
  subject <- dkge_trial_subject_chunks(source, X)
  expect_equal(dim(subject$beta), c(2L, 9L))
  expect_null(subject$y)
  expect_null(subject$trialwise_y)
  expect_equal(subject$error_model$storage, "feature_chunked")
  expect_equal(object.size(subject$design), object.size(X))
})

test_that("chunked run splits combine halves before reliability estimation", {
  set.seed(19405)
  rows <- expand.grid(repeat_id = 1:2, effect = 1:2, run = 1:4)
  X <- model.matrix(~ 0 + factor(rows$effect, levels = 1:2))
  colnames(X) <- c("e1", "e2")
  Y <- matrix(rnorm(nrow(X) * 8), nrow(X), 8)
  run <- paste0("run", rows$run)
  dense <- dkge_trial_subject(
    Y, X, split = "run", run_labels = run,
    effect_precision = "split_half"
  )
  chunked <- dkge_trial_subject_chunks(
    list(Y[, 1:3], Y[, 4:8]), X, split = "run", run_labels = run,
    effect_precision = "split_half"
  )
  expect_equal(chunked$split_betas, dense$split_betas,
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(chunked$effect_precision, dense$effect_precision,
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_true(chunked$split_provenance$run_disjoint)
})

test_that("chunk sources fail clearly on malformed or unterminated input", {
  X <- matrix(1, 4, 1, dimnames = list(NULL, "mean"))
  expect_error(dkge_trial_subject_chunks(list(), X), "no feature blocks")
  expect_error(
    dkge_trial_subject_chunks(function(i) matrix(1, 4, 1), X,
                              max_chunks = 2),
    "max_chunks"
  )
  expect_error(
    dkge_trial_subject_chunks(list(matrix(1, 3, 2)), X),
    "T x P_block"
  )
})
