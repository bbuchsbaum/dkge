beta_fake_fit <- function(n_subjects = 4L) {
  structure(
    list(
      Btil = replicate(n_subjects, matrix(0, 1, 1), simplify = FALSE),
      subject_ids = paste0("s", seq_len(n_subjects))
    ),
    class = c("dkge", "list")
  )
}

test_that("fold voxel weights never recycle across parcel counts", {
  Bts <- matrix(1:6, nrow = 2)
  expect_error(
    .dkge_subject_loader_weights(c(1, 4), Bts),
    "length 2.*3 columns|must.*parcel",
    class = "dkge_weight_dimension_error"
  )
  expect_equal(.dkge_subject_loader_weights(2, Bts), rep(2, 3))
  expect_equal(.dkge_subject_loader_weights(c(2, 2), Bts), rep(2, 3))
  expect_equal(.dkge_subject_loader_weights(c(1, 2, 3), Bts), c(1, 2, 3))
})

test_that("custom folds are an exact partition unless advanced mode is explicit", {
  fit <- beta_fake_fit()
  expect_identical(eval(formals(dkge_define_folds)$type),
                   c("subject", "custom"))
  expect_error(
    dkge_define_folds(
      fit, type = "custom",
      assignments = list(c(1L, 2L), c(2L, 3L))
    ),
    "nonoverlapping partition",
    class = "dkge_fold_partition_error"
  )
  expect_error(
    dkge_define_folds(
      fit, type = "custom",
      assignments = list(c(1L, 2L), c(3L))
    ),
    "cover every subject",
    class = "dkge_fold_partition_error"
  )

  repeated <- dkge_define_folds(
    fit, type = "custom",
    assignments = list(c(1L, 2L), c(2L, 3L, 4L)),
    partition = "repeated"
  )
  expect_true(repeated$metadata$overlap)
  expect_identical(repeated$metadata$partition, "repeated")
  expect_error(
    .dkge_normalize_folds(repeated, fit, consumer = "K-fold contrasts"),
    "K-fold contrasts does not support repeated assessment sets",
    class = "dkge_fold_partition_error"
  )
  expect_error(
    .dkge_prepare_folds(fit, repeated),
    "DKGE classification does not support repeated assessment sets",
    class = "dkge_fold_partition_error"
  )
  repeated_reversed <- repeated
  repeated_reversed$assignments <- rev(repeated$assignments)
  expect_error(
    .dkge_normalize_folds(
      repeated_reversed,
      fit,
      consumer = "K-fold contrasts"
    ),
    "K-fold contrasts does not support repeated assessment sets",
    class = "dkge_fold_partition_error"
  )

  latent_fit <- fit
  latent_fit$U <- matrix(1, nrow = 1L, ncol = 1L)
  latent_features <- replicate(
    length(latent_fit$Btil), matrix(0, nrow = 1L, ncol = 1L),
    simplify = FALSE
  )
  expect_error(
    dkge_cv_train_latent_classifier(
      latent_fit,
      factor(c("a", "b", "a", "b")),
      Z_by_subject = latent_features,
      folds = repeated
    ),
    "DKGE latent classification does not support repeated assessment sets",
    class = "dkge_fold_partition_error"
  )

  partial <- dkge_define_folds(
    fit, type = "custom",
    assignments = list(c(1L, 2L), c(3L)),
    partition = "partial"
  )
  expect_equal(partial$metadata$coverage, 3L)
  expect_identical(partial$metadata$partition, "partial")

  expect_error(
    dkge_define_folds(
      fit, type = "custom",
      assignments = list(c(1L, 2L), c(2L, 3L)),
      partition = "partial"
    ),
    "nonoverlapping partition",
    class = "dkge_fold_partition_error"
  )
  expect_error(
    dkge_define_folds(
      fit, type = "custom",
      assignments = list(c(1L, 2L), c(2L, 3L)),
      partition = "repeated"
    ),
    "cover every subject",
    class = "dkge_fold_partition_error"
  )
})

test_that("fold construction restores the caller RNG state", {
  fit <- beta_fake_fit()
  set.seed(77)
  before <- .Random.seed
  invisible(dkge_define_folds(fit, type = "subject", k = 2, seed = 9))
  expect_identical(.Random.seed, before)

  old_seed <- .Random.seed
  on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  rm(".Random.seed", envir = .GlobalEnv)
  invisible(dkge_define_folds(fit, type = "subject", k = 2, seed = 9))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("shared scalar validators reject before coercion", {
  expect_error(
    dkge_inference_spec(B = 1.9),
    "`B`.*1.9.*strictly positive integer",
    class = "dkge_validation_error"
  )
  expect_error(
    kernel_roots(diag(2), jitter = c(1e-10, 1e-8)),
    "`jitter`.*length 2.*non-negative scalar",
    class = "dkge_validation_error"
  )
  expect_error(
    dkge_transport_spec(list(matrix(0, 1, 1)), max_iter = 2.5),
    "`max_iter`.*2.5.*strictly positive integer",
    class = "dkge_validation_error"
  )
  expect_error(
    dkge_infer(list(), 1, correction = "none", alpha = 1.1),
    "`alpha`.*1.1.*probability",
    class = "dkge_validation_error"
  )
})

test_that("Sinkhorn solves on positive support and re-expands zero masses", {
  skip_if_not(exists("sinkhorn_plan_cpp", envir = asNamespace("dkge"),
                     inherits = FALSE))
  C <- matrix(c(0, 1, 1, 0), 2)
  plan <- .dkge_sinkhorn_plan(
    C, mu = c(1, 0), nu = c(0.5, 0.5),
    epsilon = 0.1, max_iter = 100L, tol = 1e-8,
    warm_start = FALSE
  )
  expect_equal(rowSums(plan), c(1, 0), tolerance = 1e-7)
  expect_equal(colSums(plan), c(0.5, 0.5), tolerance = 1e-7)
  expect_equal(plan[2, ], c(0, 0), tolerance = 0)

  supported <- .dkge_sinkhorn_plan(
    matrix(c(0, 1, 2, 1, 0, 1), nrow = 2),
    mu = c(0.4, 0.6), nu = c(0.4, 0, 0.6),
    epsilon = 0.1, max_iter = 100L, tol = 1e-8,
    warm_start = FALSE, return_diagnostics = TRUE
  )
  expect_equal(colSums(supported$plan), c(0.4, 0, 0.6), tolerance = 1e-7)
  expect_equal(supported$plan[, 2], c(0, 0), tolerance = 0)
  expect_identical(supported$diagnostics$positive_column_support, c(1L, 3L))
})

test_that("native Sinkhorn validates warm-start lengths", {
  skip_if_not(exists("sinkhorn_plan_cpp", envir = asNamespace("dkge"),
                     inherits = FALSE))
  C <- matrix(c(0, 1, 1, 0), 2)
  expect_error(
    sinkhorn_plan_cpp(
      C, c(0.5, 0.5), c(0.5, 0.5),
      epsilon = 0.1, max_iter = 10L, tol = 1e-6,
      log_u_init = numeric(0), log_v_init = NULL,
      keep_duals = TRUE
    ),
    "log_u_init.*length 2"
  )
  expect_error(
    sinkhorn_plan_cpp(
      C, c(0.5, 0.5), c(0.5, 0.5),
      epsilon = 0.1, max_iter = 10L, tol = 1e-6,
      log_u_init = NULL, log_v_init = 0,
      keep_duals = TRUE
    ),
    "log_v_init.*length 2"
  )
})

test_that("zero-signal aggregate fits return honest rank zero", {
  target <- matrix(
    0, nrow = 2, ncol = 3,
    dimnames = list(c("r1", "r2"), paste0("f", 1:3))
  )
  fit <- dkge_aggregate_fit(target, K = diag(c(1, 0)), rank = 2)

  expect_equal(fit$rank, 0L)
  expect_equal(dim(fit$U), c(2L, 0L))
  expect_length(fit$singular_values, 0L)
  expect_equal(fit$kernel_rank, 1L)
  expect_equal(fit$moment_rank, 0L)
  expect_equal(fit$effective_rank, 0L)
})

test_that("between-subject inference requires an explicit qualified method", {
  expect_null(eval(formals(dkge_between_permute)$method))
  fake <- structure(list(), class = c("dkge_between_rrr", "list"))
  expect_error(
    dkge_between_permute(fake, terms = "x", B = 5),
    "method.*chosen explicitly",
    class = "dkge_inference_compatibility_error"
  )
})
