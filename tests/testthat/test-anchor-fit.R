set.seed(1001)

test_that("dkge_fit_from_kernels reproduces pooled kernel", {
  q <- 4
  S <- 3
  K_list <- replicate(S, {
    X <- matrix(rnorm(q * 6), q)
    K <- tcrossprod(X)
    K <- K + 1e-4 * diag(q)
    K / max(1, sum(diag(K)) / q)
  }, simplify = FALSE)

  fit <- dkge_fit_from_kernels(K_list, effect_ids = paste0("e", seq_len(q)))

  expect_s3_class(fit, "dkge")
  expect_equal(nrow(fit$U), q)
  expect_equal(length(fit$subject_ids), S)
  weighted_sum <- Reduce(`+`, Map(function(K, w) w * K, K_list, fit$weights))
  expect_equal(fit$Chat, weighted_sum, tolerance = 1e-6)
  for (s in seq_len(S)) {
    expect_equal(fit$contribs[[s]], K_list[[s]], tolerance = 1e-8)
  }
  expect_equal(fit$provenance$kernel_fit$sqrt_scale, sqrt(S))
})

test_that("dkpp_select_anchors returns valid indices", {
  X <- matrix(rnorm(80 * 5), 80, 5)
  sel <- dkpp_select_anchors(X, L = 16, rho = 0.5, seed = 11)
  expect_equal(length(sel$indices), 16)
  expect_true(all(sel$indices >= 1 & sel$indices <= nrow(X)))
  expect_true(is.finite(sel$sigma) && sel$sigma > 0)
})

test_that("dkge_build_anchor_kernels outputs PSD aligned kernels", {
  features_list <- list(
    s1 = matrix(rnorm(20 * 6), 20, 6),
    s2 = matrix(rnorm(25 * 6), 25, 6)
  )
  K_item_list <- lapply(features_list, function(Fs) {
    Y <- matrix(rnorm(nrow(Fs) * 5), nrow(Fs), 5)
    K <- Y %*% t(Y)
    K + 1e-4 * diag(nrow(K))
  })

  res <- dkge_build_anchor_kernels(features_list, K_item_list, L = 12, method = "random", seed = 3)
  expect_equal(length(res), 1)
  ctx <- res[[1]]
  expect_equal(nrow(ctx$anchors), 12)
  expect_equal(length(ctx$K_aligned), length(features_list))

  for (K in ctx$K_aligned) {
    expect_equal(dim(K), c(12, 12))
    eig <- eigen((K + t(K)) / 2, symmetric = TRUE, only.values = TRUE)$values
    expect_gt(min(eig), -1e-6)
    expect_equal(sum(diag(K)), 12, tolerance = 1e-4)
  }
})

test_that("dkge_anchor_fit stores provenance and diagnostics", {
  features_list <- list(
    s1 = matrix(rnorm(18 * 5), 18, 5),
    s2 = matrix(rnorm(22 * 5), 22, 5),
    s3 = matrix(rnorm(16 * 5), 16, 5)
  )
  K_item_list <- lapply(features_list, function(Fs) {
    Y <- matrix(rnorm(nrow(Fs) * 4), nrow(Fs), 4)
    K <- Y %*% t(Y)
    K + 1e-3 * diag(nrow(K))
  })

  fit <- dkge_anchor_fit(
    features_list,
    K_item_list,
    anchors = list(L = 10, method = "dkpp", seed = 7)
  )

  expect_s3_class(fit, "dkge")
  info <- fit$provenance$anchors
  expect_true(is.list(info))
  expect_equal(info$L, 10)
  expect_equal(nrow(fit$U), 10)
  expect_true(is.data.frame(info$coverage))
  expect_equal(ncol(info$coverage), 4)

  diag <- dkge_anchor_diagnostics(fit)
  expect_true(all(c("summary", "coverage", "leverage") %in% names(diag)))
  expect_equal(nrow(diag$leverage), 10)

  contrast_dir <- dkge_anchor_contrast_from_direction(info$anchors, rep(1, ncol(info$anchors)))
  expect_equal(length(contrast_dir), 10)
  expect_equal(sqrt(sum(contrast_dir^2)), 1, tolerance = 1e-6)

  contrast_proto <- dkge_anchor_contrast_from_prototypes(info$anchors, positives = info$anchors[1:2, , drop = FALSE])
  expect_equal(length(contrast_proto), 10)
  expect_equal(sqrt(sum(contrast_proto^2)), 1, tolerance = 1e-6)
})

test_that("dkge_pipeline accepts anchor input", {
  features_list <- list(
    s1 = matrix(rnorm(15 * 4), 15, 4),
    s2 = matrix(rnorm(18 * 4), 18, 4)
  )
  K_item_list <- lapply(features_list, function(Fs) {
    Y <- matrix(rnorm(nrow(Fs) * 3), nrow(Fs), 3)
    K <- Y %*% t(Y)
    K + 1e-4 * diag(nrow(K))
  })
  anchor_spec <- dkge_input_anchor(features_list, K_item_list,
                                   anchors = list(L = 8, seed = 5L))
  contrasts <- rep(1 / sqrt(8), 8)
  pipeline_res <- dkge_pipeline(input = anchor_spec,
                                contrasts = contrasts,
                                method = "analytic",
                                inference = NULL)
  expect_s3_class(pipeline_res$fit, "dkge")
  expect_equal(pipeline_res$fit$provenance$anchors$L, 8)
  expect_true(is.list(pipeline_res$contrasts))
})

test_that("anchor fit handles fully shared item sets", {
  set.seed(42)
  n_items <- 20
  d <- 6
  features <- matrix(rnorm(n_items * d), n_items, d)
  K_item <- crossprod(matrix(rnorm(n_items * 5), 5, n_items))
  K_item <- K_item + 1e-4 * diag(n_items)

  features_list <- list(a = features, b = features)
  K_item_list <- list(a = K_item, b = K_item)

  anchor_input <- dkge_input_anchor(features_list, K_item_list,
                                    anchors = list(L = 12, seed = 9L))
  fit_anchor <- dkge_fit_from_input(anchor_input)
  # Subjects share identical projections and coverage collapses to zero
  expect_equal(fit_anchor$contribs[["a"]], fit_anchor$contribs[["b"]], tolerance = 1e-6)
  coverage <- fit_anchor$provenance$anchors$coverage
  expect_true(all(abs(coverage$p50) < 1e-6))
  expect_true(length(unique(coverage$p90)) == 1)
  expect_true(length(unique(coverage$p95)) == 1)
  # Leverage matches unit trace normalisation
  lev <- fit_anchor$provenance$anchors$leverage$per_anchor
  expect_equal(sum(lev), fit_anchor$provenance$anchors$L, tolerance = 1e-6)
})

test_that("dkge_build_anchor_kernels accepts dkge_folds assignments", {
  set.seed(123)
  S <- 4L
  subject_ids <- paste0("sub", seq_len(S))
  features_list <- lapply(seq_len(S), function(s) {
    matrix(rnorm(15), nrow = 5, ncol = 3)
  })
  names(features_list) <- subject_ids
  K_item_list <- lapply(features_list, function(Fs) {
    Y <- matrix(rnorm(nrow(Fs) * 3), nrow(Fs), 3)
    tcrossprod(Y) + 1e-6 * diag(nrow(Fs))
  })

  fit_stub <- structure(list(subject_ids = subject_ids), class = "dkge")
  folds <- as_dkge_folds(list(foldA = c(1L, 3L), foldB = c(2L, 4L)), fit_or_data = fit_stub)

  res <- dkge_build_anchor_kernels(features_list,
                                   K_item_list,
                                   folds = folds,
                                   L = 4,
                                   method = "random",
                                   seed = 5L,
                                   center = FALSE,
                                   whiten = FALSE,
                                   unit_trace = FALSE)

  expect_equal(sort(names(res)), c("foldA", "foldB"))
  expect_equal(res$foldA$test_idx, c(1L, 3L))
  expect_equal(sort(res$foldA$train_idx), c(2L, 4L))
  expect_equal(res$foldB$test_idx, c(2L, 4L))
  expect_equal(sort(res$foldB$train_idx), c(1L, 3L))
  for (ctx in res) {
    expect_equal(length(ctx$train_idx) + length(ctx$test_idx), S)
  }
})

test_that("anchor fit preserves PSD with subject subgroups sharing items", {
  set.seed(123)
  d <- 5
  features_shared <- matrix(rnorm(30 * d), 30, d)
  K_shared <- crossprod(matrix(rnorm(30 * 4), 4, 30)) + 1e-4 * diag(30)

  features_unique <- matrix(rnorm(25 * d), 25, d)
  K_unique <- crossprod(matrix(rnorm(25 * 4), 4, 25)) + 1e-4 * diag(25)

  features_list <- list(
    s1 = features_shared,
    s2 = features_shared,
    s3 = features_unique
  )
  K_item_list <- list(
    s1 = K_shared,
    s2 = K_shared,
    s3 = K_unique
  )

  anchor_input <- dkge_input_anchor(features_list, K_item_list,
                                    anchors = list(L = 14, seed = 15L))
  fit_anchor <- dkge_fit_from_input(anchor_input)

  # Subjects that share items should yield identical projections
  expect_equal(fit_anchor$contribs[["s1"]], fit_anchor$contribs[["s2"]], tolerance = 1e-6)

  # The aligned kernels remain PSD and unit-trace by construction
  eig_values <- eigen(fit_anchor$Chat, symmetric = TRUE, only.values = TRUE)$values
  expect_gt(min(eig_values), -1e-6)
  expect_equal(sum(diag(fit_anchor$K)), nrow(fit_anchor$K), tolerance = 1e-6)
})

test_that("anchor classification requires explicit weight targets", {
  set.seed(10)
  features_list <- list(
    a = matrix(rnorm(12), 4, 3),
    b = matrix(rnorm(12), 4, 3)
  )
  K_item_list <- lapply(features_list, function(Fs) {
    Y <- matrix(rnorm(nrow(Fs) * 3), nrow(Fs), 3)
    tcrossprod(Y) + 1e-4 * diag(nrow(Fs))
  })
  fit_anchor <- suppressWarnings(
    dkge_anchor_fit(features_list, K_item_list, anchors = list(L = 4, seed = 3L))
  )
  expect_error(
    dkge_classify(fit_anchor, targets = ~ 1),
    "fit\\$kernel_info\\$map"
  )
})
