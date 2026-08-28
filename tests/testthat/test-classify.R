skip_on_cran()

make_classification_fit <- function(S = 3, P = NULL, seed = 1) {
  set.seed(seed)
  factors <- list(A = list(L = 2), B = list(L = 2), time = list(L = 4))
  dk <- design_kernel(factors, basis = "effect")
  q <- nrow(dk$K)
  P <- P %||% (q + 5)  # ensure P >= q to avoid rank warnings
  betas <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
  designs <- replicate(S, diag(q), simplify = FALSE)
  fit <- dkge_fit(betas, designs, K = dk, rank = min(3, q))
  list(fit = fit, kernel = dk)
}

test_that("dkge_targets constructs effect targets", {
  fixture <- make_classification_fit()
  fit <- fixture$fit
  tg <- dkge_targets(fit, ~ A + B + A:B,
                     collapse = list(time = list(method = "mean", window = 2:3)))
  expect_gt(length(tg), 0)
  expect_true(all(vapply(tg, inherits, logical(1), "dkge_target")))
  expect_true(all(vapply(tg, function(x) nrow(x$weight_matrix) >= 1, logical(1))))
})

test_that("dkge_targets rematches both named map axes and rejects ambiguity", {
  fixture <- make_classification_fit(S = 4, seed = 19)
  fit <- fixture$fit
  reference <- dkge_targets(fit, ~ A + B + A:B, residualize = FALSE)

  permuted <- fit
  map <- fit$kernel_info$map
  row_perm <- rev(seq_len(nrow(map)))
  col_perm <- c(seq_len(ncol(map))[-1L], 1L)
  permuted$kernel_info$map <- map[row_perm, col_perm, drop = FALSE]
  observed <- dkge_targets(permuted, ~ A + B + A:B, residualize = FALSE)

  expect_identical(vapply(observed, `[[`, character(1), "name"),
                   vapply(reference, `[[`, character(1), "name"))
  for (i in seq_along(reference)) {
    expect_equal(observed[[i]]$weight_matrix,
                 reference[[i]]$weight_matrix,
                 tolerance = 1e-12)
  }

  unnamed <- fit
  dimnames(unnamed$kernel_info$map) <- NULL
  expect_error(
    dkge_targets(unnamed, ~ A),
    "map.*named cell rows and effect columns"
  )

  duplicated <- fit
  colnames(duplicated$kernel_info$map)[2L] <-
    colnames(duplicated$kernel_info$map)[1L]
  expect_error(
    dkge_targets(duplicated, ~ A),
    "map column names.*unique"
  )
})

test_that("dkge_classify returns metrics", {
  fixture <- make_classification_fit()
  fit <- fixture$fit
  cls <- dkge_classify(fit, targets = ~ A + B, n_perm = 0, seed = 11)
  expect_s3_class(cls, "dkge_classification")
  df <- as.data.frame(cls)
  expect_s3_class(df, "data.frame")
  expect_true(all(df$metric %in% cls$metric))
  expect_true(all(df$n_perm == 0))
})

test_that("dkge_classify supports logit backend", {
  fixture <- make_classification_fit()
  fit <- fixture$fit
  cls <- dkge_classify(fit, targets = ~ A + B, method = "logit", n_perm = 0)
  expect_s3_class(cls, "dkge_classification")
  df <- as.data.frame(cls)
  expect_true(all(df$metric %in% cls$metric))
})

test_that("dkge_classify lambda control grid works", {
  fixture <- make_classification_fit()
  fit <- fixture$fit
  grid <- c(1e-6, 1e-3, 1e-1)
  cls <- dkge_classify(fit, targets = ~ A, n_perm = 0, control = list(lambda_grid = grid))
  expect_true(cls$results[["A"]]$lambda %in% grid)
})

test_that("dkge_classify lambda control function works", {
  fixture <- make_classification_fit()
  fit <- fixture$fit
  lambda_fun <- function(target, fold, method, default) 5e-2
  expect_silent(
    dkge_classify(fit, targets = ~ A, n_perm = 0, control = list(lambda_fun = lambda_fun))
  )
})

test_that("classification permutations require fixed selection and recomputation", {
  expect_error(
    dkge:::.dkge_validate_classification_selection(
      n_perm = 9L, lambda = NULL, lambda_grid = NULL, lambda_fun = NULL
    ),
    "preselected scalar `lambda`",
    class = "dkge_classification_inference_error"
  )
  expect_error(
    dkge:::.dkge_validate_classification_selection(
      n_perm = 9L, lambda = 0.1,
      lambda_grid = c(0.1, 1), lambda_fun = NULL
    ),
    "preselected scalar `lambda`",
    class = "dkge_classification_inference_error"
  )

  fixture <- make_classification_fit()
  expect_error(
    dkge_classify(
      fixture$fit, targets = ~ A, n_perm = 2L, lambda = 0.1
    ),
    "randomization_recompute.*complete.*representation",
    class = "dkge_classification_inference_error"
  )

  calls <- 0L
  recompute <- function(metric, ...) {
    calls <<- calls + 1L
    c(extra_diagnostic = 99, brier = 0.25, accuracy = 0.75)
  }
  classified <- dkge_classify(
    fixture$fit, targets = ~ A, metric = c("accuracy", "brier"),
    n_perm = 2L, lambda = 0.1,
    control = list(randomization_recompute = recompute), seed = 12
  )
  expect_equal(calls, 2L)
  expect_equal(classified$results[[1]]$permutations[, "accuracy"],
               rep(0.75, 2L))
  expect_equal(classified$results[[1]]$permutations[, "brier"],
               rep(0.25, 2L))
  expect_identical(classified$lambda_selection, "preselected_external")
  expect_identical(classified$randomization_exactness,
                   "user_supplied_pipeline_recompute")
})

test_that("cell-cross classification fails closed when a fold loses rank", {
  effects <- c("e1", "e2")
  B1 <- matrix(c(1, 0, 0, 0), nrow = 2,
               dimnames = list(effects, c("p1", "p2")))
  B2 <- matrix(c(2, 0, 0, 0), nrow = 2,
               dimnames = list(effects, c("p1", "p2")))
  B3 <- matrix(c(0, 1, 0, 0), nrow = 2,
               dimnames = list(effects, c("p1", "p2")))
  X <- diag(2)
  colnames(X) <- effects
  K <- diag(2)
  dimnames(K) <- list(effects, effects)
  fit <- suppressWarnings(dkge_fit(
    list(s1 = B1, s2 = B2, s3 = B3),
    list(s1 = X, s2 = X, s3 = X),
    K = K, rank = 2L, w_method = "none", effect_scaling = "none"
  ))
  target <- diag(2)
  dimnames(target) <- list(c("class1", "class2"), effects)
  folds <- dkge_define_folds(
    fit, type = "custom", assignments = list(1L, 2L, 3L)
  )

  expect_error(
    dkge_classify(
      fit, targets = target, mode = "cell_cross", folds = folds,
      n_perm = 0L
    ),
    "Training fold 3 has effective rank 1, below fitted rank 2",
    class = "dkge_fold_rank_error"
  )
})

test_that("dkge_classify delta mode handles rank-1 targets", {
  fixture <- make_classification_fit(S = 5)  # S > rank to avoid singular covariance
  fit <- fixture$fit
  w <- matrix(c(1, -1, rep(0, nrow(fit$U) - 2)), nrow = 1)
  target <- list(
    name = "contrast",
    factors = character(0),
    labels = data.frame(),
    class_labels = c("pos", "neg"),
    weight_matrix = rbind(w, -w),
    indicator = NULL,
    residualized = FALSE,
    collapse = NULL,
    scope = "within_subject"
  )
  class(target) <- c("dkge_target", "list")
  cls <- dkge_classify(
    fit, targets = list(target), mode = "delta", lambda = 1e-3,
    n_perm = 10, scope = "signflip", seed = 5
  )
  expect_s3_class(cls, "dkge_classification")
  expect_true(all(names(cls$results[[1]]$metrics) == cls$metric))
})

test_that("delta mode computes label-aware metrics", {
  fixture <- make_classification_fit(S = 4)
  fit <- fixture$fit
  w <- matrix(c(1, -1, rep(0, nrow(fit$U) - 2)), nrow = 1)
  target <- list(
    name = "delta_target",
    factors = character(0),
    labels = data.frame(),
    class_labels = c("pos", "neg"),
    weight_matrix = rbind(w, -w),
    indicator = NULL,
    residualized = FALSE,
    collapse = NULL,
    scope = "within_subject"
  )
  class(target) <- c("dkge_target", "list")
  y <- factor(rep(c("pos", "neg"), length.out = length(fit$Btil)), levels = c("pos", "neg"))
  cls <- dkge_classify(fit,
                       targets = list(target),
                       mode = "delta",
                       metric = c("accuracy", "logloss", "brier", "auroc", "ece"),
                       y = y,
                       scope = "signflip",
                       n_perm = 0)
  res <- cls$results[[1]]
  expect_false(any(is.na(res$metrics)))
  expect_equal(res$positive_class, "pos")
  expect_equal(unname(res$subject_labels), y)
  expected_ids <- fit$subject_ids %||% paste0("subject", seq_along(y))
  expect_equal(names(res$subject_labels), expected_ids)
})

test_that("delta mode errors on mismatched labels", {
  fixture <- make_classification_fit(S = 3)
  fit <- fixture$fit
  w <- matrix(c(1, -1, rep(0, nrow(fit$U) - 2)), nrow = 1)
  target <- list(
    name = "delta_target",
    factors = character(0),
    labels = data.frame(),
    class_labels = c("pos", "neg"),
    weight_matrix = rbind(w, -w),
    indicator = NULL,
    residualized = FALSE,
    collapse = NULL,
    scope = "within_subject"
  )
  class(target) <- c("dkge_target", "list")
  y_bad <- c("pos", "neg")
  expect_error(
    dkge_classify(fit, targets = list(target), mode = "delta", y = y_bad),
    "subject labels must have length"
  )
  y_unknown <- c("pos", "pos", "maybe")
  expect_error(
    dkge_classify(fit, targets = list(target), mode = "delta", y = y_unknown),
    "unknown levels"
  )
})

test_that("dkge_confusion aggregates per-target matrices", {
  fixture <- make_classification_fit()
  fit <- fixture$fit
  cls <- dkge_classify(fit, targets = ~ A + B, n_perm = 0)
  conf_all <- dkge_confusion(cls)
  conf_matrix <- if (is.list(conf_all)) conf_all[[1]] else conf_all
  expect_true(is.matrix(conf_matrix))
  expect_equal(nrow(conf_matrix), length(cls$results[[1]]$target$class_labels))
  conf_fold <- dkge_confusion(cls, fold = 1)
  conf_matrix_fold <- if (is.list(conf_fold)) conf_fold[[1]] else conf_fold
  expect_true(is.matrix(conf_matrix_fold))
})

test_that("classification diagnostics expose per-fold data frames", {
  fixture <- make_classification_fit()
  fit <- fixture$fit
  cls <- dkge_classify(fit, targets = ~ A, n_perm = 0)
  fold_counts <- as.data.frame(cls, what = "fold_counts")
  expect_true(all(c("target", "fold", "class", "train", "test") %in% names(fold_counts)))
  lambda_df <- as.data.frame(cls, what = "lambda")
  expect_true(all(c("target", "fold", "lambda") %in% names(lambda_df)))
})

test_that("logit predictions renormalize underflowed rows", {
  model <- list(
    type = "logit",
    classes = c("case", "control"),
    beta = matrix(-1e6, nrow = 3, ncol = 2)
  )
  X <- matrix(0, nrow = 2, ncol = 2)
  probs <- .dkge_predict_logit(model, X, class_levels = c("case", "control"))
  expect_equal(rowSums(probs), rep(1, nrow(X)))
  expect_true(all(probs >= 0 & probs <= 1))
})

test_that("dkge_pipeline integrates classification", {
  fixture <- make_classification_fit()
  fit <- fixture$fit
  q <- nrow(fit$U)
  contrast <- rep(0, q)
  contrast[1] <- 1
  pipeline <- dkge_pipeline(fit = fit,
                            contrasts = contrast,
                            classification = list(targets = ~ A, n_perm = 0, seed = 2),
                            inference = NULL)
  expect_true("classification" %in% names(pipeline))
  expect_s3_class(pipeline$classification, "dkge_classification")
})
