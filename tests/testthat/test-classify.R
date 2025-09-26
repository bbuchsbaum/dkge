skip_on_cran()

make_classification_fit <- function(S = 3, P = 4, seed = 1) {
  set.seed(seed)
  factors <- list(A = list(L = 2), B = list(L = 2), time = list(L = 4))
  dk <- design_kernel(factors, basis = "effect")
  q <- nrow(dk$K)
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

test_that("dkge_classify returns metrics", {
  fixture <- make_classification_fit()
  fit <- fixture$fit
  cls <- dkge_classify(fit, targets = ~ A + B, n_perm = 5, seed = 11)
  expect_s3_class(cls, "dkge_classification")
  df <- as.data.frame(cls)
  expect_s3_class(df, "data.frame")
  expect_true(all(df$metric %in% cls$metric))
  expect_true(all(df$n_perm == 5))
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

test_that("dkge_classify delta mode handles rank-1 targets", {
  fixture <- make_classification_fit()
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
  cls <- dkge_classify(fit, targets = list(target), mode = "delta", n_perm = 10,
                       scope = "signflip", seed = 5)
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
                            classification = list(targets = ~ A, n_perm = 3, seed = 2),
                            inference = NULL)
  expect_true("classification" %in% names(pipeline))
  expect_s3_class(pipeline$classification, "dkge_classification")
})
