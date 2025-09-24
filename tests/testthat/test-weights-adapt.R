test_that("adaptive weights: kenergy, precision, kenergy_prec", {
  Q <- 6; V <- 5
  set.seed(42)
  B1 <- matrix(rnorm(Q * V), Q, V)
  B2 <- matrix(rnorm(Q * V), Q, V)
  K  <- diag(Q)

  raw_ke <- (colSums(B1^2) + colSums(B2^2)) / 2
  w_ke <- dkge:::.dkge_adapt_weights(
    B_list = list(B1, B2), adapt = "kenergy", K = K,
    winsor = 0.9999
  )
  exp_ke <- dkge:::.dkge_winsor(raw_ke, upper = 0.9999)
  expect_equal(w_ke, exp_ke, tolerance = 1e-6)

  raw_pr <- 0.5 * (1 / (colMeans(B1^2) + 1e-8) + 1 / (colMeans(B2^2) + 1e-8))
  w_pr <- dkge:::.dkge_adapt_weights(
    B_list = list(B1, B2), adapt = "precision", K = NULL,
    winsor = 0.9999
  )
  exp_pr <- dkge:::.dkge_winsor(raw_pr, upper = 0.9999)
  expect_equal(w_pr, exp_pr, tolerance = 1e-6)

  w_kp <- dkge:::.dkge_adapt_weights(
    B_list = list(B1, B2), adapt = "kenergy_prec", K = K,
    winsor = 0.9999
  )
  exp_kp <- dkge:::.dkge_winsor(raw_ke * raw_pr, upper = 0.9999)
  expect_equal(w_kp, exp_kp, tolerance = 1e-6)
})

set.seed(2024)

test_that("adaptive weighting homes in on informative voxels across rules", {
  S <- 16L; P <- 10L; q <- 2L
  informative <- 1:3
  classes <- rep(c("A", "B"), each = S / 2)

  subjects <- vector("list", S)
  B_mats <- vector("list", S)
  for (s in seq_len(S)) {
    X <- diag(q)
    B <- matrix(rnorm(q * P, sd = 0.5), nrow = q)
    if (classes[s] == "A") {
      B[1, informative] <- B[1, informative] + 4.0
    } else {
      B[1, informative] <- B[1, informative] - 4.0
    }
    B_mats[[s]] <- B
    subjects[[s]] <- dkge_subject(B, design = X, id = paste0("sub", s))
  }

  K <- diag(q)
  fit_uniform <- dkge(subjects, kernel = K, rank = 2, w_method = "none")
  y <- factor(classes, levels = c("A", "B"))

  predict_acc <- function(fit, clf, y) {
    Z_list <- dkge_project_clusters_to_latent(fit)
    Z_bar <- do.call(rbind, lapply(Z_list, colMeans))
    levels_y <- levels(y)
    preds <- character(length(y))
    for (s in seq_along(y)) {
      fold <- clf$fold_assignment[s]
      beta <- clf$beta_by_subject[[s]]
      intercept <- clf$models_by_fold[[fold]]$intercept
      if (is.null(intercept)) intercept <- 0
      score <- sum(beta * Z_bar[s, ]) + intercept
      preds[s] <- if (score >= 0) levels_y[1] else levels_y[2]
    }
    mean(preds == y)
  }

  clf_uniform <- dkge_cv_train_latent_classifier(fit_uniform, y, folds = S,
                                                 model = "lda", level = "subject",
                                                 standardize = FALSE)
  acc_uniform <- predict_acc(fit_uniform, clf_uniform, y)

  for (rule in c("kenergy", "precision", "kenergy_prec")) {
    weights_rule <- dkge:::.dkge_adapt_weights(
      B_list = B_mats,
      adapt = rule,
      K = if (rule == "precision") NULL else K,
      winsor = 0.999
    )
    if (rule == "kenergy") {
      expect_true(mean(weights_rule[informative]) >
                    mean(weights_rule[-informative]) * 5,
                  info = paste("rule", rule))
    }

    w_spec <- dkge_weights(adapt = rule, scope = "fold", mix = 1,
                           shrink = list(alpha = 1, winsor = 0.999, normalize = "mean"))
    fit_adapt <- dkge(subjects, kernel = K, rank = 2, weights = w_spec)
    clf_adapt <- dkge_cv_train_latent_classifier(fit_adapt, y, folds = S,
                                                 model = "lda", level = "subject",
                                                 standardize = FALSE)
    acc_adapt <- predict_acc(fit_adapt, clf_adapt, y)
    expect_true(acc_adapt >= acc_uniform, info = paste("rule", rule))
  }
})
