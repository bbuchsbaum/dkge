library(testthat)

make_advanced_qspace_fit <- function(S = 4L, P = 5L, seed = 7401) {
  set.seed(seed)
  effects <- c("e1", "e2", "e3")
  K <- matrix(c(
    1.4, 0.3, 0.1,
    0.3, 1.2, 0.2,
    0.1, 0.2, 1.1
  ), 3, 3, dimnames = list(effects, effects))
  subjects <- lapply(seq_len(S), function(s) {
    B <- matrix(rnorm(length(effects) * P), length(effects), P,
                dimnames = list(effects, paste0("v", seq_len(P))))
    X <- matrix(c(
      1, 0, 1,
      1, 1, 0,
      0, 1, 1,
      2, 1, 0,
      0, 2, 1,
      1, 0, 2
    ), 6, 3, byrow = TRUE)
    X[, 1] <- X[, 1] + 0.05 * (s - 1L)
    colnames(X) <- effects
    counts <- c(e1 = 2 + s, e2 = 11 - s, e3 = 3 + 2 * s)
    suppressWarnings(dkge_subject(
      B, design = X, id = paste0("s", s), effect_n = counts
    ))
  })
  data <- dkge_data(subjects)
  fit <- dkge_fit(
    data, K = K, rank = 2, w_method = "none",
    effect_weights = dkge_effect_weights("count")
  )
  list(fit = fit, data = data, subjects = subjects, K = K)
}

test_that("audit mismatch probes no longer advertise a false block factor", {
  effects <- c("e1", "e2")
  make_count_subject <- function(B, counts, id) {
    rownames(B) <- effects
    X <- diag(2)
    colnames(X) <- effects
    suppressWarnings(dkge_subject(B, X, id = id, effect_n = counts))
  }
  count_subjects <- list(
    make_count_subject(matrix(c(1, 2, 3, 4), 2, 2), c(e1 = 1, e2 = 4), "s1"),
    make_count_subject(matrix(c(4, 1, 2, 3), 2, 2), c(e1 = 9, e2 = 1), "s2")
  )
  count_fit <- dkge_fit(
    dkge_data(count_subjects), K = diag(2), rank = 1,
    w_method = "none", effect_scaling = "none",
    effect_weights = dkge_effect_weights("count")
  )
  naive_count <- do.call(cbind, lapply(count_fit$Braw, function(B) {
    count_fit$Khalf %*% t(count_fit$R) %*% B
  }))
  count_error <- norm(count_fit$Chat - tcrossprod(naive_count), "F") /
    norm(count_fit$Chat, "F")

  B1 <- matrix(c(1, 2, 3, 4), 2, 2,
               dimnames = list(c("e1", "e2"), c("v1", "v2")))
  B2 <- matrix(c(5, 6, 7, 8), 2, 2,
               dimnames = list(c("e2", "e3"), c("v1", "v2")))
  X1 <- diag(2)
  X2 <- diag(2)
  colnames(X1) <- c("e1", "e2")
  colnames(X2) <- c("e2", "e3")
  partial_data <- suppressWarnings(dkge_data(
    list(B1, B2), designs = list(X1, X2), subject_ids = c("p1", "p2")
  ))
  partial_fit <- dkge_fit(
    partial_data, K = diag(3), rank = 2, w_method = "none",
    effect_scaling = "none", missingness = "rescale"
  )
  naive_partial <- do.call(cbind, partial_fit$Braw)
  partial_error <- norm(partial_fit$Chat - tcrossprod(naive_partial), "F") /
    norm(partial_fit$Chat, "F")

  expect_gt(count_error, 0.1)
  expect_gt(partial_error, 0.2)
  for (fit in list(count_fit, partial_fit)) {
    expect_s3_class(fit, "dkge_qspace")
    expect_equal(fit$representation, "qspace_moment")
    expect_null(fit$X_concat)
    expect_null(fit$v)
  }
})

test_that("q-space loadings and contrasts follow an independent raw-beta oracle", {
  fixture <- make_advanced_qspace_fit()
  fit <- fixture$fit
  contrast <- setNames(c(1, -0.5, 0.25), fit$effects)

  manual_A <- lapply(fit$Braw, function(B) {
    t(B) %*% fit$R %*% fit$K %*% fit$U
  })
  predicted_A <- dkge_predict_loadings(fit, fit$Braw)
  expect_equal(predicted_A, manual_A, tolerance = 1e-11,
               ignore_attr = TRUE)
  expect_equal(dkge_project_clusters_to_latent(fit), manual_A,
               tolerance = 1e-11, ignore_attr = TRUE)
  expect_equal(dkge_cluster_loadings(fit), manual_A,
               tolerance = 1e-11, ignore_attr = TRUE)

  ctil <- backsolve(fit$R, contrast, transpose = FALSE)
  alpha <- t(fit$U) %*% fit$K %*% ctil
  manual_values <- lapply(manual_A, function(A) as.numeric(A %*% alpha))
  predicted <- dkge_predict(fit, fit$Braw, list(effect = contrast))
  expect_equal(unname(predicted$values), manual_values, tolerance = 1e-11,
               ignore_attr = TRUE)
})

test_that("q-space fits support audited downstream consumers or fail closed", {
  skip_if_not_installed("ggplot2")
  fixture <- make_advanced_qspace_fit(S = 6L, P = 6L)
  fit <- fixture$fit
  contrast <- setNames(c(1, -1, 0), fit$effects)

  expect_s3_class(dkge_plot_scree(fit), "ggplot")
  expect_s3_class(dkge_plot_effect_loadings(fit, comps = 1:2), "ggplot")
  contribution_plot <- dkge_plot_subject_contrib(fit, comps = 1:2)
  expect_s3_class(contribution_plot$weights, "ggplot")
  expect_s3_class(contribution_plot$energy, "ggplot")

  identity_cache <- list(
    operators = replicate(length(fit$Braw), diag(ncol(fit$Braw[[1]])),
                          simplify = FALSE)
  )
  boot <- dkge_bootstrap_qspace(
    fit, contrasts = contrast, B = 3, seed = 7402,
    transport_cache = identity_cache, medoid = 1, align = FALSE
  )
  expect_equal(boot$B, 3)
  expect_true(all(is.finite(boot$summary[[1]]$boot)))

  expect_error(
    multivarious::project(fit, diag(nrow(fit$K))),
    "requires representation='block_biprojector'"
  )
  expect_error(
    multivarious::project_block(
      fit, new_data = diag(nrow(fit$K)), block = 1,
      least_squares = TRUE
    ),
    "requires representation='block_biprojector'"
  )
  expect_error(
    multivarious::transfer(
      fit, new_data = diag(nrow(fit$K)), from = 1, to = 2
    ),
    "requires representation='block_biprojector'"
  )
  expect_error(
    dkge_project_clusters(fit, fit$Braw[[1]]),
    "requires representation='block_biprojector'"
  )
})

test_that("cell classification consumes q-space representations", {
  set.seed(7411)
  kernel <- design_kernel(
    list(A = list(L = 2), B = list(L = 2)), basis = "effect"
  )
  q <- nrow(kernel$K)
  effects <- rownames(kernel$K)
  subjects <- lapply(seq_len(6), function(s) {
    B <- matrix(rnorm(q * (q + 4)), q, q + 4,
                dimnames = list(effects, NULL))
    X <- diag(q)
    colnames(X) <- effects
    counts <- setNames(seq_len(q) + s, effects)
    suppressWarnings(dkge_subject(
      B, X, id = paste0("c", s), effect_n = counts
    ))
  })
  fit <- dkge_fit(
    dkge_data(subjects), K = kernel, rank = min(3, q),
    w_method = "none", effect_weights = dkge_effect_weights("count")
  )
  classified <- dkge_classify(
    fit, targets = ~ A, method = "lda", mode = "cell",
    n_perm = 0, verbose = FALSE
  )

  expect_s3_class(fit, "dkge_qspace")
  expect_s3_class(classified, "dkge_classification")
  expect_true(all(is.finite(unlist(classified$results[[1]]$metrics))))
})
