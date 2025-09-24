# test-info-standardize.R

library(testthat)


test_that("Classifier standardization is fold-wise (training-only)", {
  skip_on_cran()
  set.seed(123)

  # Minimal synthetic fit (reuse small sizes for speed)
  S <- 12L; q <- 4L; r <- 2L
  U <- qr.Q(qr(matrix(rnorm(q * r), q, r)))
  K <- diag(q); Khalf <- diag(q); Kihalf <- diag(q); R <- diag(q)
  P_list <- rep(60L, S)
  y <- factor(rep(c("A","B"), length.out = S), levels = c("A","B"))

  # Latent features A_s; small mean shift on dim 1
  A_list <- lapply(seq_len(S), function(s) {
    mu <- if (y[s] == "A") c(+0.4, 0) else c(-0.4, 0)
    matrix(rep(mu, each = P_list[s]), nrow = P_list[s], byrow = FALSE) +
      matrix(rnorm(P_list[s] * r, sd = 0.4), P_list[s], r)
  })
  Btil <- lapply(A_list, function(A) U %*% t(A))
  fit <- structure(list(Btil = Btil, U = U, K = K, Khalf = Khalf, Kihalf = Kihalf, R = R), class = "dkge")

  # Mock cluster->anchor pipeline (not used directly here)
  testthat::local_mocked_bindings(
    dkge_cluster_loadings = function(fit) {
      lapply(fit$Btil, function(Bs) t(Bs) %*% fit$K %*% fit$U)
    },
    .package = "dkge"
  )

  # Train with 3 folds (deterministic)
  set.seed(999)
  clf <- dkge_cv_train_latent_classifier(fit, y, folds = 3, model = "lda", level = "subject", standardize = TRUE)

  # Reconstruct per-subject mean latent vector (S×r)
  Z_bar <- do.call(rbind, lapply(dkge_project_clusters_to_latent(fit), colMeans))

  # Pick fold 1; compute the training subjects' mean/sd; compare to stored
  k <- 1L
  hold  <- clf$folds$assignments[[k]]
  train <- setdiff(seq_len(S), hold)

  mu_expected <- colMeans(Z_bar[train, , drop = FALSE])
  sd_expected <- pmax(apply(Z_bar[train, , drop = FALSE], 2, sd), 1e-8)

  expect_equal(clf$models_by_fold[[k]]$standardize$mu, as.numeric(mu_expected), tolerance = 1e-8)
  expect_equal(clf$models_by_fold[[k]]$standardize$sd, as.numeric(sd_expected), tolerance = 1e-8)
})
