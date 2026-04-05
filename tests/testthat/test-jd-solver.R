# Tests for the joint diagonalisation solver integration.

test_that("dkge_jd_solve matches pooled eigenvectors when matrices commute", {
  set.seed(1)
  q <- 4
  S <- 3
  base <- qr.Q(qr(matrix(rnorm(q^2), q)))
  eigs <- matrix(runif(q * S, 0.5, 2), nrow = q)
  A_list <- lapply(seq_len(S), function(s) {
    base %*% diag(eigs[, s]) %*% t(base)
  })
  weights <- rep(1, S)
  Chat <- Reduce(`+`, Map(function(w, A) w * A, weights, A_list))
  Chat <- (Chat + t(Chat)) / 2

  jd <- dkge_jd_solve(
    A_list = A_list,
    weights = weights,
    rank = q,
    Chat = Chat,
    control = dkge_jd_control(maxit = 50, tol = 1e-10, linesearch = FALSE)
  )
  eig <- eigen(Chat, symmetric = TRUE)

  overlap <- abs(diag(crossprod(eig$vectors, jd$Q)))
  expect_true(all(overlap > 1 - 1e-6))
})

test_that("dkge_fit solver='jd' returns K-orthonormal components and reduces off-diagonal energy", {
  set.seed(123)
  S <- 4
  q <- 5
  P <- 6
  betas <- vector("list", S)
  designs <- vector("list", S)
  for (s in seq_len(S)) {
    betas[[s]] <- matrix(rnorm(q * P), q, P)
    designs[[s]] <- diag(q)
  }
  K <- diag(q)
  fit_jd <- dkge_fit(
    betas,
    designs = designs,
    K = K,
    rank = 3,
    solver = "jd",
    jd_control = dkge_jd_control(maxit = 150, tol = 1e-8)
  )

  expect_equal(fit_jd$solver, "jd")
  expect_false(is.null(fit_jd$jd))
  expect_equal(nrow(fit_jd$U), q)
  expect_equal(ncol(fit_jd$U), ncol(fit_jd$scores_matrix))

  gram_K <- crossprod(fit_jd$U, fit_jd$K %*% fit_jd$U)
  expect_true(max(abs(gram_K - diag(diag(gram_K)))) < 1e-8)
  expect_true(max(abs(diag(gram_K) - 1)) < 1e-8)

  mask_list <- replicate(length(fit_jd$contribs), dkge:::.dkge_jd_default_mask(q), simplify = FALSE)
  identity_eval <- dkge:::.dkge_jd_eval(diag(q), fit_jd$contribs, fit_jd$weights, mask_list, compute_grad = FALSE)
  final_eval <- dkge:::.dkge_jd_eval(fit_jd$jd$Q, fit_jd$contribs, fit_jd$weights, mask_list, compute_grad = FALSE)
  expect_lt(final_eval$offdiag, identity_eval$offdiag)
})

test_that("JD solver integrates with CPCA splits", {
  set.seed(321)
  S <- 3
  q <- 4
  P <- 5
  betas <- vector("list", S)
  designs <- vector("list", S)
  for (s in seq_len(S)) {
    betas[[s]] <- matrix(rnorm(q * P), q, P)
    designs[[s]] <- diag(q)
  }
  K <- diag(q)

  fit_cpca <- dkge_fit(
    betas,
    designs = designs,
    K = K,
    rank = 2,
    solver = "jd",
    cpca_blocks = 1:2,
    cpca_part = "both",
    jd_control = dkge_jd_control(maxit = 120, tol = 1e-8)
  )

  expect_true(!is.null(fit_cpca$cpca))
  expect_true(!is.null(fit_cpca$cpca$jd_design))
  expect_true(!is.null(fit_cpca$cpca$jd_resid))
  expect_equal(nrow(fit_cpca$cpca$U_design), q)
  expect_equal(ncol(fit_cpca$cpca$U_design), 2)
  expect_equal(nrow(fit_cpca$cpca$U_resid), q)
  expect_equal(ncol(fit_cpca$cpca$U_resid), 2)
})

test_that("JD objective decreases with backtracking line search", {
  set.seed(2024)
  q <- 5
  S <- 4
  A_list <- lapply(seq_len(S), function(s) {
    M <- matrix(rnorm(q * q), q, q)
    M <- crossprod(M) + (s / 10) * diag(q)
    (M + t(M)) / 2
  })
  weights <- runif(S, 0.5, 1.5)
  ctrl <- dkge_jd_control(maxit = 40, tol = 1e-10, record = TRUE, verbose = FALSE)
  sol <- dkge_jd_solve(A_list, weights = weights, control = ctrl)
  expect_false(is.null(sol$history))
  values <- sol$history$value
  expect_gt(length(values), 1)
  diffs <- diff(values)
  expect_true(all(diffs <= 1e-10))
  expect_lt(sol$value, values[1])
})

test_that("JD gradient matches finite-difference directional derivative", {
  set.seed(99)
  q <- 3
  S <- 2
  A_list <- lapply(seq_len(S), function(s) {
    M <- matrix(rnorm(q * q), q, q)
    M <- crossprod(M) + diag(0.1, q)
    (M + t(M)) / 2
  })
  weights <- runif(S, 0.8, 1.2)
  mask_list <- replicate(S, dkge:::.dkge_jd_default_mask(q), simplify = FALSE)
  Q <- diag(q)
  eval <- dkge:::.dkge_jd_eval(Q, A_list, weights, mask_list, compute_grad = TRUE)
  H_raw <- matrix(rnorm(q * q), q, q)
  H <- 0.5 * (H_raw - t(H_raw))
  direction <- Q %*% H
  eps <- 1e-6
  Q_plus <- dkge:::.dkge_jd_retract(Q + eps * direction)
  Q_minus <- dkge:::.dkge_jd_retract(Q - eps * direction)
  val_plus <- dkge:::.dkge_jd_eval(Q_plus, A_list, weights, mask_list, compute_grad = FALSE)$value
  val_minus <- dkge:::.dkge_jd_eval(Q_minus, A_list, weights, mask_list, compute_grad = FALSE)$value
  fd <- (val_plus - val_minus) / (2 * eps)
  grad_dot <- sum(eval$grad * direction)
  expect_equal(fd, grad_dot, tolerance = 1e-5)
})

test_that("JD fit supports LOSO contrast, prediction, and transport pipelines", {
  set.seed(314)
  S <- 4
  q <- 5
  P <- 6
  betas <- vector("list", S)
  designs <- vector("list", S)
  centroids <- vector("list", S)
  for (s in seq_len(S)) {
    betas[[s]] <- matrix(rnorm(q * P), q, P)
    rownames(betas[[s]]) <- paste0("e", seq_len(q))
    designs[[s]] <- diag(q)
    colnames(designs[[s]]) <- paste0("e", seq_len(q))
    centroids[[s]] <- matrix(runif(P * 3), P, 3)
  }
  data_obj <- dkge_data(betas, designs = designs)
  fit_jd <- dkge_fit(data_obj, K = diag(q), rank = 3, solver = "jd", keep_X = TRUE)

  cvec <- rnorm(q)
  loso <- dkge_loso_contrast(fit_jd, s = 1, contrasts = cvec)
  expect_type(loso$v, "double")
  expect_length(loso$v, P)

  contrast <- dkge_contrast(fit_jd, cvec, method = "loso")
  expect_s3_class(contrast, "dkge_contrasts")

  transport <- dkge_transport_contrasts_to_medoid(fit_jd, contrast,
                                                  medoid = 1,
                                                  centroids = centroids,
                                                  betas = betas)
  expect_equal(length(transport), length(contrast$values))
  expect_true(is.matrix(transport[[1]]$subj_values))

  cv_scores <- numeric(S)
  for (holdout in seq_len(S)) {
    train_ids <- setdiff(seq_len(S), holdout)
    fit_train <- dkge_fit(dkge_data(betas[train_ids], designs = designs[train_ids]),
                          K = diag(q), rank = 3, solver = "jd")
    loadings <- dkge_project_clusters(fit_train, betas[[holdout]])
    recon <- fit_train$U %*% t(loadings)
    cv_scores[holdout] <- sum(recon * recon) / (sum(betas[[holdout]] * betas[[holdout]]) + 1e-12)
  }
  expect_true(all(is.finite(cv_scores)))
  expect_true(all(cv_scores >= 0))
  expect_true(all(cv_scores <= 1 + 1e-6))
})
