test_that("dkge_subject_model builds formula metadata", {
  dat <- data.frame(
    subject_id = paste0("s", 1:8),
    group = factor(rep(c("A", "B"), each = 4)),
    trait = c(-1.2, -0.4, 0.3, 1.1, -0.8, -0.1, 0.6, 1.4),
    age = c(22, 29, 25, 31, 24, 35, 28, 33)
  )

  design <- dkge_subject_model(~ group * trait + age, data = dat)

  expect_s3_class(design, "dkge_subject_model")
  expect_equal(design$subject_ids, dat$subject_id)
  expect_true("group:trait" %in% design$term_names)
  expect_true(length(design$term_columns[["group:trait"]]) >= 1)
  expect_equal(nrow(design$X), nrow(dat))
})

test_that("dkge_make_target wraps matrix targets with metadata", {
  Y <- matrix(rnorm(6 * 5), 6, 5)
  rownames(Y) <- paste0("s", 1:6)
  colnames(Y) <- paste0("v", 1:5)
  precision <- seq(0.5, 1.5, length.out = 5)

  target <- dkge_make_target(Y = Y, precision = precision)

  expect_s3_class(target, "dkge_target")
  expect_equal(target$Y, Y)
  expect_equal(target$feature_ids, colnames(Y))
  expect_equal(target$precision, precision)
})

test_that("dkge_between_rrr recovers rank-one formula-structured signal", {
  set.seed(211)
  dat <- data.frame(
    subject_id = paste0("s", 1:24),
    group = factor(rep(c("A", "B"), each = 12)),
    trait = c(seq(-1, 1, length.out = 12), seq(-0.7, 1.3, length.out = 12)),
    age = scale(c(30, 21, 28, 35, 25, 40, 32, 37, 24, 29, 34, 38,
                  23, 27, 36, 31, 42, 26, 33, 39, 22, 30, 35, 41),
                center = TRUE, scale = FALSE)[, 1]
  )
  design <- dkge_subject_model(~ group * trait + age, data = dat)
  X <- design$X

  a <- rep(0, ncol(X))
  names(a) <- colnames(X)
  a["groupB"] <- 1.25
  a["groupB:trait"] <- -0.7
  v <- seq(-1, 1, length.out = 9)
  B <- a %*% t(v)
  rownames(B) <- colnames(X)
  colnames(B) <- paste0("feature", seq_len(ncol(B)))
  Y <- X %*% B
  colnames(Y) <- paste0("feature", seq_len(ncol(Y)))

  target <- dkge_make_target(Y = Y, subject_ids = dat$subject_id)
  fit <- dkge_between_rrr(target, design, rank = 1)

  expect_s3_class(fit, "dkge_between_rrr")
  expect_lt(max(abs(fit$fitted - Y)), 1e-8)
  expect_equal(coef(fit), B, tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(dkge_term_map(fit, "groupB"), B["groupB", ], tolerance = 1e-8)
  expect_equal(dkge_term_map(fit, "groupB:trait"),
               B["groupB:trait", ], tolerance = 1e-8)
})

test_that("dkge_between_rrr supports target-derived feature weights", {
  set.seed(212)
  dat <- data.frame(
    subject_id = paste0("s", 1:10),
    x = seq(-1, 1, length.out = 10)
  )
  design <- dkge_subject_model(~ x, data = dat)
  X <- design$X
  B <- matrix(c(0, 2, 0, -1, 0, 0.5), nrow = ncol(X), byrow = TRUE)
  Y <- X %*% B
  target <- dkge_make_target(Y = Y,
                             subject_ids = dat$subject_id,
                             feature_weights = c(1, 2, 3))

  fit <- dkge_between_rrr(target, design, rank = 1, weights = "target")

  expect_equal(dim(fit$coef), dim(B))
  expect_true(all(is.finite(fit$coef)))
  expect_equal(fit$feature_weights, c(1, 2, 3))
})

test_that("full-rank dkge_between_rrr matches weighted least squares oracle", {
  set.seed(216)
  n <- 14L
  X <- cbind(1, x = seq(-1, 1, length.out = n), z = rnorm(n))
  rownames(X) <- paste0("s", seq_len(n))
  beta <- matrix(rnorm(ncol(X) * 5), nrow = ncol(X),
                 dimnames = list(colnames(X), paste0("f", 1:5)))
  Y <- X %*% beta + matrix(rnorm(n * 5, sd = 0.03), n, 5)
  colnames(Y) <- colnames(beta)
  rownames(Y) <- rownames(X)
  subject_w <- seq(0.5, 2, length.out = n)

  target <- dkge_make_target(Y = Y, subject_weights = subject_w)
  fit <- dkge_between_rrr(target, X, rank = ncol(X), weights = "target")

  sw <- sqrt(subject_w)
  beta_ref <- qr.solve(X * sw, Y * sw)

  expect_equal(coef(fit), beta_ref, tolerance = 1e-8, ignore_attr = TRUE)
  expect_lt(max(abs(fit$fitted - X %*% beta_ref)), 1e-8)
})

test_that("RRR residual SSE decreases monotonically with rank", {
  set.seed(217)
  n <- 16L
  X <- cbind(1, x = rnorm(n), z = rnorm(n))
  rownames(X) <- paste0("s", seq_len(n))
  Y <- X %*% matrix(rnorm(3 * 8), 3, 8) + matrix(rnorm(n * 8), n, 8)
  rownames(Y) <- rownames(X)

  target <- dkge_make_target(Y = Y)
  fit1 <- dkge_between_rrr(target, X, rank = 1)
  fit2 <- dkge_between_rrr(target, X, rank = 2)
  fit3 <- dkge_between_rrr(target, X, rank = 3)
  sse <- function(obj) sum(obj$residuals^2)

  expect_lte(sse(fit2), sse(fit1) + 1e-8)
  expect_lte(sse(fit3), sse(fit2) + 1e-8)
  expect_true(all(is.finite(fit3$coef)))
})

test_that("between RRR is invariant to subject order when IDs align", {
  set.seed(218)
  dat <- data.frame(subject_id = paste0("s", 1:15),
                    x = rnorm(15),
                    z = rnorm(15))
  design <- dkge_subject_model(~ x + z, dat)
  Y <- design$X %*% matrix(rnorm(ncol(design$X) * 6), ncol(design$X), 6)
  rownames(Y) <- dat$subject_id
  colnames(Y) <- paste0("f", 1:6)

  target <- dkge_make_target(Y = Y)
  fit_a <- dkge_between_rrr(target, design, rank = 2)

  perm <- sample(seq_len(nrow(Y)))
  target_perm <- dkge_make_target(Y = Y[perm, , drop = FALSE])
  fit_b <- dkge_between_rrr(target_perm, design, rank = 2)

  expect_equal(coef(fit_b), coef(fit_a), tolerance = 1e-8)
})

test_that("between RRR exposes fitted, residual, and prediction methods", {
  set.seed(2181)
  dat <- data.frame(
    subject_id = paste0("s", 1:14),
    group = factor(rep(c("A", "B"), each = 7)),
    trait = seq(-1, 1, length.out = 14)
  )
  design <- dkge_subject_model(~ group * trait, dat)
  Y <- design$X %*% matrix(rnorm(ncol(design$X) * 4), ncol(design$X), 4)
  rownames(Y) <- dat$subject_id
  colnames(Y) <- paste0("f", 1:4)
  fit <- dkge_between_rrr(dkge_make_target(Y = Y), design, rank = 2)

  expect_equal(fitted(fit), fit$fitted)
  expect_equal(residuals(fit), fit$residuals)
  expect_equal(predict(fit), fit$fitted)
  expect_equal(dim(predict(fit, type = "scores")), c(nrow(Y), fit$rank))

  new_dat <- data.frame(
    group = factor(c("A", "B"), levels = levels(dat$group)),
    trait = c(-0.25, 0.75)
  )
  pred <- predict(fit, newdata = new_dat)
  expect_equal(dim(pred), c(2L, ncol(Y)))
  expect_true(all(is.finite(pred)))
})

test_that("zero target feature weights mask features before weighted fitting", {
  set.seed(219)
  X <- cbind(1, x = seq(-1, 1, length.out = 10))
  rownames(X) <- paste0("s", seq_len(nrow(X)))
  Y <- X %*% matrix(rnorm(2 * 4), 2, 4)
  rownames(Y) <- rownames(X)
  colnames(Y) <- paste0("f", 1:4)
  target <- dkge_make_target(Y = Y, feature_weights = c(1, 0, 2, 0))

  fit <- dkge_between_rrr(target, X, rank = 1, weights = "target")

  expect_equal(colnames(fit$coef), c("f1", "f3"))
  expect_equal(fit$feature_weights, c(1, 2))
})

test_that("matrix design inputs without column names get stable term metadata", {
  set.seed(220)
  X <- cbind(1, rnorm(9))
  rownames(X) <- paste0("s", 1:9)
  colnames(X) <- NULL
  Y <- X %*% matrix(c(0, 1, 2, 1), 2, 2)
  rownames(Y) <- rownames(X)
  target <- dkge_make_target(Y = Y)

  fit <- dkge_between_rrr(target, X, rank = 1)

  expect_equal(rownames(coef(fit)), c("x1", "x2"))
  expect_equal(length(dkge_term_map(fit, "x2")), ncol(fit$coef))
})

test_that("between layer rejects adversarial invalid inputs clearly", {
  dat <- data.frame(subject_id = paste0("s", 1:6),
                    x = c(1, 1, 1, 2, 2, 2))
  expect_error(dkge_subject_model(~ x + I(2 * x), dat), "rank deficient")

  Y_bad <- matrix(rnorm(6 * 3), 6, 3)
  Y_bad[1, 1] <- Inf
  expect_error(dkge_make_target(Y = Y_bad), "non-finite")

  X <- cbind(1, x = seq_len(6))
  rownames(X) <- paste0("s", 1:6)
  Y <- matrix(rnorm(6 * 3), 6, 3, dimnames = list(paste0("q", 1:6), NULL))
  target <- dkge_make_target(Y = Y)
  expect_error(dkge_between_rrr(target, X, rank = 1), "present")
})

test_that("randomized RRR fuzz checks finite outputs and rank monotonicity", {
  for (seed in 301:306) {
    set.seed(seed)
    n <- sample(10:16, 1)
    k <- sample(2:4, 1)
    p <- sample(5:12, 1)
    X <- cbind(1, matrix(rnorm(n * (k - 1)), n, k - 1))
    rownames(X) <- paste0("s", seq_len(n))
    colnames(X) <- paste0("x", seq_len(k))
    Y <- X %*% matrix(rnorm(k * p), k, p) + matrix(rnorm(n * p, sd = 0.2), n, p)
    rownames(Y) <- rownames(X)
    target <- dkge_make_target(Y = Y)

    fit1 <- dkge_between_rrr(target, X, rank = 1)
    fitk <- dkge_between_rrr(target, X, rank = min(k, p))

    expect_true(all(is.finite(fit1$coef)))
    expect_true(all(is.finite(fitk$loadings_brain)))
    expect_lte(sum(fitk$residuals^2), sum(fit1$residuals^2) + 1e-8)
  }
})

test_that("between RRR handles representative S much smaller than feature count", {
  set.seed(221)
  n <- 12L
  p <- 250L
  X <- cbind(1, x = rnorm(n), z = rnorm(n))
  rownames(X) <- paste0("s", 1:n)
  Y <- X %*% matrix(rnorm(3 * p), 3, p) + matrix(rnorm(n * p, sd = 0.1), n, p)
  rownames(Y) <- rownames(X)
  target <- dkge_make_target(Y = Y)

  fit <- dkge_between_rrr(target, X, rank = 2)

  expect_equal(dim(fit$coef), c(ncol(X), p))
  expect_equal(dim(fit$loadings_brain), c(p, 2L))
  expect_true(all(is.finite(fit$singular_values)))
})

test_that("near-collinear but full-rank between-subject designs remain numerically stable", {
  set.seed(222)
  n <- 28L
  x1 <- scale(seq(-1, 1, length.out = n) + rnorm(n, sd = 1e-4),
              center = TRUE, scale = FALSE)[, 1]
  x2 <- x1 + rnorm(n, sd = 2e-6)
  z <- rnorm(n)
  X <- cbind("(Intercept)" = 1, x1 = x1, x2 = x2, z = z)
  rownames(X) <- paste0("s", seq_len(n))
  expect_equal(qr(X, tol = 1e-12)$rank, ncol(X))
  expect_gt(kappa(X), 1e5)

  beta <- matrix(rnorm(ncol(X) * 5), ncol(X), 5,
                 dimnames = list(colnames(X), paste0("f", seq_len(5))))
  Y <- X %*% beta + matrix(rnorm(n * 5, sd = 0.01), n, 5)
  rownames(Y) <- rownames(X)
  colnames(Y) <- colnames(beta)
  target <- dkge_make_target(Y = Y)

  fit <- dkge_between_rrr(target, X, rank = ncol(X), tol = 1e-12)
  beta_ref <- qr.solve(X, Y, tol = 1e-12)

  expect_equal(coef(fit), beta_ref, tolerance = 2e-4, ignore_attr = TRUE)
  expect_lt(max(abs(fit$fitted - X %*% beta_ref)), 2e-4)
})

test_that("component score target derives subject-by-component matrix from DKGE fit", {
  set.seed(213)
  q <- 3L
  S <- 4L
  P <- 5L
  effects <- paste0("e", seq_len(q))
  betas <- replicate(S, matrix(rnorm(q * P), q, P,
                               dimnames = list(effects, NULL)), simplify = FALSE)
  designs <- replicate(S, diag(q), simplify = FALSE)
  designs <- lapply(designs, function(X) {
    colnames(X) <- effects
    X
  })
  fit <- dkge_fit(betas, designs = designs, K = diag(q), rank = 2)

  target <- dkge_make_target(fit, type = "component_scores")

  expect_s3_class(target, "dkge_target")
  expect_equal(dim(target$Y), c(S, 2L))
  expect_equal(target$subject_ids, fit$subject_ids)
})

test_that("transported map targets integrate with between-subject RRR and permutation inference", {
  set.seed(223)
  data <- create_mismatched_data(P_vec = c(4, 5, 3, 6, 4), q = 3, seed = 223)
  ids <- paste0("s", seq_len(data$S))
  effects <- paste0("e", seq_len(data$q))
  data$betas <- lapply(data$betas, function(B) {
    rownames(B) <- effects
    B
  })
  data$designs <- lapply(data$designs, function(X) {
    colnames(X) <- effects
    X
  })
  fit <- dkge_fit(
    dkge_data(data$betas, designs = data$designs, subject_ids = ids),
    K = data$K,
    rank = 2
  )
  transport <- dkge_transport_spec(
    centroids = data$centroids,
    medoid = 1L,
    method = "sinkhorn",
    epsilon = 0.1
  )
  contrast <- matrix(c(1, -1, 0), ncol = 1,
                     dimnames = list(NULL, "c1"))

  expect_no_warning(
    target <- dkge_make_target(
      fit,
      type = "transported_maps",
      contrast = contrast,
      transport = transport,
      crossfit = "analytic"
    )
  )
  meta <- data.frame(
    subject_id = ids,
    group = factor(c("A", "A", "B", "B", "B")),
    trait = seq(-1, 1, length.out = data$S)
  )
  design <- dkge_subject_model(~ group * trait, meta)
  between <- dkge_between_rrr(target, design, rank = 1)
  perm <- dkge_between_permute(
    between,
    terms = c("group", "trait", "group:trait"),
    B = 9,
    seed = 77
  )

  expect_s3_class(target, "dkge_target")
  expect_equal(target$type, "transported_maps")
  expect_equal(dim(target$Y), c(data$S, data$P[1]))
  expect_equal(target$geometry$medoid, 1L)
  expect_equal(dim(between$coef), c(ncol(design$X), data$P[1]))
  expect_true(all(is.finite(between$fitted)))
  expect_true(all(perm$summary$p >= 0 & perm$summary$p <= 1))
})

test_that("dkge_between_permute returns formula-aware term tests", {
  set.seed(214)
  n <- 18L
  dat <- data.frame(
    subject_id = paste0("s", seq_len(n)),
    group = factor(rep(c("A", "B"), each = n / 2)),
    trait = c(seq(-1, 1, length.out = n / 2),
              seq(-0.8, 1.2, length.out = n / 2))
  )
  design <- dkge_subject_model(~ group * trait, dat)
  X <- design$X
  beta <- matrix(0, nrow = ncol(X), ncol = 6,
                 dimnames = list(colnames(X), paste0("f", 1:6)))
  beta["groupB:trait", ] <- c(2, 1, 0, -1, -2, -1)
  Y <- X %*% beta + matrix(rnorm(n * 6, sd = 0.05), n, 6)
  target <- dkge_make_target(Y = Y, subject_ids = dat$subject_id)
  fit <- dkge_between_rrr(target, design, rank = 1)

  perm <- dkge_between_permute(fit,
                               terms = c("group", "trait", "group:trait"),
                               B = 25,
                               seed = 99,
                               adjust = "fdr")

  expect_s3_class(perm, "dkge_between_permutation")
  expect_equal(perm$summary$term, c("group", "trait", "group:trait"))
  expect_true(all(perm$summary$p >= 0 & perm$summary$p <= 1))
  expect_equal(length(perm$tests[["group:trait"]]$null), 25L)
  expect_true(perm$summary$statistic[perm$summary$term == "group:trait"] > 0)
})

test_that("dkge_between_permute returns featurewise maxT inference and term maps", {
  set.seed(2141)
  n <- 16L
  dat <- data.frame(
    subject_id = paste0("s", seq_len(n)),
    group = factor(rep(c("A", "B"), each = n / 2)),
    trait = seq(-1, 1, length.out = n)
  )
  design <- dkge_subject_model(~ group * trait, dat)
  X <- design$X
  beta <- matrix(0, nrow = ncol(X), ncol = 5,
                 dimnames = list(colnames(X), paste0("f", 1:5)))
  beta["groupB:trait", ] <- c(2, 1, 0, -1, -2)
  Y <- X %*% beta + matrix(rnorm(n * 5, sd = 0.08), n, 5)
  target <- dkge_make_target(Y = Y, subject_ids = dat$subject_id)
  fit <- dkge_between_rrr(target, design, rank = 1)

  perm <- dkge_between_permute(
    fit,
    terms = "group:trait",
    B = 19,
    seed = 321,
    scope = "both",
    feature_adjust = "maxT"
  )

  ft <- perm$feature_tests[["group:trait"]]
  expect_equal(dim(perm$term_maps[["group:trait"]]), c(1L, ncol(Y)))
  expect_equal(length(ft$statistic), ncol(Y))
  expect_true(all(ft$p >= 0 & ft$p <= 1))
  expect_true(all(ft$p_adjusted >= 0 & ft$p_adjusted <= 1))
  expect_equal(length(ft$null_max), 19L)
})

test_that("dkge_between_permute is reproducible with a seed", {
  set.seed(215)
  dat <- data.frame(subject_id = paste0("s", 1:12),
                    x = seq(-1, 1, length.out = 12))
  design <- dkge_subject_model(~ x, dat)
  Y <- design$X %*% matrix(c(0, 1, 2, 0, -1, 0), nrow = 2, byrow = TRUE)
  Y <- Y + matrix(rnorm(length(Y), sd = 0.1), nrow(Y), ncol(Y))
  target <- dkge_make_target(Y = Y, subject_ids = dat$subject_id)
  fit <- dkge_between_rrr(target, design, rank = 1)

  a <- dkge_between_permute(fit, terms = "x", B = 20, seed = 123)
  b <- dkge_between_permute(fit, terms = "x", B = 20, seed = 123)

  expect_equal(a$tests$x$null, b$tests$x$null)
  expect_equal(a$summary$p, b$summary$p)
})

test_that("dkge_between_permute restores caller RNG state after seeded run", {
  set.seed(216)
  dat <- data.frame(subject_id = paste0("s", 1:10),
                    x = seq(-1, 1, length.out = 10))
  design <- dkge_subject_model(~ x, dat)
  Y <- design$X %*% matrix(c(0, 1, -1, 0), 2, 2) +
    matrix(rnorm(20, sd = 0.1), 10, 2)
  target <- dkge_make_target(Y = Y, subject_ids = dat$subject_id)
  fit <- dkge_between_rrr(target, design, rank = 1)

  set.seed(999)
  expected <- runif(4)
  set.seed(999)
  invisible(dkge_between_permute(fit, terms = "x", B = 5, seed = 123))
  observed <- runif(4)

  expect_equal(observed, expected)
})

test_that("singleton exchangeability blocks make Freedman-Lane null deterministic", {
  set.seed(217)
  dat <- data.frame(subject_id = paste0("s", 1:9),
                    x = seq(-1, 1, length.out = 9))
  design <- dkge_subject_model(~ x, dat)
  Y <- design$X %*% matrix(c(0, 1, 0, -1), 2, 2) +
    matrix(rnorm(18, sd = 0.05), 9, 2)
  target <- dkge_make_target(Y = Y, subject_ids = dat$subject_id)
  fit <- dkge_between_rrr(target, design, rank = 1)

  perm <- dkge_between_permute(fit,
                               terms = "x",
                               B = 6,
                               blocks = seq_len(nrow(Y)),
                               seed = 42)

  expect_equal(perm$tests$x$null, rep(perm$tests$x$statistic, 6),
               tolerance = 1e-10)
  expect_equal(perm$summary$p, 1)
})
