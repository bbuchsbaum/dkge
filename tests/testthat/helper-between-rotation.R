# Deterministic fixtures for residual-space rotation tests and calibration.

dkge_test_rotation_simulate <- function(
    seed,
    n = 20L,
    p = 24L,
    rank = 2L,
    error = c("gaussian", "t5", "skew"),
    null_term = NULL,
    power_term = NULL,
    effect_size = 1) {
  error <- match.arg(error)
  set.seed(seed)
  dat <- data.frame(
    subject_id = paste0("s", seq_len(n)),
    group = factor(rep(c("A", "B"), length.out = n)),
    trait = scale(rnorm(n), center = TRUE, scale = FALSE)[, 1],
    age = scale(rnorm(n), center = TRUE, scale = FALSE)[, 1]
  )
  design <- dkge_subject_model(~ group * trait + age, dat)
  X <- design$X

  feature_pattern <- matrix(rnorm(2L * p), 2L, p)
  feature_pattern <- feature_pattern /
    sqrt(rowMeans(feature_pattern^2))
  coefficient <- matrix(0, ncol(X), 2L,
                        dimnames = list(colnames(X), c("latent1", "latent2")))

  if (!is.null(null_term)) {
    tested <- dkge:::.dkge_between_term_columns(design, null_term)
    eligible <- setdiff(seq_len(ncol(X)), tested)
    eligible <- setdiff(eligible, match("(Intercept)", colnames(X)))
    if (!length(eligible)) stop("No nuisance rows available for the fixture.")
    coefficient[eligible, 1L] <- seq(0.35, 0.8, length.out = length(eligible))
    coefficient[eligible, 2L] <- seq(-0.55, 0.25, length.out = length(eligible))
  }
  if (!is.null(power_term)) {
    tested <- dkge:::.dkge_between_term_columns(design, power_term)
    coefficient[tested, 1L] <- effect_size
  }
  signal <- X %*% coefficient %*% feature_pattern

  gaussian_error <- function() {
    latent <- matrix(rnorm(n * 2L), n, 2L)
    pattern <- matrix(rnorm(2L * p), 2L, p)
    latent %*% pattern + matrix(rnorm(n * p, sd = 0.35), n, p)
  }
  E <- switch(
    error,
    gaussian = gaussian_error(),
    t5 = {
      base <- gaussian_error()
      row_scale <- sqrt(3 / stats::rchisq(n, df = 5))
      base * row_scale
    },
    skew = {
      latent <- matrix(rexp(n * 2L) - 1, n, 2L)
      pattern <- matrix(rnorm(2L * p), 2L, p)
      innovation <- matrix(rexp(n * p) - 1, n, p) * 0.35
      latent %*% pattern + innovation
    }
  )
  Y <- signal + E
  rownames(Y) <- dat$subject_id
  colnames(Y) <- paste0("f", seq_len(p))
  target <- dkge_make_target(Y = Y, subject_ids = dat$subject_id)
  fit <- dkge_between_rrr(target, design, rank = rank)

  list(data = dat, design = design, X = X, Y = Y, fit = fit,
       coefficient = coefficient, feature_pattern = feature_pattern)
}

dkge_test_rotation_dense_null <- function(fixture, term, B, seed) {
  dkge_test_rotation_dense_result(fixture, term, B, seed)$null
}

dkge_test_rotation_dense_result <- function(fixture, term, B, seed) {
  X <- fixture$X
  Y <- fixture$Y
  fit <- fixture$fit
  Xred <- dkge:::.dkge_between_reduced_X(fixture$design, X, term)
  full_setup <- dkge:::.dkge_between_rrr_fast_setup(X, tol = fit$tol)
  red_setup <- dkge:::.dkge_between_rrr_fast_setup(Xred, tol = fit$tol)
  rotation_setup <- dkge:::.dkge_between_rotation_setup(
    full_setup, red_setup, Y, tol = fit$tol
  )

  set.seed(seed)
  descriptors <- sample.int(.Machine$integer.max, B, replace = TRUE)
  observed_map <- dkge:::.dkge_between_term_map_from_coef(
    fit$coef, fixture$design, term, drop = FALSE
  )
  observed_feature <- dkge:::.dkge_between_feature_stat(observed_map)
  observed_red <- dkge:::.dkge_rrr_fit_core(
    Xred, Y, rank = fit$rank, tol = fit$tol, warn_rank = FALSE
  )
  observed_stat <- dkge:::.dkge_between_stat(fit, observed_red)
  draws <- lapply(descriptors, function(rotation_seed) {
    G <- dkge:::.dkge_seeded_haar_rotation(
      rotation_seed, rotation_setup$residual_dimension
    )
    Yb <- rotation_setup$Q0 %*% rotation_setup$fitted_coordinates +
      rotation_setup$Qe %*% G %*% rotation_setup$residual_coordinates
    full_b <- dkge:::.dkge_rrr_fit_core(
      X, Yb, rank = fit$rank, tol = fit$tol, warn_rank = FALSE
    )
    red_b <- dkge:::.dkge_rrr_fit_core(
      Xred, Yb, rank = fit$rank, tol = fit$tol, warn_rank = FALSE
    )
    term_map <- dkge:::.dkge_between_term_map_from_coef(
      full_b$coef, fixture$design, term, drop = FALSE
    )
    list(
      stat = dkge:::.dkge_between_stat(full_b, red_b),
      feature = dkge:::.dkge_between_feature_stat(term_map)
    )
  })
  null <- vapply(draws, `[[`, numeric(1), "stat")
  feature_null <- do.call(rbind, lapply(draws, `[[`, "feature"))
  stat_tol <- sqrt(.Machine$double.eps) * max(1, abs(observed_stat))
  feature_tol <- sqrt(.Machine$double.eps) *
    pmax(1, abs(observed_feature))
  feature_p <- vapply(seq_along(observed_feature), function(index) {
    (1 + sum(feature_null[, index] >=
               observed_feature[[index]] - feature_tol[[index]])) / (B + 1)
  }, numeric(1))
  null_max <- apply(feature_null, 1L, max)
  maxT_p <- vapply(observed_feature, function(value) {
    tolerance <- sqrt(.Machine$double.eps) * max(1, abs(value))
    (1 + sum(null_max >= value - tolerance)) / (B + 1)
  }, numeric(1))
  list(
    statistic = observed_stat,
    null = null,
    p = (1 + sum(null >= observed_stat - stat_tol)) / (B + 1),
    feature_statistic = observed_feature,
    feature_p = feature_p,
    null_max = null_max,
    maxT_p = maxT_p
  )
}
