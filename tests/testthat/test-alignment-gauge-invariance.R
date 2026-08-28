library(testthat)

court_path <- dkge_validation_path("functional-alignment-court", "court.R")
if (!is.na(court_path)) source(court_path, local = TRUE)

skip_without_court <- function() {
  skip_if(is.na(court_path), "validation sources are not available in this layout")
}

# The latent basis U is only identified up to a K-orthogonal rotation. Loadings
# rotate with it (A = B' K U -> A Q), so a transport cost built from loadings
# must be a function of the *subspace*, not of the coordinates chosen inside
# it. If it is not, every reported alignment depends on an arbitrary gauge.
#
# This is the boundary that makes K-Procrustes load-bearing: within one basis
# the gauge cancels and no alignment is needed; across fold bases the subspaces
# genuinely differ and the cancellation no longer applies.

fa_gauge_random_rotation <- function(r, seed) {
  set.seed(seed)
  qr.Q(qr(matrix(stats::rnorm(r * r), r, r)))
}

fa_gauge_operators <- function(dat, cfg, features) {
  # fa_court_map() takes one entry per contrast, each a list of per-subject
  # vectors. Only the operators are of interest here, so a single zero field
  # is enough to drive the mapper.
  zero_field <- stats::setNames(
    lapply(seq_along(features), function(s) rep(0, nrow(features[[s]]))),
    names(features)
  )
  fa_court_map(list(gauge = zero_field), features, dat, cfg, "gauge")$operators
}

test_that("transport is invariant to the arbitrary rotation of a shared basis", {
  skip_on_cran()
  skip_without_court()
  cfg <- fa_court_base_config()
  cfg$S <- 8L
  dat <- fa_court_generate(cfg, seed = 7301L)
  fit <- fa_court_fit(dat, cfg)
  loadings <- dkge:::.dkge_fit_subject_loadings(fit)

  Q <- fa_gauge_random_rotation(ncol(fit$U), seed = 11L)
  expect_equal(crossprod(Q), diag(ncol(Q)), tolerance = 1e-12)
  rotated <- lapply(loadings, function(A) A %*% Q)

  base_ops <- fa_gauge_operators(dat, cfg, loadings)
  rotated_ops <- fa_gauge_operators(dat, cfg, rotated)

  # Not "close": the cost is a function of pairwise distances between
  # row-normalised loadings, and an orthogonal Q preserves those exactly.
  expect_lt(fa_court_operator_distance(rotated_ops, base_ops), 1e-10)
})

test_that("a rotated basis is the same subspace, and a fold basis is not", {
  skip_on_cran()
  skip_without_court()
  cfg <- fa_court_base_config()
  cfg$S <- 8L
  dat <- fa_court_generate(cfg, seed = 7302L)
  fit <- fa_court_fit(dat, cfg)
  K <- fit$K
  U <- fit$U
  Q <- fa_gauge_random_rotation(ncol(U), seed = 13L)

  # A gauge change leaves the K-projector untouched...
  proj <- function(V) V %*% t(V) %*% K
  expect_equal(proj(U %*% Q), proj(U), tolerance = 1e-10)

  # ...whereas a held-out fold basis is a genuinely different subspace, so the
  # projector moves. This is why Procrustes cannot be dismissed as cosmetic
  # once features from different folds are compared or averaged.
  fold <- dkge_loso_contrast(fit, s = 1L, contrasts = dat$contrasts[, 1])
  drift <- norm(proj(fold$basis) - proj(U), "F")
  expect_gt(drift, 1e-6)

  # And the fold basis is still K-orthonormal, so the drift is subspace
  # movement rather than a broken normalisation.
  expect_equal(
    crossprod(fold$basis, K %*% fold$basis),
    diag(ncol(fold$basis)),
    tolerance = 1e-8
  )
})
