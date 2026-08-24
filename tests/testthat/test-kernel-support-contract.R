library(testthat)

make_kernel_support_fixture <- function(gamma = 1, S = 5L, P = 12L,
                                        seed = 824L) {
  q <- 7L
  K <- matrix(0, q, q)
  K[1L, 1L] <- 1
  for (pair in list(2:3, 4:5, 6:7)) {
    K[pair, pair] <- matrix(c(1, gamma, gamma, 1), 2L, 2L)
  }
  set.seed(seed)
  B_list <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
  X_list <- replicate(S, diag(q), simplify = FALSE)
  list(K = K, B_list = B_list, X_list = X_list, q = q, S = S)
}

capture_dkge_warnings <- function(expr) {
  classes <- character(0)
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      classes <<- c(classes, class(w)[[1L]])
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, classes = classes)
}

test_that("PSD roots preserve the exact support at every positive scale", {
  fixture <- make_kernel_support_fixture()
  reference <- dkge:::.dkge_kernel_geometry(fixture$K)

  expect_equal(reference$rank, 4L)
  expect_equal(reference$nullity, 3L)
  for (scale in c(1e-12, 1, 1e12)) {
    roots <- dkge:::.dkge_kernel_geometry(scale * fixture$K)
    expect_equal(roots$rank, 4L)
    expect_equal(roots$support_projector, reference$support_projector,
                 tolerance = 1e-8)
    expect_equal(roots$Khalf %*% roots$Khalf, scale * fixture$K,
                 tolerance = 1e-8)
    expect_equal(roots$Kihalf %*% (scale * fixture$K) %*% roots$Kihalf,
                 roots$support_projector, tolerance = 1e-8)
  }

  exact <- kernel_roots(fixture$K)
  regularized <- kernel_roots(fixture$K, jitter = 1e-8)
  expect_equal(exact$rank, 4L)
  expect_false(exact$regularized)
  expect_equal(regularized$rank, fixture$q)
  expect_true(regularized$regularized)
  expect_equal(regularized$original_rank, 4L)
})

test_that("fit rank is capped by image(K) and remains exactly K-orthonormal", {
  fixture <- make_kernel_support_fixture()
  warned <- capture_dkge_warnings(
    dkge_fit(
      fixture$B_list, fixture$X_list, fixture$K,
      rank = fixture$q, ridge = 0.2, w_method = "none",
      effect_scaling = "none"
    )
  )
  fit <- warned$value

  expect_contains(warned$classes, "dkge_kernel_rank_warning")
  expect_equal(fit$rank_requested, fixture$q)
  expect_equal(fit$rank, 4L)
  expect_equal(fit$kernel_rank, 4L)
  expect_equal(fit$kernel_nullity, 3L)
  expect_true(fit$rank_reduced)
  expect_equal(fit$kernel_diagnostics$status, "singular")
  expect_equal(crossprod(fit$U, fit$K %*% fit$U), diag(4), tolerance = 1e-8)

  null_energy <- fit$kernel_support_projector - diag(fixture$q)
  expect_equal(unname(null_energy %*% fit$Chat), matrix(0, fixture$q, fixture$q),
               tolerance = 1e-8)
})

test_that("fixed-geometry CV cannot reward a kernel for deleting target energy", {
  fixture <- make_kernel_support_fixture(S = 4L)
  q <- fixture$q
  isotropic_B <- replicate(fixture$S, diag(q), simplify = FALSE)
  isotropic_X <- replicate(fixture$S, diag(q), simplify = FALSE)
  near <- make_kernel_support_fixture(gamma = 0.999999)$K

  cv <- suppressWarnings(dkge_cv_kernel_grid(
    isotropic_B, isotropic_X,
    K_grid = list(identity = diag(q), near = near, singular = fixture$K),
    rank = 4L, w_method = "none",
    kernel_rank_policy = "allow_singular"
  ))

  expect_equal(cv$validation$kind, "effect_space_identity")
  expect_true(all(cv$table$admissible))
  expect_equal(cv$table$mean, rep(4 / q, 3L), tolerance = 1e-8)
  expect_lt(max(cv$table$mean), 0.999)
  expect_equal(cv$table$kernel_rank, c(q, q, 4L))

  scaled <- suppressWarnings(dkge_cv_kernel_grid(
    isotropic_B, isotropic_X,
    K_grid = list(small = 1e-12 * near, large = 1e12 * near),
    rank = 4L, w_method = "none"
  ))
  expect_equal(scaled$table$mean[1L], scaled$table$mean[2L], tolerance = 1e-8)
})

test_that("CV exposes and excludes impossible kernel-rank pairs", {
  fixture <- make_kernel_support_fixture(S = 4L)
  q <- fixture$q
  isotropic_B <- replicate(fixture$S, diag(q), simplify = FALSE)
  isotropic_X <- replicate(fixture$S, diag(q), simplify = FALSE)

  grid_warnings <- capture_dkge_warnings(dkge_cv_kernel_grid(
    isotropic_B, isotropic_X,
    K_grid = list(identity = diag(q), singular = fixture$K),
    rank = q, w_method = "none"
  ))
  grid <- grid_warnings$value
  singular_row <- grid$table[grid$table$kernel == "singular", , drop = FALSE]
  expect_false(singular_row$admissible)
  expect_match(singular_row$reason, "full-rank policy")
  expect_equal(grid$pick, "identity")
  expect_contains(grid_warnings$classes, "dkge_cv_kernel_rank_warning")

  expect_error(
    dkge_cv_rank_loso(
      isotropic_B, isotropic_X, fixture$K,
      ranks = 1:4, w_method = "none"
    ),
    class = "dkge_cv_kernel_rank_error"
  )

  rank_warnings <- capture_dkge_warnings(dkge_cv_rank_loso(
    isotropic_B, isotropic_X, fixture$K,
    ranks = 1:q, w_method = "none",
    kernel_rank_policy = "allow_singular"
  ))
  rank_cv <- rank_warnings$value
  expect_equal(rank_cv$table$param, 1:4)
  expect_equal(rank_cv$inadmissible_ranks, 5:7)
  expect_contains(rank_warnings$classes, "dkge_cv_rank_warning")

  joint <- suppressWarnings(dkge_cv_kernel_rank(
    isotropic_B, isotropic_X,
    K_grid = list(identity = diag(q), singular = fixture$K),
    ranks = 1:q, top_k = 2L, w_method = "none"
  ))
  expect_equal(joint$pick$kernel, "identity")
  expect_equal(joint$pick$rank, q)

  # A rank-deficient candidate can win the cheap alignment screen. It must not
  # consume the only top-k slot under the default full-rank policy.
  support_B <- replicate(fixture$S, fixture$K, simplify = FALSE)
  guarded_screen <- suppressWarnings(dkge_cv_kernel_rank(
    support_B, isotropic_X,
    K_grid = list(identity = diag(q), singular = fixture$K),
    ranks = 1:3, top_k = 1L, w_method = "none"
  ))
  expect_equal(guarded_screen$tables$alignment$kernel[[1L]], "singular")
  expect_equal(attr(guarded_screen$tables$alignment, "top"), "identity")
  expect_equal(guarded_screen$pick$kernel, "identity")

  saturation <- capture_dkge_warnings(dkge_cv_kernel_grid(
    isotropic_B, isotropic_X,
    K_grid = list(identity = diag(q), scaled_identity = 2 * diag(q)),
    rank = q, w_method = "none"
  ))
  expect_true(saturation$value$saturated)
  expect_contains(saturation$classes, "dkge_cv_saturation_warning")
})

test_that("contrast diagnostics fail on null targets and reveal query collisions", {
  fixture <- make_kernel_support_fixture()
  fit <- suppressWarnings(dkge_fit(
    fixture$B_list, fixture$X_list, fixture$K,
    rank = 4L, w_method = "none", effect_scaling = "none"
  ))

  pure_null <- c(0, 1, -1, 0, 0, 0, 0)
  expect_error(
    dkge_contrast(fit, pure_null, method = "loso", align = FALSE),
    class = "dkge_kernel_contrast_error"
  )

  supported <- c(0, 1, 1, 0, 0, 0, 0)
  supported_result <- expect_no_warning(
    dkge_contrast(fit, supported, method = "loso", align = FALSE)
  )
  expect_equal(supported_result$metadata$kernel_estimability$status, "estimable")
  expect_equal(supported_result$metadata$kernel_estimability$null_fraction, 0,
               tolerance = 1e-10)

  collapsed <- capture_dkge_warnings(dkge_contrast(
    fit,
    list(amplitude = c(0, 1, 0, 0, 0, 0, 0),
         pm_slope = c(0, 0, 1, 0, 0, 0, 0)),
    method = "loso", align = FALSE
  ))
  expect_contains(collapsed$classes, "dkge_kernel_contrast_warning")
  expect_contains(collapsed$classes, "dkge_kernel_contrast_collision_warning")
  expect_true(collapsed$value$metadata$kernel_query_pairs$collision)
  expect_equal(
    collapsed$value$metadata$kernel_query_pairs$query_correlation,
    1, tolerance = 1e-10
  )
  expect_equal(
    collapsed$value$values$amplitude,
    collapsed$value$values$pm_slope,
    tolerance = 1e-10
  )

  preflight <- dkge_contrast_diagnostics(
    fit,
    list(amplitude = c(0, 1, 0, 0, 0, 0, 0),
         pm_slope = c(0, 0, 1, 0, 0, 0, 0),
         difference = pure_null)
  )
  expect_equal(preflight$kernel$rank, 4L)
  expect_equal(preflight$estimability$status,
               c("partially_estimable", "partially_estimable", "null"))
  expect_true(preflight$pairs$collision[1L])
})

test_that("a zero-rank kernel fails before fitting", {
  fixture <- make_kernel_support_fixture()
  expect_error(
    dkge_fit(fixture$B_list, fixture$X_list, matrix(0, fixture$q, fixture$q),
             rank = 1L, w_method = "none"),
    class = "dkge_kernel_rank_error"
  )
})
