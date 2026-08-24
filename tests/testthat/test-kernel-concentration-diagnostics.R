library(testthat)

make_concentrated_kernel_fit <- function(gamma = 0.98, seed = 12L) {
  set.seed(seed)
  effects <- c("amplitude", "pm_slope")
  K <- matrix(c(1, gamma, gamma, 1), 2L, 2L,
              dimnames = list(effects, effects))
  B_list <- replicate(5L, {
    B <- matrix(rnorm(2L * 20L), 2L, 20L)
    rownames(B) <- effects
    B
  }, simplify = FALSE)
  X <- diag(2L)
  dimnames(X) <- list(effects, effects)
  X_list <- replicate(5L, X, simplify = FALSE)
  fit <- suppressWarnings(dkge_fit(
    B_list, X_list, K,
    rank = 2L, w_method = "none", effect_scaling = "none"
  ))
  list(fit = fit, K = K, B_list = B_list, X_list = X_list)
}

capture_warning_classes <- function(expr) {
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

test_that("kernel diagnostics expose scale-invariant spectral concentration", {
  K <- diag(c(1, rep(0.01, 5L)))
  expected_pr <- sum(diag(K))^2 / sum(diag(K)^2)

  small <- kernel_roots(1e-12 * K)
  large <- kernel_roots(1e12 * K)

  expect_equal(small$rank, 6L)
  expect_equal(small$effective_rank_pr, expected_pr, tolerance = 1e-10)
  expect_equal(small$effective_rank_fraction, expected_pr / 6,
               tolerance = 1e-10)
  expect_equal(small$leading_eigenvalue_share, 1 / 1.05,
               tolerance = 1e-10)
  expect_equal(large$effective_rank_pr, small$effective_rank_pr,
               tolerance = 1e-10)
  expect_equal(large$effective_rank_fraction, small$effective_rank_fraction,
               tolerance = 1e-10)
  expect_equal(large$leading_eigenvalue_share,
               small$leading_eigenvalue_share, tolerance = 1e-10)
})

test_that("practical query collinearity has its own default tolerance", {
  fixture <- make_concentrated_kernel_fit()
  contrasts <- list(
    amplitude = c(1, 0),
    pm_slope = c(0, 1)
  )

  diagnostics <- dkge_contrast_diagnostics(fixture$fit, contrasts)
  strict <- dkge_contrast_diagnostics(
    fixture$fit, contrasts, collinearity_tol = 1e-4
  )

  expect_equal(diagnostics$estimability$status, rep("estimable", 2L))
  expect_equal(strict$estimability, diagnostics$estimability)
  expect_equal(diagnostics$pairs$input_correlation, 0, tolerance = 1e-12)
  expect_equal(diagnostics$pairs$query_correlation, 0.98, tolerance = 1e-10)
  expect_true(diagnostics$pairs$collision)
  expect_false(strict$pairs$collision)
  expect_equal(diagnostics$summary$n_collisions, 1L)
  expect_equal(diagnostics$summary$max_abs_query_correlation, 0.98,
               tolerance = 1e-10)
  expect_equal(unname(diagnostics$summary$max_query_pair),
               c("amplitude", "pm_slope"))
  expect_equal(diagnostics$kernel$effective_rank_pr,
               sum(eigen(fixture$K, symmetric = TRUE)$values)^2 /
                 sum(eigen(fixture$K, symmetric = TRUE)$values^2),
               tolerance = 1e-10)

  warned <- capture_warning_classes(dkge_contrast(
    fixture$fit, contrasts, method = "loso", align = FALSE
  ))
  expect_contains(warned$classes, "dkge_kernel_contrast_collision_warning")
  expect_equal(warned$value$metadata$kernel_query_summary$n_collisions, 1L)

  disabled <- capture_warning_classes(dkge_contrast(
    fixture$fit, contrasts, method = "loso", align = FALSE,
    collinearity_tol = NULL
  ))
  expect_false("dkge_kernel_contrast_collision_warning" %in% disabled$classes)
})

test_that("CV surfaces a predictive winner with concentrated kernel spectrum", {
  q <- 6L
  S <- 5L
  B_list <- lapply(seq_len(S), function(s) {
    B <- matrix(0, q, 2L)
    B[1L, 1L] <- 1
    B[s + 1L, 2L] <- 3
    B
  })
  X_list <- replicate(S, diag(q), simplify = FALSE)
  K_concentrated <- diag(c(1, rep(0.01, q - 1L)))

  warned <- capture_warning_classes(dkge_cv_kernel_grid(
    B_list, X_list,
    K_grid = list(identity = diag(q), concentrated = K_concentrated),
    rank = 2L, w_method = "none"
  ))
  cv <- warned$value
  selected <- cv$table[cv$table$kernel == "concentrated", , drop = FALSE]

  expect_equal(cv$pick, "concentrated")
  expect_equal(cv$table$mean, c(0, 0.1), tolerance = 1e-12)
  expect_equal(selected$kernel_rank, q)
  expect_lt(selected$kernel_effective_rank_pr, 2)
  expect_lt(selected$kernel_effective_rank_to_rank, 0.75)
  expect_gt(selected$kernel_leading_eigenvalue_share, 0.95)
  expect_true(cv$kernel_concentration$concentrated)
  expect_contains(warned$classes, "dkge_cv_kernel_concentration_warning")

  disabled <- capture_warning_classes(dkge_cv_kernel_grid(
    B_list, X_list,
    K_grid = list(identity = diag(q), concentrated = K_concentrated),
    rank = 2L, w_method = "none",
    kernel_concentration_threshold = NULL
  ))
  expect_false("dkge_cv_kernel_concentration_warning" %in% disabled$classes)
  expect_false(disabled$value$kernel_concentration$concentrated)
})

test_that("component contrasts isolate coordinates while dual saliences do not", {
  fixture <- make_concentrated_kernel_fit()
  fit <- fixture$fit

  component_contrasts <- dkge_component_contrasts(fit)
  alpha <- crossprod(
    fit$U,
    fit$K %*% backsolve(fit$R, component_contrasts, transpose = FALSE)
  )
  dual_alpha <- crossprod(
    fit$U,
    fit$K %*% backsolve(fit$R, fit$K %*% fit$U, transpose = FALSE)
  )

  expect_equal(dim(component_contrasts), c(2L, 2L))
  expect_equal(rownames(component_contrasts), fit$effects)
  expect_equal(colnames(component_contrasts), c("LV1", "LV2"))
  expect_equal(alpha, diag(2L), tolerance = 1e-10)
  expect_gt(max(abs(dual_alpha - diag(2L))), 0.9)

  first <- dkge_component_contrasts(fit, comps = 1L)
  expect_equal(dim(first), c(2L, 1L))
  expect_equal(colnames(first), "LV1")

  rank_seven <- structure(
    list(U = diag(7L), R = diag(7L), effects = paste0("effect", 1:7)),
    class = "dkge"
  )
  expect_equal(dim(dkge_component_contrasts(rank_seven)), c(7L, 7L))
})

test_that("new concentration thresholds fail early and readably", {
  fixture <- make_concentrated_kernel_fit()
  contrasts <- list(amplitude = c(1, 0), pm_slope = c(0, 1))
  expect_error(
    dkge_contrast_diagnostics(fixture$fit, contrasts, collinearity_tol = 1),
    "collinearity_tol"
  )

  expect_error(
    dkge_cv_kernel_grid(
      fixture$B_list, fixture$X_list,
      K_grid = list(identity = diag(2L)), rank = 1L,
      kernel_concentration_threshold = 0
    ),
    "kernel_concentration_threshold"
  )
})
