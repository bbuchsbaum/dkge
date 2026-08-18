library(testthat)

# test-target.R
# Subject-by-feature target construction (dkge_make_target).

make_target_fit <- function(S = 3, q = 3, P = 4, T = 12, seed = 42) {
  set.seed(seed)
  betas <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
  designs <- replicate(S, {
    X <- matrix(rnorm(T * q), T, q)
    qr.Q(qr(X))
  }, simplify = FALSE)
  dkge_fit(betas, designs, K = diag(q), rank = 2)
}

test_that("transported_maps keeps every contrast when there are no centroids", {
  fit <- make_target_fit()
  con <- dkge_contrast(fit, list(c1 = c(1, -1, 0), c2 = c(0, 1, -1)),
                       method = "loso", align = FALSE)

  target <- dkge_make_target(fit, type = "transported_maps", contrast_obj = con)

  # Contrasts 2..k used to be silently discarded on this path.
  expect_equal(ncol(target$Y), 2L * length(con$values[[1]][[1]]))
  expect_equal(target$feature_ids,
               c(paste0("c1:", 1:4), paste0("c2:", 1:4)))
  expect_equal(target$subject_ids, fit$subject_ids)

  # Each block holds that contrast's per-subject values, subjects in fit order.
  for (nm in c("c1", "c2")) {
    block <- target$Y[, startsWith(colnames(target$Y), paste0(nm, ":")),
                      drop = FALSE]
    expected <- do.call(rbind, con$values[[nm]])
    expect_equal(unname(block), unname(expected), tolerance = 1e-12,
                 info = nm)
  }
})

test_that("transported_maps feature ids are stable across calls", {
  fit <- make_target_fit()
  con <- dkge_contrast(fit, list(c1 = c(1, -1, 0), c2 = c(0, 1, -1)),
                       method = "loso", align = FALSE)
  first <- dkge_make_target(fit, type = "transported_maps", contrast_obj = con)
  second <- dkge_make_target(fit, type = "transported_maps", contrast_obj = con)
  expect_identical(first$feature_ids, second$feature_ids)
  expect_identical(first$subject_ids, second$subject_ids)
  expect_equal(first$Y, second$Y)

  # Unnamed contrasts fall back to positional contrast names, still stable.
  con_unnamed <- dkge_contrast(fit, list(c(1, -1, 0), c(0, 1, -1)),
                               method = "loso", align = FALSE)
  target <- dkge_make_target(fit, type = "transported_maps",
                             contrast_obj = con_unnamed)
  expect_equal(target$feature_ids,
               c(paste0("contrast1:", 1:4), paste0("contrast2:", 1:4)))
})

test_that("a single contrast keeps its own prefix", {
  fit <- make_target_fit()
  con <- dkge_contrast(fit, list(only = c(1, 0, -1)), method = "loso",
                       align = FALSE)
  target <- dkge_make_target(fit, type = "transported_maps", contrast_obj = con)
  expect_equal(target$feature_ids, paste0("only:", 1:4))
  expect_equal(target$type, "transported_maps")
  expect_equal(target$feature_space, "transported")
})
