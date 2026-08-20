# Executable contracts for transport semantics, diagnostics, and extension APIs.

library(testthat)

test_that("joint plans and application operators obey distinct conservation laws", {
  mu <- c(0.6, 0.4)
  nu <- c(0.2, 0.3, 0.5)
  plan <- rbind(c(0.1, 0.2, 0.3),
                c(0.1, 0.1, 0.2))

  intensive <- dkge:::.dkge_transport_operator(plan, mu, nu, "intensive")
  extensive <- dkge:::.dkge_transport_operator(plan, mu, nu, "extensive")

  expect_equal(rowSums(plan), mu, tolerance = 1e-15)
  expect_equal(colSums(plan), nu, tolerance = 1e-15)
  expect_equal(colSums(intensive), rep(1, 3), tolerance = 1e-15)
  expect_equal(as.numeric(t(intensive) %*% rep(7, 2)), rep(7, 3),
               tolerance = 1e-14)
  expect_equal(rowSums(extensive), rep(1, 2), tolerance = 1e-15)
  expect_equal(sum(t(extensive) %*% c(2, 5)), 7, tolerance = 1e-14)
})

test_that("medoid transport preserves constants or total mass end to end", {
  A_list <- list(
    diag(2),
    rbind(c(1, 0), c(0, 1), c(1, 1))
  )
  centroids <- list(
    rbind(c(0, 0, 0), c(2, 0, 0)),
    rbind(c(0, 0, 0), c(2, 0, 0), c(1, 1, 0))
  )
  sizes <- list(c(2, 1), c(1, 2, 3))

  field <- dkge_transport_to_medoid_sinkhorn(
    list(rep(7, 2), rep(7, 3)), A_list, centroids, sizes = sizes,
    medoid = 2, value_type = "intensive", warm_start = FALSE
  )
  expect_equal(field$subj_values, matrix(7, 2, 3), tolerance = 2e-3)
  expect_equal(field$operators[[2]], diag(3))
  expect_equal(diag(field$plans[[2]]), sizes[[2]] / sum(sizes[[2]]))

  totals <- dkge_transport_to_medoid_sinkhorn(
    list(c(2, 5), c(1, 2, 4)), A_list, centroids, sizes = sizes,
    medoid = 2, value_type = "extensive", warm_start = FALSE
  )
  expect_equal(sum(totals$subj_values[1, ]), 7, tolerance = 1e-10)
  expect_equal(rowSums(totals$operators[[1]]), rep(1, 2), tolerance = 1e-10)
  expect_equal(totals$subj_values[2, ], c(1, 2, 4))
})

test_that("Sinkhorn cache keys include every cost entry and reuse exact solves", {
  dkge_clear_sinkhorn_cache()
  on.exit(dkge_clear_sinkhorn_cache(), add = TRUE)
  env <- dkge:::.dkge_sinkhorn_cache
  C1 <- matrix(seq_len(64) / 64, 8, 8)
  C2 <- matrix(rev(seq_len(64)) / 64, 8, 8)
  mu <- nu <- rep(1 / 8, 8)

  first <- dkge:::.dkge_sinkhorn_plan(
    C1, mu, nu, epsilon = 1, max_iter = 2000, tol = 1e-8,
    return_diagnostics = TRUE
  )
  second <- dkge:::.dkge_sinkhorn_plan(
    C2, mu, nu, epsilon = 1, max_iter = 2000, tol = 1e-8,
    return_diagnostics = TRUE
  )
  expect_false(first$diagnostics$cache_hit)
  expect_false(second$diagnostics$cache_hit)
  expect_equal(length(setdiff(ls(env, all.names = TRUE), ".order")), 2L)

  repeated <- dkge:::.dkge_sinkhorn_plan(
    C1, mu, nu, epsilon = 1, max_iter = 2000, tol = 1e-8,
    return_diagnostics = TRUE
  )
  expect_true(repeated$diagnostics$cache_hit)
  expect_equal(repeated$diagnostics$iterations, 0L)
  expect_identical(repeated$plan, first$plan)
})

test_that("non-converged Sinkhorn state is diagnosed and never cached", {
  dkge_clear_sinkhorn_cache()
  on.exit(dkge_clear_sinkhorn_cache(), add = TRUE)
  set.seed(817)
  C <- matrix(sample(seq(0, 100, length.out = 64)), 8, 8)
  mu <- nu <- rep(1 / 8, 8)

  expect_warning(
    solved <- dkge:::.dkge_sinkhorn_plan(
      C, mu, nu, epsilon = 1e-4, max_iter = 1L, tol = 1e-12,
      return_diagnostics = TRUE
    ),
    "did not converge"
  )
  expect_false(solved$diagnostics$converged)
  expect_gt(solved$diagnostics$marginal_error, solved$diagnostics$tolerance)
  expect_equal(length(setdiff(ls(dkge:::.dkge_sinkhorn_cache,
                               all.names = TRUE), ".order")), 0L)
})

test_that("default Sinkhorn budget converges on a deterministic representative problem", {
  set.seed(20260820)
  unit_rows <- function(X) X / sqrt(rowSums(X^2))
  source_feat <- unit_rows(matrix(rnorm(12 * 4), 12, 4))
  target_feat <- unit_rows(matrix(rnorm(15 * 4), 15, 4))
  source_xyz <- matrix(rnorm(12 * 3, sd = 20), 12, 3)
  target_xyz <- matrix(rnorm(15 * 3, sd = 20), 15, 3)

  mapping <- fit_mapper(
    dkge_mapper_spec("sinkhorn"),
    source_feat = source_feat,
    target_feat = target_feat,
    source_xyz = source_xyz,
    target_xyz = target_xyz
  )
  diagnostics <- mapping$fit_info$diagnostics
  expect_true(diagnostics$converged)
  expect_lte(diagnostics$marginal_error, diagnostics$tolerance)
  expect_lte(diagnostics$iterations, 5000L)
})

test_that("transport method aliases are honest and unsupported kNN fails early", {
  centroids <- list(matrix(0, 2, 3))
  expect_warning(alias <- dkge_transport_spec(centroids, method = "sinkhorn_cpp"),
                 "deprecated")
  expect_equal(alias$method, "sinkhorn")
  expect_equal(dkge_transport_spec(centroids, method = "ridge")$method, "ridge")
  expect_error(dkge_transport_spec(centroids, method = "knn"), "arg")
  expect_equal(dkge:::.dkge_resolve_mapper_spec(NULL, "ols", list())$strategy,
               "ols")
})

test_that("transported-map targets honor the method stored in their spec", {
  loadings <- list(diag(2), matrix(c(1, 0, 0, 1, 1, 1), 3, 2, byrow = TRUE))
  centroids <- list(matrix(0, 2, 3), matrix(0, 3, 3))
  values <- list(c(2, 4), c(1, 3, 5))
  fit <- structure(list(subject_ids = c("s1", "s2")), class = "dkge")
  transport <- dkge_transport_spec(centroids, medoid = 1, method = "ridge")

  target <- dkge_make_target(
    fit = fit,
    type = "transported_maps",
    values = values,
    loadings = loadings,
    transport = transport
  )
  expect_equal(target$provenance$transport$mapper$strategy, "ridge")
  expect_equal(dim(target$Y), c(2, 2))
})

test_that("dense mapper constructor supports extension-defined S3 backends", {
  fit_name <- "fit_mapper.dkge_mapper_contractprobe"
  apply_name <- "apply_mapper.dkge_mapper_fit_contractprobe"
  assign(fit_name, function(spec, ...) {
    scale <- if (is.null(spec$pars$scale)) 1 else spec$pars$scale
    structure(list(scale = scale),
              class = "dkge_mapper_fit_contractprobe")
  }, envir = .GlobalEnv)
  assign(apply_name, function(fitted_mapper, values, ...) {
    fitted_mapper$scale * values
  }, envir = .GlobalEnv)
  on.exit(rm(list = c(fit_name, apply_name), envir = .GlobalEnv), add = TRUE)

  spec <- dkge_mapper("contractprobe", scale = 3)
  expect_s3_class(spec, "dkge_mapper_contractprobe")
  fitted <- fit_mapper(spec)
  expect_equal(apply_mapper(fitted, c(2, -1)), c(6, -3))
  expect_error(dkge_mapper("bad-name"), "identifier")
})
