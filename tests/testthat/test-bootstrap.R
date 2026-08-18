library(testthat)

set.seed(100)
S <- 3
q <- 3
P <- 4
T <- 40

betas <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
designs <- replicate(S, qr.Q(qr(matrix(rnorm(T * q), T, q))), simplify = FALSE)
centroids <- replicate(S, matrix(runif(P * 3), P, 3), simplify = FALSE)

fit <- dkge(betas, designs = designs, K = diag(q), rank = 2)
fit$centroids <- centroids

transport_loadings <- dkge_transport_loadings_to_medoid(fit,
                                                        medoid = 1,
                                                        centroids = centroids,
                                                        mapper = dkge_mapper_spec("sinkhorn", epsilon = 0.05))
subject_maps <- transport_loadings$subjects[[1]]
values_medoid <- lapply(seq_len(nrow(subject_maps)), function(i) subject_maps[i, ])

cache <- transport_loadings$cache

contrasts <- c(1, -1, 0)

vox_map <- diag(ncol(subject_maps))



test_that("projection bootstrap returns expected shapes", {
  boot <- dkge_bootstrap_projected(values_medoid, B = 20, voxel_operator = vox_map)
  expect_equal(length(boot$medoid$mean), ncol(subject_maps))
  expect_equal(dim(boot$medoid$boot), c(20, ncol(subject_maps)))
  expect_equal(dim(boot$voxel$boot), c(20, ncol(subject_maps)))
})


test_that("q-space bootstrap runs with cached transport", {
  boot_q <- dkge_bootstrap_qspace(fit,
                                  contrasts = contrasts,
                                  B = 10,
                                  seed = 123,
                                  transport_cache = cache,
                                  medoid = 1,
                                  voxel_operator = vox_map,
                                  scheme = "poisson")
  expect_equal(boot_q$B, 10)
  expect_equal(length(boot_q$summary), 1)
  expect_equal(ncol(boot_q$summary[[1]]$boot_medoid), ncol(subject_maps))
  expect_true(boot_q$summary[[1]]$medoid$sd[1] >= 0)
})


test_that("analytic bootstrap falls back gracefully", {
  boot_a <- dkge_bootstrap_analytic(fit,
                                    contrasts = contrasts,
                                    B = 8,
                                    seed = 321,
                                    transport_cache = cache,
                                    medoid = 1,
                                    scheme = "bayes",
                                    voxel_operator = vox_map,
                                    perturb_tol = 0.5)
  expect_equal(boot_a$B, 8)
  expect_true(boot_a$fallbacks >= 0)
  expect_equal(ncol(boot_a$summary[[1]]$boot_medoid), ncol(subject_maps))
})


test_that("linear pooling reproduces the re-pooled Chat exactly", {
  set.seed(11)
  Sl <- 4; ql <- 5; Pl <- 8; Tl <- 30
  betas_l <- replicate(Sl, matrix(rnorm(ql * Pl), ql, Pl), simplify = FALSE)
  designs_l <- replicate(Sl, qr.Q(qr(matrix(rnorm(Tl * ql), Tl, ql))), simplify = FALSE)
  fit_l <- dkge(betas_l, designs_l, K = diag(ql), rank = 2)

  # Default fits pool linearly, so the O(1) matvec branch is the live one.
  expect_false(dkge:::.dkge_fit_pool_is_nonlinear(fit_l))

  contrib_matrix <- vapply(fit_l$contribs, as.numeric, numeric(ql * ql))
  expect_equal(matrix(contrib_matrix %*% as.numeric(fit_l$weights), ql, ql),
               unname(fit_l$Chat), tolerance = 1e-12)

  set.seed(99)
  for (i in seq_len(5)) {
    xi <- stats::rpois(Sl, lambda = 1)
    linear <- matrix(contrib_matrix %*% (as.numeric(fit_l$weights) * xi), ql, ql)
    repooled <- dkge:::.dkge_repool_fit(fit_l, sample_weights = xi)
    expect_equal(linear, unname(repooled$Chat), tolerance = 1e-12)
  }
})


test_that("non-linear pooling without stored moments is an error, not a silent fallback", {
  # `.dkge_repool_fit()` returns NULL when the fit carries no q-space effect
  # moments. Falling through to the linear matvec there would report a
  # differently-normalised bootstrap as if it matched the fit.
  fit_bad <- fit
  fit_bad$missingness <- "rescale"
  fit_bad$effect_moments <- NULL
  expect_true(dkge:::.dkge_fit_pool_is_nonlinear(fit_bad))
  expect_null(dkge:::.dkge_repool_fit(fit_bad))

  expect_error(
    dkge_bootstrap_qspace(fit_bad,
                          contrasts = contrasts,
                          B = 2,
                          seed = 7,
                          transport_cache = cache,
                          medoid = 1,
                          scheme = "poisson"),
    "does not carry the q-space effect moments"
  )
})


test_that("q-space bootstrap fast path matches the re-pooling path", {
  boot_fast <- dkge_bootstrap_qspace(fit,
                                     contrasts = contrasts,
                                     B = 6,
                                     seed = 4242,
                                     transport_cache = cache,
                                     medoid = 1,
                                     scheme = "poisson")

  # Force the re-pooling branch and confirm the results are identical.
  boot_repool <- with_mocked_bindings(
    dkge_bootstrap_qspace(fit,
                          contrasts = contrasts,
                          B = 6,
                          seed = 4242,
                          transport_cache = cache,
                          medoid = 1,
                          scheme = "poisson"),
    .dkge_fit_pool_is_nonlinear = function(fit) TRUE
  )

  expect_equal(boot_fast$summary[[1]]$boot, boot_repool$summary[[1]]$boot,
               tolerance = 1e-12)
  expect_equal(boot_fast$summary[[1]]$mean, boot_repool$summary[[1]]$mean,
               tolerance = 1e-12)
  expect_equal(boot_fast$summary[[1]]$sd, boot_repool$summary[[1]]$sd,
               tolerance = 1e-12)
})
