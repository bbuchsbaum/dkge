library(testthat)

set.seed(100)
S <- 3
q <- 3
P <- 4
T <- 40

betas <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
designs <- replicate(S, qr.Q(qr(matrix(rnorm(T * q), T, q))), simplify = FALSE)
centroids <- replicate(S, matrix(runif(P * 3), P, 3), simplify = FALSE)

fit <- dkge(betas, designs = designs, kernel = diag(q), rank = 2)
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
