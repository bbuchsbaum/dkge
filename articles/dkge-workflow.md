# DKGE Workflow

DKGE compresses per-subject GLM beta matrices into a shared latent space
that respects the experimental design. This vignette walks through the
core pipeline: data prep → fit → component inspection → inference →
bootstrap.

## Data preparation

Three subjects, three design effects, four spatial clusters each:

``` r
S <- 3; q <- 3; P <- 4; T <- 30

betas   <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
designs <- replicate(S, {
  X <- matrix(rnorm(T * q), T, q)
  qr.Q(qr(X))   # orthonormalise for a clean compressed covariance
}, simplify = FALSE)
centroids <- replicate(S, matrix(runif(P * 3, -15, 15), P, 3), simplify = FALSE)

subjects    <- lapply(seq_len(S), function(s)
  dkge_subject(betas[[s]], design = designs[[s]], id = paste0("sub", s)))
data_bundle <- dkge_data(subjects)
```

[`dkge_subject()`](https://bbuchsbaum.github.io/dkge/reference/dkge_subject.md)
pairs each beta matrix with its design matrix.
[`dkge_data()`](https://bbuchsbaum.github.io/dkge/reference/dkge_data.md)
harmonises effect ordering and subject IDs across the list.

## Fitting

``` r
fit <- dkge(data_bundle, K = diag(data_bundle$q), rank = 2)
fit$centroids <- centroids  # convenience: attach for transport helpers; pass explicitly in production

round(fit$sdev, 3)   # component singular values
#> [1] 7.576 6.439
fit$weights          # per-subject MFA block weights
#> [1] 1.329986 0.648805 1.021209
```

The identity kernel `diag(q)` applies no design smoothing — use
[`design_kernel()`](https://bbuchsbaum.github.io/dkge/reference/design_kernel.md)
for structured factorial designs.

## Component structure

Cluster-wise loadings on the first component for each subject:

``` r
loadings <- lapply(fit$Btil, function(Bts) t(Bts) %*% fit$K %*% fit$U)
lapply(loadings, function(A) round(A[, 1, drop = FALSE], 3))
#> [[1]]
#>             [,1]
#> cluster_1 -1.473
#> cluster_2 -0.483
#> cluster_3 -1.128
#> cluster_4  3.494
#> 
#> [[2]]
#>             [,1]
#> cluster_1 -3.427
#> cluster_2  1.721
#> cluster_3  4.820
#> cluster_4 -1.828
#> 
#> [[3]]
#>            [,1]
#> cluster_1 1.728
#> cluster_2 2.217
#> cluster_3 1.179
#> cluster_4 0.048
```

Project new subjects onto the fixed basis with
[`dkge_project_block()`](https://bbuchsbaum.github.io/dkge/reference/dkge_project_block.md)
or
[`dkge_project_clusters()`](https://bbuchsbaum.github.io/dkge/reference/dkge_project_clusters.md).

## Component-level inference

[`dkge_component_stats()`](https://bbuchsbaum.github.io/dkge/reference/dkge_component_stats.md)
transports loadings to a shared medoid parcellation and returns test
statistics with FDR adjustment:

``` r
comp <- dkge_component_stats(
  fit,
  mapper     = list(strategy = "sinkhorn", epsilon = 0.05, lambda_spa = 0.5),
  centroids  = centroids,
  inference  = list(type = "parametric"),
  medoid     = 1L
)
head(comp$summary)
#>   component cluster       stat         p     p_adj significant
#> 1         1       1 -1.4591757 0.2819165 0.5634129       FALSE
#> 2         1       2  0.6948130 0.5590391 0.5634129       FALSE
#> 3         2       1 -0.6862882 0.5634129 0.5634129       FALSE
#> 4         2       2  1.0197970 0.4151049 0.5634129       FALSE
```

`comp$transport` holds subject × medoid matrices after alignment;
`comp$statistics` holds raw component vectors.

## Mapping to voxel space

``` r
voxels          <- replicate(S, matrix(runif(10 * 3, -20, 20), 10, 3), simplify = FALSE)
component_values <- lapply(loadings, function(A) A[, 1])

voxel_maps <- dkge_transport_to_voxels(fit,
  values = component_values, voxels = voxels, mapper = "ridge")
round(voxel_maps$value, 3)
```

`voxel_maps$value` is the group-average consensus map;
`voxel_maps$subj_values` preserves per-subject interpolations.

## Bootstrap inference

Cache the transport operators once; reuse them across all bootstrap
draws:

``` r
cache <- dkge_prepare_transport(fit, centroids = centroids, medoid = 1)
values_medoid <- lapply(seq_len(nrow(comp$transport[[1]])),
                        function(i) comp$transport[[1]][i, ])

# Resample already-transported subject vectors
boot_proj <- dkge_bootstrap_projected(values_medoid, B = 200, seed = 99)

# Include subspace re-estimation uncertainty
boot_q <- dkge_bootstrap_qspace(fit, contrasts = c(1, -1, 0), B = 200,
                                transport_cache = cache, medoid = 1, seed = 99)

boot_proj$medoid$sd[1]
#> [1] 0.3662273
boot_q$summary[[1]]$medoid$sd[1]
#> NULL
```

[`dkge_bootstrap_analytic()`](https://bbuchsbaum.github.io/dkge/reference/dkge_bootstrap_analytic.md)
offers a cheaper alternative via first-order perturbation theory when
thousands of draws are needed.

## Where next?

| Task                      | Function                                                                                                                                                                                       |
|---------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Hypothesis contrasts      | [`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md), [`dkge_pipeline()`](https://bbuchsbaum.github.io/dkge/reference/dkge_pipeline.md)                           |
| Out-of-sample scoring     | [`dkge_project_blocks()`](https://bbuchsbaum.github.io/dkge/reference/dkge_project_blocks.md)                                                                                                  |
| Structured design kernels | [`design_kernel()`](https://bbuchsbaum.github.io/dkge/reference/design_kernel.md)                                                                                                              |
| Classification / decoding | [`dkge_classify()`](https://bbuchsbaum.github.io/dkge/reference/dkge_classify.md) — see [`vignette("dkge-classification")`](https://bbuchsbaum.github.io/dkge/articles/dkge-classification.md) |
| Inference details         | [`vignette("dkge-contrasts-inference")`](https://bbuchsbaum.github.io/dkge/articles/dkge-contrasts-inference.md)                                                                               |
