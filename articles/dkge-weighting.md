# Weighting Strategies in DKGE

DKGE has separate controls for locations within a subject, model-level
spatial regularization, whole subjects, and spatial transport. Each
changes a different estimand. This vignette shows where the controls
enter the computation and what must be reported when they are used.

Weighting is needed when “one unit, one vote” is scientifically
wrong—for example, when parcels have unequal area or some measurements
are demonstrably less reliable. It is not a generic way to make a result
cleaner. Every weight answers two questions: *what is the unit being
reweighted, and at what stage?*

``` text
subject beta block --location scaling--> Laplacian solve --> subject moment --subject weights--> group basis --transport masses--> common map
```

### Adaptive voxel weighting (preview)

DKGE can also rescale locations adaptively before fold bases are
constructed.
[`dkge_weights()`](https://bbuchsbaum.github.io/dkge/reference/dkge_weights.md)
combines a prespecified prior, such as an external reliability map, with
training-fold statistics such as design-kernel energy or precision.
Training-fold scope closes one leakage path but does not establish
calibration or improved sensitivity.
[`vignette("dkge-adaptive-weighting")`](https://bbuchsbaum.github.io/dkge/articles/dkge-adaptive-weighting.md)
shows the specification, structural checks, and post-hoc updates with
[`dkge_update_weights()`](https://bbuchsbaum.github.io/dkge/reference/dkge_update_weights.md).

## Which weighting control changes which object?

| Layer | Unit receiving weight or coupling | Argument | Enters before | Typical rationale |
|----|----|----|----|----|
| Location metric | parcel, cluster, or voxel within subject | `Omega_list` / `omega` | subject moment | area, precision, or a declared spatial metric |
| Spatial regularizer | graph edge joining neighboring locations | `spatial` | subject moment and every reconstructed field | penalize spatial roughness inside the model |
| Subject | entire subject block | `w_method`, `w_tau`, or `weights` | pooled eigensolve | prevent one high-energy block from dominating, or apply external subject mass |
| Transport | source and target spatial mass | `sizes` | common-map summary | preserve parcel area or voxel density during alignment |

Effect-by-subject precision weights for unequal trial counts are a
fourth, separate layer covered in
[`vignette("dkge-partial-effect-spaces")`](https://bbuchsbaum.github.io/dkge/articles/dkge-partial-effect-spaces.md),
with the precision-weighted mean derived in
[`vignette("dkge-unbalanced-trialwise")`](https://bbuchsbaum.github.io/dkge/articles/dkge-unbalanced-trialwise.md).
Adaptive training-fold weights are covered in
[`vignette("dkge-adaptive-weighting")`](https://bbuchsbaum.github.io/dkge/articles/dkge-adaptive-weighting.md).

`spatial = dkge_spatial_regularizer(...)` is an operator rather than
another diagonal weight. It couples neighboring columns through
`(I + lambda * L)^{-1}`, changes the fitted group basis, and is reused
for components, held-out contrasts, and prediction. DKGE applies
adaptive column scaling first, this sparse spatial solve second, and
`Omega_list` as the final metric in the subject moment. See
[`vignette("dkge-spatial-regularization")`](https://bbuchsbaum.github.io/dkge/articles/dkge-spatial-regularization.md)
for graph construction and lambda selection. By contrast, smoothing in
[`dkge_anchor_aggregate()`](https://bbuchsbaum.github.io/dkge/reference/dkge_anchor_aggregate.md)
is a post-fit display/aggregation operation and cannot regularize the
learned solution.

The example below changes one layer at a time so that a changed output
can be attributed to a named assumption.

## Example Dataset

The dataset below has three subjects with different cluster counts,
sizes, and reliability scores. That heterogeneity is what the weighting
controls act on.

``` r

library(dkge)
library(purrr)
set.seed(42)

S <- 3
q <- 4
P <- c(6, 5, 7)
T <- 40

betas <- map(P, ~ matrix(rnorm(q * .x), q, .x))
designs <- map(seq_len(S), function(s) {
  X <- matrix(rnorm(T * q), T, q)
  qr.Q(qr(X))
})

centroids <- map(P, ~ matrix(runif(.x * 3, min = -20, max = 20), .x, 3))
cluster_sizes <- map(P, ~ sample(10:40, .x, replace = TRUE))
reliability_scores <- list(
  rep(0.95, P[1]),
  rep(0.20, P[2]),  # subject 2 is the deliberately low-reliability block
  rep(0.80, P[3])
)
```

`cluster_sizes` stands in for surface area, and `reliability_scores` for
a per-cluster quality metric such as split-half correlation. Both are
supplied by you; DKGE does not compute either.

## Baseline Fit (No Spatial Weights)

``` r

bundle <- dkge_data(betas, designs = designs)
k_identity <- diag(bundle$q)

fit_equal <- dkge(bundle, K = k_identity, rank = 2, w_method = "none")
fit_equal$weights
#> [1] 1 1 1
```

With equal weights every subject contributes uniformly to the shared
basis, so cluster-level variation reflects the planted signal rather
than any weighting choice.

## Size-Weighted Blocks

Per-cluster mass values emphasise larger parcels while keeping a
diagonal information structure, so clusters stay independent of one
another.

``` r

omega_size <- cluster_sizes

fit_size <- dkge(bundle,
                 K = k_identity,
                 Omega_list = omega_size,
                 rank = 2,
                 w_method = "mfa_sigma1",
                 w_tau = 0.25)
fit_size$weights
#> [1] 0.732 1.218 1.049
```

The subject blocks now differ in total weighted energy, because size
weighting scales each block by its own mass. `w_method` decides whether
that difference survives into the basis.

``` r

size_loadings <- lapply(fit_size$Btil, function(Bts) t(Bts) %*% fit_size$K %*% fit_size$U)
map(size_loadings, ~ round(.x[, 1, drop = FALSE], 3))
#> [[1]]
#>             [,1]
#> cluster_1  0.661
#> cluster_2  2.439
#> cluster_3  1.978
#> cluster_4 -1.359
#> cluster_5 -5.871
#> cluster_6 -1.961
#> 
#> [[2]]
#>             [,1]
#> cluster_1  1.120
#> cluster_2  0.248
#> cluster_3  1.681
#> cluster_4 -4.505
#> cluster_5  1.341
#> 
#> [[3]]
#>             [,1]
#> cluster_1 -2.295
#> cluster_2  0.966
#> cluster_3  1.208
#> cluster_4 -4.167
#> cluster_5  0.178
#> cluster_6  0.420
#> cluster_7 -0.617
```

## Reliability Weighting

Reliability scores enter the same `omega` slot as size weights: supply a
per- cluster vector and the fit downweights clusters you have flagged as
noisy.

``` r

omega_reliability <- map2(cluster_sizes, reliability_scores, `*`)

fit_reliability <- dkge(bundle,
                        K = k_identity,
                        Omega_list = omega_reliability,
                        rank = 2,
                        w_method = "mfa_sigma1",
                        w_tau = 0.25)

weight_comparison <- rbind(
  size = fit_size$weights,
  size_x_reliability = fit_reliability$weights
)
round(weight_comparison, 3)
#>                     [,1] [,2]  [,3]
#> size               0.732 1.22 1.049
#> size_x_reliability 0.430 1.97 0.604

weighted_contribution_norm <- function(fit) {
  fit$weights * vapply(
    fit$contribs,
    function(x) sqrt(sum(x * x)),
    numeric(1)
  )
}
leverage_comparison <- rbind(
  size = weighted_contribution_norm(fit_size),
  size_x_reliability = weighted_contribution_norm(fit_reliability)
)
round(leverage_comparison, 1)
#>                    [,1] [,2] [,3]
#> size               1538 1198 1233
#> size_x_reliability  858  387  568
```

![Grouped bar chart showing how size-only and size-times-reliability
weighting change each subject's contribution
scale.](dkge-weighting_files/figure-html/reliability-leverage-plot-1.png)

The two visible rows are deliberately different. Subject 2 has the
lowest reliability scores. Its scalar MFA weight increases because
inverse-leading- singular-value scaling partly compensates for a
lower-energy block; the scalar weight alone is therefore not a measure
of leverage. The second table reports the direct contribution scale
entering `Chat` (subject weight times the Frobenius norm of that
subject’s unweighted contribution). On this controlled example, adding
reliability reduces subject 2’s contribution relative to the size-only
fit.

In real data, a reliability vector might come from inverse first-level
variance, split-half or test-retest stability, a prespecified
quality-control measure, or resampling stability. The source matters: a
weight estimated from the same held-out outcome can leak information,
while an external or training-fold estimate has a clearer
interpretation. When combining area and reliability, multiply only when
the resulting estimand is intended, then rescale to a transparent
reference such as mean one.

Keep weight values positive, and add a small ridge such as `+ 1e-6` when
a weight construction can produce zeros.

## Spatial Covariance (Smoothness) Weighting

With cluster coordinates you can replace diagonal weights with a smooth
covariance matrix. The example below builds a Gaussian affinity matrix,
which rewards neighboring parcels for activating together.

``` r

gaussian_cov <- function(coords, bandwidth = 15) {
  D <- as.matrix(dist(coords))
  K <- exp(-(D / bandwidth)^2)
  (K + t(K)) / 2 + diag(1e-6, nrow(K))
}

omega_smooth <- lapply(centroids, gaussian_cov)

fit_smooth <- dkge(bundle,
                   K = k_identity,
                   Omega_list = omega_smooth,
                   rank = 2,
                   w_method = "none")
fit_smooth$weights
#> [1] 1 1 1
```

These matrices define a spatial metric: DKGE right-multiplies each beta
block by the symmetric square root of `Omega`. That is weighting or
smoothing, not whitening (which would use an inverse square root). A
Gaussian affinity couples nearby parcels in the fitted energy
calculation. Whether that coupling is scientifically helpful is an
assumption to check against an uncoupled fit, not an automatic
correction for sharp boundaries.

## Subject-Level Scaling

`w_method` controls block-level reweighting, applied after the spatial
weights have acted on individual clusters. The default Multiple Factor
Analysis scaling shrinks dominant subjects by the inverse squared
leading singular value of their contribution.

``` r

fit_mfa <- dkge(bundle,
                K = k_identity,
                Omega_list = omega_reliability,
                rank = 2,
                w_method = "mfa_sigma1",
                w_tau = 0.25)

fit_energy <- dkge(bundle,
                   K = k_identity,
                   Omega_list = omega_reliability,
                   rank = 2,
                   w_method = "energy",
                   w_tau = 0)

rbind(mfa = fit_mfa$weights,
      energy = fit_energy$weights,
      size_mfa = fit_size$weights)
#>           [,1] [,2]  [,3]
#> mfa      0.430 1.97 0.604
#> energy   0.223 2.34 0.439
#> size_mfa 0.732 1.22 1.049
```

`w_tau` shrinks extreme weights toward equal subject contributions,
which matters most when one subject would otherwise dominate.

These are **pooled-moment weights**: they determine how strongly each
beta block contributes to the DKGE eigensolve. Cross-fitted folds
re-normalize the stored raw MFA/energy scores over the training subjects
and then reapply `w_tau`; a held-out subject therefore cannot change the
training fold’s weight scale. They are not automatically reused as
inverse-variance weights in a one-sample group test. Any group-level
subject weighting is a separate estimand. Current
[`dkge_infer_aligned()`](https://bbuchsbaum.github.io/dkge/reference/dkge_infer_aligned.md)
inference supports equal-subject weighting only. Explicit fixed weights
may be recorded for descriptive rendering and conditional
projected-bootstrap summaries, but they are not a calibrated weighted
group test.

## Weighting During Transport

Transport takes its own `sizes` argument. It mirrors `omega` but acts at
a different stage, *after* the components are estimated. Supplying
consistent sizes there keeps the medoid representation faithful to
parcel surface area.

``` r

fit_mfa$centroids <- centroids  # attach for transport helpers; pass explicitly in production

comp_equal <- suppressWarnings(dkge_component_stats(
  fit_mfa,
  mapper = dkge_mapper_spec("sinkhorn", epsilon = 0.05, max_iter = 2000),
  centroids = centroids,
  inference = NULL,
  medoid = 1L
))

comp_weighted <- suppressWarnings(dkge_component_stats(
  fit_mfa,
  mapper = dkge_mapper_spec("sinkhorn", epsilon = 0.05, max_iter = 2000),
  centroids = centroids,
  sizes = cluster_sizes,
  inference = NULL,
  medoid = 1L
))

data.frame(equal = head(comp_equal$summary$mean),
           size_weighted = head(comp_weighted$summary$mean))
#>    equal size_weighted
#> 1 -0.712        -0.723
#> 2 -1.662        -1.655
#> 3 -0.714        -0.742
#> 4  0.494         0.544
#> 5  4.677         4.478
#> 6  0.837         1.603
```

This deprecated component helper is descriptive only: it has no p-values
or significance fields. Supplying `sizes` changes the displayed
reference-support means because parcel mass changes the transport plan.
Group subject weighting remains equal unless fixed weights are declared
in a typed aligned-map object.

## What to check before reporting

Construct subjects with
[`dkge_subject()`](https://bbuchsbaum.github.io/dkge/reference/dkge_subject.md):
weights stay attached to the data through pipelines and streaming
loaders. In the legacy descriptive voxel mapper, pass cluster sizes to
[`dkge_transport_to_voxels()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_to_voxels.md)
via `sizes`, or volumetric weighting is lost. That helper learns from
full-fit loadings and does not produce typed inferential alignment; use
the functional-alignment workflow for group statistics.

Keep vectors positive and add a mild ridge such as `+ 1e-6` to
covariance matrices. Revisit each fit’s `weights`, singular values, and
eigenspectra whenever you change a weighting assumption; they show how
contributions are balanced across subjects.

These controls tailor an analysis to study-specific quality metrics or
prior spatial beliefs without changing the core fitting pipeline.

## Where to go next

- [`vignette("dkge-adaptive-weighting")`](https://bbuchsbaum.github.io/dkge/articles/dkge-adaptive-weighting.md)
  derives weights from the data instead of from a prior.
- [`vignette("dkge-partial-effect-spaces")`](https://bbuchsbaum.github.io/dkge/articles/dkge-partial-effect-spaces.md)
  covers effect-by-subject precision, which is a separate layer from
  these four.
- [`vignette("dkge-dense-rendering")`](https://bbuchsbaum.github.io/dkge/articles/dkge-dense-rendering.md)
  shows where transport weighting acts.
