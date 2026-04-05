# DKGE Workflow

Design-Kernel Group Embedding (DKGE) provides a principled approach for
summarizing subject-level general linear models (GLMs) into a shared
latent space that respects the underlying structure in the design
matrix. This vignette provides a comprehensive guided tour that takes
you from raw subject beta coefficients through to component-level
statistical inference and voxelwise visualization of group-level
patterns.

## Why DKGE?

Multi-subject neuroimaging studies typically generate thousands of
per-subject beta maps that are all tied to the same experimental design
structure. Traditional group analysis approaches rely on coordinate-wise
averaging or mass-univariate statistical models, both of which
unfortunately discard the rich smoothness and dependency patterns that
are encoded within the experimental design itself.

DKGE addresses these limitations through a three-pronged approach.
First, it harmonizes subject-level GLM outputs so that effects across
different subjects share a common measurement scale. Second, it rotates
the multiblock beta arrays into orthogonal components by leveraging a
design kernel that explicitly encodes similarity relationships between
different experimental effects. Finally, it exposes a comprehensive
suite of utilities for statistical inference, contrast testing, and
optimal-transport-based map alignment that can operate across both
subjects and spatial locations.

## Pipeline At A Glance

The DKGE workflow follows a logical sequence of steps that build upon
each other to produce interpretable group-level results:

1.  **Data preparation**: Begin by organizing per-subject beta
    coefficients, design matrices, and optional spatial weights using
    [`dkge_subject()`](https://bbuchsbaum.github.io/dkge/reference/dkge_subject.md)
    and
    [`dkge_data()`](https://bbuchsbaum.github.io/dkge/reference/dkge_data.md)
    functions.

2.  **Kernel specification**: Define a design kernel `K` that captures
    the similarity relationships among your experimental design effects.

3.  **Model fitting**: Execute the core DKGE algorithm by calling
    [`dkge()`](https://bbuchsbaum.github.io/dkge/reference/dkge.md) (or
    the lower-level
    [`dkge_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit.md))
    to estimate the shared latent basis that spans across subjects.

4.  **Component inspection**: Examine the fitted model by inspecting
    subject-specific weights, component eigenvalues, and the resulting
    loadings structure.

5.  **Consensus mapping and inference**: Transport individual component
    maps into a consensus space and perform statistical inference using
    [`dkge_component_stats()`](https://bbuchsbaum.github.io/dkge/reference/dkge_component_stats.md)
    or specialized contrast utilities.

6.  **Voxelwise visualization** (optional): Interpolate the consensus
    maps back to voxel space using
    [`dkge_transport_to_voxels()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_to_voxels.md)
    to enable detailed spatial visualization.

## Required Ingredients

To successfully apply DKGE, you’ll need to prepare several key data
components:

- **Beta coefficient matrices**: These are q × P_s matrices containing
  GLM coefficients, where q represents the number of experimental
  effects and P_s represents the number of spatial clusters or voxels
  for subject s.

- **Design matrices**: Each subject contributes a T_s × q matrix of
  regressors, where T_s is the number of time points for that subject.
  The column names in these matrices serve to define the canonical
  labels for each experimental effect.

- **Design kernel `K`**: This is a positive semi-definite matrix that
  operates over the experimental effects. You can use `diag(q)` as a
  neutral starting choice, or construct more sophisticated structured
  kernels using the
  [`design_kernel()`](https://bbuchsbaum.github.io/dkge/reference/design_kernel.md)
  function.

- **Spatial weights** (optional): These can be provided as vectors or
  matrices that describe the reliability or effective size of each
  spatial cluster or voxel, allowing the algorithm to appropriately
  weight different brain regions.

- **Spatial centroids** (optional but recommended): These are 3-D
  coordinate locations for each cluster, which become essential when
  transporting statistical values across different subjects or when
  mapping results back to voxel space.

## Walkthrough: Simulated Study

To illustrate the DKGE workflow in practice, we’ll work through a
complete example using simulated data. The following code generates a
minimal dataset consisting of three subjects, three experimental design
effects, and four spatial clusters per subject. In your own analyses,
you would naturally replace these simulated components with real beta
coefficients, design matrices, and spatial centroids derived from your
actual neuroimaging study.

``` r
library(dkge)
set.seed(2024)

S <- 3        # subjects
q <- 3        # effects per design
P <- 4        # clusters per subject
T <- 30       # time points per subject design

betas <- replicate(S, matrix(rnorm(q * P), q, P), simplify = FALSE)
designs <- replicate(S, {
  X <- matrix(rnorm(T * q), T, q)
  qr.Q(qr(X))  # orthonormalise columns to stabilise the fit
}, simplify = FALSE)

centroids <- replicate(S, matrix(runif(P * 3, min = -15, max = 15), P, 3), simplify = FALSE)

subjects <- lapply(seq_len(S), function(s) {
  dkge_subject(betas[[s]], design = designs[[s]], id = paste0("sub", s))
})

data_bundle <- dkge_data(subjects)
identity_kernel <- diag(data_bundle$q)
```

The
[`dkge_subject()`](https://bbuchsbaum.github.io/dkge/reference/dkge_subject.md)
function serves to package each beta matrix together with its
corresponding design matrix and any optional weighting information.
Meanwhile,
[`dkge_data()`](https://bbuchsbaum.github.io/dkge/reference/dkge_data.md)
takes care of harmonizing the effect ordering across subjects and
standardizing subject identifiers throughout the dataset. When working
with real neuroimaging studies, you would typically pass parcellated
cluster sizes through the `omega` argument, which allows the algorithm
to appropriately weight different subjects based on the spatial extent
of their data during the fitting process.

### Fitting the DKGE model

Now we can proceed to fit a rank-2 DKGE model using an identity kernel,
which means we’re not imposing any smoothing relationships between
different experimental effects. The object returned by the fitting
process inherits from the `multiblock_biprojector` class, which ensures
that it integrates seamlessly with the broader ecosystem of
`multivarious` projection helper functions.

``` r
fit <- dkge(data_bundle, K = identity_kernel, rank = 2)
fit$centroids <- centroids  # convenience: attach for transport helpers; pass explicitly in production

round(fit$sdev, 3)
#> [1] 7.58 6.44
fit$weights
#> [1] 1.330 0.649 1.021
```

The fitted model object contains several key components that summarize
the results. The `fit$sdev` element reports the singular values
associated with each of the retained components, providing insight into
their relative importance. The `fit$weights` element displays the
per-subject block weights that result from any MFA-style scaling that
was applied during fitting. Finally, the `fit$U` and `fit$v` elements
contain the component loadings in effect space and cluster space,
respectively, which encode how the original variables relate to the
discovered latent components.

### Inspecting component structure

The eigenvalues stored in `fit$evals` provide crucial information about
how much variance each discovered component explains within the K-metric
space defined by your design kernel. To gain deeper insight into the
structure of these components, the following code demonstrates how to
extract and examine the leading component loadings for each individual
subject.

``` r
loadings <- lapply(fit$Btil, function(Bts) {
  t(Bts) %*% fit$K %*% fit$U
})

lapply(loadings, function(A) round(A[, 1, drop = FALSE], 3))
#> [[1]]
#>             [,1]
#> cluster_1 -1.473
#> cluster_2 -0.483
#> cluster_3 -1.128
#> cluster_4  3.494
#> 
#> [[2]]
#>            [,1]
#> cluster_1 -3.43
#> cluster_2  1.72
#> cluster_3  4.82
#> cluster_4 -1.83
#> 
#> [[3]]
#>            [,1]
#> cluster_1 1.728
#> cluster_2 2.217
#> cluster_3 1.179
#> cluster_4 0.048
```

The resulting `loadings` matrices contain the cluster-wise contributions
for each discovered component, revealing how different spatial regions
contribute to the overall component structure. When you need to apply
this fitted model to new data, you can leverage either
[`dkge_project_block()`](https://bbuchsbaum.github.io/dkge/reference/dkge_project_block.md)
or
[`dkge_project_clusters()`](https://bbuchsbaum.github.io/dkge/reference/dkge_project_clusters.md)
to score additional datasets against the established basis.

### Component-level inference and medoid transport

The
[`dkge_component_stats()`](https://bbuchsbaum.github.io/dkge/reference/dkge_component_stats.md)
function performs a crucial step in the analysis pipeline by
transporting component loadings from individual subject spaces to a
shared reference parcellation (known as the medoid). Once this transport
is complete, the function computes statistical summaries for each
cluster in the reference space. Different mapper strategies offer
various ways to balance spatial proximity against embedding distance
relationships. In this example, we employ Sinkhorn optimal transport
with tight entropic regularization to achieve a robust alignment.

``` r
comp <- dkge_component_stats(
  fit,
  mapper = list(strategy = "sinkhorn", epsilon = 0.05, lambda_spa = 0.5),
  centroids = centroids,
  inference = list(type = "parametric"),
  medoid = 1L
)

head(comp$summary)
#>   component cluster   stat     p p_adj significant
#> 1         1       1 -1.459 0.282 0.563       FALSE
#> 2         1       2  0.695 0.559 0.563       FALSE
#> 3         2       1 -0.686 0.563 0.563       FALSE
#> 4         2       2  1.020 0.415 0.563       FALSE
```

The returned object is a comprehensive list that provides multiple
perspectives on the statistical results:

- `summary`: A tidy data frame containing test statistics, raw p-values,
  and false discovery rate (FDR) adjustments for easy interpretation and
  downstream analysis.
- `statistics`: The raw component-wise statistic vectors in their
  original form before any tidying operations were applied.
- `transport`: Subject-level matrices with dimensions (subjects × medoid
  clusters) that capture the aligned data after the optimal transport
  procedure.

### Mapping components to voxel space

When your analysis includes voxel-level features such as spatial
centroids or PCA-reduced spatial signatures, the
[`dkge_transport_to_voxels()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_to_voxels.md)
function provides a powerful mechanism for interpolating the transported
component values onto a full voxel grid. This capability enables
high-resolution visualization of your results without forcing you to
commit to any single parcellation scheme, thereby preserving flexibility
in how you present and interpret your findings.

``` r
voxels <- replicate(S, matrix(runif(10 * 3, min = -20, max = 20), 10, 3), simplify = FALSE)
component_values <- lapply(loadings, function(A) A[, 1])

voxel_maps <- dkge_transport_to_voxels(
  fit,
  values = component_values,
  voxels = voxels,
  mapper = "ridge"
)

round(voxel_maps$value, 3)
```

The resulting object provides both group-level and individual-level
perspectives on the interpolated results. The `voxel_maps$value` element
contains the averaged consensus map that synthesizes information across
all subjects, while `voxel_maps$subj_values` preserves the
subject-specific interpolations for cases where individual variation is
of particular interest.

## Bootstrap Shortcuts

When conducting bootstrap-based statistical inference, it’s crucial to
ensure that every resampled statistic exists within the same reference
coordinate system. By caching and re-using a fixed optimal transport
mapping, we can guarantee this consistency across bootstrap iterations.
The following approach demonstrates how to cache the subject-to-medoid
transport operators once and then efficiently feed them to the
specialized bootstrap helper functions:

``` r
cache <- dkge_prepare_transport(fit, centroids = centroids, medoid = 1)
values_medoid <- lapply(seq_len(nrow(comp$transport[[1]])),
                        function(i) comp$transport[[1]][i, ])
boot_proj <- dkge_bootstrap_projected(values_medoid, B = 200, seed = 99)
boot_q <- dkge_bootstrap_qspace(fit, contrasts = c(1, -1, 0), B = 200,
                                transport_cache = cache, medoid = 1, seed = 99)
boot_proj$medoid$sd[1]
#> [1] 0.366
boot_q$summary[[1]]$medoid$sd[1]
#> NULL
```

These two bootstrap approaches capture different sources of variability
in your statistical inference. The
[`dkge_bootstrap_projected()`](https://bbuchsbaum.github.io/dkge/reference/dkge_bootstrap_projected.md)
function resamples the already-transported subject vectors directly,
focusing on the variability that remains after the group-level latent
space has been fixed. In contrast,
[`dkge_bootstrap_qspace()`](https://bbuchsbaum.github.io/dkge/reference/dkge_bootstrap_qspace.md)
incorporates additional uncertainty by including variability that arises
from re-estimating the group-level subspace itself on each bootstrap
iteration.

For computationally demanding scenarios where thousands of bootstrap
replicates are required, the
[`dkge_bootstrap_analytic()`](https://bbuchsbaum.github.io/dkge/reference/dkge_bootstrap_analytic.md)
function offers an efficient alternative. This approach applies the same
cached transport mapping but updates the subspace using fast first-order
perturbation theory, only falling back to the more expensive exact
eigendecomposition when the perturbation becomes too large to trust the
linear approximation.

## Where to next?

Having completed this foundational walkthrough, you’re now equipped to
explore more advanced DKGE capabilities:

- **Hypothesis testing**: Leverage
  [`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md),
  [`dkge_loso_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_loso_contrast.md),
  or the comprehensive
  [`dkge_pipeline()`](https://bbuchsbaum.github.io/dkge/reference/dkge_pipeline.md)
  function to perform sophisticated hypothesis-driven contrasts and
  cross-validation procedures.

- **Out-of-sample prediction**: Apply
  [`dkge_project_blocks()`](https://bbuchsbaum.github.io/dkge/reference/dkge_project_blocks.md)
  to score held-out subjects or evaluate the model performance on
  entirely new parcellation schemes.

- **Advanced visualization**: Combine the voxel transport capabilities
  with
  [`dkge_paint_medoid_map()`](https://bbuchsbaum.github.io/dkge/reference/dkge_paint_medoid_map.md)
  (available in the companion `dkge.neuroim2` package) or integrate with
  external visualization frameworks to create publication-ready
  group-level activation maps.

- **Structured kernels**: Explore the rich possibilities offered by
  [`design_kernel()`](https://bbuchsbaum.github.io/dkge/reference/design_kernel.md)
  to encode smoothness relationships across ordinal or circular
  experimental factors, or to capture complex multi-factor interaction
  patterns in your design.

For comprehensive details on function arguments and additional utility
functions, consult the complete function reference documentation
available in the `man/` directory.
