# Model-Level Spatial Regularization in DKGE

Suppose each subject contributes a $`q \times P_s`$ GLM beta matrix: the
rows are the same $`q`$ effects, while the columns are that subject’s
voxels or parcels. DKGE solves its eigensystem in the small
$`q`$-dimensional effect space. That does **not** mean spatial
regularization must be a cosmetic step after fitting. A spatial operator
applied to the beta columns changes the effect-space moment and
therefore can change the learned basis itself.

This article follows one estimand:

> Learn the rank-one effect direction that best explains spatial beta
> fields after a prespecified graph-Laplacian regularization, then
> obtain held-out subject fields under that same regularized model.

The expected outputs are a $`q \times 1`$ group basis, one regularized
loading field per subject, a held-out score for each candidate penalty,
and LOSO contrast fields that inherit the selected operator.

``` text
q x P_s beta fields -> sparse spatial solve -> q x q pooled moment -> group basis -> held-out fields
```

## Where does the regularizer enter?

Let $`L_s`$ be a combinatorial graph Laplacian on subject $`s`$’s
spatial units. DKGE uses the Tikhonov resolvent

``` math
H_{s,\lambda} = (I + \lambda L_s)^{-1}, \qquad \lambda \ge 0,
```

and replaces the design-scaled beta block $`\widetilde B_s`$ by
$`\widetilde B_s H_{s,\lambda}^{\mathsf T}`$. With a spatial metric
$`\Omega_s`$, the subject contribution becomes

``` math
M_s(\lambda) =
\widetilde B_s H_{s,\lambda}^{\mathsf T}\,
\Omega_s\,
H_{s,\lambda}\widetilde B_s^{\mathsf T}.
```

Note that $`H_{s,\lambda}`$ appears **twice**: the resolvent smooths the
beta block, and the moment is then formed from the smoothed block, so
the effective moment-level smoothing is the squared resolvent rather
than $`H_{s,\lambda}`$ itself. Under an identity spatial metric this is
$`\widetilde B_s H_{s,\lambda}^2 \widetilde B_s^{\mathsf T}`$.

If you also supply voxel weights via
[`dkge_weights()`](https://bbuchsbaum.github.io/dkge/reference/dkge_weights.md),
be aware that the two per-unit weightings sit on opposite sides of the
smoother: voxel weights multiply the beta columns *before* smoothing (so
a zero-weight unit becomes a hole that its neighbours diffuse into),
while $`\Omega_s`$ acts *after* it, as a metric on the smoothed field.
To exclude a unit completely, remove the matching beta column, graph
node, domain label, and every aligned per-unit input together. Changing
only `coords` would violate the fitted spatial-domain contract.

Thus $`\lambda`$ changes the pooled $`q \times q`$ moment that DKGE
diagonalizes. The same $`H_{s,\lambda}`$ is then used for components,
cross-fitted contrasts, bootstrap inputs, transport inputs, and
prediction. It is never necessary to form a dense inverse: DKGE stores a
sparse Cholesky factorization and solves against beta fields as needed.
At $`\lambda=0`$, the operation is exactly the identity.

The graph and the design kernel play different roles. $`L_s`$ says which
**spatial units** should vary together. $`K`$ says which **experimental
effects** should vary together. Neither can substitute for the other.

## Build an example with a rough distractor

The example has two effects at 40 ordered locations. The first effect
varies smoothly; the second alternates at every location and has larger
amplitude. An unsmoothed rank-one fit therefore prefers the rough field.

``` r

S <- 6L
P <- 40L
effects <- c("smooth", "rough")
locations <- sprintf("location%02d", seq_len(P))

coords <- cbind(x = seq_len(P) - 1, y = 0, z = 0)
rownames(coords) <- locations

smooth_truth <- sin(seq(0, 2 * pi, length.out = P))
rough_truth <- 1.5 * (-1)^seq_len(P)

B_list <- setNames(lapply(seq_len(S), function(s) {
  B <- rbind(
    smooth = smooth_truth + rnorm(P, sd = 0.08),
    rough = rough_truth + rnorm(P, sd = 0.08)
  )
  colnames(B) <- locations
  B
}), sprintf("sub%02d", seq_len(S)))

X0 <- diag(2)
colnames(X0) <- effects
X_list <- setNames(replicate(S, X0, simplify = FALSE), names(B_list))

K <- diag(2)
dimnames(K) <- list(effects, effects)

c(beta_rows = nrow(B_list[[1]]), spatial_units = ncol(B_list[[1]]))
#>     beta_rows spatial_units 
#>             2            40
```

## Construct the graph with `adjoin`

[`dkge_spatial_regularizer()`](https://bbuchsbaum.github.io/dkge/reference/dkge_spatial_regularizer.md)
delegates coordinate-based graph construction to
[`adjoin::spatial_laplacian()`](https://rdrr.io/pkg/adjoin/man/spatial_laplacian.html).
DKGE then validates the Laplacian, binds it to the beta-column domain,
and records the graph recipe and `adjoin` version. Here a distance just
above one connects adjacent locations on the line.

``` r

spatial_demo <- dkge_spatial_regularizer(
  coords = coords,
  lambda = 2,
  dthresh = 1.01,
  nnk = 3,
  weight_mode = "binary",
  normalized = FALSE,
  handle_isolates = "keep_zero"
)
spatial_demo
#> <dkge_spatial_regularizer>
#>   source   : adjoin 
#>   domains  : shared 
#>   lambda   : 2 
#>   status   : active 
#>   units    : 40 
#>   edges    : 39
```

`dthresh` is a distance cutoff **in the units of `coords`**, and the
defaults (`dthresh = 1.42`, `nnk = 27`) assume unit-spaced voxel
indices, where 1.42 picks up the face and edge neighbours of a 3x3x3
neighbourhood. Real coordinates are usually in millimetres, where that
threshold connects nothing: every unit is isolated, `L` is the zero
matrix, and the resolvent is the identity, so the fit is bit-identical
to an unregularized one even though `lambda` is positive. DKGE warns
when a graph ends up with no edges, and the printed object reports the
edge count, which is the fastest way to confirm the graph is what you
meant:

``` r

mm_coords <- cbind(x = seq(0, 33, by = 3), y = 0, z = 0)

# Default threshold on millimetre coordinates: warns, and smooths nothing.
dkge_spatial_regularizer(coords = mm_coords, lambda = 2)
#> Warning: The spatial graph has no edges, so the resolvent is the identity and
#> positive lambda = 2 will have no effect. `dthresh` is measured in the units of
#> `coords`; the default assumes unit-spaced indices. With this graph dthresh =
#> 1.42, so millimetre coordinates would need a threshold on the millimetre scale
#> (e.g. slightly above the voxel spacing).
#> <dkge_spatial_regularizer>
#>   source   : adjoin 
#>   domains  : shared 
#>   lambda   : 2 
#>   status   : inert 
#>   units    : 12 
#>   edges    : 0 
#>   NOTE     : edgeless graphs give H = I; lambda has no effect

# A threshold on the millimetre scale connects the chain.
dkge_spatial_regularizer(coords = mm_coords, lambda = 2,
                         dthresh = 3.1, nnk = 6)
#> <dkge_spatial_regularizer>
#>   source   : adjoin 
#>   domains  : shared 
#>   lambda   : 2 
#>   status   : active 
#>   units    : 12 
#>   edges    : 11
```

DKGE uses four spatial-effectiveness states:

| Status     | Meaning                                                   |
|------------|-----------------------------------------------------------|
| `inactive` | `lambda = 0`; the explicit identity model                 |
| `inert`    | positive `lambda`, but every graph is edgeless            |
| `partial`  | positive `lambda`, with only some subject graphs edgeless |
| `active`   | positive `lambda`, and every graph can smooth             |

Direct construction and fitting warn for `inert` or `partial`
specifications. Spatial CV is stricter: if any positive candidate is
requested but every graph is edgeless, the candidate scores cannot
identify `lambda`, so
[`dkge_cv_spatial_grid()`](https://bbuchsbaum.github.io/dkge/reference/dkge_cv_spatial_grid.md)
stops before fitting folds. A partially edgeless CV warns once with the
affected domain identifiers and continues.

Because `coords` is one matrix, this is a shared spatial domain: every
subject must have the same number and ordering of columns. The matching
location names make that contract explicit. If subjects have different
parcellations, supply a named list of coordinate matrices or precomputed
Laplacians instead; DKGE will bind one operator to each subject and fail
on missing or mismatched domains. Names on that list are matched to the
fitted subject IDs. For raw lists, declare them explicitly with
`dkge_data(B_list, X_list, subject_ids = names(B_list))` before fitting.

The combinatorial Laplacian is intentional. It preserves constant
fields, so DKGE currently requires `normalized = FALSE`. For a surface
mesh or a scientifically defined adjacency, construct the Laplacian in
`adjoin` and pass it through `laplacian =` rather than approximating the
topology with Euclidean distance.

## Verify that smoothing changes the solution

We first compare an unsmoothed fit with a deliberately visible penalty.
Both use the same identity design kernel and disable subject weighting,
so the only changed assumption is the spatial regularizer.

``` r

fit_raw <- dkge_fit(
  B_list, X_list, K,
  rank = 1,
  w_method = "none",
  effect_scaling = "none"
)

fit_demo <- dkge_fit(
  B_list, X_list, K,
  rank = 1,
  w_method = "none",
  effect_scaling = "none",
  spatial = spatial_demo
)

effect_comparison <- data.frame(
  effect = effects,
  raw = abs(fit_raw$U[, 1]),
  regularized = abs(fit_demo$U[, 1])
)
transform(effect_comparison,
          raw = round(raw, 3), regularized = round(regularized, 3))
#>   effect   raw regularized
#> 1 smooth 0.002       0.999
#> 2  rough 1.000       0.036
```

Component signs are arbitrary, hence the absolute values. The raw
component is almost entirely the high-amplitude rough effect. Under the
regularized moment, the dominant direction rotates toward the smooth
effect. This is the key distinction from post-hoc image smoothing: the
fitted $`U`$ changed.

The reconstructed subject field uses the same operator. For a compact
spatial diagnostic, define graph-path roughness as the sum of squared
adjacent differences.

``` r

roughness <- function(x) sum(diff(x)^2)

map_raw <- dkge_project_btil(fit_raw, fit_raw$Btil[[1]])[, 1]
map_demo <- dkge_project_btil(
  fit_demo, fit_demo$Btil[[1]], subject = 1
)[, 1]

map_raw <- map_raw * sign(stats::cor(map_raw, rough_truth))
map_demo <- map_demo * sign(stats::cor(map_demo, smooth_truth))

round(c(raw = roughness(map_raw), regularized = roughness(map_demo)), 3)
#>         raw regularized 
#>     353.960       0.411
```

![Two line plots: the unsmoothed rank-one loading alternates sharply at
adjacent locations, while the regularized loading varies
smoothly.](dkge-spatial-regularization_files/figure-html/plot-spatial-fields-1.png)

A deliberately strong spatial penalty rotates the rank-one solution away
from the alternating field and yields a smooth loading field. Component
signs were aligned only for display.

The fitted object retains auditable domain and graph diagnostics.
Identical fingerprints here are expected because all subjects share one
ordered graph.

``` r

spatial_diagnostics <- dkge_diagnostics(fit_demo)$spatial
spatial_diagnostics[c("status", "requested", "effective", "fully_effective")]
#> $status
#> [1] "active"
#> 
#> $requested
#> [1] TRUE
#> 
#> $effective
#> [1] TRUE
#> 
#> $fully_effective
#> [1] TRUE
spatial_diagnostics$diagnostics[
  , c("subject", "n_units", "n_edges", "lambda", "effective", "domain_mode")
]
#>   subject n_units n_edges lambda effective   domain_mode
#> 1   sub01      40      39      2      TRUE label_matched
#> 2   sub02      40      39      2      TRUE label_matched
#> 3   sub03      40      39      2      TRUE label_matched
#> 4   sub04      40      39      2      TRUE label_matched
#> 5   sub05      40      39      2      TRUE label_matched
#> 6   sub06      40      39      2      TRUE label_matched
```

## Select lambda on raw held-out fields

The visible $`\lambda=2`$ comparison was chosen to explain the
mechanism, not to justify that penalty. Use
[`dkge_cv_spatial_grid()`](https://bbuchsbaum.github.io/dkge/reference/dkge_cv_spatial_grid.md)
when $`\lambda`$ is to be selected from the data. Keep the graph, design
kernel, rank, and validation geometry fixed, and include zero as the
unsmoothed candidate.

``` r

cv_spatial <- dkge_cv_spatial_grid(
  B_list, X_list, K,
  spatial = spatial_demo,
  lambdas = c(0, 0.1, 0.5, 2),
  rank = 1,
  w_method = "none",
  effect_scaling = "none"
)

transform(
  cv_spatial$table[, c("lambda", "mean", "se", "admissible")],
  mean = round(mean, 4),
  se = round(se, 4)
)
#>   lambda   mean     se admissible
#> 1    0.0 0.8249 0.0014       TRUE
#> 2    0.1 0.8248 0.0015       TRUE
#> 3    0.5 0.1766 0.0014       TRUE
#> 4    2.0 0.1760 0.0014       TRUE
cv_spatial$pick
#> [1] 0.1
```

Each candidate basis is trained without the held-out subject, but all
candidates are scored against the **unsmoothed** held-out beta block in
one fixed effect-space geometry. Prespecified effect scaling, location
weights, and `Omega_list` remain in that fixed geometry; only the
candidate Laplacian solve is omitted. A candidate therefore cannot
improve its own score by shrinking the field in both the numerator and
denominator. `best` is the largest mean score; `pick` is the largest,
smoothest $`\lambda`$ within one standard error of that best score. In
this example the selected value is 0.1, much weaker than the explanatory
$`\lambda=2`$.

As with kernel selection, a saturated warning means the held-out
criterion cannot discriminate the candidates well. Also remember that
$`\lambda`$’s scale depends on the graph weights: compare penalties only
while holding the graph construction fixed. If graph topology or edge
weights are themselves learned from the outcome, they require their own
nested validation rather than being tuned on the same folds.

Refit the complete cohort with the selected specification returned by
CV:

``` r

fit_selected <- dkge_fit(
  B_list, X_list, K,
  rank = 1,
  w_method = "none",
  effect_scaling = "none",
  spatial = cv_spatial$spatial
)
dkge_diagnostics(fit_selected)$spatial[c("active", "lambda", "shared")]
#> $active
#> [1] TRUE
#> 
#> $lambda
#> [1] 0.1
#> 
#> $shared
#> [1] TRUE
```

## Carry the operator into LOSO inference

Do not smooth a final $`t`$- or $`p`$-map. Instead, construct held-out
subject fields under the fitted spatial model and perform inference on
those fields.
[`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md)
does this automatically because the operator is part of the fit.

``` r

rough_contrast <- c(smooth = 0, rough = 1)
loso_raw <- dkge_contrast(fit_raw, rough_contrast, method = "loso")
loso_selected <- dkge_contrast(
  fit_selected, rough_contrast, method = "loso"
)

loso_roughness <- data.frame(
  model = c("raw", "selected"),
  median_roughness = c(
    median(vapply(loso_raw$values[[1]], roughness, numeric(1))),
    median(vapply(loso_selected$values[[1]], roughness, numeric(1)))
  )
)
transform(loso_roughness,
          median_roughness = round(median_roughness, 2))
#>      model median_roughness
#> 1      raw           357.70
#> 2 selected           184.03
```

Transport to an identified subject support, anchor field, or voxel grid
happens **after** these regularized subject fields are constructed. It
aligns spatial domains; it does not redefine the fitted effect-space
solution.

## Keep the spatial controls distinct

Several DKGE arguments mention space, but they do not estimate the same
thing.

| Control | What it changes | Changes fitted basis? |
|----|----|---:|
| `spatial = dkge_spatial_regularizer(...)` | Correlation across neighboring beta columns through $`H_{s,\lambda}`$ | Yes |
| `Omega_list` or adaptive location weights | Spatial mass, precision, or reliability in the subject moment | Yes |
| transport mapper and masses | Alignment of already estimated subject fields to a common domain | No |
| `dkge_anchor_aggregate(L=, lambda=)` | Optional smoothing of an already transported anchor average | No |

The operations can be combined, but order matters: DKGE applies adaptive
column scaling first, the Laplacian solve second, and `Omega_list` as
the final spatial metric in the moment. See
[`vignette("dkge-weighting")`](https://bbuchsbaum.github.io/dkge/articles/dkge-weighting.md)
for the complete weighting pipeline.

Two further boundaries are deliberate:

- Active spatial regularization is not available with
  `debias = "analytic"`. The current analytic contract stores diagonal
  residual variances, while $`H`$ induces spatial noise covariance that
  those variances do not identify. DKGE fails closed; use
  `debias = "split_half"` or `"none"`.
- A shared graph is automatically reused for prediction. A fit with
  subject-specific graphs requires a new `spatial` specification aligned
  to the prediction subjects; smoothing cannot be added only at
  prediction time.
- [`dkge_fit_from_kernels()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit_from_kernels.md)
  and anchor-input fits do not expose physical beta columns: their
  synthetic factor columns are algebraic, not spatial. They reject
  model-level `spatial` regularization rather than treating those
  factors as neighboring locations.

## Reporting checklist

For a reproducible analysis, report:

1.  the spatial domain and whether the graph was shared or
    subject-specific;
2.  how the graph was built (`dthresh`, `nnk`, edge weights, isolate
    policy, or the provenance of a precomputed Laplacian);
3.  the candidate $`\lambda`$ grid, validation geometry, selection rule,
    and whether scores saturated;
4.  the selected $`\lambda`$, latent rank, design kernel, `Omega_list`,
    and other weighting choices;
5.  that inference used regularized held-out subject fields, with no
    direct smoothing of the final inferential statistic.

Continue with
[`vignette("dkge-contrasts-inference")`](https://bbuchsbaum.github.io/dkge/articles/dkge-contrasts-inference.md)
for the inferential estimand,
[`vignette("dkge-weighting")`](https://bbuchsbaum.github.io/dkge/articles/dkge-weighting.md)
for the neighboring weighting controls, or
[`vignette("dkge-dense-rendering")`](https://bbuchsbaum.github.io/dkge/articles/dkge-dense-rendering.md)
for transport and rendering after the regularized subject fields have
been estimated.
