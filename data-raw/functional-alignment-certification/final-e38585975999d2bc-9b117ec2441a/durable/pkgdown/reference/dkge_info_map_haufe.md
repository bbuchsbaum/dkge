# Haufe-style encoding maps from latent classifiers

Converts discriminative weights into encoding (activation) patterns in
latent space using the Haufe transform \\a_z = \Sigma\_{zz} \beta\\,
then projects them to subject cluster space and through the renderer
pipeline.

## Usage

``` r
dkge_info_map_haufe(
  fit,
  clf,
  renderer,
  Z_by_subject = NULL,
  lambda = 0,
  to_vox = TRUE,
  inference = c("none", "signflip", "parametric")
)
```

## Arguments

- fit:

  Fitted `dkge` object.

- clf:

  Cross-fitted classifier returned by
  [`dkge_cv_train_latent_classifier()`](https://bbuchsbaum.github.io/dkge/reference/dkge_cv_train_latent_classifier.md).

- renderer:

  Renderer produced by
  [`dkge_build_renderer()`](https://bbuchsbaum.github.io/dkge/reference/dkge_build_renderer.md).

- Z_by_subject:

  Optional list of latent cluster features used to estimate fold
  covariances when they are not stored in `clf`.

- lambda:

  Non-negative smoothing parameter applied via
  [`dkge_anchor_aggregate()`](https://bbuchsbaum.github.io/dkge/reference/dkge_anchor_aggregate.md).

- to_vox:

  Logical; when `TRUE` and a decoder is available in `renderer`, a dense
  voxel map is produced.

- inference:

  Must be `"none"`. Legacy `"signflip"` and `"parametric"` modes are
  retired; use operator-bound aligned maps from
  [`dkge_transport_contrasts_to_reference()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_contrasts_to_reference.md)
  with
  [`dkge_infer_aligned()`](https://bbuchsbaum.github.io/dkge/reference/dkge_infer_aligned.md).

## Value

List mirroring
[`dkge_info_map_from_classifier()`](https://bbuchsbaum.github.io/dkge/reference/dkge_info_map_from_classifier.md)
with `meta$kind = "haufe"`.

## Examples

``` r
# \donttest{
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 6, P = 20, snr = 5
)
fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 2)
#> Warning: Argument 'kernel' is deprecated; use 'K' instead.
centroids <- lapply(toy$B_list, function(B) matrix(rnorm(ncol(B) * 3), ncol(B), 3))
renderer <- suppressWarnings(dkge_build_renderer(fit,
                                centroids = centroids,
                                anchor_xyz = matrix(rnorm(20 * 3), 20, 3),
                                anchor_n = 20,
                                anchor_method = "sample"))
y <- factor(rep(c("class1", "class2"), length.out = length(fit$Btil)))
folds <- dkge_define_folds(fit, type = "custom",
                           assignments = list(c(1, 2), c(3, 4), c(5, 6)))
clf <- dkge_cv_train_latent_classifier(fit, y, folds = folds)
enc <- dkge_info_map_haufe(fit, clf, renderer, inference = "none", to_vox = FALSE)
length(enc$mean_anchor)
#> [1] 20
# }
```
