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

  Fitted \`dkge\` object.

- clf:

  Cross-fitted classifier returned by
  \[dkge_cv_train_latent_classifier()\].

- renderer:

  Renderer produced by \[dkge_build_renderer()\].

- Z_by_subject:

  Optional list of latent cluster features used to estimate fold
  covariances when they are not stored in \`clf\`.

- lambda:

  Non-negative smoothing parameter applied via
  \[dkge_anchor_aggregate()\].

- to_vox:

  Logical; when \`TRUE\` and a decoder is available in \`renderer\`, a
  dense voxel map is produced.

- inference:

  One of \`"none"\`, \`"signflip"\`, or \`"parametric"\` specifying the
  group inference routine applied to anchor maps.

## Value

List mirroring \[dkge_info_map_from_classifier()\] with \`meta\$kind =
"haufe"\`.

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
renderer <- dkge_build_renderer(fit,
                                centroids = centroids,
                                anchor_xyz = matrix(rnorm(20 * 3), 20, 3),
                                anchor_n = 20,
                                anchor_method = "sample")
y <- factor(rep(c("class1", "class2"), length.out = length(fit$Btil)))
folds <- dkge_define_folds(fit, type = "custom",
                           assignments = list(c(1, 2), c(3, 4), c(5, 6)))
clf <- dkge_cv_train_latent_classifier(fit, y, folds = folds)
enc <- dkge_info_map_haufe(fit, clf, renderer, inference = "none", to_vox = FALSE)
length(enc$mean_anchor)
#> [1] 20
# }
```
