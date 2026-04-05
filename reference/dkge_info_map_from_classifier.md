# Decoder-style information map from latent classifier weights

Pulls a cross-fitted latent weight vector \\\beta^{(-s)}\\ back to each
subject's cluster space via \\m_s = A_s \beta^{(-s)}\\, transports those
maps to anchors, aggregates across subjects, and optionally performs
inference.

## Usage

``` r
dkge_info_map_from_classifier(
  fit,
  betas,
  renderer,
  lambda = 0,
  to_vox = TRUE,
  inference = c("none", "signflip", "parametric")
)
```

## Arguments

- fit:

  Fitted \`dkge\` object.

- betas:

  Either a numeric vector of length \`r\` or a list of length \`S\`
  (number of subjects) containing one latent weight vector per subject.

- renderer:

  Renderer produced by \[dkge_build_renderer()\].

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

A list containing the mean anchor field, optional voxel map, per-subject
anchor maps, and inference outputs (t- and p-values when requested).

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
info <- dkge_info_map_from_classifier(fit, clf$beta_by_subject, renderer, to_vox = FALSE)
length(info$mean_anchor)
#> [1] 20
# }
```
