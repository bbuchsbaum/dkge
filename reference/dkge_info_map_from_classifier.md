# Decoder-style information map from latent classifier weights

Pulls a cross-fitted latent weight vector \\\beta^{(-s)}\\ back to each
subject's cluster space via \\m_s = A_s \beta^{(-s)}\\, transports those
maps to anchors, and aggregates across subjects for descriptive display.
Legacy renderer-based inference is refused because the renderer can
learn correspondence from the same full-fit features.

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

  Fitted `dkge` object.

- betas:

  Either a numeric vector of length `r` or a list of length `S` (number
  of subjects) containing one latent weight vector per subject.

- renderer:

  Renderer produced by
  [`dkge_build_renderer()`](https://bbuchsbaum.github.io/dkge/reference/dkge_build_renderer.md).

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

A list containing the mean anchor field, optional voxel map, per-subject
anchor maps, and a descriptive provenance receipt. Inferential fields
are always `NULL`.

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
info <- dkge_info_map_from_classifier(fit, clf$beta_by_subject, renderer, to_vox = FALSE)
length(info$mean_anchor)
#> [1] 20
# }
```
