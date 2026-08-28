# Group-LOCO anchor importance (zeroing proxy)

Approximates the influence of each anchor by zeroing its neighborhood in
the subject maps derived from classifier weights and accumulating the
lost margin magnitude across subjects.

## Usage

``` r
dkge_info_map_loco(
  fit,
  clf,
  renderer,
  neighborhoods = NULL,
  k_nn = 32,
  aggregate = c("mean", "sum")
)
```

## Arguments

- fit:

  Fitted `dkge` object.

- clf:

  Cross-fitted classifier from
  [`dkge_cv_train_latent_classifier()`](https://bbuchsbaum.github.io/dkge/reference/dkge_cv_train_latent_classifier.md).

- renderer:

  Renderer produced by
  [`dkge_build_renderer()`](https://bbuchsbaum.github.io/dkge/reference/dkge_build_renderer.md).

- neighborhoods:

  Optional list of integer vectors defining anchor neighborhoods. When
  `NULL`, the function derives them from the renderer's anchor graph.

- k_nn:

  When neighborhoods are not supplied, limits the automatically derived
  neighborhood size.

- aggregate:

  Aggregation method across subjects: `"mean"` or `"sum"`.

## Value

A list with LOCO scores per anchor, per-subject anchor maps, and the
neighborhood definition.
