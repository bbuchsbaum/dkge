# Render aligned maps or an existing support-level statistic

Render aligned maps or an existing support-level statistic

## Usage

``` r
dkge_render_aligned(
  renderer,
  x,
  statistic = c("weighted_mean", "mean", "median"),
  decode = TRUE
)
```

## Arguments

- renderer:

  A
  [`dkge_renderer()`](https://bbuchsbaum.github.io/dkge/reference/dkge_renderer.md)
  object.

- x:

  A
  [`dkge_aligned_maps()`](https://bbuchsbaum.github.io/dkge/reference/dkge_aligned_maps.md)
  object or an already computed support-level vector, matrix, or named
  list.

- statistic:

  Aggregation used only when `x` contains aligned subject rows.
  `"weighted_mean"` uses the explicit weights recorded on `x`.

- decode:

  Logical; apply the support decoder when one exists.

## Value

A `dkge_rendered_alignment` with support values and optional decoded
values. No correspondence is learned by this function.
