# Aggregate anchor fields with optional Laplacian smoothing

This is a post-fit rendering operation: it smooths an already
transported anchor average and does not change the DKGE basis,
components, or held-out subject fields. To place a spatial penalty
inside the fitted solution, use
[`dkge_spatial_regularizer()`](https://bbuchsbaum.github.io/dkge/reference/dkge_spatial_regularizer.md)
through the `spatial` argument of
[`dkge_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit.md).

## Usage

``` r
dkge_anchor_aggregate(anchor_list, subj_weights = NULL, L = NULL, lambda = 0)
```

## Arguments

- anchor_list:

  List of anchor-valued vectors (length `Q`).

- subj_weights:

  Optional subject weights applied during the average.

- L:

  Optional graph Laplacian for Tikhonov regularisation.

- lambda:

  Non-negative smoothing parameter. `0` disables smoothing.

## Value

A list with the smoothed field `y`, the raw weighted mean `ybar`,
`coverage` (weighted contribution counts), and placeholder `ess` values.

## Examples

``` r
anchor_list <- list(rnorm(10), rnorm(10), rnorm(10))
agg <- dkge_anchor_aggregate(anchor_list)
length(agg$y)
#> [1] 10
```
