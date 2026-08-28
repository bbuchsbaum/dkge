# Render values with a legacy descriptive renderer (deprecated)

This compatibility helper never performs inference. It consumes the
correspondence already stored in a legacy renderer and labels the result
as descriptive. New workflows should use
[`dkge_render_aligned()`](https://bbuchsbaum.github.io/dkge/reference/dkge_render_aligned.md).

## Usage

``` r
dkge_render_subject_values(renderer, values_list, lambda = 0, to_vox = TRUE)
```

## Arguments

- renderer:

  Object produced by
  [`dkge_build_renderer()`](https://bbuchsbaum.github.io/dkge/reference/dkge_build_renderer.md).

- values_list:

  List of per-subject value vectors (aligned with centroids).

- lambda:

  Optional Laplacian smoothing strength applied in anchor space.

- to_vox:

  Logical; when `TRUE` (default) and a decoder is available, voxel maps
  are produced.

## Value

A list with `anchor` (dense anchor field), optional `voxel` map, the
aggregation diagnostics, and the intermediate per-subject anchor maps.

## Examples

``` r
# \donttest{
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 3, P = 15, snr = 5
)
fit <- dkge(toy$B_list, toy$X_list, K = toy$K, rank = 2)
centroids <- lapply(toy$B_list, function(B) matrix(rnorm(ncol(B) * 3), ncol(B), 3))
renderer <- suppressWarnings(dkge_build_renderer(fit,
                                centroids = centroids,
                                anchor_xyz = matrix(rnorm(20 * 3), 20, 3),
                                anchor_n = 20,
                                anchor_method = "sample"))
values_list <- lapply(centroids, function(C) rnorm(nrow(C)))
out <- suppressWarnings(dkge_render_subject_values(renderer, values_list))
length(out$anchor)
#> [1] 20
# }
```
