# Build or validate anchor coordinates in MNI space

Build or validate anchor coordinates in MNI space

## Usage

``` r
dkge_make_anchors(
  xyz = NULL,
  anchors = NULL,
  n_anchor = 20000L,
  method = c("kmeans", "sample"),
  seed = NULL
)
```

## Arguments

- xyz:

  Optional \\N \times 3\\ matrix of grey-matter coordinates (in mm).
  Used when anchors need to be derived from voxel locations.

- anchors:

  Optional precomputed \\Q \times 3\\ anchor matrix. When supplied it is
  validated and returned unchanged.

- n_anchor:

  Target number of anchors when deriving them from \`xyz\`.

- method:

  Selection strategy when deriving anchors. \`"kmeans"\` (default) uses
  k-means centroids, \`"sample"\` performs a uniform subsample.

- seed:

  Optional random seed for reproducibility.

## Value

A numeric \\Q \times 3\\ matrix of anchor coordinates.

## Examples

``` r
xyz <- matrix(rnorm(200 * 3), 200, 3)
anchors <- dkge_make_anchors(xyz = xyz, n_anchor = 10, method = "sample", seed = 1)
dim(anchors)
#> [1] 10  3
```
