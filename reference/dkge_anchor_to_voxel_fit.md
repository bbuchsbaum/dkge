# Fit a sparse anchor-to-voxel decoder

Each voxel is represented as a convex combination of its \`k\` nearest
anchors using Gaussian weights. The result can be reused to decode any
anchor-level map.

## Usage

``` r
dkge_anchor_to_voxel_fit(anchors, vox_xyz, k = 8, sigma = NULL)
```

## Arguments

- anchors:

  Anchor coordinate matrix (\`Q x 3\`).

- vox_xyz:

  Voxel coordinate matrix (\`V x 3\`).

- k:

  Number of anchors per voxel (default 8).

- sigma:

  Optional Gaussian length-scale. If \`NULL\`, it is set to the square
  root of the median squared distance between voxels and their nearest
  anchors.

## Value

A decoder object storing neighbour indices, weights, and a sparse matrix
implementing the transformation.

## Examples

``` r
anchors <- matrix(rnorm(20 * 3), 20, 3)
vox_xyz <- matrix(rnorm(50 * 3), 50, 3)
decoder <- dkge_anchor_to_voxel_fit(anchors, vox_xyz, k = 4)
decoder$params$k
#> [1] 4
```
