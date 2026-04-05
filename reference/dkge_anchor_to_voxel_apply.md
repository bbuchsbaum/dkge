# Decode anchor values to voxel space

Decode anchor values to voxel space

## Usage

``` r
dkge_anchor_to_voxel_apply(decoder, anchor_values)
```

## Arguments

- decoder:

  Object returned by \[dkge_anchor_to_voxel_fit()\].

- anchor_values:

  Numeric vector of length \`decoder\$n_anchors\`.

## Value

Numeric vector of voxel values.

## Examples

``` r
anchors <- matrix(rnorm(20 * 3), 20, 3)
vox_xyz <- matrix(rnorm(50 * 3), 50, 3)
decoder <- dkge_anchor_to_voxel_fit(anchors, vox_xyz, k = 4)
vox_vals <- dkge_anchor_to_voxel_apply(decoder, rnorm(nrow(anchors)))
length(vox_vals)
#> [1] 50
```
