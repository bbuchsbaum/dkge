# Project multiple cluster/voxel vectors

Project multiple cluster/voxel vectors

## Usage

``` r
dkge_project_clusters(fit, B, omega_vec = NULL, w = 1)
```

## Arguments

- fit:

  A \`dkge\` object.

- B:

  qxP matrix of cluster betas.

- omega_vec:

  Optional vector of per-cluster weights.

- w:

  Optional subject weight.

## Value

Pxrank matrix of projected coordinates.
