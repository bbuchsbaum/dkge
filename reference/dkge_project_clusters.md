# Project multiple cluster/voxel vectors

Project multiple cluster/voxel vectors

## Usage

``` r
dkge_project_clusters(fit, B, omega_vec = NULL, w = 1, subject = NULL)
```

## Arguments

- fit:

  A `dkge` object.

- B:

  qxP matrix of cluster betas.

- omega_vec:

  Optional vector of per-cluster weights.

- w:

  Optional subject weight.

- subject:

  Optional training-subject index or id used to select a
  subject-specific spatial operator.

## Value

Pxrank matrix of projected coordinates.
