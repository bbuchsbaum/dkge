# Project a new cluster/voxel vector onto DKGE components

Project a new cluster/voxel vector onto DKGE components

## Usage

``` r
dkge_project_cluster(fit, b, omega = 1, w = 1)
```

## Arguments

- fit:

  A \`dkge\` object.

- b:

  Numeric vector of length q (effects).

- omega:

  Optional scalar or matrix weight.

- w:

  Optional scalar subject weight.

## Value

Numeric vector of length \`rank\` (component scores).
