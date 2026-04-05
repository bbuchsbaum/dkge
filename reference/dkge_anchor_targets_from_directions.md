# Assemble anchor targets from feature-space directions

Converts named direction vectors into anchor weight rows using
\[dkge_anchor_contrast_from_direction()\].

## Usage

``` r
dkge_anchor_targets_from_directions(
  anchors,
  directions,
  sigma = NULL,
  normalize = TRUE
)
```

## Arguments

- anchors:

  Matrix of anchor coordinates (\`L x d\`).

- directions:

  Either a named list of numeric vectors (length \`d\`) or a matrix
  whose rows are named directions in the same feature space as
  \`anchors\`.

- sigma:

  Optional bandwidth forwarded to
  \[dkge_anchor_contrast_from_direction()\].

- normalize:

  Logical; when \`TRUE\` (default) the resulting weights are
  L2-normalised.

## Value

Matrix with one row per supplied direction and \`nrow(anchors)\`
columns.

## Examples

``` r
anchors <- matrix(rnorm(10 * 4), 10, 4)
dirs <- list(classA = rnorm(4), classB = rnorm(4))
W <- dkge_anchor_targets_from_directions(anchors, dirs)
dim(W)
#> [1]  2 10
```
