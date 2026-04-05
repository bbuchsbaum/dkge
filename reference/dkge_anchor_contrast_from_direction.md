# Build an anchor contrast from a feature-space direction

Build an anchor contrast from a feature-space direction

## Usage

``` r
dkge_anchor_contrast_from_direction(
  anchors,
  direction,
  sigma = NULL,
  normalize = TRUE
)
```

## Arguments

- anchors:

  Matrix of anchor coordinates (\`L x d\`).

- direction:

  Numeric vector (length \`d\`) describing a linear probe in the shared
  feature space.

- sigma:

  Optional bandwidth; defaults to the anchor median heuristic.

- normalize:

  Logical; L2-normalise the resulting contrast.

## Value

Numeric vector of length \`L\` suitable for \[dkge_contrast()\].

## Examples

``` r
anchors <- matrix(rnorm(10 * 4), 10, 4)
direction <- rnorm(4)
w <- dkge_anchor_contrast_from_direction(anchors, direction)
length(w)
#> [1] 10
```
