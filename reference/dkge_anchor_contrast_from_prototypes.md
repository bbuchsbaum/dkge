# Build an anchor contrast from prototype feature sets

Build an anchor contrast from prototype feature sets

## Usage

``` r
dkge_anchor_contrast_from_prototypes(
  anchors,
  positives,
  negatives = NULL,
  sigma = NULL,
  normalize = TRUE
)
```

## Arguments

- anchors:

  Matrix of anchor coordinates (\`L x d\`).

- positives:

  Matrix (or vector) of positive prototypes in feature space.

- negatives:

  Optional matrix (or vector) of negative prototypes.

- sigma:

  Optional bandwidth for the prototype kernel. Defaults to the median
  distance between anchors and prototypes.

- normalize:

  Logical; L2-normalise the resulting contrast.

## Value

Numeric vector of length \`L\` suitable for \[dkge_contrast()\].

## Examples

``` r
anchors <- matrix(rnorm(10 * 4), 10, 4)
pos <- anchors[1:2, , drop = FALSE]
neg <- anchors[3:4, , drop = FALSE]
w <- dkge_anchor_contrast_from_prototypes(anchors, positives = pos, negatives = neg)
length(w)
#> [1] 10
```
