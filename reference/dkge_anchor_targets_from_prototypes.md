# Assemble anchor targets from prototype sets

Convenience helper that stacks calls to
\[dkge_anchor_contrast_from_prototypes()\] so the resulting matrix can
be fed directly to \[dkge_classify()\] as a weight specification.

## Usage

``` r
dkge_anchor_targets_from_prototypes(
  anchors,
  prototypes,
  negatives = NULL,
  sigma = NULL,
  normalize = TRUE
)
```

## Arguments

- anchors:

  Matrix of anchor coordinates (\`L x d\`).

- prototypes:

  Named list whose elements are matrices (rows = prototypes in the same
  feature space as \`anchors\`). List names become class labels.

- negatives:

  Optional named list of matrices providing negative prototypes per
  class (matched by name). When \`NULL\`, classes are contrasted against
  the origin.

- sigma:

  Optional bandwidth passed to
  \[dkge_anchor_contrast_from_prototypes()\]. Defaults to the per-class
  median heuristic.

- normalize:

  Logical indicator forwarded to
  \[dkge_anchor_contrast_from_prototypes()\].

## Value

Matrix with one row per class and \`nrow(anchors)\` columns.

## Examples

``` r
anchors <- matrix(rnorm(10 * 4), 10, 4)
proto <- list(
  classA = anchors[1:2, , drop = FALSE],
  classB = anchors[3:4, , drop = FALSE]
)
W <- dkge_anchor_targets_from_prototypes(anchors, proto)
dim(W)
#> [1]  2 10
```
