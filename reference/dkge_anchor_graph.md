# Construct a kNN anchor graph and Laplacian

Construct a kNN anchor graph and Laplacian

## Usage

``` r
dkge_anchor_graph(
  anchors,
  k = 10,
  sigma = NULL,
  weight_mode = c("heat", "normalized", "binary"),
  type = c("mutual", "normal", "asym")
)
```

## Arguments

- anchors:

  \\Q \times 3\\ matrix of anchor coordinates.

- k:

  Number of neighbours used in the graph (default 10).

- sigma:

  Optional Gaussian length-scale (mm). If \`NULL\`, the \`adjoin\`
  defaults are used.

- weight_mode:

  Edge weighting scheme passed to \[adjoin::graph_weights()\]. Defaults
  to \`"heat"\`.

- type:

  Graph symmetrisation strategy (see \[adjoin::graph_weights()\]);
  \`"mutual"\` helps enforce symmetry.

## Value

A list containing the neighbour graph, sparse adjacency \`W\`, degree
matrix \`D\`, and Laplacian \`L\`.

## Examples

``` r
anchors <- matrix(rnorm(40 * 3), 40, 3)
g <- dkge_anchor_graph(anchors, k = 3)
names(g)
#> [1] "graph" "W"     "D"     "L"    
```
