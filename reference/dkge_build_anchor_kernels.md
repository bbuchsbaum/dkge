# Fold-aware anchor kernel construction

Projects per-subject item kernels onto a shared anchor basis derived
from a pooled training feature set inside each cross-validation fold.
The resulting aligned kernels all have dimension \`L x L\` with
identical ordering.

## Usage

``` r
dkge_build_anchor_kernels(
  features_list,
  K_item_list,
  folds = NULL,
  L = 128,
  method = c("kcenter", "kmeanspp", "random", "dkpp"),
  seed = 1,
  sigma = NULL,
  rho = 0.5,
  fill = c("kcenter", "kmeanspp"),
  center = TRUE,
  whiten = TRUE,
  eps = 1e-06,
  unit_trace = TRUE,
  item_weights = NULL
)
```

## Arguments

- features_list:

  List of length \`S\`; element \`s\` is an \`n_s x d\` feature matrix
  for subject \`s\`.

- K_item_list:

  List of length \`S\`; element \`s\` is an \`n_s x n_s\` PSD item
  kernel aligned with \`features_list\[\[s\]\]\`.

- folds:

  Optional fold specification. Accepts \`NULL\` (single context), a list
  of held-out subject indices, or a \`dkge_folds\` object.

- L:

  Number of anchors.

- method:

  Anchor selection strategy (\`"kcenter"\`, \`"kmeanspp"\`,
  \`"random"\`, or \`"dkpp"\`).

- seed:

  Integer seed applied per fold (incremented by fold index).

- sigma:

  Optional bandwidth; when \`NULL\`, fold-specific heuristics are used.

- rho:

  Fraction used by the DPP stage when \`method = "dkpp"\`.

- fill:

  Completion strategy for \`method = "dkpp"\`.

- center:

  Logical; when \`TRUE\`, subtract column means of the feature response
  matrix before forming the projection.

- whiten:

  Logical; apply whitening with the anchor Gram inverse square root.

- eps:

  Diagonal jitter used during whitening.

- unit_trace:

  Logical; trace-normalise each subject kernel to maintain comparable
  scale.

- item_weights:

  Optional list of numeric vectors providing per-item weights (\`length
  == n_s\`).

## Value

A list indexed by folds; each element contains the anchors, bandwidth,
aligned kernels, anchor identifiers, and the fold's train/test indices.

## Examples

``` r
set.seed(1)
features_list <- list(
  s1 = matrix(rnorm(30 * 5), 30, 5),
  s2 = matrix(rnorm(40 * 5), 40, 5),
  s3 = matrix(rnorm(35 * 5), 35, 5)
)
K_item_list <- lapply(features_list, function(X) {
  Z <- matrix(rnorm(nrow(X) * 4), nrow(X), 4)
  tcrossprod(Z)
})
built <- dkge_build_anchor_kernels(features_list, K_item_list, L = 8, method = "dkpp")
dim(built[[1]]$K_aligned[[1]])
#> [1] 8 8
```
