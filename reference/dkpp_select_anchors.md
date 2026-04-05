# Determinantal k-means++ (d-kpp) anchor selection

Combines a greedy k-DPP seeding phase with either farthest-first or D^2
sampling to reach \`L\` anchors. The kernel bandwidth defaults to the
median distance heuristic computed on the training pool.

## Usage

``` r
dkpp_select_anchors(
  X_train,
  L = 128,
  rho = 0.5,
  sigma = NULL,
  fill = c("kcenter", "kmeanspp"),
  seed = 1,
  min_gain = 1e-12
)
```

## Arguments

- X_train:

  Matrix of training feature vectors pooled across subjects.

- L:

  Target number of anchors.

- rho:

  Fraction of anchors chosen by the DPP stage (0 \< rho \<= 1).

- sigma:

  Optional RBF bandwidth. When \`NULL\`, the median heuristic is applied
  on \`X_train\`.

- fill:

  Strategy for completing the anchor set after the DPP stage:
  \`"kcenter"\` (farthest-first) or \`"kmeanspp"\`.

- seed:

  Integer seed for reproducibility.

- min_gain:

  Minimum marginal gain tolerated during the DPP phase.

## Value

List with \`indices\` (row indices in \`X_train\`), \`sigma\`, and
selection diagnostics.

## Examples

``` r
X_train <- matrix(rnorm(20 * 3), 20, 3)
sel <- dkpp_select_anchors(X_train, L = 6, seed = 1)
length(sel$indices)
#> [1] 6
```
