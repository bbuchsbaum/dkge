# Build a plotting contrast basis for DKGE effects

Returns a q-by-m matrix whose columns define interpretable coordinates
for plotting DKGE saliences. With cell-space design metadata, the
default is a factorial basis with a grand-mean column plus main
effects/interactions. With declared effect-space metadata or no
recoverable design metadata, the default is the identity basis over
effect rows. Declared cell metadata is matched by unique names;
incomplete or ambiguous mappings are errors.

## Usage

``` r
dkge_design_basis(
  fit,
  basis = NULL,
  include_intercept = TRUE,
  coding = c("sum", "helmert", "poly"),
  normalize = c("unit_K", "unit_l2", "none")
)
```

## Arguments

- fit:

  Fitted \`dkge\` object.

- basis:

  Optional custom q-by-m matrix. Row names, when present, are matched to
  \`fit\$effects\`.

- include_intercept:

  Include a grand-mean column for automatically built factorial bases.

- coding:

  Contrast coding for automatically built factorial bases: \`"sum"\`,
  \`"helmert"\`, or \`"poly"\`.

- normalize:

  Column scaling applied to automatically built bases, including the
  identity fallback used when the fit carries no design-cell metadata:
  \`"unit_K"\` (default) gives every column unit K-norm, \\\sqrt{c^\top
  K c} = 1\\; \`"unit_l2"\` gives unit Euclidean norm; and \`"none"\`
  leaves the coding matrix unscaled. User-supplied \`basis\` matrices
  are never rescaled.

## Value

A numeric matrix with attributes \`term\` and \`source\`.

## Details

Basis columns are scored against saliences through \\C^\top K U\\, so
the comparability of coordinates across columns is governed by the K
metric, not by Euclidean length. Normalizing to unit Euclidean norm
(\`"unit_l2"\`) leaves columns with different K-norms — a grand-mean
column is typically inflated relative to interaction columns purely by
the metric — which is why \`"unit_K"\` is the default.

## Examples

``` r
dk <- design_kernel(list(task = list(L = 2), measure = list(L = 3)),
                    basis = "cell", normalize = "none")
labels <- rownames(dk$K)
q <- length(labels)
set.seed(1)
B_list <- replicate(3, matrix(rnorm(q * 6), q, 6,
                              dimnames = list(labels, NULL)), simplify = FALSE)
X_list <- replicate(3, matrix(diag(q), q, q,
                              dimnames = list(labels, labels)), simplify = FALSE)
fit <- dkge(B_list, X_list, K = dk, rank = 2)
C <- dkge_design_basis(fit)
attr(C, "term")
#> [1] "grand_mean"   "task"         "measure"      "measure"      "task:measure"
#> [6] "task:measure"
# Columns have unit K-norm by default.
round(diag(t(C) %*% fit$K %*% C), 12)
#>     grand_mean           task       measure1       measure2 task1:measure1 
#>              1              1              1              1              1 
#> task1:measure2 
#>              1 
```
