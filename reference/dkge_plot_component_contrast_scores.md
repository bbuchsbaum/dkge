# Plot DKGE component saliences in a contrast basis

Plot DKGE component saliences in a contrast basis

## Usage

``` r
dkge_plot_component_contrast_scores(
  fit,
  basis = NULL,
  comps = NULL,
  include_intercept = TRUE,
  coding = c("sum", "helmert", "poly"),
  normalize = c("unit_K", "unit_l2", "none"),
  type = c("heatmap", "bars")
)
```

## Arguments

- fit:

  Fitted \`dkge\` object.

- basis:

  Optional custom q-by-m matrix. Row names, when present, are matched to
  \`fit\$effects\`.

- comps:

  Components to include (defaults to first min(rank, 6)).

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

- type:

  Plot type: \`"heatmap"\` for compact comparison, \`"bars"\` for a
  component-by-component coordinate plot.

## Value

A ggplot object.

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
dkge_plot_component_contrast_scores(fit, comps = 1:2)
```
