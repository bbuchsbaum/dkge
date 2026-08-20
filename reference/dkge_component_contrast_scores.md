# Express DKGE component saliences in a contrast basis

Computes \\C' K U\\, where columns of \`C\` are user-supplied or
inferred design contrasts. This is a coordinate summary of the raw
saliences, not a replacement for inspecting \`K

## Usage

``` r
dkge_component_contrast_scores(
  fit,
  basis = NULL,
  comps = NULL,
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

## Value

Tidy data frame with one row per component and contrast-basis column.

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
head(dkge_component_contrast_scores(fit, comps = 1:2))
#>           term       contrast component component_id       score basis_source
#> 1   grand_mean     grand_mean       LV1            1 -0.74762851    factorial
#> 2         task           task       LV1            1 -0.58378655    factorial
#> 3      measure       measure1       LV1            1 -0.13232831    factorial
#> 4      measure       measure2       LV1            1 -0.31308419    factorial
#> 5 task:measure task1:measure1       LV1            1 -0.02067990    factorial
#> 6 task:measure task1:measure2       LV1            1  0.01723399    factorial
```
