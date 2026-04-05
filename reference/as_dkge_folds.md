# Convert to DKGE fold assignments

Coerce various fold specifications into the \`dkge_folds\` structure
used by dkge cross-fitting helpers. Methods may use \`fit_or_data\` to
resolve subject identifiers when necessary.

## Usage

``` r
as_dkge_folds(x, fit_or_data = NULL, ...)
```

## Arguments

- x:

  Object describing fold assignments

- fit_or_data:

  Optional \`dkge\` or \`dkge_data\` object used to resolve subject
  identifiers

- ...:

  Additional arguments passed to methods

## Value

Object with class \`dkge_folds\`

## See also

\`vignette("dkge-classification", package = "dkge")\` for an example of
fold conversion in practice.

## Examples

``` r
folds <- as_dkge_folds(list(fold1 = 1:2, fold2 = 3:4))
folds$k
#> [1] 2
```
