# Convert to a DKGE design kernel

Coerce arbitrary objects into a DKGE kernel bundle. The default method
preserves existing behaviour by accepting matrices or list objects with
a \`\$K\` component and optional metadata.

## Usage

``` r
as_dkge_kernel(x, ...)
```

## Arguments

- x:

  Object containing kernel information

- ...:

  Additional arguments passed to methods

## Value

List with entries \`K\` (q x q matrix) and optional \`info\`

## See also

\`vignette("dkge-classification", package = "dkge")\` for a worked
example that uses these generics to integrate hyperdesign inputs.

## Examples

``` r
K <- diag(3)
as_dkge_kernel(K)$K
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
```
