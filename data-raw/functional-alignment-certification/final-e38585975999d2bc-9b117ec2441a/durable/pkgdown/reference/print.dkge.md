# Print a fitted DKGE model

Reports the fitted dimensions, spatial-smoothing effectiveness, and the
subject-level block weighting used in the pooled moment. When selected,
MFA/energy normalization weights are not automatically inverse-variance
weights for group inference; `w_method = "none"` reports equal moment
weights.

## Usage

``` r
# S3 method for class 'dkge'
print(x, ...)
```

## Arguments

- x:

  A fitted `dkge` object.

- ...:

  Unused.

## Value

`x`, invisibly.
