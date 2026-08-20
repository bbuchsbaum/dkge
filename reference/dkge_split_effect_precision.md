# Export split-half reliability as effect precision

Converts the across-spatial-unit agreement between two stored effect
estimates into a finite precision in \`\[0, 1\]\`. Negative or undefined
correlations map to zero. A deterministic count factor shrinks effects
with few observations in either half.

## Usage

``` r
dkge_split_effect_precision(
  subject,
  method = c("spearman_brown", "positive_r2"),
  min_features = 3L,
  count_prior = 2
)
```

## Arguments

- subject:

  A \`dkge_subject\` containing two \`split_betas\` matrices.

- method:

  Reliability transform. \`"spearman_brown"\` applies the Spearman-Brown
  correction to the non-negative half correlation; \`"positive_r2"\`
  squares the non-negative correlation.

- min_features:

  Minimum number of spatial units required to estimate a correlation.

- count_prior:

  Non-negative pseudo-count controlling low-count shrinkage. Zero
  disables count shrinkage.

## Value

A named numeric effect-precision vector in \`\[0, 1\]\`. Diagnostics are
stored in the \`"diagnostics"\` attribute.
