# Refit a DKGE object with a new voxel-weight specification

Convenience helper that rebuilds a DKGE fit using the original inputs
but a different \[dkge_weights()\] specification. The original fit must
have been constructed with \`keep_inputs = TRUE\` (the default for
\[dkge()\]), so the underlying \`dkge_data\` bundle is available.

## Usage

``` r
dkge_update_weights(fit, weights = NULL)
```

## Arguments

- fit:

  A fitted \`dkge\` object.

- weights:

  A \`dkge_weights\` specification to apply. When \`NULL\`, the existing
  weight spec stored in \`fit\$weight_spec\` is reused.

## Value

A new \`dkge\` object with updated voxel weights.
