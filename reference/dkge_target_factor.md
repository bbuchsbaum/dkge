# Target helper for a single factor

Target helper for a single factor

## Usage

``` r
dkge_target_factor(
  fit,
  factor,
  residualize = TRUE,
  collapse = NULL,
  scope = "within_subject"
)
```

## Arguments

- fit:

  dkge object containing \`kernel_info\$map\` metadata.

- factor:

  Factor name.

- residualize:

  Logical; if \`TRUE\` (default) residualise higher-order targets
  against previously constructed lower-order targets.

- collapse:

  Optional named list describing how to collapse factors that do not
  appear in a target. Each entry may be \`"mean"\`, \`list(method =
  "mean", window = 3:8)\`, or a numeric vector of length equal to the
  number of levels providing custom weights (automatically normalised).

- scope:

  Permutation/exchangeability scope stored with each target (default
  "within_subject").
