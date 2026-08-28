# Apply a fitted functional template to subject-level values

Apply a fitted functional template to subject-level values

## Usage

``` r
dkge_align_to_template(
  template,
  support,
  values,
  contrast_ids = NULL,
  subject_weights = NULL,
  allow_nonconverged = FALSE
)
```

## Arguments

- template:

  A template returned by
  [`dkge_fit_functional_template()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit_functional_template.md).

- support:

  The same identified
  [`dkge_reference_support()`](https://bbuchsbaum.github.io/dkge/reference/dkge_reference_support.md)
  used to fit the template.

- values:

  Either a
  [`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md)
  result or one list of subject vectors (or a named list of such lists)
  per contrast. A typed contrast result binds estimator, family,
  operator, and output provenance. Raw values can be transported for
  description, but their result is always inferentially ineligible.
  Typed input must match the exact contrast result and family bound when
  the template features were built. Same-data residualized templates
  additionally reject raw values.

- contrast_ids:

  Optional contrast identifiers.

- subject_weights:

  Optional fixed subject aggregation weights. Named weights are
  reordered by exact subject ID; unnamed weights are positional in
  fitted subject order. Equal subject weighting is the default and is
  distinct from DKGE fit-level MFA weights.

- allow_nonconverged:

  Logical; permit construction of explicitly ineligible aligned maps
  from a non-converged template.

## Value

A `dkge_aligned_maps` object. Typed contrast input produces an
operator-bound object with its fitted alignment attached; raw input
produces a descriptive/ineligible object without an inferential
application receipt.
