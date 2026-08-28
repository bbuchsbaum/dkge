# Store descriptive subject-by-support maps

This constructor accepts values that are already on a reference support,
so it cannot prove that a fitted correspondence produced them. Its
result is always descriptive and inferentially ineligible, even when
`fitted_alignment` itself is eligible. Use
[`dkge_transport_contrasts_to_reference()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_contrasts_to_reference.md)
or
[`dkge_align_to_template()`](https://bbuchsbaum.github.io/dkge/reference/dkge_align_to_template.md)
to apply fitted operators and obtain operator-bound aligned maps.

## Usage

``` r
dkge_aligned_maps(
  values,
  fitted_alignment,
  subject_ids = NULL,
  contrast_ids = NULL,
  subject_weights = NULL,
  estimand = NULL
)
```

## Arguments

- values:

  Named list of matrices, one per contrast. Every matrix must be finite
  with one row per subject and one column per reference location.

- fitted_alignment:

  A typed `dkge_fitted_alignment`; used to identify the support, not to
  certify that it produced `values`.

- subject_ids:

  Optional subject identifiers.

- contrast_ids:

  Optional contrast identifiers.

- subject_weights:

  Optional fixed group weights. Named weights are reordered by exact
  subject ID; unnamed weights are positional in fitted subject order.
  Equal subject weighting is the default; DKGE fit-level MFA weights are
  never inherited silently.

- estimand:

  Optional explicit estimand description.

## Value

An immutable, descriptive `dkge_aligned_maps` object. It can be rendered
but is refused by group-inference functions.
