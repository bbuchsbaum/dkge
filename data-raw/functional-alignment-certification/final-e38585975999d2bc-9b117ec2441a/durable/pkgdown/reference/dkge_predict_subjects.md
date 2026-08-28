# Convenience prediction for subject collections

Harmonises a variety of beta inputs (matrices, `dkge_subject` objects,
or `dkge_data` bundles) before forwarding to
[`dkge_predict()`](https://bbuchsbaum.github.io/dkge/reference/dkge_predict.md).
This allows callers to work with tidy inputs without manually assembling
`B_list` structures.

## Usage

``` r
dkge_predict_subjects(
  object,
  betas,
  contrasts,
  ids = NULL,
  return_loadings = TRUE,
  spatial = NULL
)
```

## Arguments

- object:

  dkge \| dkge_stream \| dkge_model.

- betas:

  Subject data. Accepts a matrix, list of matrices, `dkge_subject`
  objects, or a `dkge_data` bundle.

- contrasts:

  List or matrix accepted by
  [`dkge_predict()`](https://bbuchsbaum.github.io/dkge/reference/dkge_predict.md).

- ids:

  Optional subject identifiers overriding those inferred from `betas`.

- return_loadings:

  Logical; when TRUE, include projected loadings in the result bundle.

- spatial:

  Optional
  [`dkge_spatial_regularizer()`](https://bbuchsbaum.github.io/dkge/reference/dkge_spatial_regularizer.md)
  for prediction subjects. A fit with a shared spatial domain reuses its
  stored specification automatically. A fit with subject-specific
  domains requires this argument.

## Value

Output from
[`dkge_predict()`](https://bbuchsbaum.github.io/dkge/reference/dkge_predict.md)
with harmonised subject names.
