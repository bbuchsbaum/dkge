# Convenience prediction for subject collections

Harmonises a variety of beta inputs (matrices, \`dkge_subject\` objects,
or \`dkge_data\` bundles) before forwarding to \[dkge_predict()\]. This
allows callers to work with tidy inputs without manually assembling
\`B_list\` structures.

## Usage

``` r
dkge_predict_subjects(
  object,
  betas,
  contrasts,
  ids = NULL,
  return_loadings = TRUE
)
```

## Arguments

- object:

  dkge \| dkge_stream \| dkge_model.

- betas:

  Subject data. Accepts a matrix, list of matrices, \`dkge_subject\`
  objects, or a \`dkge_data\` bundle.

- contrasts:

  List or matrix accepted by \[dkge_predict()\].

- ids:

  Optional subject identifiers overriding those inferred from \`betas\`.

- return_loadings:

  Logical; when TRUE, include projected loadings in the result bundle.

## Value

Output from \[dkge_predict()\] with harmonised subject names.
