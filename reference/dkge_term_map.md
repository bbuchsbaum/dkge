# Extract a term-specific coefficient map

Extract a term-specific coefficient map

## Usage

``` r
dkge_term_map(object, term, contrast = NULL, drop = TRUE)
```

## Arguments

- object:

  A \`dkge_between_rrr\` object.

- term:

  Formula term or model-matrix column name.

- contrast:

  Optional numeric contrast over model-matrix columns. When supplied,
  \`term\` is used only as a label.

- drop:

  Logical; drop single-row results to a vector.

## Value

Numeric vector or matrix of coefficients in target feature space.
