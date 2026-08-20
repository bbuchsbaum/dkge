# Build a subject-by-feature target for between-subject DKGE models

Creates the explicit boundary object between the DKGE representation
stage and between-subject multivariate models. The returned target
stores a subject-by-feature matrix plus optional coverage, precision,
and weighting metadata.

## Usage

``` r
dkge_make_target(
  fit = NULL,
  type = c("matrix", "component_scores", "transported_maps"),
  Y = NULL,
  contrast = NULL,
  contrast_obj = NULL,
  transport = NULL,
  values = NULL,
  loadings = NULL,
  centroids = NULL,
  sizes = NULL,
  medoid = NULL,
  mapper = NULL,
  crossfit = c("none", "analytic", "loso", "kfold"),
  feature_ids = NULL,
  subject_ids = NULL,
  coverage = NULL,
  precision = NULL,
  subject_weights = NULL,
  feature_weights = NULL,
  geometry = NULL,
  provenance = NULL,
  ...
)
```

## Arguments

- fit:

  Optional \`dkge\` object used to derive DKGE targets.

- type:

  Target type. \`"matrix"\` wraps a supplied matrix,
  \`"component_scores"\` averages subject component expressions, and
  \`"transported_maps"\` builds a common-space matrix from contrast
  values.

- Y:

  Optional subject-by-feature matrix.

- contrast:

  Optional DKGE contrast used when \`type = "transported_maps"\`.

- contrast_obj:

  Optional precomputed \`dkge_contrasts\` object.

- transport:

  Optional \`dkge_transport_spec\` used for transported maps.

- values:

  Optional subject values. May be a matrix or a list of vectors.

- loadings:

  Optional subject loading matrices used by transport.

- centroids, sizes, medoid, mapper:

  Transport inputs overriding \`transport\`.

- crossfit:

  Cross-fitting method used when computing \`contrast\`.

- feature_ids, subject_ids:

  Optional identifiers.

- coverage:

  Optional coverage matrix/vector aligned with \`Y\`.

- precision:

  Optional precision matrix/vector aligned with \`Y\`.

- subject_weights:

  Optional subject reliability weights.

- feature_weights:

  Optional feature reliability weights.

- geometry:

  Optional feature geometry metadata.

- provenance:

  Optional provenance metadata.

- ...:

  Additional arguments forwarded to contrast/transport helpers.

## Value

Object of class \`dkge_target\`.
