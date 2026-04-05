# Prepare subject-to-medoid transport operators for reuse

Computes and caches subject-to-medoid transport matrices so downstream
routines (e.g. bootstraps) can reuse a fixed consensus mapping without
re-solving the transport problem on every call.

## Usage

``` r
dkge_prepare_transport(
  fit,
  centroids = NULL,
  loadings = NULL,
  betas = NULL,
  sizes = NULL,
  mapper = "sinkhorn",
  medoid = 1L,
  ...
)
```

## Arguments

- fit:

  A \`dkge\` object.

- centroids:

  List of subject centroid matrices. Defaults to the centroids stored on
  \`fit\` or \`fit\$input\`.

- loadings:

  Optional list of subject loadings (\`P_s x r\`). When omitted, they
  are recomputed from \`fit\$Btil\` or the supplied \`betas\`.

- betas:

  Optional list of subject betas used to recompute loadings when
  \`loadings\` is \`NULL\`.

- sizes:

  Optional list of cluster masses.

- mapper:

  Mapper specification or shorthand passed to \[dkge_mapper_spec()\].

- medoid:

  Index (1-based) of the reference subject.

- ...:

  Additional mapper arguments such as \`epsilon\` or \`lambda_spa\`.

## Value

A list containing cached transport objects: \`operators\`,
\`mapper_spec\`, \`feature_list\`, \`size_list\`, \`feature_ref\`,
\`size_ref\`, \`centroids\`, and \`medoid\`.
