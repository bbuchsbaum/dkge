# Transport component loadings to a medoid parcellation

Transport component loadings to a medoid parcellation

## Usage

``` r
dkge_transport_loadings_to_medoid(
  fit,
  medoid,
  centroids,
  loadings = NULL,
  betas = NULL,
  sizes = NULL,
  mapper = NULL,
  method = c("sinkhorn", "sinkhorn_cpp"),
  transport_cache = NULL,
  ...
)
```

## Arguments

- fit:

  A \`dkge\` object used to compute the loadings.

- medoid:

  Integer index of the reference subject (1-based).

- centroids:

  List of subject cluster centroids (each P_s x 3 matrix).

- loadings:

  Optional list of subject loadings (P_s x r). When omitted, they are
  recomputed from \`betas\`.

- betas:

  Optional list of subject betas used to recompute loadings when
  \`loadings\` is \`NULL\`.

- sizes:

  Optional list of cluster masses (defaults to uniform weights).

- mapper:

  Optional mapper specification created by \[dkge_mapper_spec()\]. When
  \`NULL\`, defaults to Sinkhorn with the supplied parameters.

- method:

  Backwards-compatible alias; currently only "sinkhorn" is supported.

- transport_cache:

  Optional cache from \[dkge_prepare_transport()\]. When supplied,
  cached operators are reused for all components.

- ...:

  Additional parameters passed when building the default mapper
  specification (e.g. \`epsilon\`, \`lambda_emb\`).

## Value

List with \`group\` (medoid cluster vectors per component), \`subjects\`
(per-subject transported values), and \`cache\` (transport cache reused
for future calls).
