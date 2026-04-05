# Transport subject contrasts to a medoid parcellation

Transport subject contrasts to a medoid parcellation

## Usage

``` r
dkge_transport_contrasts_to_medoid(
  fit,
  contrast_obj,
  medoid,
  centroids = NULL,
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

  A \`dkge\` object used to compute the contrasts.

- contrast_obj:

  A \`dkge_contrasts\` result.

- medoid:

  Integer index of the reference subject (1-based).

- centroids:

  List of subject cluster centroids (each P_s x 3 matrix).

- loadings:

  Optional list of subject loadings (P_s x r).

- betas:

  Optional list of subject betas used to recompute loadings.

- sizes:

  Optional list of cluster masses.

- mapper:

  Optional mapper specification created by \[dkge_mapper_spec()\].

- method:

  Backwards-compatible alias; currently only "sinkhorn" is supported.

- transport_cache:

  Optional cache from \[dkge_prepare_transport()\]. When supplied,
  cached operators are reused for every contrast.

- ...:

  Additional parameters passed when building the default mapper
  specification.

## Value

Named list of transport results (one per contrast) with an attached
\`cache\` element for reuse.
