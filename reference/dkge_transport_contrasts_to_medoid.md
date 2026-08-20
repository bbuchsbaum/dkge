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
  method = c("sinkhorn", "ridge", "ols", "sinkhorn_cpp"),
  transport_cache = NULL,
  provenance = NULL,
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

  Mapper strategy (\`"sinkhorn"\`, \`"ridge"\`, or \`"ols"\`). The
  legacy \`"sinkhorn_cpp"\` name is a deprecated alias for
  \`"sinkhorn"\`.

- transport_cache:

  Optional cache from \[dkge_prepare_transport()\]. When supplied,
  cached operators are reused for every contrast.

- provenance:

  Optional transport provenance declaration. Descriptive loading-derived
  transport is recorded when omitted.

- ...:

  Additional parameters passed when building the default mapper
  specification.

## Value

Named list of transport results (one per contrast) with an attached
\`cache\` element for reuse.
