# Transport component loadings for legacy descriptive display

This compatibility helper learns correspondence from full-fit component
loadings and is therefore descriptive/ineligible. It cannot supply an
inferential alignment. New workflows should build typed alignment
features and call
[`dkge_prepare_alignment()`](https://bbuchsbaum.github.io/dkge/reference/dkge_prepare_alignment.md).

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
  method = c("sinkhorn", "ridge", "ols", "sinkhorn_cpp"),
  transport_cache = NULL,
  ...
)
```

## Arguments

- fit:

  A `dkge` object used to compute the loadings.

- medoid:

  Integer index of the reference subject (1-based).

- centroids:

  List of subject cluster centroids (each P_s x 3 matrix).

- loadings:

  Optional list of subject loadings (P_s x r). When omitted, they are
  recomputed from `betas`.

- betas:

  Optional list of subject betas used to recompute loadings when
  `loadings` is `NULL`.

- sizes:

  Optional list of cluster masses (defaults to uniform weights).

- mapper:

  Optional mapper specification created by
  [`dkge_mapper_spec()`](https://bbuchsbaum.github.io/dkge/reference/dkge_mapper_spec.md).
  When `NULL`, defaults to Sinkhorn with the supplied parameters.

- method:

  Mapper strategy (`"sinkhorn"`, `"ridge"`, or `"ols"`). The legacy
  `"sinkhorn_cpp"` name is a deprecated alias for `"sinkhorn"`.

- transport_cache:

  Optional fitted alignment from
  [`dkge_prepare_transport()`](https://bbuchsbaum.github.io/dkge/reference/dkge_prepare_transport.md).
  Reuse requires exact structural fingerprints.

- ...:

  Additional parameters passed when building the default mapper
  specification (e.g. `epsilon`, `lambda_emb`).

## Value

List with `group` (medoid cluster vectors per component), `subjects`
(per-subject transported values), and `cache` (transport cache reused
for future calls).
