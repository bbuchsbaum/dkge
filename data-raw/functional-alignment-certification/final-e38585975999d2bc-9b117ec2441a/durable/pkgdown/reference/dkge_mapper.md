# Create a pluggable DKGE anchor mapper for dense rendering

Constructs a mapper descriptor for the legacy **dense rendering / anchor
pipeline** — used when projecting subject-space voxel/parcel values onto
a set of 3-D spatial anchor points (e.g., medoid centroids). Pass the
result to deprecated
[`dkge_build_renderer()`](https://bbuchsbaum.github.io/dkge/reference/dkge_build_renderer.md)
or
[`dkge_render_subject_values()`](https://bbuchsbaum.github.io/dkge/reference/dkge_render_subject_values.md),
or directly to
[`fit_mapper()`](https://bbuchsbaum.github.io/dkge/reference/fit_mapper.md)
together with `subj_points` / `anchor_points` matrices. These operations
are descriptive and do not make an alignment inferentially eligible. Use
[`dkge_prepare_alignment()`](https://bbuchsbaum.github.io/dkge/reference/dkge_prepare_alignment.md)
for that contract.

Use
[`dkge_mapper_spec()`](https://bbuchsbaum.github.io/dkge/reference/dkge_mapper_spec.md)
instead when you need a **transport pipeline mapper** that operates in
an abstract feature space (ridge regression, Sinkhorn OT over
embeddings) for functions such as
[`dkge_prepare_transport()`](https://bbuchsbaum.github.io/dkge/reference/dkge_prepare_transport.md)
or
[`dkge_transport_spec()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_spec.md).

## Usage

``` r
dkge_mapper(type = c("knn", "sinkhorn", "ridge", "gw"), ...)
```

## Arguments

- type:

  Mapper backend identifier: `"knn"` (barycentric kNN), `"sinkhorn"` (OT
  over point clouds), `"ridge"`, `"gw"`, or a custom identifier whose
  `fit_mapper.dkge_mapper_<type>()` method is supplied by an extension.
  The ridge and Gromov-Wasserstein backends require plugins.

- ...:

  Backend-specific parameters stored within the mapper object (e.g. `k`,
  `sigx`, `sigz` for kNN; `epsilon` for Sinkhorn).

## Value

A `dkge_mapper` S3 descriptor consumed by
[`fit_mapper()`](https://bbuchsbaum.github.io/dkge/reference/fit_mapper.md)
and
[`apply_mapper()`](https://bbuchsbaum.github.io/dkge/reference/apply_mapper.md).

## Examples

``` r
spec <- dkge_mapper(type = "knn", k = 3, sigx = 1)
subj_points <- matrix(rnorm(12 * 3), 12, 3)
anchor_points <- matrix(rnorm(6 * 3), 6, 3)
fit <- fit_mapper(spec, subj_points = subj_points, anchor_points = anchor_points)
y_anchor <- apply_mapper(fit, rnorm(nrow(subj_points)))
length(y_anchor)
#> [1] 6
```
