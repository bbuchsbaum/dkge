# Specify a DKGE mapper strategy for the transport pipeline

Creates a mapper specification for the **transport pipeline** — used
when mapping subject-level contrast/loading vectors from a source
feature space (e.g., parcel embeddings) to a common reference space.
Pass the result to
[`dkge_transport_spec()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_spec.md),
[`dkge_prepare_transport()`](https://bbuchsbaum.github.io/dkge/reference/dkge_prepare_transport.md),
or directly to
[`fit_mapper()`](https://bbuchsbaum.github.io/dkge/reference/fit_mapper.md)
together with `source_feat` / `target_feat` matrices.

Use
[`dkge_mapper()`](https://bbuchsbaum.github.io/dkge/reference/dkge_mapper.md)
for the legacy descriptive dense-rendering pipeline.
[`dkge_build_renderer()`](https://bbuchsbaum.github.io/dkge/reference/dkge_build_renderer.md)
and
[`dkge_render_subject_values()`](https://bbuchsbaum.github.io/dkge/reference/dkge_render_subject_values.md)
are deprecated display helpers and do not establish inferential
correspondence. For group functional alignment, use
[`dkge_prepare_alignment()`](https://bbuchsbaum.github.io/dkge/reference/dkge_prepare_alignment.md)
and
[`dkge_render_aligned()`](https://bbuchsbaum.github.io/dkge/reference/dkge_render_aligned.md).

## Usage

``` r
dkge_mapper_spec(type = c("sinkhorn", "ridge", "ols"), ..., name = NULL)
```

## Arguments

- type:

  Mapping strategy identifier: `"sinkhorn"` (optimal transport),
  `"ridge"` (ridge regression), or `"ols"` (ordinary least squares).

- ...:

  Strategy-specific hyperparameters stored in the specification (e.g.
  `lambda`, `epsilon`, `lambda_emb`, `lambda_spa`).

- name:

  Optional user-facing name for diagnostics.

## Value

A `dkge_mapper_spec` object consumed by
[`fit_mapper()`](https://bbuchsbaum.github.io/dkge/reference/fit_mapper.md)
and
[`predict_mapper()`](https://bbuchsbaum.github.io/dkge/reference/predict_mapper.md).

## Examples

``` r
spec <- dkge_mapper_spec(type = "ridge", lambda = 1e-2)
source_feat <- matrix(rnorm(10 * 2), 10, 2)
target_feat <- matrix(rnorm(8 * 2), 8, 2)
mapping <- fit_mapper(spec, source_feat = source_feat, target_feat = target_feat)
mapped <- predict_mapper(mapping, rnorm(10))
length(mapped)
#> [1] 8
```
