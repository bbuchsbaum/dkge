# Specify a DKGE mapper strategy for the transport pipeline

Creates a mapper specification for the \*\*transport pipeline\*\* — used
when mapping subject-level contrast/loading vectors from a source
feature space (e.g., parcel embeddings) to a common reference space.
Pass the result to \[dkge_transport_spec()\],
\[dkge_prepare_transport()\], or directly to \[fit_mapper()\] together
with \`source_feat\` / \`target_feat\` matrices.

Use \[dkge_mapper()\] instead when you need a \*\*dense rendering /
anchor mapper\*\* that works in 3-D spatial coordinates (kNN
barycentric, Sinkhorn OT over point clouds) for functions such as
\[dkge_build_renderer()\] or \[dkge_render_subject_values()\].

## Usage

``` r
dkge_mapper_spec(type = c("sinkhorn", "ridge", "ols"), ..., name = NULL)
```

## Arguments

- type:

  Mapping strategy identifier: \`"sinkhorn"\` (optimal transport),
  \`"ridge"\` (ridge regression), or \`"ols"\` (ordinary least squares).

- ...:

  Strategy-specific hyperparameters stored in the specification (e.g.
  \`lambda\`, \`epsilon\`, \`lambda_emb\`, \`lambda_spa\`).

- name:

  Optional user-facing name for diagnostics.

## Value

A \`dkge_mapper_spec\` object consumed by \[fit_mapper()\] and
\[predict_mapper()\].

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
