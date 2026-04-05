# Create a pluggable DKGE anchor mapper for dense rendering

Constructs a mapper descriptor for the \*\*dense rendering / anchor
pipeline\*\* — used when projecting subject-space voxel/parcel values
onto a set of 3-D spatial anchor points (e.g., medoid centroids). Pass
the result to \[dkge_build_renderer()\],
\[dkge_render_subject_values()\], or directly to \[fit_mapper()\]
together with \`subj_points\` / \`anchor_points\` matrices.

Use \[dkge_mapper_spec()\] instead when you need a \*\*transport
pipeline mapper\*\* that operates in an abstract feature space (ridge
regression, Sinkhorn OT over embeddings) for functions such as
\[dkge_prepare_transport()\] or \[dkge_transport_spec()\].

## Usage

``` r
dkge_mapper(type = c("knn", "sinkhorn", "ridge", "gw"), ...)
```

## Arguments

- type:

  Mapper backend identifier: \`"knn"\` (barycentric kNN, fully
  implemented), \`"sinkhorn"\` (OT over point clouds), \`"ridge"\`, or
  \`"gw"\` (Gromov-Wasserstein; latter three require external plugins).

- ...:

  Backend-specific parameters stored within the mapper object (e.g.
  \`k\`, \`sigx\`, \`sigz\` for kNN; \`epsilon\` for Sinkhorn).

## Value

A \`dkge_mapper\` S3 descriptor consumed by \[fit_mapper()\] and
\[apply_mapper()\].

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
