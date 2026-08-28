# Prepare a legacy descriptive renderer (deprecated)

This helper fits correspondence directly from caller-supplied geometry
or features and is retained only for descriptive displays. It does not
create a typed inferential alignment. Use
[`dkge_prepare_alignment()`](https://bbuchsbaum.github.io/dkge/reference/dkge_prepare_alignment.md)
followed by
[`dkge_renderer()`](https://bbuchsbaum.github.io/dkge/reference/dkge_renderer.md)
and
[`dkge_render_aligned()`](https://bbuchsbaum.github.io/dkge/reference/dkge_render_aligned.md)
for new workflows.

## Usage

``` r
dkge_build_renderer(
  fit,
  centroids,
  anchors = NULL,
  anchor_xyz = NULL,
  anchor_n = 20000L,
  anchor_method = c("kmeans", "sample"),
  anchor_seed = NULL,
  vox_xyz = NULL,
  mapper = dkge_mapper("knn", k = 8, sigx = 3),
  graph_k = NULL,
  decoder_k = 8,
  reliabilities = NULL,
  subject_feats = NULL,
  anchor_feats = NULL,
  feat_lambda = NULL,
  feat_sigma = NULL,
  subject_weights = NULL
)
```

## Arguments

- fit:

  Fitted `dkge` object.

- centroids:

  List of per-subject centroid matrices (`P_s x 3`). Must align with
  `fit$Btil`.

- anchors:

  Optional precomputed anchor coordinate matrix. When `NULL`, anchors
  are derived from `anchor_xyz` if provided, otherwise from `vox_xyz`.

- anchor_xyz:

  Optional matrix of candidate points used to derive anchors via
  [`dkge_make_anchors()`](https://bbuchsbaum.github.io/dkge/reference/dkge_make_anchors.md).
  Ignored when `anchors` is supplied.

- anchor_n:

  Number of anchors to draw when constructing them from coordinates.

- anchor_method:

  Method passed to
  [`dkge_make_anchors()`](https://bbuchsbaum.github.io/dkge/reference/dkge_make_anchors.md)
  when anchors are derived. Defaults to `"kmeans"`.

- anchor_seed:

  Optional seed forwarded to
  [`dkge_make_anchors()`](https://bbuchsbaum.github.io/dkge/reference/dkge_make_anchors.md).

- vox_xyz:

  Optional voxel coordinates for constructing a decoder.

- mapper:

  Mapper specification created with
  [`dkge_mapper()`](https://bbuchsbaum.github.io/dkge/reference/dkge_mapper.md).
  Defaults to the barycentric kNN mapper.

- graph_k:

  Optional integer; when provided, an anchor graph of this neighborhood
  size is constructed for subsequent smoothing.

- decoder_k:

  Number of anchors per voxel when building the decoder.

- reliabilities:

  Optional list of per-subject reliability vectors passed to the mapper
  during fitting.

- subject_feats:

  Optional list of matrices supplying latent features per subject
  cluster. When provided and the mapper consumes latent information
  (e.g., Sinkhorn), they are forwarded via `subj_feats`.

- anchor_feats:

  Optional anchor-level feature matrix aligned with `anchors`. Derived
  automatically by pooling subject features when `NULL` and
  `subject_feats` are provided.

- feat_lambda:

  Feature cost weight passed to Sinkhorn mappers. Ignored by kNN.

- feat_sigma:

  Feature bandwidth used when computing feature costs.

- subject_weights:

  Optional fixed subject weights for descriptive aggregation. Equal
  subject weighting is the default; fit-level MFA weights are never
  inherited silently.

## Value

A deprecated `dkge_legacy_renderer` bundling anchors, optional
graph/decoder, fitted per-subject mappers, fixed subject weights, and an
ineligible/descriptive alignment receipt.

## Examples

``` r
# \donttest{
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 3, P = 15, snr = 5
)
fit <- dkge(toy$B_list, toy$X_list, K = toy$K, rank = 2)
centroids <- lapply(toy$B_list, function(B) matrix(rnorm(ncol(B) * 3), ncol(B), 3))
renderer <- suppressWarnings(dkge_build_renderer(fit,
                                centroids = centroids,
                                anchor_xyz = matrix(rnorm(20 * 3), 20, 3),
                                anchor_n = 20,
                                anchor_method = "sample"))
length(renderer$anchors)
#> [1] 60
# }
```
