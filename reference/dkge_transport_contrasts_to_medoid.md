# Transport subject contrasts to a medoid parcellation (deprecated)

`dkge_transport_contrasts_to_medoid()` is retained for one migration
cycle. It preserves the legacy descriptive/fold-safe behavior where
safe, while current cache-provenance checks still fail closed. New
analyses should use
[`dkge_transport_contrasts_to_reference()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_contrasts_to_reference.md),
which requires typed features and an auditable reference-selection
channel.

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
  reference_selection = NULL,
  alignment_features = NULL,
  alignment_mode = c("fold_safe", "independent", "contrast_orthogonal", "descriptive"),
  ...
)
```

## Arguments

- fit:

  A fitted `dkge` object.

- contrast_obj:

  Cross-fitted contrasts from
  [`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md).

- medoid:

  Legacy integer reference-subject index.

- centroids:

  Named subject centroid matrices, or `NULL` to use `fit`.

- loadings:

  Optional legacy loose loading matrices.

- betas:

  Optional legacy beta matrices used to derive loadings.

- sizes:

  Optional subject parcel masses.

- mapper:

  Fixed mapper specification.

- method:

  Legacy mapper strategy.

- transport_cache:

  Optional exactly matching fitted alignment.

- reference_selection:

  Optional typed selection object supported by the migration shim.

- alignment_features:

  Typed functional features.

- alignment_mode:

  Legacy feature-provenance mode.

- ...:

  Mapper parameters used only when `mapper` is shorthand.
