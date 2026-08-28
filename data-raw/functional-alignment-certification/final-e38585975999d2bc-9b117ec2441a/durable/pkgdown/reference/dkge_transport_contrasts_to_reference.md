# Align cross-fitted contrasts to an identified reference support

Strict reference-oriented counterpart to the deprecated
[`dkge_transport_contrasts_to_medoid()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_contrasts_to_medoid.md).
Typed alignment features are mandatory, preventing an inferential call
from silently falling back to full-fit loadings. Reference selection
follows the same held-out, geometry-only, or explicit rules as
[`dkge_prepare_alignment()`](https://bbuchsbaum.github.io/dkge/reference/dkge_prepare_alignment.md).

## Usage

``` r
dkge_transport_contrasts_to_reference(
  fit,
  contrast_obj,
  alignment_features,
  centroids = NULL,
  sizes = NULL,
  mapper = dkge_mapper_spec("sinkhorn"),
  reference_selection = NULL,
  validation_features = NULL,
  reference_subject = NULL,
  selection_method = c("auto", "functional_heldout", "geometry_only", "explicit",
    "descriptive_training"),
  transport_cache = NULL,
  ...
)
```

## Arguments

- fit:

  A fitted `dkge` object.

- contrast_obj:

  Cross-fitted contrasts from
  [`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md).

- alignment_features:

  Typed functional features.

- centroids:

  Named subject centroid matrices, or `NULL` to use `fit`.

- sizes:

  Optional subject parcel masses.

- mapper:

  Fixed mapper specification.

- reference_selection, validation_features, reference_subject:

  Reference selection object, held-out validation channel, or fixed
  subject; see
  [`dkge_prepare_alignment()`](https://bbuchsbaum.github.io/dkge/reference/dkge_prepare_alignment.md).

- selection_method:

  Reference-selection channel; see
  [`dkge_prepare_alignment()`](https://bbuchsbaum.github.io/dkge/reference/dkge_prepare_alignment.md).

- transport_cache:

  Optional exactly matching fitted alignment.

- ...:

  Mapper parameters used only when `mapper` is shorthand.

## Value

A named transport result with attached `dkge_fitted_alignment`,
`dkge_aligned_maps`, and `dkge_reference_selection` objects.
