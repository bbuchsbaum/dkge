# Prepare a reference-oriented fitted alignment

This is the strict high-level replacement for constructing transport
state around an unnamed or default subject index. It requires a typed
functional feature object. With a second independent feature channel it
selects a functional medoid by held-out symmetric reconstruction; with
no held-out channel it selects a geometry-only medoid. A caller may
instead supply a fixed `reference_subject` or a previously fitted
selection object.

## Usage

``` r
dkge_prepare_alignment(
  fit,
  alignment_features,
  centroids = NULL,
  sizes = NULL,
  mapper = dkge_mapper_spec("sinkhorn"),
  reference_selection = NULL,
  validation_features = NULL,
  reference_subject = NULL,
  selection_method = c("auto", "functional_heldout", "geometry_only", "explicit",
    "descriptive_training"),
  ...
)
```

## Arguments

- fit:

  A fitted `dkge` object.

- alignment_features:

  Typed functional features from
  [`dkge_alignment_features()`](https://bbuchsbaum.github.io/dkge/reference/dkge_alignment_features.md).
  Independent features are recommended.

- centroids:

  Named subject centroid matrices, or `NULL` to use geometry stored on
  `fit`.

- sizes:

  Optional positive subject parcel masses.

- mapper:

  Fixed mapper specification.

- reference_selection:

  Optional immutable result from
  [`dkge_select_reference_subject()`](https://bbuchsbaum.github.io/dkge/reference/dkge_select_reference_subject.md).

- validation_features:

  Optional second, separately identified independent channel with a
  validated provenance receipt, used for held-out functional reference
  selection.

- reference_subject:

  Optional fixed subject index or ID. This creates an explicit
  reference, not a medoid.

- selection_method:

  Reference-selection channel. `"auto"` chooses held-out functional
  selection when `validation_features` is present, explicit selection
  when `reference_subject` is present, and geometry-only selection
  otherwise. Training-feature selection is opt-in and ineligible.

- ...:

  Mapper parameters used only when `mapper` is a shorthand.

## Value

A typed `dkge_fitted_alignment` with immutable reference-selection,
feature, support, template, mapper, and numerical receipts.
