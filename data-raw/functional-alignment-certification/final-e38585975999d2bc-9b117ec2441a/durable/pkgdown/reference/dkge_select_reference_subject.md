# Select an auditable reference subject for functional alignment

A reference subject supplies the target support; it is a *medoid* only
when selected by a stated cohort criterion. The recommended functional
criterion fits each directed mapper on one provenance-declared
independent feature channel and scores it on a second, separately
identified channel. DKGE validates the receipts and nonidentity of those
channels; the study owner is responsible for the scientific independence
claim. Pair losses are symmetrized and the subject with the smallest
mean held-out normalized reconstruction error is selected. Mapper cost
weights, spatial scale, epsilon, and feature normalization are fixed
across every candidate.

## Usage

``` r
dkge_select_reference_subject(
  centroids,
  alignment_features = NULL,
  validation_features = NULL,
  sizes = NULL,
  mapper = dkge_mapper_spec("sinkhorn"),
  method = c("functional_heldout", "geometry_only", "explicit", "descriptive_training"),
  reference_subject = NULL,
  subject_ids = NULL,
  provenance = NULL,
  tie_tolerance = 1e-10
)
```

## Arguments

- centroids:

  Named list of subject coordinate matrices.

- alignment_features:

  Typed training features from
  [`dkge_alignment_features()`](https://bbuchsbaum.github.io/dkge/reference/dkge_alignment_features.md).

- validation_features:

  A second typed independent feature object for held-out functional
  selection.

- sizes:

  Optional positive parcel masses.

- mapper:

  Fixed mapper specification used for every candidate/direction.

- method:

  Selection channel.

- reference_subject:

  Required subject index or ID for explicit mode.

- subject_ids:

  Optional subject IDs, primarily for geometry/explicit selection when
  no feature object supplies them.

- provenance:

  Optional explicit-selection provenance.

- tie_tolerance:

  Relative tolerance for deterministic score ties; ties are broken by
  lexicographic subject ID, never input position.

## Value

An immutable `dkge_reference_selection` object.

## Details

Geometry-only selection is ancillary and uses the same symmetric
reconstruction criterion on coordinates with the functional cost
disabled. Explicit selection is recorded as fixed by the caller. Reusing
the training functional features for both fit and score is available
only through the clearly ineligible `"descriptive_training"` method.
