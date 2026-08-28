# Attach functional template features to a reference support

The template owns target functional features and masses. It does not own
coordinates or rendering topology. `feature_source = "independent"` is
the recommended functional source. DKGE verifies that its provenance
receipt includes a non-empty `independent_data_hash`; the caller remains
responsible for the scientific claim that the identified
acquisition/training channel is statistically independent. Missing
receipt evidence fails closed.

## Usage

``` r
dkge_functional_template(
  support,
  features,
  masses = NULL,
  feature_source = c("independent", "geometry_only", "same_data_residualized",
    "fully_recomputed", "descriptive_adaptive"),
  provenance = NULL,
  fitting = NULL,
  eligibility = NULL,
  alignment_features_hash = NULL,
  subject_ids = NULL,
  reference_selection = NULL
)
```

## Arguments

- support:

  A
  [`dkge_reference_support()`](https://bbuchsbaum.github.io/dkge/reference/dkge_reference_support.md)
  object.

- features:

  Finite target feature matrix with one row per support point.

- masses:

  Positive target masses. Defaults to equal masses.

- feature_source:

  Provenance category for correspondence features.

- provenance:

  Optional evidence, including `independent_data_hash` for an
  independently trained template.

- fitting, eligibility:

  Reserved internal receipts. Public callers must leave them `NULL`;
  fitted templates are created only by
  [`dkge_fit_functional_template()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit_functional_template.md).

- alignment_features_hash:

  Optional structural hash of the typed source feature object from which
  the template was learned.

- subject_ids:

  Optional identities of subjects used to learn the template.

- reference_selection:

  Optional immutable reference-selection receipt identifying how a
  subject-derived display support was chosen.

## Value

An immutable descriptive template/initializer. It cannot mint fitted
correspondence or inferential eligibility.
