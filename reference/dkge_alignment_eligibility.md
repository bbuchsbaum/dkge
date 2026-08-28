# Classify inferential eligibility of a fitted alignment

Eligibility is composite. Independent correspondence features are not an
exactness guarantee when contrast values still depend on a same-data,
rank-truncated latent span. The frozen v1 court therefore makes
`same_data_rank_truncated` approximate even with independent features.

## Usage

``` r
dkge_alignment_eligibility(
  feature_source = .dkge_alignment_feature_sources,
  estimator_source = .dkge_alignment_estimator_sources,
  recompute_under_null = FALSE,
  source_verified = FALSE,
  reference_selection_status = c("eligible", "approximate", "ineligible")
)
```

## Arguments

- feature_source:

  Correspondence-feature provenance category.

- estimator_source:

  Provenance of the contrast-generating latent span.

- recompute_under_null:

  Logical; whether every beta-dependent quantity and correspondence
  operator is recomputed under each valid null action.

- source_verified:

  Logical; whether the claimed feature source has the receipt evidence
  required by its contract. This validates provenance, not the
  scientific independence of an external data-generating process.

- reference_selection_status:

  Inferential status of the reference-support selection channel.

## Value

A typed `dkge_alignment_eligibility` record.
