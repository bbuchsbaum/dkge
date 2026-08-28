# Group inference on already aligned subject maps

This is the typed inference boundary for functional alignment. It
consumes one row per subject on one identified support, never fits
correspondence, and refuses ineligible alignment objects by default.
Fit-level MFA pooling weights are not inherited; the aligned object
records its group weighting.

## Usage

``` r
dkge_infer_aligned(
  aligned_maps,
  inference = c("signflip", "parametric"),
  correction = c("maxT", "fdr", "bonferroni", "none"),
  n_perm = 2000,
  alpha = 0.05,
  allow_approximate_alignment = FALSE
)
```

## Arguments

- aligned_maps:

  An operator-bound
  [`dkge_aligned_maps()`](https://bbuchsbaum.github.io/dkge/reference/dkge_aligned_maps.md)
  object returned by
  [`dkge_transport_contrasts_to_reference()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_contrasts_to_reference.md)
  or
  [`dkge_align_to_template()`](https://bbuchsbaum.github.io/dkge/reference/dkge_align_to_template.md).

- inference:

  `"signflip"` or `"parametric"`.

- correction:

  Multiple-testing correction.

- n_perm:

  Number of sign flips.

- alpha:

  Significance level.

- allow_approximate_alignment:

  Permit an explicitly labelled approximate alignment.
  Descriptive/ineligible objects remain errors.

## Value

A `dkge_inference` object carrying the aligned estimand and weighting.

## Details

The subject is the sampling and exchangeability unit. Locations and
contrasts within a subject row are not independent replicates. Sign-flip
inference requires joint row-sign invariance under the null: one sign is
applied to each subject simultaneously across every prespecified
contrast-by-location cell. With `correction = "maxT"`, each null draw
takes the maximum statistic over that complete joint family. Marginal
permutation p-values use the same shared sign matrix, and FDR or
Bonferroni correction is applied once after flattening that complete
family (then split back into the original contrast shapes).
`correction = "none"` leaves those joint-draw marginal p-values
unadjusted.

Correspondence is fixed before this boundary. Independent-feature and
same-data-residualized alignments can still depend on a same-cohort,
rank-truncated latent span. The frozen calibration court's full-pipeline
re-estimation exact null-action comparator passed, but the
negative-control promotion gate did not; these frozen-plan modes
therefore remain `"approximate"`, fail closed by default, and require an
explicit override.

Parametric inference performs cellwise one-sample t tests and assumes
independent subjects plus an adequate normal approximation; it does not
define a max-T randomization distribution. Current aligned group
inference supports equal subject weighting only. Fit-level MFA weights
and explicit descriptive/bootstrap weights are not imported into this
group test.
