# Construct over-ranked functional-alignment features

Functional correspondence is represented in compact coordinates for
`image(K)`: if `K = L_K L_K'`, subject features are
`G_s = Bmodel_s' L_K`. This retains `rank(K)` coordinates independently
of the low fitted estimation rank. With no parcel preprocessing,
`Btil_s' Khalf = G_s V_K'`, so all Euclidean feature distances are
exactly preserved by the compact representation.

## Usage

``` r
dkge_alignment_features(
  fit,
  contrast_obj,
  feature_source = c("independent", "same_data_residualized"),
  independent_betas = NULL,
  independent_data_hash = NULL,
  effect_noise_cov = NULL,
  control = NULL
)
```

## Arguments

- fit:

  A fitted `dkge` object.

- contrast_obj:

  A `dkge_contrasts` result to which the features are immutably bound.

- feature_source:

  Either independent features (recommended/default) or
  covariance-residualized same-data features.

- independent_betas:

  For independent mode, one raw q-by-P beta matrix per fitted subject,
  on the same effect and parcel ordering.

- independent_data_hash:

  Stable non-empty identifier for the independent acquisition/training
  data. Supplying the primary fitted betas is rejected.

- effect_noise_cov:

  Optional list of q-by-q subject effect covariances for same-data
  residualization. By default these are read from `fit$subjects`.

- control:

  Numerical gates from
  [`dkge_alignment_feature_control()`](https://bbuchsbaum.github.io/dkge/reference/dkge_alignment_feature_control.md).

## Value

An immutable `dkge_alignment_features` object. Its `$features` list can
be supplied directly to the transport fitter through the typed object.

## Details

`feature_source = "independent"` is the default and requires separate
beta maps plus an explicit independent-data identifier. The independent
maps are transformed by the design-only fitted ruler but never by
beta-adaptive voxel weights or fitted spatial smoothing.

`feature_source = "same_data_residualized"` is an opt-in approximate
mode. It reconstructs each LOSO contrast as `G_s gamma_s`, then
conditions the whole feature matrix on the joint contrast family using
`Sigma_G,s = L_K' R' Lambda_s R L_K`. Its covariance oracle is exact
under the stated separable model conditional on fixed `gamma_s`; because
`gamma_s` depends on the estimated rank-truncated held-out basis,
frozen-plan sign-flip inference remains approximate. Full re-estimation
under every valid sign action is the exact same-data reference. The
preprocessing receipt verifies the implemented effect/parcel
transformations; it records separability as a model assumption and does
not empirically establish it. Both feature modes require contrasts from
a pooled, non-CPCA fit. CPCA and JD estimators fail closed until their
exact training-fold estimator can be replayed.
