# Numerical gates for functional-alignment features

These controls deliberately fail closed when a low-dimensional contrast
family consumes the useful functional feature space or when the
conditional residualization problem is not identified. No generalized
inverse is used to conceal a singular contrast covariance.

## Usage

``` r
dkge_alignment_feature_control(
  kernel_tolerance = 1e-10,
  rank_tolerance = 1e-09,
  max_condition = 1e+08,
  min_feature_rank = 2L,
  min_residual_rank = 2L,
  min_effective_rank = 1.25,
  min_retained_energy = 0.05,
  max_row_collapse_fraction = 0.25,
  row_norm_tolerance = 1e-07,
  reconstruction_tolerance = 1e-08,
  covariance_orthogonality_tolerance = 1e-08
)
```

## Arguments

- kernel_tolerance:

  Relative eigentolerance defining `image(K)`.

- rank_tolerance:

  Relative singular/eigenvalue tolerance for numerical ranks.

- max_condition:

  Largest permitted condition number for the contrast covariance.

- min_feature_rank:

  Minimum observed rank for independent features.

- min_residual_rank:

  Minimum numerical rank after residualization.

- min_effective_rank:

  Minimum participation-ratio rank after residualization.

- min_retained_energy:

  Minimum observed residual-to-original squared Frobenius energy ratio.

- max_row_collapse_fraction:

  Largest permitted fraction of initially nonzero rows that collapse
  numerically after residualization.

- row_norm_tolerance:

  Relative row-norm threshold used to diagnose collapse.

- reconstruction_tolerance:

  Relative tolerance for reconstructing each LOSO contrast as
  `G_s gamma_s`.

- covariance_orthogonality_tolerance:

  Relative tolerance for the analytic conditional-covariance oracle.

## Value

A typed alignment-feature control object.
