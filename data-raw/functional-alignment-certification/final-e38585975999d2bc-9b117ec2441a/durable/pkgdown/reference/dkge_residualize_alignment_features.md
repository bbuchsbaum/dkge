# Remove a contrast family from functional features conditionally

For feature row `g`, generating directions `Gamma`, and
effect-coordinate covariance `Sigma_G`, this computes
`f = g - Sigma_G Gamma (Gamma' Sigma_G Gamma)^-1 Gamma' g`. With row
matrices this is `F = G - (G Gamma) C^-1 Gamma' Sigma_G`. Under a
separable spatial/effect covariance, the returned analytic oracle
verifies `Cov(F_p, v_q) = 0` for every parcel pair, up to the arbitrary
spatial scalar `rho[p,q]`.

## Usage

``` r
dkge_residualize_alignment_features(G, gamma, Sigma_G, control = NULL)
```

## Arguments

- G:

  `P` by `k` kernel-image feature matrix.

- gamma:

  `k` by `m` matrix of contrast-generating directions.

- Sigma_G:

  `k` by `k` effect-coordinate covariance.

- control:

  Numerical gates from
  [`dkge_alignment_feature_control()`](https://bbuchsbaum.github.io/dkge/reference/dkge_alignment_feature_control.md).

## Value

A typed result containing residual features, reconstructed contrast
values, the checked conditional coefficient, covariance oracle, and
diagnostics.
