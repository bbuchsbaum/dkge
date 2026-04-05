# Fit DKGE bases on CPCA-filtered components

Computes separate DKGE bases on the CPCA design-aligned or residual
parts of the compressed covariance. Helpful for decomposing the latent
space into interpretable sub-structures.

## Usage

``` r
dkge_fit_cpca(
  fit,
  blocks = NULL,
  T = NULL,
  part = c("design", "resid", "both"),
  rank = NULL,
  ridge = 0
)
```

## Arguments

- fit:

  A \`dkge\` object from \[dkge_fit()\] or \[dkge()\].

- blocks:

  Optional integer vector of effect indices defining the subspace (used
  when \`T\` is not supplied).

- T:

  Optional explicit qxq0 basis matrix; overrides \`blocks\` when given.

- part:

  Which portion to return: "design", "resid", or "both".

- rank:

  Target rank for the returned bases (defaults to \`ncol(fit\$U)\`).

- ridge:

  Optional ridge term applied before the eigen decompositions.

## Value

List containing the requested bases/eigenvalues and the projector.
