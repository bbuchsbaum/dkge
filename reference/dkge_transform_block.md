# Preprocess a subject block into DKGE training space

Applies the same transformations used during fitting (pooled design
Cholesky factor, kernel whitening, optional spatial weights, and subject
weighting) to a subject's beta matrix.

## Usage

``` r
dkge_transform_block(fit, B_s, Omega_s = NULL, w_s = NULL)
```

## Arguments

- fit:

  A \`dkge\` object.

- B_s:

  qxP matrix of subject betas.

- Omega_s:

  Optional weights (vector length P or PxP matrix) matching the columns
  of \`B_s\`.

- w_s:

  Optional subject-level weight (defaults to 1 when omitted).

## Value

qxP matrix in the DKGE training space.
