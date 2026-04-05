# Consensus K-orthonormal basis (K-Procrustes mean)

Iteratively aligns bases, takes a weighted average, and retraction via
K-orthonormalization until convergence.

## Usage

``` r
dkge_consensus_basis_K(
  U_list,
  K,
  weights = NULL,
  Kroots = NULL,
  max_iter = 50,
  tol = 1e-06,
  allow_reflection = TRUE
)
```

## Arguments

- U_list:

  list of qxr K-orthonormal bases

- K:

  qxq design kernel

- weights:

  optional numeric weights (default equal)

- Kroots:

  optional precomputed kernel roots

- max_iter:

  maximum iterations

- tol:

  convergence tolerance on 1 - mean principal cosines

- allow_reflection:

  passed to alignment step

## Value

list(U, iters, converged, gaps, scores)
