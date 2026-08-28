# Positive-semidefinite kernel roots

Positive-semidefinite kernel roots

## Usage

``` r
kernel_roots(K, jitter = 0, tol = 1e-10)

dkge_kernel_roots(K, jitter = 0, tol = 1e-10)
```

## Arguments

- K:

  Positive semi-definite kernel matrix.

- jitter:

  Optional non-negative diagonal regularization added explicitly to `K`
  before computing roots. The default, zero, preserves an exact null
  space. Use a positive value only when changing the kernel geometry is
  scientifically intended.

- tol:

  Relative eigentolerance used to define numerical support.

## Value

List with the exact square root, Moore–Penrose inverse square root,
eigenstructure, support projectors, numerical rank diagnostics,
participation-ratio effective rank, and leading-eigenvalue share.
