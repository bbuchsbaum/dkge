# Robust kernel roots

Robust kernel roots

## Usage

``` r
kernel_roots(K, jitter = 1e-10)

dkge_kernel_roots(K, jitter = 1e-10)
```

## Arguments

- K:

  Positive semi-definite kernel matrix.

- jitter:

  Small diagonal jitter added before inversion.

## Value

List with \`Khalf\`, \`Kihalf\`, eigenvalues, eigenvectors, and basic
diagnostics.
