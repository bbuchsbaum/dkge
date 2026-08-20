# Range-space roots for a positive-semidefinite kernel

Computes the true symmetric square root and the Moore-Penrose inverse
square root. Eigenvalues at or below the applied absolute-plus-relative
tolerance remain exactly zero; null directions are never jittered into
the metric.

## Usage

``` r
kernel_roots(
  K,
  jitter = NULL,
  absolute_tolerance = .Machine$double.xmin,
  relative_tolerance = 1e-08
)

dkge_kernel_roots(
  K,
  jitter = NULL,
  absolute_tolerance = .Machine$double.xmin,
  relative_tolerance = 1e-08
)
```

## Arguments

- K:

  Positive semi-definite kernel matrix.

- jitter:

  Deprecated absolute truncation threshold. When supplied it is used as
  \`absolute_tolerance\`; eigenvalues below it are discarded, not raised
  to it.

- absolute_tolerance:

  Non-negative absolute spectral tolerance.

- relative_tolerance:

  Non-negative tolerance relative to the largest kernel eigenvalue.

## Value

List with \`Khalf\`, \`Kihalf\`, eigenvalues, eigenvectors, numerical
rank/nullity, retained range, and the applied tolerances.
