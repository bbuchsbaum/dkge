# K-orthogonal projector onto span(T) in effect space

Builds both the K-selfadjoint projector operating in effect coordinates
and the Euclidean projector in the \\K^{1/2}\\ metric. These projectors
are useful for CPCA-style splits of the compressed covariance.

## Usage

``` r
dkge_projector_K(T, K)
```

## Arguments

- T:

  qxq0 matrix whose columns span the subspace of interest.

- K:

  qxq positive semi-definite design kernel.

## Value

List with entries \`P_K\` (effect-space projector) and \`P_hat\`
(projector in the \\K^{1/2}\\ metric).

## Details

All computations take place in the metric induced by the design kernel
\`K\`. Choosing a non-identity kernel therefore changes the projector:
smooth kernels diffuse energy across correlated effects, whereas a
diagonal kernel keeps the split tied to the raw effect indices. When the
design-aligned subspace is not well represented by coordinate axes,
supply \`cpca_T\` with columns that already capture the kernel geometry
so the CPCA basis respects those relationships.

## Examples

``` r
K <- diag(5)
T <- matrix(rnorm(10), 5, 2)
P <- dkge_projector_K(T, K)
dim(P$P_K)
#> [1] 5 5
```
