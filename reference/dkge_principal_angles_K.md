# Principal angles between two DKGE bases in the K metric

Computes subspace principal angles between two bases using the DKGE
design kernel as the inner-product metric. This is useful for comparing
stratified or resampled DKGE fits after K-Procrustes alignment
diagnostics.

## Usage

``` r
dkge_principal_angles_K(U1, U2, K, orthonormalize = TRUE)
```

## Arguments

- U1, U2:

  Basis matrices with the same row dimension as \`K\`. Columns need not
  be K-orthonormal.

- K:

  Design kernel defining the inner product.

- orthonormalize:

  Logical; K-orthonormalize \`U1\` and \`U2\` first (default \`TRUE\`).

## Value

Numeric vector of principal angles in radians.

## Details

Principal angles are a property of the \*subspaces\* spanned by \`U1\`
and \`U2\`, not of the particular bases used to represent them. Both
inputs are therefore reduced to an ordinary Euclidean basis of \\K^{1/2}
U\\ via `svd(K^{1/2} U)` (the internal helper \`.dkge_k_span_basis()\`)
before the cosines \\\sigma_i(Q_1^\top Q_2)\\ are formed, so any two
bases of the same subspace (for example \`U\` and \`2 \* U\`) return
angles of zero. Set \`orthonormalize = FALSE\` only when the inputs are
already known to satisfy \\U^\top K U = I\\; the singular values are
otherwise not cosines and the clamp to \\\[0, 1\]\\ would silently
report perfect alignment.

The returned vector has one entry per shared dimension, ordered from the
smallest angle (most aligned direction) to the largest. Rank-deficient
inputs are truncated to the subspace they actually span, so the length
is \`min(rank_K(U1), rank_K(U2))\`, which can be smaller than
\`min(ncol(U1), ncol(U2))\`; a basis of numerical rank 0 is an error
rather than a fabricated subspace. With \`orthonormalize = FALSE\` the
inputs are taken at face value and all \`min(ncol(U1), ncol(U2))\`
singular values are returned.

## Examples

``` r
K <- diag(c(2, 1, 1))
U1 <- cbind(c(1, 0, 0), c(0, 1, 0))
# A rescaled basis of the same subspace has zero angles.
dkge_principal_angles_K(U1, 3 * U1, K)
#> [1] 0 0
# A K-orthogonal direction sits at 90 degrees.
dkge_principal_angles_K(U1[, 1, drop = FALSE], U1[, 2, drop = FALSE], K)
#> [1] 1.570796
```
