# Robust K-orthonormalization

Ensures the columns of \`W\` are orthonormal with respect to the design
kernel metric: U^T K U = I.

## Usage

``` r
dkge_k_orthonormalize(W, K, Kroots = NULL)
```

## Arguments

- W:

  qxr matrix (columns = basis vectors)

- K:

  qxq design kernel (PSD)

- Kroots:

  optional precomputed kernel roots from \`.dkge_kernel_roots\`

## Value

qxr matrix with K-orthonormal columns

## Examples

``` r
K <- diag(5)
W <- matrix(rnorm(10), 5, 2)
U <- dkge_k_orthonormalize(W, K)
# Verify K-orthonormality
round(t(U) %*% K %*% U, 10)
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
```
