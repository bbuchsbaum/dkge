# K-orthogonal Procrustes alignment

Aligns basis \`U\` to reference \`Uref\` by solving max_R tr((U_ref^T K
U) R) subject to R^T R = I.

## Usage

``` r
dkge_procrustes_K(Uref, U, K, allow_reflection = TRUE)
```

## Arguments

- Uref:

  reference basis (qxr)

- U:

  basis to align (qxr)

- K:

  qxq design kernel

- allow_reflection:

  logical; if FALSE, forces det(R)=+1

## Value

list(U_aligned, R, d=sum(singular values), cosines, det)
