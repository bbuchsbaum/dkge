# K-orthogonal Procrustes alignment

Aligns basis \`U\` to reference \`Uref\` by solving max_R tr((U_ref^T K
U) R) subject to R^T R = I.

## Usage

``` r
dkge_procrustes_K(Uref, U, K, allow_reflection = TRUE)
```

## Arguments

- Uref:

  K-orthonormal reference basis (qxr)

- U:

  K-orthonormal basis to align (qxr)

- K:

  qxq finite, symmetric, positive-semidefinite design kernel

- allow_reflection:

  logical; if FALSE, forces det(R)=+1

## Value

A list with the aligned basis \`U_aligned\`, orthogonal rotation \`R\`,
achieved objective \`d\`, unconstrained objective \`unconstrained_d\`,
principal \`cosines\`, and rotation \`determinant\`. When reflection is
forbidden, \`d\` can be smaller than \`unconstrained_d\`.
