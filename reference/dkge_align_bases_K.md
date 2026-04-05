# Align multiple bases to a reference

Align multiple bases to a reference

## Usage

``` r
dkge_align_bases_K(
  U_list,
  K,
  ref = 1L,
  allow_reflection = TRUE,
  weights = NULL
)
```

## Arguments

- U_list:

  list of qxr bases

- K:

  qxq design kernel

- ref:

  reference basis or index (default 1)

- allow_reflection:

  logical passed to \`dkge_procrustes_K\`

- weights:

  optional numeric weights (stored for convenience)

## Value

list(U_aligned, R, Uref, score, weights)
