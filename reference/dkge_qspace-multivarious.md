# Fail-closed multivarious methods for q-space DKGE fits

Pair-normalized, debiased, ridged, CPCA, and JD fits do not inherit from
\`multiblock_biprojector\`. These methods reject physical block
projection instead of inventing incompatible loadings.

## Usage

``` r
# S3 method for class 'dkge_qspace'
project(x, new_data, ...)

# S3 method for class 'dkge_qspace'
project_block(x, new_data, block, least_squares = TRUE, ...)

# S3 method for class 'dkge_qspace'
transfer(x, new_data, from, to, opts = list(), ...)
```

## Arguments

- x:

  A \`dkge_qspace\` fit.

- new_data:

  Unused; required by the multivarious generics.

- ...:

  Unused.

- block, least_squares, from, to, opts:

  Unused arguments of the corresponding multivarious generics.
