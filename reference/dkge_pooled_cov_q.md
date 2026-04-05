# Pooled design-space covariance and Cholesky factor

Computes the qxq pooled design covariance and the corresponding Cholesky
factor needed for fast kernel alignment screening.

## Usage

``` r
dkge_pooled_cov_q(B_list, X_list, Omega_list = NULL)
```

## Arguments

- B_list:

  List of qxP subject beta matrices.

- X_list:

  List of Txq subject design matrices.

- Omega_list:

  Optional list of spatial weights.

## Value

List containing \`C\` (pooled covariance in the ruler metric), \`R\`
(upper-triangular Cholesky factor), and \`G\` (pooled design Gram
matrix).

## Examples

``` r
toy <- dkge_sim_toy(
  factors = list(cond = list(L = 3)),
  active_terms = "cond", S = 3, P = 10, snr = 5
)
pooled <- dkge_pooled_cov_q(toy$B_list, toy$X_list)
dim(pooled$C)
#> [1] 4 4
```
