# Freeze a DKGE fit into a compact model for prediction

Freeze a DKGE fit into a compact model for prediction

## Usage

``` r
dkge_freeze(fit)
```

## Arguments

- fit:

  a dkge or dkge_stream object

## Value

list with U, K, R and class 'dkge_model'

## Examples

``` r
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 3, P = 20, snr = 5
)
fit <- dkge_fit(toy$B_list, toy$X_list, toy$K, rank = 2)
model <- dkge_freeze(fit)
print(model)
#> $U
#>             [,1]         [,2]
#> [1,]  1.11632046 -0.500964364
#> [2,] -0.05817352  0.004677653
#> [3,]  0.05618580  0.015020009
#> [4,] -0.08978873 -0.045960778
#> [5,]  1.73165688  3.871022992
#> 
#> $K
#>               [,1]          [,2]          [,3]          [,4]          [,5]
#> [1,]  6.666667e-01 -3.021212e-18  4.399644e-18 -2.400376e-18  5.693294e-17
#> [2,] -3.021212e-18  3.333333e-01 -4.744909e-18  1.242446e-18  2.788273e-18
#> [3,]  4.399644e-18 -4.744909e-18  3.333333e-01 -2.005616e-18  1.164681e-17
#> [4,] -2.400376e-18  1.242446e-18 -2.005616e-18  5.555556e-02 -6.170064e-18
#> [5,]  5.693294e-17  2.788273e-18  1.164681e-17 -6.170064e-18  5.555556e-02
#> 
#> $R
#>          effect1  effect2  effect3  effect4  effect5
#> effect1 1.732051 0.000000 0.000000 0.000000 0.000000
#> effect2 0.000000 1.732051 0.000000 0.000000 0.000000
#> effect3 0.000000 0.000000 1.732051 0.000000 0.000000
#> effect4 0.000000 0.000000 0.000000 1.732051 0.000000
#> effect5 0.000000 0.000000 0.000000 0.000000 1.732051
#> 
#> $effects
#> [1] "effect1" "effect2" "effect3" "effect4" "effect5"
#> 
#> attr(,"class")
#> [1] "dkge_model"
```
