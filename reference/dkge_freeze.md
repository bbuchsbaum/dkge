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
#>              [,1]        [,2]
#> [1,]  0.826322086 0.903668379
#> [2,]  0.031522415 0.009575469
#> [3,] -1.277794194 1.169044492
#> [4,]  0.061535214 0.005058690
#> [5,]  0.002274033 0.005972575
#> 
#> $K
#>                  A            B1            B2          A:B1          A:B2
#> A     6.666667e-01 -3.021212e-18  4.399644e-18 -2.400376e-18  5.693294e-17
#> B1   -3.021212e-18  3.333333e-01 -4.744909e-18  1.242446e-18  2.788273e-18
#> B2    4.399644e-18 -4.744909e-18  3.333333e-01 -2.005616e-18  1.164681e-17
#> A:B1 -2.400376e-18  1.242446e-18 -2.005616e-18  5.555556e-02 -6.170064e-18
#> A:B2  5.693294e-17  2.788273e-18  1.164681e-17 -6.170064e-18  5.555556e-02
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
#> [1] "A"    "B1"   "B2"   "A:B1" "A:B2"
#> 
#> attr(,"class")
#> [1] "dkge_model"
```
