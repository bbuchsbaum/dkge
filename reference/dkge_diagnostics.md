# Summarise DKGE diagnostics

Provides a compact list of variance explained, subject weights, and rank
metadata for quick inspection.

## Usage

``` r
dkge_diagnostics(fit)
```

## Arguments

- fit:

  A \`dkge\` object.

## Value

List with variance table, subject weights, and rank info.

## Examples

``` r
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 3, P = 15, snr = 5
)
fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 2)
#> Warning: Argument 'kernel' is deprecated; use 'K' instead.
diag <- dkge_diagnostics(fit)
names(diag)
#> [1] "variance"      "weights"       "rank"          "q"            
#> [5] "n_subjects"    "voxel_weights" "weight_spec"  
```
