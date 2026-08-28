# Summarize DKGE diagnostics

Provides a compact list of variance explained, subject weights, rank,
and kernel-support metadata for quick inspection.

## Usage

``` r
dkge_diagnostics(fit)
```

## Arguments

- fit:

  A `dkge` object.

## Value

List with variance table, subject weights, rank info, and scalar kernel
diagnostics, including numerical rank/nullity, condition, status,
participation-ratio effective rank, leading-eigenvalue share, and any
model-level spatial-regularization provenance.

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
#> [5] "kernel"        "n_subjects"    "voxel_weights" "weight_spec"  
#> [9] "spatial"      
```
