# Construct a DKGE global effect grid

Construct a DKGE global effect grid

## Usage

``` r
dkge_effect_grid(factors, scope = NULL, block_factors = NULL, sep = ":")
```

## Arguments

- factors:

  Named list of factor specifications. Each entry may be a design-kernel
  factor spec, a vector of level labels, or a scalar level count. List
  specs accept \`type\`, \`L\`, \`levels\`, \`values\`, and \`l\`
  (positive finite length-scale for ordinal/circular/continuous
  factors).

- scope:

  Named character vector or list assigning each factor to \`"within"\`
  or \`"between"\`. Defaults to \`"within"\` for all factors.

- block_factors:

  Optional factor names that should define independent blocks when
  passed to \[design_kernel()\].

- sep:

  Separator used when constructing cell labels.

## Value

A \`dkge_effect_grid\` object with factor specs, scopes, cells, labels,
and the default kernel terms implied by the factor order. A one-factor
grid has one default main-effect term because its full interaction is
the same term.

## Examples

``` r
grid <- dkge_effect_grid(
  factors = list(
    cue = c("valid", "invalid"),
    load = list(L = 3, type = "ordinal"),
    group = 2
  ),
  scope = c(group = "between"),
  block_factors = "group"
)
grid$cell_labels
#>  [1] "valid:load1:group1"   "valid:load1:group2"   "valid:load2:group1"  
#>  [4] "valid:load2:group2"   "valid:load3:group1"   "valid:load3:group2"  
#>  [7] "invalid:load1:group1" "invalid:load1:group2" "invalid:load2:group1"
#> [10] "invalid:load2:group2" "invalid:load3:group1" "invalid:load3:group2"
grid$scope
#>       cue      load     group 
#>  "within"  "within" "between" 

# Feed the grid to design_kernel() so the kernel indexes the same cells
kern <- design_kernel(grid, basis = "cell")
identical(rownames(kern$K), grid$cell_labels)
#> [1] TRUE
```
