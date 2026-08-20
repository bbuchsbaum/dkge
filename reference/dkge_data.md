# Bundle subject-level inputs for DKGE

Bundle subject-level inputs for DKGE

## Usage

``` r
dkge_data(
  betas,
  designs = NULL,
  omega = NULL,
  subject_ids = NULL,
  effects = NULL
)
```

## Arguments

- betas:

  List of subject records (matrices or \`dkge_subject\` objects)

- designs:

  Optional list of design matrices (ignored when \`betas\` already
  contain subjects)

- omega:

  Optional list of cluster weights

- subject_ids:

  Optional subject identifiers

- effects:

  Optional character vector pinning the global effect order. It must be
  a permutation of the effects observed across subjects (the union when
  coverage is partial). Supply \`dkge_effect_grid()\$cell_labels\` here
  so that the bundle, the design kernel, and the grid all index effects
  identically; without it the union is ordered by first appearance
  across subjects, which depends on the order of \`betas\`.

## Value

An object of class \`dkge_data\`

## Examples

``` r
betas <- replicate(3, matrix(rnorm(5 * 80), 5, 80), simplify = FALSE)
designs <- replicate(
  3,
  matrix(
    rnorm(150 * 5), 150, 5,
    dimnames = list(NULL, paste0("eff", 1:5))
  ),
  simplify = FALSE
)
data <- dkge_data(betas, designs)
data$effects
#> [1] "eff1" "eff2" "eff3" "eff4" "eff5"

# Pin a specific effect order
dkge_data(betas, designs, effects = paste0("eff", 5:1))$effects
#> [1] "eff5" "eff4" "eff3" "eff2" "eff1"
```
