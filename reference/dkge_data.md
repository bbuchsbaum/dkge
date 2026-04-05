# Bundle subject-level inputs for DKGE

Bundle subject-level inputs for DKGE

## Usage

``` r
dkge_data(betas, designs = NULL, omega = NULL, subject_ids = NULL)
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
```
