# Control parameters for DKGE joint diagonalisation

Control parameters for DKGE joint diagonalisation

## Usage

``` r
dkge_jd_control(
  maxit = 200L,
  step = 0.5,
  step_decay = 0.5,
  step_min = 1e-06,
  tol = 1e-07,
  linesearch = TRUE,
  armijo = 1e-04,
  verbose = FALSE,
  record = FALSE
)
```

## Arguments

- maxit:

  Maximum number of optimization iterations.

- step:

  Initial step size for the projected gradient descent.

- step_decay:

  Multiplicative decay applied during backtracking.

- step_min:

  Minimum step size allowed before terminating.

- tol:

  Convergence tolerance applied to both gradient norm and off-diagonal
  energy.

- linesearch:

  Logical; when \`TRUE\` perform Armijo backtracking.

- armijo:

  Armijo condition constant used during backtracking.

- verbose:

  Logical; emit per-iteration progress when \`TRUE\`.

- record:

  Logical; store per-iteration diagnostics.

## Value

A list of control parameters.

## Examples

``` r
ctrl <- dkge_jd_control(maxit = 10, verbose = FALSE)
names(ctrl)
#> [1] "maxit"      "step"       "step_decay" "step_min"   "tol"       
#> [6] "linesearch" "armijo"     "verbose"    "record"    
```
