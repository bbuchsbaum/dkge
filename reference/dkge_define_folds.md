# Define folds for K-fold cross-validation

Create fold assignments for subjects or time points in DKGE analysis.
Supports multiple strategies including random subject-level folds,
time-based splits, run-based partitions, and custom user-defined
assignments.

## Usage

``` r
dkge_define_folds(
  fit,
  type = c("subject", "time", "run", "custom"),
  k = 5,
  runs = NULL,
  assignments = NULL,
  seed = NULL,
  align = FALSE,
  ...
)
```

## Arguments

- fit:

  A \`dkge\` object or \`dkge_data\` bundle

- type:

  Type of fold definition: - \`"subject"\`: Random assignment of
  subjects to folds (default) - \`"time"\`: Split each subject's time
  series into temporal blocks - \`"run"\`: Use experimental runs as
  natural folds - \`"custom"\`: User provides fold assignments directly

- k:

  Number of folds (ignored for type="custom")

- runs:

  For type="run", a list of run indicators per subject

- assignments:

  For type="custom", a list of fold assignments

- seed:

  Random seed for reproducible fold assignment

- align:

  Logical; if TRUE (default) compute Procrustes alignment/consensus when
  folds are evaluated.

- ...:

  Additional arguments for specific fold types

## Value

A \`dkge_folds\` object containing: - \`type\`: The fold type used -
\`k\`: Number of folds - \`assignments\`: List specifying which data
belongs to each fold - \`metadata\`: Additional information about fold
creation

## Details

Different fold types serve different purposes:

\*\*Subject-level folds\*\* (\`type = "subject"\`): Assigns entire
subjects to folds. This maintains subject independence and is
appropriate when subjects are exchangeable. Each fold will have
approximately S/k subjects.

\*\*Time-based folds\*\* (\`type = "time"\`): Splits each subject's time
series into k temporal blocks. Useful for assessing temporal stability
or when early vs late responses differ. Requires access to original time
series dimensions.

\*\*Run-based folds\*\* (\`type = "run"\`): Uses experimental runs as
natural folds. Common in fMRI where runs provide natural breaks.
Requires run indicators.

\*\*Custom folds\*\* (\`type = "custom"\`): Full control over fold
assignments. Supply a list where each element specifies indices for that
fold.

## Examples

``` r
# \donttest{
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 6, P = 20, snr = 5
)
fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 2)
#> Warning: Argument 'kernel' is deprecated; use 'K' instead.

# Random subject-level 5-fold
folds <- dkge_define_folds(fit, type = "subject", k = 5)

# Custom fold specification (explicit holdouts)
folds_custom <- dkge_define_folds(
  fit,
  type = "custom",
  assignments = list(fold1 = 1:2, fold2 = 3:4, fold3 = 5:6)
)
folds_custom$k
#> [1] 3
# }
```
