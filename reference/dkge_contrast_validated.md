# Dual-path DKGE contrasts with coverage diagnostics

Computes DKGE contrasts twice on the same fold structure: once with an
observed-only coverage policy (typically \`missingness = "rescale"\`)
and once with a completed/penalised policy (e.g., \`missingness =
"shrink"\`). Coverage metadata from the fit and folds is returned
together with simple sensitivity summaries.

## Usage

``` r
dkge_contrast_validated(
  fit,
  contrasts,
  folds = NULL,
  ridge = 0,
  parallel = FALSE,
  verbose = FALSE,
  align = FALSE,
  observed_missingness = c("rescale", "mask", "none"),
  completed_missingness = c("shrink", "rescale", "none"),
  observed_args = list(),
  completed_args = list(),
  ...
)
```

## Arguments

- fit:

  Fitted \[dkge()\] object.

- contrasts:

  Contrast specification accepted by \[dkge_contrast()\].

- folds:

  Fold definition (integer \`k\`, \`dkge_folds\`, data frame, or list).

- ridge:

  Ridge added to held-out compressed matrices.

- parallel:

  Logical; whether to parallelise held-out projections.

- verbose:

  Logical; emit progress messages.

- align:

  Logical; align held-out bases across folds.

- observed_missingness:

  Coverage policy for the observed-only path.

- completed_missingness:

  Coverage policy for the completed path.

- observed_args:

  Optional list of arguments for the observed policy (e.g.,
  \`list(min_pairs = 2)\` for masking).

- completed_args:

  Optional list of arguments for the completed policy.

- ...:

  Additional arguments forwarded to the underlying cross-fitting helper.

## Value

A list with class \`dkge_contrast_validated\` containing: -
\`observed\`, \`completed\`: outputs from the respective paths. -
\`summary\`: data frame with weighted means and sensitivity indices. -
\`provenance\`: coverage metadata (effect IDs, subject masks, per-fold
pair-count matrices).

## Examples

``` r
# \donttest{
toy <- dkge_sim_toy(
  factors = list(cond = list(L = 3)),
  active_terms = "cond", S = 4, P = 15, snr = 5
)
fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 2)
#> Warning: Argument 'kernel' is deprecated; use 'K' instead.
q <- nrow(fit$U)
c_vec <- rep(0, q)
c_vec[2] <- 1
c_vec[3] <- -1
res <- dkge_contrast_validated(fit,
                               contrasts = list(cond = c_vec),
                               folds = 3)
res$summary
#>   contrast estimate_observed estimate_completed delta sensitivity
#> 1     cond        0.06764641         0.06764641     0           0
# }
```
