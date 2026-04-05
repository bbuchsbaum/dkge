# Cross-validated classification on DKGE effect patterns

Cross-validated classification on DKGE effect patterns

## Usage

``` r
dkge_classify(
  fit,
  targets,
  y = NULL,
  method = c("lda", "logit"),
  folds = NULL,
  lambda = NULL,
  metric = c("accuracy", "logloss"),
  mode = c("auto", "cell", "cell_cross", "delta"),
  standardize_within_fold = NULL,
  residualize = TRUE,
  collapse = NULL,
  restrict_factors = NULL,
  n_perm = 0L,
  scope = NULL,
  class_weights = c("none", "balanced", "inverse"),
  ridge = 0,
  control = NULL,
  blocks = NULL,
  parallel = FALSE,
  verbose = FALSE,
  seed = NULL
)
```

## Arguments

- fit:

  dkge object.

- targets:

  Target specification consumed by \[dkge_targets()\] or a list of
  \`dkge_target\` objects.

- y:

  Optional subject-level labels for delta-mode targets. Can be a vector
  (recycled across all delta targets) or a list matching \`targets\`;
  values are coerced to factors using each target's \`class_labels\`.

- method:

  Classifier backend: "lda" (default) or "logit" (one-vs-rest ridge
  logistic regression).

- folds:

  Cross-fitting specification. \`NULL\` (default) performs LOSO. Integer
  values request subject-level K-fold. A \`dkge_folds\` object is also
  accepted. Other coercible inputs are routed through
  \[as_dkge_folds()\].

- lambda:

  Optional ridge parameter passed to the classifier backend.

- metric:

  Evaluation metric(s); defaults to c("accuracy", "logloss").

- mode:

  Decoding mode: "auto" (default), "cell", "cell_cross", or "delta".
  Cell-cross trains on held-in subjects and tests generalisation to the
  held-out subject.

- standardize_within_fold:

  Logical indicating whether to z-score features using training data
  inside each fold. When \`NULL\` (default), standardisation is enabled
  automatically for \`mode = "cell_cross"\` and disabled otherwise.

- residualize:

  Forwarded to \[dkge_targets()\].

- collapse:

  Forwarded to \[dkge_targets()\] for factor collapsing.

- restrict_factors:

  Optional factor subset for \`spec = "fullcell"\`.

- n_perm:

  Number of permutations for empirical p-values (default 0).

- scope:

  Override permutation scope (otherwise target scope used).

- class_weights:

  Class weighting scheme ("none", "balanced", "inverse").

- ridge:

  Optional ridge added when recomputing held-out bases (default 0).

- control:

  Optional list of advanced controls (power users). Recognised entries:
  \`lambda_grid\` (numeric vector of candidate penalties) and
  \`lambda_fun\` (function returning a lambda per target/fold with
  signature \`function(target, fold, method, default)\`). Defaults to
  \`NULL\`, leaving the standard \`lambda\` behaviour unchanged.

- blocks:

  Optional vector identifying within-subject blocks (e.g., runs or
  sessions) used to constrain permutations. Length must match the number
  of subjects in \`fit\` when supplied.

- parallel:

  Logical; reserved for future parallelism hooks.

- verbose:

  Logical; print progress messages.

- seed:

  Optional random seed applied before permutations.

## Details

Anchor-based fits produced by \[dkge_anchor_fit()\] do not retain the
design-factor metadata that \`dkge_targets()\` expects. In that setting
you must supply explicit weight matrices (rows = classes, columns =
effects) or ready-made \[\`dkge_target\`\] objects—helpers such as
\[dkge_anchor_targets_from_prototypes()\] and
\[dkge_anchor_targets_from_directions()\] can be used to construct them.

## Examples

``` r
# Simulate toy data with design kernel (includes kernel_info for targets)
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 6, P = 20, snr = 5
)
# Use design_kernel() result to retain factor metadata
kern <- design_kernel(
  factors = list(A = list(L = 2), B = list(L = 3)),
  basis = "effect"
)
fit <- dkge(toy$B_list, toy$X_list, kernel = kern, rank = 2)
#> Warning: Argument 'kernel' is deprecated; use 'K' instead.

# Decode factor A using cell-mode classification
# \donttest{
clf <- dkge_classify(fit, targets = ~A, method = "lda")
clf
#> DKGE Classification
#> --------------------
#> Targets: 1
#> Classifier: lda
#> Metrics: accuracy, logloss
#> Permutations: 0
#>   A: accuracy=1.000, logloss=0.000
# }
```
