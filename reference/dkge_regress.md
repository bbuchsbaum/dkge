# Multi-output regression on DKGE effects

Fits cross-validated regressors that map DKGE effect embeddings to
continuous targets supplied per effect (e.g., semantic vectors, ratings,
reference maps). The default engine uses ridge-regularised
multi-response regression via glmnet.

## Usage

``` r
dkge_regress(
  fit,
  Y,
  effect_ids = rownames(Y),
  folds = NULL,
  k = 5,
  features = "effects",
  engine = "glmnet",
  alpha = 0,
  standardize = TRUE,
  seed = 1,
  return_models = FALSE
)
```

## Arguments

- fit:

  A \`dkge\` object containing an effect-space basis (\`fit\$U\`).

- Y:

  Numeric matrix of shape \`(#effects) x (#outputs)\` with row names
  that identify the effects/items to be modelled.

- effect_ids:

  Optional character vector of effect identifiers to include (defaults
  to \`rownames(Y)\`).

- folds:

  Optional list giving held-out effect IDs for each fold. When \`NULL\`,
  random K-fold splits are generated.

- k:

  Number of folds to create when \`folds\` is \`NULL\` (default \`5\`).

- features:

  Either the string "effects" (use the fitted DKGE effect embeddings) or
  a custom function \`function(fit, ids, ctx)\` returning a matrix of
  predictors for the supplied effect IDs. The function may use the
  optional \`ctx\` argument to receive fold-specific context objects.

- engine:

  Regression engine. One of "glmnet" (default ridge multi-output
  regression), "lm" (ordinary least squares solved jointly across
  outputs), or a list with \`fit\`/\`predict\` closures: \`list(fit =
  function(X, Y) {...}, predict = function(model, X) {...})\`.

- alpha:

  Elastic-net mixing parameter passed to glmnet (0 for ridge, 1 for
  lasso).

- standardize:

  Logical; when \`TRUE\` (default) the glmnet engine standardises
  predictors internally.

- seed:

  Random seed used when generating folds.

- return_models:

  Logical; when \`TRUE\`, attach the fitted model for each fold to the
  returned object.

## Value

An object of class \`"dkge_regress"\` containing predictions (\`pred\`),
observed targets (\`truth\`), fold assignments (\`folds\`), summary
metrics (\`metrics\`), and optionally the per-fold models (\`models\`).

## Examples

``` r
# \donttest{
set.seed(1)
Q <- 4L
effects <- paste0("e", seq_len(Q))
betas <- replicate(3, matrix(rnorm(Q * 6), Q, 6, dimnames = list(effects, NULL)), simplify = FALSE)
designs <- replicate(3, diag(Q), simplify = FALSE)
designs <- lapply(designs, function(X) { colnames(X) <- effects; X })
fit <- dkge_fit(betas, designs = designs, K = diag(Q))
W <- matrix(rnorm(ncol(fit$U) * 2L), ncol = 2L)
rownames(W) <- colnames(fit$U)
Y <- fit$U %*% W
rownames(Y) <- effects
res <- dkge_regress(fit, Y, k = 2, engine = "lm")
print(res)
#> dkge_regress results
#>   RMSE (micro): 1.6879
#>   R^2  (micro): -6.3174
#>   R^2  (macro): -6.2059
#>   Cosine mean : -0.2194
# }
```
