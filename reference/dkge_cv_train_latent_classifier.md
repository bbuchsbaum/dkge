# Cross-fitted linear classifiers in the DKGE latent space

Trains a binary linear classifier on DKGE latent features using
subject-level folds (LOSO or K-fold). The result stores one weight
vector \\\beta^{(-s)}\\ per held-out subject so downstream
information-map routines can construct bias-aware decoder or Haufe maps.

## Usage

``` r
dkge_cv_train_latent_classifier(
  fit,
  y,
  Z_by_subject = NULL,
  folds = 5,
  model = c("lda", "ridge_logit", "lsvm"),
  ridge = 0.001,
  level = c("subject", "sample"),
  standardize = TRUE
)
```

## Arguments

- fit:

  Fitted \`dkge\` object.

- y:

  Either a factor of length \`S\` (subject-level labels) or a list of
  length \`S\` supplying per-sample labels for each subject.

- Z_by_subject:

  Optional list of latent feature matrices (\`P_s x r\`). When \`NULL\`,
  \[dkge_project_clusters_to_latent()\] is used.

- folds:

  Either an integer \`K\` or a \`dkge_folds\` object created with
  \[dkge_define_folds()\].

- model:

  Binary classifier: \`"lda"\` (pooled covariance linear discriminant,
  default), \`"ridge_logit"\` (ridge-penalised logistic regression via
  glmnet when available), or \`"lsvm"\` (placeholder falling back to
  LDA).

- ridge:

  Ridge penalty added to the pooled covariance matrix when using LDA.
  Ensures numerical stability when the latent dimension is high.

- level:

  Training granularity: \`"subject"\` averages each subject's clusters,
  while \`"sample"\` stacks all cluster samples from training subjects.

- standardize:

  Logical; when \`TRUE\` (default) latent features are standardised
  within each training fold before fitting. Stored weight vectors are
  converted back to the original latent scale.

## Value

An object of class \`dkge_clf\` containing fold models, per-subject
weight vectors, and metadata required by downstream mapping utilities.

## Examples

``` r
# \donttest{
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 6, P = 20, snr = 5
)
fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 2)
#> Warning: Argument 'kernel' is deprecated; use 'K' instead.
y <- factor(rep(c("class1", "class2"), length.out = length(fit$Btil)))
folds <- dkge_define_folds(fit, type = "custom",
                           assignments = list(c(1, 2), c(3, 4), c(5, 6)))
clf <- dkge_cv_train_latent_classifier(fit, y, folds = folds)
length(clf$beta_by_subject)
#> [1] 6
# }
```
