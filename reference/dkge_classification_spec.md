# Classification specification helper

Classification specification helper

## Usage

``` r
dkge_classification_spec(
  targets,
  method = c("lda", "logit"),
  folds = NULL,
  lambda = NULL,
  metric = c("accuracy", "logloss"),
  mode = c("auto", "cell", "cell_cross", "delta"),
  ...
)
```

## Arguments

- targets:

  Target specification accepted by \[dkge_classify()\].

- method:

  Classifier backend ("lda" or "logit").

- folds:

  Optional fold specification.

- lambda:

  Optional ridge parameter.

- metric:

  Classification metrics to report.

- mode:

  Decoding mode passed to \[dkge_classify()\]: "auto" selects
  automatically, "cell" uses per-cell embeddings, "cell_cross" uses
  cross-validated per-cell embeddings, and "delta" uses pairwise deltas.

- ...:

  Additional fields stored on the spec (e.g., \`n_perm\`, \`scope\`).

## Value

Object with class \`dkge_classification_spec\`.

## Examples

``` r
cls <- dkge_classification_spec(targets = ~ condition, method = "lda")
```
