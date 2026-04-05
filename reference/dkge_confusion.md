# Fold-wise confusion matrices for DKGE classification

Fold-wise confusion matrices for DKGE classification

## Usage

``` r
dkge_confusion(x, target = NULL, fold = NULL)
```

## Arguments

- x:

  A `dkge_classification` object.

- target:

  Optional target name or index. When `NULL`, all targets are returned.

- fold:

  Optional fold identifier (index or numeric label). When `NULL`,
  confusion matrices are summed across folds.

## Value

A confusion matrix (when a single target is requested) or a named list
of matrices.
