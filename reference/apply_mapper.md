# Apply a fitted mapper to values

Alias to \[predict_mapper()\] that mirrors the terminology used by the
dense rendering helpers.

## Usage

``` r
apply_mapper(fitted_mapper, values, ...)
```

## Arguments

- fitted_mapper:

  Mapping object returned by \[fit_mapper()\].

- values:

  Numeric vector or matrix of source-space values.

- ...:

  Optional backend-specific arguments (e.g. reliabilities).

## Value

Numeric vector or matrix of mapped values.

## Details

For a fitted Sinkhorn dense mapper, the default treats \`values\` as an
intensive field and applies the target-conditional operator, so a
constant source field remains constant. Reliability supplied while
fitting defines the source mass and is not multiplied into values a
second time. Supplying \`reliab\` explicitly here requests an additional
post-fit reweighting. Setting \`normalize_by_reliab = FALSE\` in that
method returns the raw joint plan action rather than a
target-conditional field.

## Examples

``` r
spec <- dkge_mapper("knn", k = 3, sigx = 1)
subj_points <- matrix(rnorm(12 * 3), 12, 3)
anchor_points <- matrix(rnorm(6 * 3), 6, 3)
fitted <- fit_mapper(spec, subj_points = subj_points, anchor_points = anchor_points)
out <- apply_mapper(fitted, rnorm(nrow(subj_points)))
length(out)
#> [1] 6
```
