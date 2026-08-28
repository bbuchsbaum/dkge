# Inference specification helper

Inference specification helper

## Usage

``` r
dkge_inference_spec(
  B = 2000L,
  tail = c("two.sided", "greater", "less"),
  center = c("mean", "median", "none"),
  allow_approximate_alignment = FALSE
)
```

## Arguments

- B:

  Number of permutations for sign-flip inference.

- tail:

  Tail of the test: "two.sided", "greater", or "less".

- center:

  Location statistic. Only `"mean"` is supported. Legacy `"median"` and
  `"none"` values now fail at construction because the downstream max-T
  statistic is a studentized mean.

- allow_approximate_alignment:

  Logical; explicitly permit inference from an estimator or fitted
  alignment labelled `"approximate"`. The default is fail-closed;
  ineligible states are never permitted.

## Value

Object with class `dkge_inference_spec`.

## Examples

``` r
infer <- dkge_inference_spec(B = 1000, tail = "two.sided")
```
