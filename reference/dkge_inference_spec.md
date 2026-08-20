# Inference specification helper

Inference specification helper

## Usage

``` r
dkge_inference_spec(
  B = 2000L,
  tail = c("two.sided", "greater", "less"),
  center = "mean"
)
```

## Arguments

- B:

  Number of permutations for sign-flip inference.

- tail:

  Tail of the test: "two.sided", "greater", or "less".

- center:

  Location statistic for permutations. Only \`"mean"\` is implemented by
  the beta inference service.

## Value

Object with class \`dkge_inference_spec\`.

## Examples

``` r
infer <- dkge_inference_spec(B = 1000, tail = "two.sided")
```
