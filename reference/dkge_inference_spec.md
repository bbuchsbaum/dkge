# Inference specification helper

Inference specification helper

## Usage

``` r
dkge_inference_spec(
  B = 2000L,
  tail = c("two.sided", "greater", "less"),
  center = c("mean", "median", "none")
)
```

## Arguments

- B:

  Number of permutations for sign-flip inference.

- tail:

  Tail of the test: "two.sided", "greater", or "less".

- center:

  Centering method for permutations: "mean", "median", or "none".

## Value

Object with class \`dkge_inference_spec\`.

## Examples

``` r
infer <- dkge_inference_spec(B = 1000, tail = "two.sided")
```
