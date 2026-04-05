# Construct an inference service

Construct an inference service

## Usage

``` r
dkge_inference_service(spec = NULL, ...)
```

## Arguments

- spec:

  Inference specification (list or \`dkge_inference_spec\`).

- ...:

  Additional parameters merged into the specification.

## Value

Object of class \`dkge_inference_service\`.

## Examples

``` r
inference_srv <- dkge_inference_service(dkge_inference_spec(B = 1000))
```
