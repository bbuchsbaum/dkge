# Construct a legacy descriptive transport service

Pipeline transport services are retained for descriptive migration
output. They cannot be composed with inference; use the typed
reference-oriented alignment workflow for inferential maps.

## Usage

``` r
dkge_transport_service(spec = NULL, ...)
```

## Arguments

- spec:

  Transport specification (list or `dkge_transport_spec`).

- ...:

  Additional key-value pairs merged into the specification.

## Value

Object of class `dkge_transport_service`.

## Examples

``` r
transport_srv <- dkge_transport_service(dkge_transport_spec(centroids = list(matrix(0, 2, 3))))
```
