# Construct a transport service

Construct a transport service

## Usage

``` r
dkge_transport_service(spec = NULL, ...)
```

## Arguments

- spec:

  Transport specification (list or \`dkge_transport_spec\`).

- ...:

  Additional key-value pairs merged into the specification.

## Value

Object of class \`dkge_transport_service\`.

## Examples

``` r
transport_srv <- dkge_transport_service(dkge_transport_spec(centroids = list(matrix(0, 2, 3))))
```
