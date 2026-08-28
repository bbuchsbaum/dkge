# Run a dkge service object

Generic dispatcher for executing pre-configured service objects.
Dispatches to the appropriate typed runner based on the service class.

## Usage

``` r
dkge_run_service(service, ...)

# S3 method for class 'dkge_contrast_service'
dkge_run_service(service, ...)

# S3 method for class 'dkge_inference_service'
dkge_run_service(service, ...)

# S3 method for class 'dkge_transport_service'
dkge_run_service(service, ...)

# Default S3 method
dkge_run_service(service, ...)
```

## Arguments

- service:

  A service object created by
  [`dkge_contrast_service()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast_service.md),
  [`dkge_inference_service()`](https://bbuchsbaum.github.io/dkge/reference/dkge_inference_service.md),
  or
  [`dkge_transport_service()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_service.md)

- ...:

  Additional arguments passed to the typed runner

## Value

Result from the typed service runner
