# Construct a contrast service

Packages the arguments required by
[`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md)
into a reusable service object. The service can be passed to
[`dkge_pipeline()`](https://bbuchsbaum.github.io/dkge/reference/dkge_pipeline.md)
(via the `service` argument when it gains support) or executed manually
with `.dkge_run_contrast_service()`.

## Usage

``` r
dkge_contrast_service(method = c("loso", "kfold", "analytic"), ridge = 0, ...)
```

## Arguments

- method:

  Cross-fitting strategy ("loso", "kfold", or "analytic").

- ridge:

  Ridge penalty forwarded to
  [`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md).

- ...:

  Additional arguments stored with the service (e.g. `folds`,
  `parallel`).

## Value

Object of class `dkge_contrast_service`.

## Examples

``` r
contrast_srv <- dkge_contrast_service(method = "loso", ridge = 0)
```
