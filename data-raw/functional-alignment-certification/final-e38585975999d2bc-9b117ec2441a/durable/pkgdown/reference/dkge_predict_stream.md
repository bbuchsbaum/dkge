# Streaming prediction for new subjects via a loader

Streaming prediction for new subjects via a loader

## Usage

``` r
dkge_predict_stream(object, loader, contrasts, spatial = NULL)
```

## Arguments

- object:

  dkge \| dkge_stream \| dkge_model

- loader:

  object with n(), B(s) methods (and optional X(s))

- contrasts:

  list or matrix as in dkge_predict()

- spatial:

  Optional
  [`dkge_spatial_regularizer()`](https://bbuchsbaum.github.io/dkge/reference/dkge_spatial_regularizer.md)
  for prediction subjects. A fit with a shared spatial domain reuses its
  stored specification automatically. A fit with subject-specific
  domains requires this argument.

## Value

list(values=list per subject, A_list=list of loadings)
