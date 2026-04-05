# Streaming prediction for new subjects via a loader

Streaming prediction for new subjects via a loader

## Usage

``` r
dkge_predict_stream(object, loader, contrasts)
```

## Arguments

- object:

  dkge \| dkge_stream \| dkge_model

- loader:

  object with n(), B(s) methods (and optional X(s))

- contrasts:

  list or matrix as in dkge_predict()

## Value

list(values=list per subject, A_list=list of loadings)
