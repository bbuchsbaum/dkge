# Predict DKGE loadings for new subjects (out-of-sample)

Predict DKGE loadings for new subjects (out-of-sample)

## Usage

``` r
dkge_predict_loadings(object, B_list, spatial = NULL)
```

## Arguments

- object:

  dkge \| dkge_stream \| dkge_model

- B_list:

  list of qxP_s beta matrices for new subjects

- spatial:

  Optional
  [`dkge_spatial_regularizer()`](https://bbuchsbaum.github.io/dkge/reference/dkge_spatial_regularizer.md)
  for prediction subjects. A fit with a shared spatial domain reuses its
  stored specification automatically. A fit with subject-specific
  domains requires this argument.

## Value

list of P_sxr loadings (A_s) for each subject
