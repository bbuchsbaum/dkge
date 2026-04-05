# Predict DKGE loadings for new subjects (out-of-sample)

Predict DKGE loadings for new subjects (out-of-sample)

## Usage

``` r
dkge_predict_loadings(object, B_list)
```

## Arguments

- object:

  dkge \| dkge_stream \| dkge_model

- B_list:

  list of qxP_s beta matrices for new subjects

## Value

list of P_sxr loadings (A_s) for each subject
