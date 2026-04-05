# Predict DKGE contrasts for new subjects (out-of-sample)

Predict DKGE contrasts for new subjects (out-of-sample)

## Usage

``` r
dkge_predict(object, B_list, contrasts, return_loadings = TRUE)
```

## Arguments

- object:

  dkge \| dkge_stream \| dkge_model

- B_list:

  list of qxP_s betas

- contrasts:

  list of named q-vectors or a qxk matrix (columns are contrasts)

- return_loadings:

  logical; if TRUE also return A_list

## Value

list(A_list=..., values = list of per-contrast subject vectors)
