# Build prior weights from a mask

Build prior weights from a mask

## Usage

``` r
dkge_weights_prior_mask(mask, value_in = 1, value_out = 0)
```

## Arguments

- mask:

  Logical or numeric vector identifying voxels in the ROI.

- value_in, value_out:

  Weights inside/outside the ROI (defaults 1/0).

## Value

Numeric vector of prior weights.
