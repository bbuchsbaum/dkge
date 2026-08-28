# Build prior weights from ROI labels

Build prior weights from ROI labels

## Usage

``` r
dkge_weights_prior_roi(labels, roi_values = NULL)
```

## Arguments

- labels:

  Integer/factor ROI labels per voxel.

- roi_values:

  Optional named numeric multipliers per ROI.

## Value

Numeric vector of prior weights.
