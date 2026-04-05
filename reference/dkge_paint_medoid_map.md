# Paint medoid cluster values back to a label volume

Paint medoid cluster values back to a label volume

## Usage

``` r
dkge_paint_medoid_map(values, labels, out_file = NULL)
```

## Arguments

- values:

  Vector or matrix of values defined on the medoid parcellation.

- labels:

  Medoid labels describing the target parcellation.

- out_file:

  Optional path to save the rendered map.

## Value

A \`BrainVolume\` or file path, depending on \`out_file\`.
