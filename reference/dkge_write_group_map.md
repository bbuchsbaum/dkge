# Write a group map as NIfTI using a medoid label image

Paints medoid-level values onto the reference parcellation and either
returns a \`neuroim2::BrainVolume\` or writes it to disk.

## Usage

``` r
dkge_write_group_map(
  group_values,
  medoid_labels,
  label_table = NULL,
  out_file = NULL
)
```

## Arguments

- group_values:

  Numeric vector of medoid-cluster values (length Q). When named,
  entries are matched to label IDs before fallback to positional order.

- medoid_labels:

  A \`neuroim2::BrainVolume\` containing integer medoid labels.

- label_table:

  Optional data frame with cluster metadata (currently unused).

- out_file:

  Optional output path (\`.nii\` or \`.nii.gz\`). When \`NULL\`, the
  painted volume is returned without writing to disk.

## Value

Either the output path (when \`out_file\` is supplied) or a
\`neuroim2::BrainVolume\`.

## Examples

``` r
# \donttest{
if (requireNamespace("neuroim2", quietly = TRUE)) {
  labels <- neuroim2::read_vol(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2"))
  vol <- dkge_write_group_map(group_values = 1, medoid_labels = labels)
  class(vol)
}
#> [1] "DenseNeuroVol"
#> attr(,"package")
#> [1] "neuroim2"
# }
```
