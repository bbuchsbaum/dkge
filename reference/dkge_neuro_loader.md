# Build a streaming loader backed by neuroim2 objects

The returned loader exposes \`n()\`, \`X(s)\`, \`B(s)\`, and
\`Omega(s)\` methods compatible with streaming DKGE fits.

## Usage

``` r
dkge_neuro_loader(design_objs, bv_list, labels_list = NULL, omega_fun = NULL)
```

## Arguments

- design_objs:

  List of design matrices or \`fmridesign\` objects.

- bv_list:

  List of \`neuroim2::NeuroVec\` objects, \`ClusteredNeuroVec\` objects,
  or file paths readable by \`neuroim2::read_vec()\`.

- labels_list:

  List of \`neuroim2::NeuroVol\` label volumes (same length as
  \`bv_list\`). Can be NULL if bv_list contains ClusteredNeuroVec
  objects.

- omega_fun:

  Optional function mapping a label volume or ClusteredNeuroVec to
  cluster weights. Default uses cluster sizes as weights.

## Value

Loader list suitable for streaming DKGE fits.

## Examples

``` r
# \donttest{
if (requireNamespace("neuroim2", quietly = TRUE)) {
  labels <- neuroim2::read_vol(system.file("extdata", "global_mask2.nii.gz", package = "neuroim2"))
  vols <- lapply(1:6, function(i) labels * i)
  bv <- neuroim2::vec_from_vols(vols)
  x_mat <- cbind(intercept = 1, trend = seq_len(length(vols)))
  loader <- dkge_neuro_loader(
    design_objs = list(x_mat, x_mat),
    bv_list = list(bv, bv),
    labels_list = list(labels, labels)
  )
  loader$n()
}
#> [1] 2
# }
```
