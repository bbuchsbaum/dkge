# Unified inference for DKGE contrasts

High-level interface for statistical inference on DKGE contrasts with
integrated cross-fitting and multiple testing correction.

## Usage

``` r
dkge_infer(
  fit,
  contrasts,
  method = c("loso", "kfold", "analytic"),
  inference = c("signflip", "freedman-lane", "parametric"),
  correction = c("maxT", "fdr", "bonferroni", "none"),
  n_perm = 2000,
  alpha = 0.05,
  transported = FALSE,
  transport = NULL,
  allow_approximate_alignment = FALSE,
  ...
)
```

## Arguments

- fit:

  A `dkge` object from
  [`dkge_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit.md)
  or [`dkge()`](https://bbuchsbaum.github.io/dkge/reference/dkge.md)

- contrasts:

  Contrast specification (see
  [`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md))

- method:

  Cross-fitting method: "loso", "kfold", or "analytic"

- inference:

  Inference type:

  - `"signflip"`: Sign-flip permutation test (default)

  - `"freedman-lane"`: Freedman-Lane permutation (requires adapters)

  - `"parametric"`: Parametric t-test (assumes normality)

- correction:

  Multiple testing correction:

  - `"maxT"`: Family-wise error rate via max-T (default)

  - `"fdr"`: False discovery rate (Benjamini-Hochberg)

  - `"bonferroni"`: Bonferroni correction

  - `"none"`: No correction

- n_perm:

  Number of permutations for non-parametric tests

- alpha:

  Significance level for corrections

- transported:

  Deprecated. Must be `FALSE`. Rendering/alignment is no longer learned
  inside this inference helper.

- transport:

  Deprecated. Must be `NULL`. Use the typed two-stage workflow
  [`dkge_transport_contrasts_to_reference()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_contrasts_to_reference.md)
  then
  [`dkge_infer_aligned()`](https://bbuchsbaum.github.io/dkge/reference/dkge_infer_aligned.md).

- allow_approximate_alignment:

  Logical; permit an explicitly labelled `"approximate"` alignment or
  same-data rank-truncated LOSO, K-fold, or analytic estimator.
  Ineligible/descriptive states are always refused. The default is
  fail-closed.

- ...:

  Additional arguments passed to
  [`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md)
  and inference functions

## Value

An object of class `dkge_inference` containing:

- `contrasts`: The contrast results from cross-fitting

- `statistics`: Test statistics per cluster/voxel

- `p_values`: Raw p-values

- `p_adjusted`: Adjusted p-values based on correction method

- `significant`: Logical indicators of significance

- `method`: Cross-fitting method used

- `inference`: Inference type used

- `correction`: Correction method applied

- `metadata`: Additional information about the analysis

## Details

This function integrates the cross-fitting machinery from
[`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md)
with various statistical inference procedures. It first computes
contrast values using the specified cross-fitting method, then applies
the chosen inference procedure to obtain p-values, and finally applies
multiple testing correction.

This helper performs inference only in the contrast result's native
support. It never learns correspondence. When subject supports differ,
first build typed alignment features, transport with
[`dkge_transport_contrasts_to_reference()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_contrasts_to_reference.md),
extract its `dkge_aligned_maps` object, and pass that object to
[`dkge_infer_aligned()`](https://bbuchsbaum.github.io/dkge/reference/dkge_infer_aligned.md).

The workflow is:

1.  Compute contrast values via cross-fitting (LOSO/K-fold/analytic)

2.  Apply inference procedure (sign-flip/Freedman-Lane/parametric)

3.  Apply multiple testing correction (maxT/FDR/Bonferroni)

Max-T uses one subject-sign action and one maximum over the entire
prespecified contrast-by-location family. Its conditional FWER guarantee
requires joint row-sign invariance of the supplied subject matrix.
Standard rank-truncated LOSO/K-fold DKGE estimates retain same-data
latent-span dependence and are therefore labelled approximate unless the
full estimator is rebuilt under every null action. The analytic method
is also approximate: it uses a first-order approximation to the LOSO
basis except where its diagnostic selects an exact-LOSO fallback.

## See also

[`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md),
[`dkge_infer_aligned()`](https://bbuchsbaum.github.io/dkge/reference/dkge_infer_aligned.md),
[`dkge_transport_contrasts_to_reference()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_contrasts_to_reference.md),
[`dkge_signflip_maxT()`](https://bbuchsbaum.github.io/dkge/reference/dkge_signflip_maxT.md),
[`dkge_freedman_lane()`](https://bbuchsbaum.github.io/dkge/reference/dkge_freedman_lane.md)

## Examples

``` r
# Simulate and fit
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 6, P = 15, snr = 5
)
fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 2)
#> Warning: Argument 'kernel' is deprecated; use 'K' instead.

# LOSO with sign-flip and maxT correction (fast with few perms for example)
# \donttest{
results <- dkge_infer(
  fit, c(1, rep(0, 4)), n_perm = 100,
  allow_approximate_alignment = TRUE
)
results
#> DKGE Inference Results
#> ----------------------
#> Cross-fitting: loso
#> Inference: signflip
#> Correction: maxT
#> Contrasts: 1
#> Alpha level: 0.05
#>   contrast1: 0/15 significant
# }
```
