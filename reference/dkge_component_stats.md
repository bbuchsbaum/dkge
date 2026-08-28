# Deprecated descriptive component consensus

This legacy helper transports full-fit component loadings with
correspondence learned from those same loadings. It is retained only for
descriptive summaries and cannot compute p-values, confidence claims, or
significance. Use
[`dkge_prepare_alignment()`](https://bbuchsbaum.github.io/dkge/reference/dkge_prepare_alignment.md),
apply correspondence with
[`dkge_transport_contrasts_to_reference()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_contrasts_to_reference.md)
or
[`dkge_align_to_template()`](https://bbuchsbaum.github.io/dkge/reference/dkge_align_to_template.md),
and then call
[`dkge_infer_aligned()`](https://bbuchsbaum.github.io/dkge/reference/dkge_infer_aligned.md)
for the typed inferential workflow.

## Usage

``` r
dkge_component_stats(
  fit,
  mapper = "sinkhorn",
  centroids = NULL,
  sizes = NULL,
  inference = NULL,
  medoid = 1L,
  components = NULL,
  adjust = "fdr",
  ...
)

dkge_write_component_stats(fit, file, ...)
```

## Arguments

- fit:

  A fitted `dkge` object.

- mapper:

  Mapper strategy (string or
  [`dkge_mapper_spec()`](https://bbuchsbaum.github.io/dkge/reference/dkge_mapper_spec.md)).
  Defaults to "sinkhorn".

- centroids:

  Optional list of subject centroid matrices; defaults to centroids
  stored in `fit` if available.

- sizes:

  Optional list of cluster masses (one vector per subject).

- inference:

  Deprecated. Must be `NULL`; legacy full-fit correspondence is
  descriptive/ineligible and cannot enter an inference boundary.

- medoid:

  Reference subject index for the descriptive display.

- components:

  Optional vector of component indices or names; default is all
  components.

- adjust:

  Deprecated and ignored because inferential p-values are no longer
  produced by this helper.

- ...:

  Additional mapper-specific parameters (e.g. `epsilon`).

- file:

  Path to the CSV file where component statistics will be written.

## Value

A list with fields:

- `summary`: tidy data frame of descriptive means and standard
  deviations.

- `statistics`: per-component mean vectors.

- `transport`: per-component transported subject matrices.

- `eligibility`: an ineligible/descriptive alignment receipt.

## Examples

``` r
# \donttest{
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 3, P = 15, snr = 5
)
fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 2)
#> Warning: Argument 'kernel' is deprecated; use 'K' instead.
centroids <- lapply(toy$B_list, function(B) matrix(rnorm(ncol(B) * 3), ncol(B), 3))
res <- suppressWarnings(dkge_component_stats(fit,
                            centroids = centroids,
                            mapper = "ridge",
                            inference = NULL,
                            components = 1))
head(res$summary)
#>   component cluster        mean        sd n_subjects
#> 1         1       1  1.22679721 0.1664723          3
#> 2         1       2  0.09321496 0.1170358          3
#> 3         1       3  1.86502739 0.2972021          3
#> 4         1       4  1.86158567 0.3286176          3
#> 5         1       5 -0.18470602 0.0899974          3
#> 6         1       6  1.80812321 0.3264447          3
# }
```
