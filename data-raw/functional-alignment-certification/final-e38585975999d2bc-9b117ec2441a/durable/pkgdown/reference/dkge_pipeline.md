# End-to-end DKGE workflow

Fits DKGE (if needed), computes cross-fitted contrasts, optionally
produces legacy descriptive transport output, and performs
native-support sign-flip inference. Pipeline transport and inference
cannot be composed. Functional alignment inference uses
[`dkge_transport_contrasts_to_reference()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_contrasts_to_reference.md)
followed by
[`dkge_infer_aligned()`](https://bbuchsbaum.github.io/dkge/reference/dkge_infer_aligned.md).

## Usage

``` r
dkge_pipeline(
  fit = NULL,
  input = NULL,
  betas = NULL,
  designs = NULL,
  kernel = NULL,
  omega = NULL,
  spatial = NULL,
  contrasts,
  transport = NULL,
  inference = NULL,
  classification = NULL,
  method = c("loso", "kfold", "analytic"),
  ridge = 0,
  ...
)
```

## Arguments

- fit:

  Optional pre-computed `dkge` object. If `NULL`, provide `betas`,
  `designs`, and `kernel` to fit inside the pipeline.

- input:

  Optional DKGE input descriptor created with
  [`dkge_input_anchor()`](https://bbuchsbaum.github.io/dkge/reference/dkge_input_anchor.md)
  or future helpers. When supplied (and `fit` is `NULL`),
  `dkge_pipeline()` will build the fit via
  [`dkge_fit_from_input()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit_from_input.md).

- betas, designs, kernel:

  Inputs passed to
  [`dkge()`](https://bbuchsbaum.github.io/dkge/reference/dkge.md) when
  neither `fit` nor `input` is supplied.

- omega:

  Optional spatial weights forwarded to
  [`dkge()`](https://bbuchsbaum.github.io/dkge/reference/dkge.md).

- spatial:

  Optional model-level
  [`dkge_spatial_regularizer()`](https://bbuchsbaum.github.io/dkge/reference/dkge_spatial_regularizer.md)
  forwarded only to the raw-beta fitting stage. Current anchor input
  descriptors do not expose a physical beta-column domain and therefore
  reject this argument.

- contrasts:

  Contrast specification as accepted by
  [`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md).

- transport:

  Either a legacy descriptive transport specification/service or `NULL`.
  It cannot be combined with `inference`.

- inference:

  Either an inference specification/service or `NULL` (the default).
  Same-data rank-truncated inference is approximate and requires an
  explicit `allow_approximate_alignment = TRUE` in the inference spec.

- classification:

  Optional specification passed to
  [`dkge_classify()`](https://bbuchsbaum.github.io/dkge/reference/dkge_classify.md).

- method:

  Cross-fitting strategy for contrasts (default "loso").

- ridge:

  Optional ridge added during held-out decompositions.

- ...:

  Additional arguments passed to
  [`dkge()`](https://bbuchsbaum.github.io/dkge/reference/dkge.md) when
  fitting inside the pipeline, or to
  [`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md).

## Value

List containing the fit, diagnostics, raw contrast values, optional
legacy descriptive maps, and optional native-support inference results.

## Examples

``` r
# Simulate toy data
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 5, P = 25, snr = 5
)

# Run pipeline with LOSO contrasts
result <- dkge_pipeline(
  betas = toy$B_list,
  designs = toy$X_list,
  kernel = toy$K,
  contrasts = c(1, rep(0, 4)),  # first effect
  method = "loso"
)
names(result)
#> [1] "fit"            "diagnostics"    "contrasts"      "transport"     
#> [5] "inference"      "classification"
```
