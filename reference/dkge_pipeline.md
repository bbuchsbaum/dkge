# End-to-end DKGE workflow

Fits DKGE (if needed), computes cross-fitted contrasts, optionally
transports them to a medoid parcellation, and performs sign-flip
inference.

## Usage

``` r
dkge_pipeline(
  fit = NULL,
  input = NULL,
  betas = NULL,
  designs = NULL,
  kernel = NULL,
  omega = NULL,
  contrasts,
  transport = NULL,
  inference = list(),
  classification = NULL,
  method = c("loso", "kfold", "analytic"),
  ridge = 0,
  ...
)
```

## Arguments

- fit:

  Optional pre-computed \`dkge\` object. If \`NULL\`, provide \`betas\`,
  \`designs\`, and \`kernel\` to fit inside the pipeline.

- input:

  Optional DKGE input descriptor created with \[dkge_input_anchor()\] or
  future helpers. When supplied (and \`fit\` is \`NULL\`),
  \`dkge_pipeline()\` will build the fit via \[dkge_fit_from_input()\].

- betas, designs, kernel:

  Inputs passed to \[dkge()\] when neither \`fit\` nor \`input\` is
  supplied.

- omega:

  Optional spatial weights forwarded to \[dkge()\].

- contrasts:

  Contrast specification as accepted by \[dkge_contrast()\].

- transport:

  Either a transport specification/service or \`NULL\`.

- inference:

  Either an inference specification/service or \`NULL\`.

- classification:

  Optional specification passed to \[dkge_classify()\].

- method:

  Cross-fitting strategy for contrasts (default "loso").

- ridge:

  Optional ridge added during held-out decompositions.

- ...:

  Additional arguments passed to \[dkge()\] when fitting inside the
  pipeline, or to \[dkge_contrast()\].

## Value

List containing the fit, diagnostics, raw contrast values, transported
maps (if requested), and inference results.

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
