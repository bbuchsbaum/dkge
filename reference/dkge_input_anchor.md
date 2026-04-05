# Anchor-based DKGE input descriptor

Builds an input specification that projects subject item kernels onto a
shared anchor basis before fitting DKGE. Use this object with
\[dkge_fit_from_input()\] or supply it to \[dkge_pipeline()\] via the
\`input\` argument. The descriptor is immutable; downstream calls may
extend its \`dkge_args\` field with additional DKGE fitting options.

## Usage

``` r
dkge_input_anchor(
  features_list,
  K_item_list,
  folds = NULL,
  anchors = list(),
  design_kernel = NULL,
  dkge_args = list()
)
```

## Arguments

- features_list:

  List of subject feature matrices (n_s x d).

- K_item_list:

  List of subject item kernels (n_s x n_s PSD).

- folds:

  Optional fold structure passed to \[dkge_build_anchor_kernels()\].

- anchors:

  Optional list overriding anchor selection defaults (see
  \[dkge_anchor_fit()\]).

- design_kernel:

  Optional design kernel supplied to \[dkge_fit_from_kernels()\].
  Defaults to the identity in effect space.

- dkge_args:

  Optional list of arguments forwarded to the DKGE fitter after anchor
  kernels are constructed (e.g. \`w_method\`, \`cpca_part\`).

## Value

Object of class \`dkge_input_anchor\` (inherits from \`dkge_input\`).

## Examples

``` r
set.seed(1)
features_list <- list(
  s1 = matrix(rnorm(20 * 4), 20, 4),
  s2 = matrix(rnorm(25 * 4), 25, 4),
  s3 = matrix(rnorm(22 * 4), 22, 4)
)
K_item_list <- lapply(features_list, function(X) tcrossprod(matrix(rnorm(nrow(X) * 3), nrow(X), 3)))
input <- dkge_input_anchor(features_list, K_item_list, anchors = list(L = 6, seed = 1))
class(input)
#> [1] "dkge_input_anchor" "dkge_input"       
```
