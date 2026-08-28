# Low-level descriptive Sinkhorn transport to a reference subject

This compatibility primitive fits correspondence from caller-supplied
loose features (`A_list`) and therefore returns an inferentially
ineligible fitted alignment. It is useful for descriptive transport
diagnostics. Inferential workflows should use typed features with
[`dkge_prepare_alignment()`](https://bbuchsbaum.github.io/dkge/reference/dkge_prepare_alignment.md)
or
[`dkge_transport_contrasts_to_reference()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_contrasts_to_reference.md).

`dkge_transport_to_medoid_sinkhorn_cpp()` is a deprecated compatibility
alias. The reference-oriented function already uses the compiled
Sinkhorn backend.

## Usage

``` r
dkge_transport_to_reference_sinkhorn(
  v_list,
  A_list,
  centroids,
  sizes = NULL,
  reference_subject,
  lambda_emb = 1,
  lambda_spa = 0.5,
  sigma_mm = 15,
  epsilon = 0.05,
  max_iter = 5000L,
  tol = 1e-04,
  value_type = c("intensive", "extensive"),
  warm_start = TRUE,
  transport_cache = NULL
)

dkge_transport_to_medoid_sinkhorn_cpp(
  v_list,
  A_list,
  centroids,
  sizes = NULL,
  medoid,
  lambda_emb = 1,
  lambda_spa = 0.5,
  sigma_mm = 15,
  epsilon = 0.05,
  max_iter = 5000L,
  tol = 1e-04,
  value_type = c("intensive", "extensive"),
  warm_start = TRUE,
  return_plans = FALSE,
  transport_cache = NULL
)
```

## Arguments

- v_list:

  List of subject-level cluster values (length P_s each).

- A_list:

  List of subject loadings (P_s x r).

- centroids:

  List of subject cluster centroids (each P_s x 3 matrix).

- sizes:

  Optional list of cluster masses (defaults to uniform weights).

- reference_subject:

  Integer index of the fixed reference subject (1-based). This argument
  does not select or certify a medoid.

- medoid:

  Deprecated name for `reference_subject`.

- lambda_emb, lambda_spa:

  Cost weights for embedding and spatial terms.

- sigma_mm:

  Spatial rescaling (millimetres).

- epsilon, max_iter, tol:

  Sinkhorn parameters.

- value_type:

  Value semantics. `"intensive"` transports field values as
  target-conditional averages and preserves constants; `"extensive"`
  distributes source totals and preserves their sum.

- warm_start:

  Logical; reuse converged dual variables for an identical cost-and-mass
  problem.

- transport_cache:

  Optional fitted alignment returned by
  [`dkge_prepare_transport()`](https://bbuchsbaum.github.io/dkge/reference/dkge_prepare_transport.md).
  Cached operators are reused only after every structural input is
  fingerprint-validated.

- return_plans:

  Logical; if TRUE, include transport plans in the output.

## Value

List containing summary statistics, transported subject maps, and
per-subject joint `plans`, application `operators`, and solver
`diagnostics`.
