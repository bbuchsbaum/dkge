# Transport cluster values to a medoid via entropic Sinkhorn OT

\`dkge_transport_to_medoid_sinkhorn_cpp()\` is a deprecated
compatibility alias. The main function already uses the compiled
Sinkhorn backend.

## Usage

``` r
dkge_transport_to_medoid_sinkhorn(
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

- medoid:

  Integer index of the reference subject (1-based).

- lambda_emb, lambda_spa:

  Cost weights for embedding and spatial terms.

- sigma_mm:

  Spatial rescaling (millimetres).

- epsilon, max_iter, tol:

  Sinkhorn parameters.

- value_type:

  Value semantics. \`"intensive"\` transports field values as
  target-conditional averages and preserves constants on positive-mass
  targets; \`"extensive"\` distributes source totals and preserves their
  sum over positive-mass sources. Null-mass target columns (intensive)
  and source rows (extensive) are represented by zeros in the
  application operator.

- warm_start:

  Logical; reuse converged dual variables for an identical cost-and-mass
  problem.

- transport_cache:

  Optional cache returned by \[dkge_prepare_transport()\]. When
  supplied, cached operators are reused and the mapper configuration
  stored in the cache takes precedence.

- return_plans:

  Logical; if TRUE, include transport plans in the output.

## Value

List containing summary statistics, transported subject maps, and
per-subject joint \`plans\`, application \`operators\`, and solver
\`diagnostics\`.
