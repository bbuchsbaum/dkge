# Transport cluster values to a medoid via entropic Sinkhorn OT

Transport cluster values to a medoid via entropic Sinkhorn OT

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
  max_iter = 200,
  tol = 1e-06,
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
  max_iter = 300,
  tol = 1e-06,
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

- transport_cache:

  Optional cache returned by \[dkge_prepare_transport()\]. When
  supplied, cached operators are reused and the mapper configuration
  stored in the cache takes precedence.

- return_plans:

  Logical; if TRUE, include transport plans in the output.

## Value

List containing summary statistics, transported subject maps, and the
per-subject transport operators.
