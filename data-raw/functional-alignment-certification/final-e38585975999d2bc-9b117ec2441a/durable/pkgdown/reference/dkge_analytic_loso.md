# Analytic LOSO contrast using eigenvalue perturbation

Approximates leave-one-subject-out contrast values using first-order
eigenvalue perturbation theory, avoiding full eigen-decomposition per
subject.

## Usage

``` r
dkge_analytic_loso(fit, s, contrasts, tol = 1e-06, fallback = TRUE, ridge = 0)
```

## Arguments

- fit:

  A `dkge` object from
  [`dkge_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit.md)

- s:

  Subject index (1-based) to leave out

- contrasts:

  Contrast vector in the original design basis (length q)

- tol:

  Numerical tolerance for determining when to fall back to full eigen

- fallback:

  If TRUE, fall back to full eigen when perturbation may be unstable

- ridge:

  Ridge regularization parameter (default 0)

## Value

List with fields:

- `v`: Cluster contrast values for subject s

- `alpha`: Contrast coordinates in latent space

- `basis`: Approximated held-out basis U^(-s)

- `method`: "analytic" or "fallback" if full eigen was used

## Details

This function implements the first-order eigenvalue perturbation
approximation described in the paper. For the held-out compressed
covariance:

Delta_s = Chat^(-s) - Chat

where `Chat^(-s)` is formed with the same fold-local subject-weight
normalization and shrinkage as exact LOSO. The eigenvalues and
eigenvectors are updated using:

- deltalambda_j = v_j^T Delta_s v_j (eigenvalue shift)

- deltav_j = Sum over k!=j of (v_k^T Delta_s v_j)/(lambda_j - lambda_k)
  v_k (eigenvector rotation)

After the exact fold moment has been constructed, this replaces its
O(q^3) eigen-decomposition with an O(q^2 r) first-order eigensystem
update, where r is the rank. The approximation is accurate when:

1.  The fold perturbation norm is small relative to the fitted
    eigensystem

2.  Eigenvalue gaps are large (well-separated components)

3.  The resulting first-order rotation coefficients remain below the
    gate

When these conditions are violated (detected via condition number or
eigenvalue gaps), the function can fall back to full
eigen-decomposition.

Fallback diagnostics use a closed reason vocabulary with this
precedence: `pair_normalized_pooling`, `covariance_aware_moment`,
`nonuniform_voxel_weights`, `missing_full_decomposition`,
`dimension_mismatch`, `eigengap`, and `perturbation_magnitude`.
Structural reasons are checked before numerical perturbation thresholds,
so a large perturbation cannot mask the more basic fact that the fitted
moment does not support this approximation. The q-by-q perturbation
itself is exact for the fold pooling contract; only its eigensystem
update is first-order. A successful approximation reports `analytic`.

## References

Golub, G. H., & Van Loan, C. F. (2013). Matrix computations (4th ed.).
Stewart, G. W., & Sun, J. (1990). Matrix perturbation theory.

## Examples

``` r
# \donttest{
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 4, P = 20, snr = 5
)
fit <- dkge(toy$B_list, toy$X_list, kernel = toy$K, rank = 2)
#> Warning: Argument 'kernel' is deprecated; use 'K' instead.
c_vec <- c(1, -1, 0, 0, 0)
result <- dkge_analytic_loso(fit, s = 1, contrasts = c_vec)
result$method
#> [1] "fallback"
# }
```
