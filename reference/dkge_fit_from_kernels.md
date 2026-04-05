# Fit DKGE from precomputed subject effect kernels

Converts a list of subject-level effect kernels \\K_s \in \mathbb{R}^{q
\times q}\\ into synthetic GLM inputs that reuse \[dkge_fit()\] without
modifying the core implementation. Each kernel is factorised into a
symmetric square root, scaled to keep the pooled design metric
unchanged, and paired with an identity design matrix so the resulting
DKGE fit matches the supplied kernels.

## Usage

``` r
dkge_fit_from_kernels(
  K_list,
  effect_ids,
  subject_ids = NULL,
  design_kernel = NULL,
  sqrt_tol = 1e-10,
  ...
)
```

## Arguments

- K_list:

  List of symmetric positive semi-definite matrices sharing the same
  effect ordering.

- effect_ids:

  Character vector of length \\q\\ naming the shared effect (anchor)
  indices.

- subject_ids:

  Optional character vector naming subjects. Defaults to the names of
  \`K_list\` or sequential identifiers.

- design_kernel:

  Optional design kernel passed to \[dkge_fit()\]. Defaults to the \\q
  \times q\\ identity matrix, which matches the whitened anchor setup.

- sqrt_tol:

  Eigenvalue tolerance used when extracting square roots.

- ...:

  Additional arguments forwarded to \[dkge_fit()\].

## Value

A \`dkge\` object identical to one obtained from \[dkge_fit()\], with
provenance annotated to record the kernel-driven construction.

## Examples

``` r
# \donttest{
q <- 5
Ks <- replicate(3, {
  X <- matrix(rnorm(q * q), q)
  S <- crossprod(X)
  S / sqrt(sum(diag(S)))
}, simplify = FALSE)
fit <- dkge_fit_from_kernels(Ks, effect_ids = paste0("z", seq_len(q)))
# }
```
