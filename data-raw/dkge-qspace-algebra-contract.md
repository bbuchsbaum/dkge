# DKGE q-space moment and fit-object contract

This note is the executable algebra specification for Mote ticket
`bd-01M03KTTKQHAE2BYW6C7RBQY2E`. It separates the common q-space eigensystem
from the stronger multiblock singular-value-decomposition contract.

## Common q-space estimator

Let `M` be the pooled effect-space moment after subject weighting, optional
effect-pair normalization, missingness handling, and optional noise correction.
With the pooled design ruler `R` and symmetric kernel root `Khalf`, DKGE solves

```
Chat = Khalf %*% t(R) %*% M %*% R %*% Khalf.
```

`Chat` is always symmetrized before solving. For `solver = "pooled"`, DKGE uses
the ordinary double-precision symmetric eigensolver. Write its descending full
spectrum as `Chat = Q diag(lambda) t(Q)`. Only eigenpairs satisfying the
documented scale-relative positive-eigenvalue tolerance are retained. Negative
eigenvalues remain in `eig_values_full` and the moment diagnostics; they are not
silently converted into positive signal components.

For retained `Q_r` and positive `lambda_r`:

```
U       = Kihalf %*% Q_r
sdev    = sqrt(lambda_r)
s       = Q_r %*% diag(sdev)
Chat_r+ = s %*% t(s)
```

The kernel inverse root is a Moore-Penrose inverse on the numerical support of
`K`. Consequently, every advertised component must satisfy
`t(U) %*% K %*% U = I`. The common reconstruction claim is only that `s %*%
t(s)` is the retained positive spectral part of `Chat`; it is not necessarily
the full moment when the requested rank is truncated or `Chat` is indefinite.

The public fields have these meanings:

- `Chat`: the exact symmetrized q by q matrix presented to the eigensolver,
  including any requested ridge or CPCA filtering;
- `eig_values_full`, `eig_vectors_full`, and `evals`: the full descending
  symmetric eigendecomposition of `Chat` for the pooled solver;
- `U`: the retained K-orthonormal effect basis;
- `sdev`: square roots of retained positive eigenvalues;
- `s` and `scores_matrix`: the q by r retained positive q-space factor `s`
  above; these are not subject scores;
- `rank`: the number of retained positive components, which can be smaller than
  `rank_requested`;
- `moment_diagnostics`: the negative count and mass before positive-subspace
  truncation.

## Representation classes

### Exact block biprojector

An exact multiblock representation is advertised only when all of the following
hold:

1. the pooled symmetric eigensolver is used;
2. no effect-pair normalization or missingness transformation is active;
3. no analytic or split-half debiasing is active;
4. no CPCA filtering is active; and
5. no ridge has been added to `Chat`.

Subject weights, voxel weights, and a positive-semidefinite spatial `Omega_s`
are factorizable right/whole-block weights and are allowed. Structurally missing
effect rows are zeroed before `R` and `Khalf` mix effect coordinates. The exact
training matrix is

```
Xstar_s = sqrt(w_s) * Khalf %*% t(R) %*% Bzero_s %*%
          sqrt(Omega_s and voxel weights)
Xstar   = cbind(Xstar_1, ..., Xstar_S).
```

For this representation, DKGE guarantees

```
Chat             = Xstar %*% t(Xstar)
v                = t(Xstar) %*% Q_r %*% diag(1 / sdev)
t(v) %*% v       = I
s %*% t(v)       = the rank-r SVD reconstruction of Xstar.
```

Only this representation inherits from `multiblock_biprojector`, advertises
`v`, and supports `dkge_project_blocks()` / `dkge_project_block()`.

### Q-space moment representation

Effect-pair reliability normalization, `rescale`, `mask`, `shrink`, analytic
debiasing, split-half cross-moments, ridge, CPCA, and joint diagonalization do
not in general share the physical block factor above. Pairwise division and
entrywise masking can destroy positive semidefiniteness; analytic subtraction
and symmetrized split-half products can be indefinite. Even when the realized
matrix happens to be PSD and therefore has a canonical spectral factor, that
factor is not a subject-by-voxel block loading matrix.

Such fits use `representation = "qspace_moment"`, do not inherit from
`multiblock_biprojector`, and set `v` and `X_concat` to `NULL`. `U`, `sdev`, `s`,
the full spectrum, q-space prediction, LOSO contrasts, and moment diagnostics
remain meaningful under the definitions above. Requests for a physical block
projection or `keep_X = TRUE` fail with an actionable message.

`solver = "jd"` is also q-space-only. Its solver-specific diagonal values are
not claimed to reconstruct `Chat`; weighted low-rank/ALS solving is Tier 3 and
is explicitly outside Tiers 0 to 2.

## Subject and downstream projections

`dkge_predict_loadings()` and contrast routines project subject effect matrices
onto `K %*% U`; they remain q-space operations and are valid for both
representations. Standardized cluster coordinates that divide by `sdev` are
defined as coordinates relative to the retained positive q-space spectrum, not
as rows of a physical `v`, unless `representation = "block_biprojector"`.

Any downstream function that requires physical block loadings must call the
representation guard and fail closed. It must not synthesize `v` from raw
blocks when `Chat` was produced by a different moment estimator.

## Executable enforcement map

`tests/testthat/test-fit-algebra-contract.R` enforces:

| Failure mode | Test family and invariant |
|---|---|
| ordinary assembly drifts from the fit matrix | differential oracle: `Chat == tcrossprod(X_concat)` |
| missing rows are zeroed after coordinate mixing | regression oracle with non-diagonal `R` and `K` |
| a non-factorizable moment advertises physical loadings | contract test: q-space class, `v = NULL`, no multiblock inheritance |
| invalid right singular vectors look plausible | contract test: `crossprod(v) == I` whenever `v` exists |
| q-space reconstruction is ambiguous | spectral oracle: `tcrossprod(s)` equals the retained positive eigensystem |
| debiasing creates negative modes | adversarial analytic-debias test retaining negative diagnostics but only positive components |
| downstream code consumes absent block geometry | fail-closed tests for block projection and `keep_X = TRUE` |
| kernel support is singular | adversarial PSD-kernel test: `crossprod(U, K %*% U) == I` |
| effect ordering changes the estimator | metamorphic simultaneous row/column permutation test |

Later tickets add randomized oracle matrices, resampling equivalence, scale and
conditioning stress tests, and performance budgets. Those later tests may
tighten tolerances but may not weaken the representation boundary defined here.
