# DKGE Concepts: Kernels, Inference, and What the Method Actually Does

Read this before your first fit. It explains what DKGE is actually
decomposing, what choosing a kernel does to that decomposition, and at
which level — component, contrast, or feature — you can claim
statistical evidence. Workflow mechanics live in
[`vignette("dkge-workflow")`](https://bbuchsbaum.github.io/dkge/articles/dkge-workflow.md);
kernel construction details in
[`vignette("dkge-design-kernels")`](https://bbuchsbaum.github.io/dkge/articles/dkge-design-kernels.md);
the PLS-by-PLS comparison in
[`vignette("dkge-vs-pls")`](https://bbuchsbaum.github.io/dkge/articles/dkge-vs-pls.md).
This vignette is the orientation that ties them together.

## What DKGE actually decomposes

A naive group decomposition would eigendecompose the raw subject
covariance $`C = \sum_s B_s B_s^\top`$. DKGE never does that. It
eigendecomposes the kernel-whitened covariance

``` math
\widehat C \;=\; K^{1/2} \, C \, K^{1/2}
```

and maps the resulting eigenvectors back through $`K^{-1/2}`$ to obtain
the group basis $`U`$. Two consequences follow that are worth
internalising before choosing a kernel:

1.  **$`K`$ is a metric, not a basis.** DKGE does not project effects
    onto pre-chosen kernel directions; it changes the inner product
    under which the covariance is diagonalised. Components are free to
    be any direction in q-space — the kernel only re-weights how “large”
    each direction looks during the eigensolve.
2.  **The data $`C`$ still does the work.** If the empirical covariance
    has a dominant direction, that direction will surface regardless of
    $`K`$ (subject to scaling). The kernel amplifies directions it
    considers smooth / structurally plausible and attenuates rough or
    “unrelated” ones. It does not manufacture structure that the data
    does not already support.

In aggregate-target mode
([`dkge_aggregate_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_aggregate_fit.md))
the same machinery applies to $`Y Y^\top`$ where $`Y`$ is a cells ×
features matrix — useful when you want a PLSC-style decomposition that
still admits a design kernel. The subject-level fit
[`dkge_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit.md)
differs only in that $`C`$ is the pooled, design-row-standardised
subject covariance, not raw cell means.

### Identity-K is the unconstrained baseline

`K = diag(q)` collapses the whitening to the identity, leaving you with
an eigendecomposition of $`C`$ itself (after the standard q-space
compression). There is no kernel-imposed similarity, smoothness, or
hierarchy among effects. For an aggregate fit this is essentially the
left-singular-vector structure of the uncentered second moment
$`Y Y^\top`$ (the default `dkge_aggregate_fit(..., center = "none")`) —
ordinary cell-mean PCA, the “let the data speak” mode. Use
`center = "grand"` or `"column"` when you want a centered SVD instead.

For a full subject-level fit, identity-$`K`$ is *not* equivalent to
ignoring design entirely: rows are still GLM effects, and DKGE still
row-standardises through the pooled design Cholesky factor $`R`$. What
identity-$`K`$ means is “do not couple effects beyond their GLM/design
scaling”.

### Structured K is prior-regularised

Anything else — `design_kernel(factors, terms, ...)`, ordinal kernels,
block factors, term-weighted compositions — biases the decomposition
toward directions that lie in high-prior-variance subspaces and away
from directions the kernel calls rough. This is a feature, not a bug; it
is how design information gets into the latent geometry. But it does
mean components from a structured fit cannot be reported as if they had
“emerged” from a neutral decomposition.

**Discipline:** when reporting structured-$`K`$ components, present an
identity-$`K`$ refit alongside as the unconstrained baseline. If the
same components dominate both, the structure is in the data. If the
kernel substantially changes the loadings or eigenspectrum, the
structure is partly prior.

## Identity vs structured K — a small worked example

A 2 × 2 aggregate cell-mean matrix with four cells and forty features.
We pass the matrix straight to
[`dkge_aggregate_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_aggregate_fit.md)
— the dedicated
[`dkge_aggregate_target()`](https://bbuchsbaum.github.io/dkge/reference/dkge_aggregate_target.md)
constructor is only needed when you want subject- or weight-aware
bootstrap/permutation later (see
[`?dkge_aggregate_target`](https://bbuchsbaum.github.io/dkge/reference/dkge_aggregate_target.md)).

``` r

P <- 40
Y <- matrix(rnorm(4 * P), 4, P)
rownames(Y) <- c("A1:B1", "A1:B2", "A2:B1", "A2:B2")
```

Fit twice — identity-$`K`$ baseline and a factor-structured kernel:

``` r

fit_I <- dkge_aggregate_fit(Y, K = diag(4), rank = 3)

K_struct <- design_kernel(
  factors = list(A = list(L = 2, type = "nominal"),
                 B = list(L = 2, type = "nominal")),
  terms   = list("A", "B", c("A", "B")),
  rho     = c(A = 1, B = 1, "A:B" = 0.3),
  basis   = "cell",
  normalize = "unit_trace"
)
fit_K <- dkge_aggregate_fit(Y, K = K_struct, rank = 3)

data.frame(
  component   = 1:3,
  identity_K  = round(fit_I$singular_values, 3),
  structured_K = round(fit_K$singular_values, 3)
)
#>     component identity_K structured_K
#> LV1         1       7.46         4.64
#> LV2         2       6.96         3.61
#> LV3         3       6.64         2.63
```

Two questions to ask of this table:

- **Did the eigenspectrum change shape?** A flat-to-steep shift means
  the kernel is concentrating energy in directions it favours.
- **Did the leading loadings change?** Compare the two bases in the
  kernel metric with
  `dkge_principal_angles_K(fit_I$U, fit_K$U, fit_K$K)`. Columns of `U`
  are K-orthonormal, so `crossprod(fit_I$U, fit_K$U)` is not a rotation
  matrix and cannot be near-identity. Small principal angles mean the
  data already aligned with the kernel; large angles mean the kernel is
  doing the talking.

Either result is acceptable; the obligation is to *know which* and to
say so.

## How DKGE conducts inference

There are three distinct levels at which you can make a claim. They use
different machinery and answer different questions. Conflating them is
the single most common reporting error.

| Level | Question | Tool | Unit |
|----|----|----|----|
| **Component** | “Is this latent direction reliable?” | singular values; [`dkge_aggregate_permute()`](https://bbuchsbaum.github.io/dkge/reference/dkge_aggregate_permute.md) / [`dkge_between_permute()`](https://bbuchsbaum.github.io/dkge/reference/dkge_between_permute.md); identity vs structured $`K`$ refit; [`dkge_plot_subspace_stability()`](https://bbuchsbaum.github.io/dkge/reference/dkge_plot_subspace_stability.md) | latent dimension |
| **Contrast** | “Is this q-space contrast non-zero?” | [`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md) (LOSO / K-fold), [`dkge_bootstrap_qspace()`](https://bbuchsbaum.github.io/dkge/reference/dkge_bootstrap_qspace.md), analytic SEs | linear combination of effects |
| **Feature** | “Is this cluster / voxel salient under this contrast or component?” | [`dkge_component_stats()`](https://bbuchsbaum.github.io/dkge/reference/dkge_component_stats.md), [`dkge_transport_contrasts_to_medoid()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_contrasts_to_medoid.md), FDR adjustment | spatial unit (cluster, anchor, voxel) |

The structure is hierarchical: component-level inference asks whether a
direction exists; contrast-level inference asks whether a chosen
direction has non-zero magnitude; feature-level inference asks where,
spatially, that magnitude is concentrated. Every published finding
should name its level.

A few practical notes:

- **Subject is the unit of inference for between-subject claims.**
  Permutation and bootstrap routines resample subjects, never cells,
  anchors, or voxels. This applies to both
  [`dkge_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit.md)
  and aggregate-target workflows.
- **Cross-fitting (LOSO / K-fold) is the bias-aware default for
  contrasts.** In-fold contrasts inherit the optimism of the same data
  used to fit $`U`$; `method = "loso"` removes that.
- **Identity vs structured $`K`$ refits double as a sanity check for
  component-level claims.** This is the natural DKGE analogue of the
  multifer-style “is this LV real?” diagnostic: a component that
  survives ablating the kernel structure is one the data is asserting,
  not one the prior is asserting.

## Where DKGE sits among related methods

DKGE is not a new family of decomposition. It is a metric-aware version
of methods you already know, plus a transport and inference layer on
top.

| Method | Decomposes | Design role | Output |
|----|----|----|----|
| **PCA on cell means** | $`Y Y^\top`$ | none | data-driven axes of cell variation |
| **PLSC** (McIntosh & Lobaugh) | cross-block $`C^\top M`$ via SVD | orthonormal contrast block | covariance-maximising LVs between design and brain |
| **`dkge_aggregate_fit(K = I)`** | $`Y Y^\top`$ | none beyond cell labels | identical in spirit to cell-mean PCA |
| **`dkge_aggregate_fit(K = K_struct)`** | $`K^{1/2} Y Y^\top K^{1/2}`$ | metric / preconditioner | LVs in the kernel-induced inner product |
| **`dkge_fit(K = I)`** | pooled $`\sum_s \widetilde B_s \widetilde B_s^\top`$ | GLM row-standardisation only | unconstrained subject-level group basis |
| **`dkge_fit(K = K_struct)`** | $`K^{1/2}`$ (pooled cov) $`K^{1/2}`$ | metric / preconditioner | design-aware subject-level group basis |

Two broader relatives, for orientation:

- **RDA / CCA** project the response onto a subspace spanned by a design
  block. DKGE never restricts the response space; it changes the metric
  used on it. An RDA fit with rank $`r`$ and identity-metric is closer
  to a *strict-K* analogue where $`K`$ is rank-deficient and supported
  on the design subspace.
- **Ridge-regularised linear discriminants** (e.g. ridge LDA) likewise
  modify the metric of the underlying scatter problem. DKGE’s
  structured-$`K`$ refit is a closer cousin to ridge LDA than to either
  PCA or PLSC: both pre- and post-multiply the empirical second-moment
  matrix by a fixed PSD weight before solving.

If you find yourself wanting to say “this LV is design-driven”, the
right comparison is identity-$`K`$ DKGE, not PCA on the raw betas —
those two only coincide for aggregate fits.

## Before you fit

A short checklist that maps the above to actions:

1.  **Decide which estimand you want.** Aggregate cell means
    ([`dkge_aggregate_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_aggregate_fit.md))
    or subject-level group structure
    ([`dkge()`](https://bbuchsbaum.github.io/dkge/reference/dkge.md) /
    [`dkge_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit.md)).
    The inference unit is subjects in both cases.
2.  **Always run identity-$`K`$ at least once.** It is the unconstrained
    baseline. Save the eigenspectrum, the leading loadings, and one or
    two contrast LOSO results.
3.  **Add structure deliberately.** Build a
    [`design_kernel()`](https://bbuchsbaum.github.io/dkge/reference/design_kernel.md)
    only for terms you have prior reason to couple. Use
    [`dkge_cv_rank_loso()`](https://bbuchsbaum.github.io/dkge/reference/dkge_cv_rank_loso.md)
    to choose rank and kernel parameters jointly (see
    [`vignette("dkge-design-kernels")`](https://bbuchsbaum.github.io/dkge/articles/dkge-design-kernels.md)).
4.  **Compare identity vs structured $`K`$ refits before reporting.**
    Diff singular values; compute
    `dkge_principal_angles_K(U_I, U_K, K_struct)`; note whether
    structure is data-asserted or prior-asserted. A raw
    `crossprod(U_I, U_K)` is not a valid rotation comparison because the
    bases use the K metric.
5.  **Name the inference level.** Component, contrast, or feature. Use
    the matching tool from the table above; do not transfer a
    permutation p-value from one level to another.

A finding that survives all five steps is one you can defend.

## See also

- [`vignette("dkge-workflow")`](https://bbuchsbaum.github.io/dkge/articles/dkge-workflow.md)
  — end-to-end pipeline
- [`vignette("dkge-design-kernels")`](https://bbuchsbaum.github.io/dkge/articles/dkge-design-kernels.md)
  — kernel construction and CV-based selection
- [`vignette("dkge-contrasts-inference")`](https://bbuchsbaum.github.io/dkge/articles/dkge-contrasts-inference.md)
  — LOSO, analytic, and bootstrap inference details
- [`vignette("dkge-components")`](https://bbuchsbaum.github.io/dkge/articles/dkge-components.md)
  — component interpretation and rotation
- [`vignette("dkge-vs-pls")`](https://bbuchsbaum.github.io/dkge/articles/dkge-vs-pls.md)
  — detailed comparison with classical PLSC
- [`vignette("dkge-between-subjects")`](https://bbuchsbaum.github.io/dkge/articles/dkge-between-subjects.md)
  — subject-level RRR and permutation tests
- [`vignette("dkge-partial-effect-spaces")`](https://bbuchsbaum.github.io/dkge/articles/dkge-partial-effect-spaces.md)
  — partial cell coverage, effect-precision weighting, and trialwise
  debiasing
