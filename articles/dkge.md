# Getting Started with DKGE

``` r
library(dkge)
set.seed(42)
```

Design-Kernel Group Embedding (DKGE) is an R package for multi-subject
fMRI analysis. It works in the compressed **q-space** of design effects
(typically 10–100 dimensions) rather than voxel-space, enabling
memory-efficient group fits, leave-one-subject-out contrasts, and
optimal-transport alignment to a common parcellation — all without ever
forming large subject-by-voxel matrices.

## Minimal example

Five subjects, three effects, twenty voxels each:

``` r
S <- 5; q <- 3; P <- 20
subjects <- lapply(seq_len(S), function(s) {
  beta   <- matrix(rnorm(q * P), q, P)
  design <- diag(q); colnames(design) <- paste0("eff", seq_len(q))
  dkge_subject(beta, design = design, id = paste0("sub", s))
})

fit <- dkge(subjects, K = diag(q), rank = 2)
fit          # print method shows eigenspectrum and subject count
#> Multiblock Bi-Projector object:
#>   Projection matrix dimensions:  100 x 2 
#>   Block indices:
#>     Block 1: 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20
#>     Block 2: 21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40
#>     Block 3: 41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60
#>     Block 4: 61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80
#>     Block 5: 81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100
```

The returned object is a `multiblock_biprojector` with q-space loadings
in `fit$U`, compressed projections in `fit$Btil`, and diagnostics in
`fit$weights` and `fit$evals`.

## Reading order

| Vignette                                                                                          | Purpose                                                                                                           |
|---------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------|
| **You are here**                                                                                  | Overview and entry point                                                                                          |
| [Workflow](https://bbuchsbaum.github.io/dkge/articles/dkge-workflow.md)                           | Complete fit → contrast → transport → visualise pipeline                                                          |
| [Design Kernels](https://bbuchsbaum.github.io/dkge/articles/dkge-design-kernels.md)               | Encoding factor structure with [`design_kernel()`](https://bbuchsbaum.github.io/dkge/reference/design_kernel.md)  |
| [Contrasts & Inference](https://bbuchsbaum.github.io/dkge/articles/dkge-contrasts-inference.md)   | LOSO contrasts, analytic and bootstrap inference                                                                  |
| [Classification](https://bbuchsbaum.github.io/dkge/articles/dkge-classification.md)               | Cross-validated decoding with [`dkge_classify()`](https://bbuchsbaum.github.io/dkge/reference/dkge_classify.md)   |
| [Weighting Strategies](https://bbuchsbaum.github.io/dkge/articles/dkge-weighting.md)              | Spatial (`Omega_list`), subject, and transport weights                                                            |
| [Adaptive Voxel Weighting](https://bbuchsbaum.github.io/dkge/articles/dkge-adaptive-weighting.md) | Fold-safe [`dkge_weights()`](https://bbuchsbaum.github.io/dkge/reference/dkge_weights.md) API                     |
| [Components](https://bbuchsbaum.github.io/dkge/articles/dkge-components.md)                       | Inspecting group-level components and loadings                                                                    |
| [Plotting](https://bbuchsbaum.github.io/dkge/articles/dkge-plotting.md)                           | [`theme_dkge()`](https://bbuchsbaum.github.io/dkge/reference/theme_dkge.md), scree plots, contribution maps       |
| [Dense Rendering](https://bbuchsbaum.github.io/dkge/articles/dkge-dense-rendering.md)             | Mapping cluster values to voxel space                                                                             |
| [Anchors](https://bbuchsbaum.github.io/dkge/articles/dkge-anchors.md)                             | [`dkge_anchor_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_anchor_fit.md) and anchor-based prediction |
| [CPCA](https://bbuchsbaum.github.io/dkge/articles/dkge-cpca.md)                                   | Contrastive PCA inside the DKGE framework                                                                         |
| [DKGE vs PLS](https://bbuchsbaum.github.io/dkge/articles/dkge-vs-pls.md)                          | Comparison with partial least squares approaches                                                                  |
| [Performance](https://bbuchsbaum.github.io/dkge/articles/dkge-performance.md)                     | Mapper selection, warm starts, custom backends                                                                    |
| [Architecture](https://bbuchsbaum.github.io/dkge/articles/dkge-architecture.md)                   | Design document: planned modular refactor (not current API)                                                       |

## Key concepts at a glance

**`dkge_subject(beta, design)`** — wraps a single subject’s `q × P` beta
matrix and `T × q` design matrix. Supply optional `omega` (cluster
weights) and `id`.

**`dkge_data(betas, designs)`** — bundles a list of subjects, aligns
partially-overlapping effect sets, and records provenance.

**`dkge(subjects, K, rank)`** — main entry point. `K` is the `q × q`
design kernel; pass `design_kernel(factors, basis = "effect")` for
structured factor designs.

**`dkge_contrast(fit, contrast_vec, method = "loso")`** —
leave-one-subject-out contrasts in q-space; pair with
[`dkge_transport_contrasts_to_medoid()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_contrasts_to_medoid.md)
to map to a common parcellation.

**`dkge_classify(fit, targets, method = "lda")`** — fully
cross-validated decoding with optional permutation tests. Supply targets
via `dkge_targets(fit, ~ factor1 + factor2)` or a plain contrast weight
matrix.

## Getting help

- [`?dkge`](https://bbuchsbaum.github.io/dkge/reference/dkge.md) — main
  fitting function
- [`?dkge_contrast`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md)
  — LOSO and k-fold contrasts
- [`?dkge_classify`](https://bbuchsbaum.github.io/dkge/reference/dkge_classify.md)
  — classification and decoding
- [`?design_kernel`](https://bbuchsbaum.github.io/dkge/reference/design_kernel.md)
  — structured factor kernels
- Issues: <https://github.com/bbuchsbaum/dkge/issues>
