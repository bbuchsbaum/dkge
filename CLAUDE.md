# CLAUDE.md

Guidance for Claude Code (and other agents — `AGENTS.md` is a symlink to this file) when working in this repository.

## Package overview

`dkge` implements **Design-Kernel Group Embedding** for design-aware cluster-level fMRI analysis. Heavy linear algebra runs in q×q effect space (q ≈ 10–100), never in P×P voxel space. The package covers batch and streaming estimators, leave-one-subject-out cross-fitting, kernel-based design alignment, optimal transport to medoid parcellations, between-subject modeling, and a classification/RRR layer on top of the latent embedding.

Whitepaper: `data-raw/algo.md`.

## Build and development commands

```bash
R CMD build .
R CMD INSTALL dkge_*.tar.gz
R CMD check --as-cran dkge_*.tar.gz

Rscript -e "devtools::load_all()"
Rscript -e "devtools::document()"     # regenerate man pages from roxygen2
Rscript -e "devtools::test()"         # ~87 testthat files in tests/testthat
Rscript -e "pkgdown::build_site()"
```

## Architecture

### R/ — main code (~60 files, grouped by concern)

- **Core fit & prediction**: `dkge-fit.R`, `dkge-fit-core.R`, `dkge-fit-from-kernels.R`, `dkge-predict.R`, `dkge-input.R`, `dkge-pipeline.R`, `dkge-services.R`
- **Cross-fitting & resampling**: `dkge-loso.R`, `dkge-folds.R`, `dkge-kfold.R`, `dkge-cv.R`, `dkge-bootstrap.R`, `dkge-permute.R`
- **Design kernels & alignment**: `design-kernel.R`, `dkge-procrustes.R`, `dkge-align-data.R`, `dkge-align-effects.R`, `dkge-jd.R`, `compat_neuralign.R`
- **Anchors / CPCA / projection**: `dkge-anchor-build.R`, `dkge-anchor-fit.R`, `dkge-anchor-targets.R`, `dkge-cpca.R`, `dkge-project.R`, `dkge-target.R`, `dkge-targets.R`
- **Between-subject layer**: `dkge-between-rrr.R`, `dkge-between-infer.R`, `dkge-subject-model.R`
- **Inference & contrasts**: `dkge-inference.R`, `dkge-analytic.R`, `dkge-contrast.R`, `dkge-contrast-validated.R`, `dkge-regress.R`
- **Classification on latent codes**: `dkge-classify.R`, `dkge-classify-core.R`, `dkge-classify-cv.R`, `dkge-classify-backends.R`, `dkge-latent-clf.R`, `dkge-latent-utils.R`
- **Components, weights, metrics**: `dkge-components.R`, `dkge-weights.R`, `dkge-metrics.R`
- **OT, mappers, neuro IO**: `dkge-transport.R`, `dkge-mapper.R`, `dkge-info-maps.R`, `dkge-neuroim2.R`, `dkge-voxel.R`
- **Rendering, plotting, output**: `dkge-render-core.R`, `dkge-plot.R`, `dkge-plot-suite.R`, `dkge-write.R`
- **Misc**: `dkge-data.R`, `dkge-sim.R`, `dkge-specs.R`, `dkge-utils.R`, `globals.R`, `hyperdesign-generics.R`, `zzz.R`, `RcppExports.R`

### src/ — C++ acceleration

- `pwdist.cpp` — pairwise distances
- `sinkhorn.cpp` — Sinkhorn optimal transport
- `RcppExports.cpp` — auto-generated

### future/ — staged work-in-progress (NOT exported)

`dkge-stream.R`, `dkge-cpca.R`, `dkge-cv.R`, `dkge-mixed.R`, `dkge-pipeline.R`, `dkge-neuroim2-stream.R`, `dkge-sinkhorn.R`, `dkge-sinkhorn-cpp.R`, `dkge-write.R`. These are not loaded by the package; treat them as design drafts. If you graduate one, move it into `R/` and reconcile any name collision with the existing file.

### Vignettes (`vignettes/*.Rmd`)

`dkge`, `dkge-workflow`, `dkge-architecture`, `dkge-design-kernels`, `dkge-anchors`, `dkge-cpca`, `dkge-components`, `dkge-contrasts-inference`, `dkge-classification`, `dkge-between-subjects`, `dkge-weighting`, `dkge-adaptive-weighting`, `dkge-dense-rendering`, `dkge-plotting`, `dkge-performance`, `dkge-vs-pls`.

## Mathematical foundation

- Compress to q×q matrices where q = number of design effects.
- Design kernel `K` (q×q PSD) encodes factor structure: nominal, ordinal, circular, interactions.
- LOSO cross-fitting yields unbiased contrasts in q-space.
- Optimal transport maps subject parcels to a medoid parcellation.

Key data shapes:

- `B_list`: list of q×P_s GLM beta matrices per subject
- `X_list`: list of T_s×q design matrices per subject
- `K`: q×q design kernel
- `U`: q×r group latent basis (K-orthonormal)

## Implementation rules

- All heavy linear algebra in q×q space. **Never** form P_s×P_s.
- Double precision; symmetrize after every quadratic form.
- Add a small ridge to eigenvalues when matrices are near-singular.
- Streaming variants are two-pass: pass 1 computes the pooled design Cholesky `R`; pass 2 accumulates compressed covariance.
- K-orthonormality in `dkge_fit()` is achieved by eigendecomposition, not Procrustes.

## Dependencies

**Imports**: `fmridesign`, `fmrireg`, `FNN`, `future.apply`, `Matrix`, `methods`, `adjoin`, `neuroim2`, `stats`, `Rcpp`, `multivarious`, `utils`, `ggplot2`, `glmnet`.

**Suggests**: `fmriAR`, `microbenchmark`, `neuralign`, `pkgdown`, `remotes`, `testthat`, `T4transport`, `knitr`, `rmarkdown`, `patchwork`, `ggrepel`, `tidyr`, `purrr`.

**LinkingTo**: `Rcpp`, `RcppArmadillo`.

**Remotes (GitHub)**: `bbuchsbaum/{neuralign, neuroim2, fmridesign, fmriAR, fmrireg, adjoin}`.

## neuralign migration

A migration replaces dkge's K-Procrustes alignment code with neuralign's equivalents. See `docs/MIGRATION_NEURALIGN.md`.

- 1:1 map: `dkge_k_orthonormalize`/`dkge_procrustes_K`/`dkge_align_bases_K`/`dkge_consensus_basis_K` → `k_orthonormalize`/`k_procrustes`/`k_align_bases`/`k_consensus_basis`.
- Core `dkge_fit()` is unaffected (eigendecomposition path).
- Only post-fit analysis (bootstrap, analytic, folds, plot, sim) calls K-Procrustes — 12 call sites across 6 files.
- `compat_neuralign.R` is the adapter seam.
- 24 golden tests gate each phase.
- Design kernels, Sinkhorn OT, spatial mappers, and effect-kernel alignment stay in dkge.
