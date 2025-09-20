# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

This is an R package implementing Design-Kernel Group Embedding (DKGE) for cluster-level fMRI analysis. It provides streaming and memory-light estimators, leave-one-subject-out cross-fitting, kernel-based design alignment, and optimal transport to medoid parcellations.

## Build and Development Commands

```bash
# Build and install the package
R CMD build .
R CMD INSTALL dkge_*.tar.gz

# Check package for CRAN submission
R CMD check --as-cran dkge_*.tar.gz

# Load package in R for testing
Rscript -e "devtools::load_all()"

# Document package (generate man pages from roxygen2)
Rscript -e "devtools::document()"

# Run tests (when implemented)
Rscript -e "devtools::test()"

# Build package documentation site
Rscript -e "pkgdown::build_site()"
```

## Architecture

### Core Components

1. **Main algorithms** (R/):
   - `dkge-fit.R`: Core DKGE batch fit and LOSO contrasts in q-space
   - `dkge-predict.R`: Out-of-sample prediction with frozen basis
   - `dkge-cpca.R`: CPCA inside-span filtering
   - `dkge-procrustes.R`: Procrustes alignment utilities
   - `dkge-pipeline.R`: High-level pipeline functions

2. **Future implementations** (future/):
   - `dkge-stream.R`: Streaming two-pass implementation
   - `dkge-cv.R`: Cross-validation and model selection
   - `dkge-inference.R`: Statistical inference (sign-flip, Freedman-Lane)
   - `dkge-sinkhorn.R`, `dkge-sinkhorn-cpp.R`: Optimal transport
   - `dkge-mixed.R`: Mixed-effects extensions
   - `dkge-neuroim2-stream.R`: Integration with neuroim2 for on-disk data
   - `dkge-write.R`: Output utilities

3. **C++ acceleration** (src/):
   - `pwdist.cpp`: Pairwise distance computations

### Mathematical Foundation

The package implements algorithms from `data-raw/algo.md`, which specifies:
- Compression to q×q matrices where q = number of design effects (typically 10-100)
- Design kernel K encoding factor structure (nominal, ordinal, circular, interactions)
- LOSO cross-fitting for unbiased estimation
- Optimal transport for mapping to medoid parcellations

Key data structures:
- `B_list`: List of q×P_s GLM beta matrices per subject
- `X_list`: List of T_s×q design matrices per subject
- `K`: q×q design kernel (positive semidefinite)
- `U`: q×r group latent basis (K-orthonormal)

## Package Dependencies

Core dependencies:
- `fmridesign`, `fmrireg`: GLM and design matrix handling
- `neuroim2`: Neuroimaging data structures
- `lme4`, `lmerTest`: Mixed-effects models
- `future.apply`: Parallel processing

Remote dependencies (GitHub):
- `bbuchsbaum/neuroim2`
- `bbuchsbaum/fmridesign`
- `bbuchsbaum/fmriAR`
- `bbuchsbaum/fmrireg`

## Key Implementation Notes

- All heavy linear algebra operates in q×q space for memory efficiency
- Never form large P_s×P_s matrices
- Use double precision and symmetrization for numerical stability
- Add small ridge to eigenvalues when needed
- Streaming variants process subjects in two passes: first to compute the pooled design Cholesky factor R, then to accumulate compressed covariance
