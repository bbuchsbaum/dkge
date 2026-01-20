# Architecture

**Analysis Date:** 2026-01-19

## Pattern Overview

**Overall:** Domain-Driven Design with Layered Service Architecture

**Key Characteristics:**
- All heavy linear algebra operates in qxq design space (q = 10-100 effects), never forming large PxP cluster matrices
- Streaming-first design with two-pass algorithms for memory efficiency
- Leave-one-subject-out (LOSO) cross-fitting as the default for unbiased estimation
- S3 generics enable pluggable components (mappers, loaders, subjects)

## Layers

**Data Layer:**
- Purpose: Subject-level data ingestion, validation, and bundling
- Location: `R/dkge-data.R`
- Contains: `dkge_subject()` constructors (matrix, NeuroVec, ClusteredNeuroVec), `dkge_data()` bundler
- Depends on: `neuroim2`, `fmrireg` for neuroimaging types
- Used by: Fit layer, prediction layer

**Kernel Layer:**
- Purpose: Build design similarity kernels encoding factorial structure
- Location: `R/design-kernel.R`
- Contains: `design_kernel()` factory, kernel root computation, kernel alignment
- Depends on: Base R stats (contr.sum, contr.helmert)
- Used by: Fit layer

**Fit Layer (Core Engine):**
- Purpose: Estimate shared latent basis U from multi-subject GLM betas
- Location: `R/dkge-fit.R`, `R/dkge-fit-core.R`
- Contains: `dkge_fit()`, internal stages (`.dkge_fit_prepare`, `.dkge_fit_accumulate`, `.dkge_fit_solve`, `.dkge_fit_assemble`)
- Depends on: Data layer, kernel layer, `multivarious` (multiblock_biprojector)
- Used by: Contrast layer, prediction layer, transport layer

**Contrast Layer:**
- Purpose: Compute cross-fitted design contrasts with multiple strategies
- Location: `R/dkge-contrast.R`, `R/dkge-analytic.R`, `R/dkge-kfold.R`, `R/dkge-folds.R`
- Contains: `dkge_contrast()` dispatcher, LOSO/K-fold/analytic implementations
- Depends on: Fit layer, fold definitions
- Used by: Inference layer, transport layer

**Transport Layer:**
- Purpose: Map subject parcellations to a common medoid reference via optimal transport
- Location: `R/dkge-transport.R`, `R/dkge-mapper.R`
- Contains: `dkge_transport_contrasts_to_medoid()`, Sinkhorn solver, mapper S3 system
- Depends on: Fit layer, contrast layer, C++ Sinkhorn implementation
- Used by: Inference layer, rendering layer

**Inference Layer:**
- Purpose: Statistical testing with permutation-based FWER control
- Location: `R/dkge-inference.R`
- Contains: `dkge_infer()`, `dkge_signflip_maxT()`, `dkge_freedman_lane()` scaffold
- Depends on: Contrast layer, transport layer
- Used by: Pipeline layer

**Prediction Layer:**
- Purpose: Out-of-sample loadings and contrasts with frozen basis
- Location: `R/dkge-predict.R`
- Contains: `dkge_freeze()`, `dkge_predict()`, `dkge_predict_loadings()`, `predict.dkge()` S3 method
- Depends on: Fit layer (frozen model components)
- Used by: End users, classification layer

**Classification Layer:**
- Purpose: Cross-validated decoding of effect patterns
- Location: `R/dkge-classify.R`, `R/dkge-targets.R`
- Contains: `dkge_classify()`, target specification system, LDA/logistic backends
- Depends on: Fit layer, fold layer
- Used by: Pipeline layer

**Rendering Layer:**
- Purpose: Dense anchor-based group maps for visualization
- Location: `R/dkge-render-core.R`, `R/dkge-anchor-build.R`, `R/dkge-anchor-fit.R`
- Contains: Anchor construction, kNN graphs, anchor-to-voxel decoders
- Depends on: `neighborweights`, `FNN`
- Used by: End users for visualization

**Pipeline/Service Layer:**
- Purpose: High-level orchestration of fit-contrast-transport-inference workflows
- Location: `R/dkge-pipeline.R`, `R/dkge-services.R`
- Contains: `dkge_pipeline()`, service constructors (contrast_service, transport_service, inference_service)
- Depends on: All other layers
- Used by: End users as primary entry point

## Data Flow

**Standard DKGE Workflow:**

1. **Input Construction:** Subject betas + designs -> `dkge_subject()` -> `dkge_data()` bundle
2. **Kernel Setup:** Factor specifications -> `design_kernel()` -> K matrix with metadata
3. **Model Fitting:** Data + K -> `dkge_fit()` -> dkge object with U, R, Chat, Btil
4. **Contrast Computation:** Fit + contrast vector -> `dkge_contrast()` -> per-subject cluster values
5. **Transport:** Contrasts + centroids -> `dkge_transport_contrasts_to_medoid()` -> aligned group map
6. **Inference:** Transported maps -> `dkge_signflip_maxT()` -> p-values with FWER control

**LOSO Cross-Fitting Flow:**

1. For each subject s:
   - Compute Chat^(-s) = Chat - w_s * S_s (remove subject contribution)
   - Eigendecompose Chat^(-s) -> U^(-s)
   - Project held-out subject: A_s = Btil_s^T K U^(-s)
   - Compute contrast: v_s = A_s * alpha^(-s)

**State Management:**
- `dkge` objects inherit from `multivarious::multiblock_biprojector`
- Key state: U (basis), R (Cholesky ruler), K (kernel), Btil (standardized betas), Chat (compressed covariance), contribs (per-subject S_s)
- Frozen models via `dkge_freeze()` retain only U, K, R for lightweight prediction

## Key Abstractions

**dkge_subject:**
- Purpose: Single subject record with beta matrix, design, and optional omega weights
- Examples: `R/dkge-data.R` lines 68-126
- Pattern: S3 dispatch on input type (matrix, NeuroVec, ClusteredNeuroVec, list)

**dkge_data:**
- Purpose: Multi-subject bundle with harmonized effects, subject IDs, cluster IDs
- Examples: `R/dkge-data.R` lines 271-363
- Pattern: Input normalization ensuring consistent effect naming across subjects

**dkge (fit object):**
- Purpose: Complete DKGE model state after fitting
- Examples: `R/dkge-fit-core.R` lines 366-497
- Pattern: Inherits multiblock_biprojector; stores all components for downstream analysis

**dkge_contrasts:**
- Purpose: Cross-fitted contrast results with method metadata
- Examples: `R/dkge-contrast.R` lines 244-250
- Pattern: Named list structure with `values`, `method`, `contrasts`, `metadata`

**dkge_mapper_spec / dkge_mapper:**
- Purpose: Pluggable transport strategy specifications
- Examples: `R/dkge-mapper.R` lines 11-35
- Pattern: S3 generics `fit_mapper()` / `apply_mapper()` for extensibility

**dkge_weights:**
- Purpose: Voxel-level weighting specification for adaptive accumulation
- Examples: `R/dkge-weights.R` lines 33-59
- Pattern: Declarative spec with prior, adapt, combine, shrink parameters

## Entry Points

**Primary User Entry:**
- Location: `R/dkge-data.R` (`dkge()` function, lines 440-465)
- Triggers: User provides betas, designs, kernel
- Responsibilities: Bundle inputs, validate kernel, call dkge_fit, optionally store inputs

**Pipeline Entry:**
- Location: `R/dkge-pipeline.R` (`dkge_pipeline()` function, lines 28-105)
- Triggers: User provides fit (or inputs to build one), contrasts, optional transport/inference specs
- Responsibilities: Orchestrate full workflow, return comprehensive results bundle

**Prediction Entry:**
- Location: `R/dkge-predict.R` (`predict.dkge()` S3 method, lines 326-328)
- Triggers: User provides fitted model, newdata with betas and contrasts
- Responsibilities: Apply frozen basis to new subjects

## Error Handling

**Strategy:** Fail-fast with informative messages; validate at layer boundaries

**Patterns:**
- `stopifnot()` for dimension/type assertions
- Named conditions for recoverable errors (e.g., `dkge_transport_needed` when cluster counts differ)
- Graceful fallback in analytic LOSO when perturbation theory breaks down
- Warning messages for deprecated arguments with migration guidance

## Cross-Cutting Concerns

**Logging:** Message-based progress via `verbose` parameters; no dedicated logging framework

**Validation:** Input validation in constructors (`dkge_subject`, `dkge_data`); dimension checks throughout

**Parallelism:** `future.apply::future_lapply()` for subject-level parallelism; OpenMP in C++ code

**Numerical Stability:**
- Small ridge added to kernel eigenvalues (jitter parameter)
- Symmetrization of covariance matrices before eigen
- Backsolve with Cholesky for R^(-1)c
- Double precision throughout

---

*Architecture analysis: 2026-01-19*
