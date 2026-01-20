# External Integrations

**Analysis Date:** 2026-01-19

## APIs & External Services

**None Detected:**
- This package does not integrate with external web APIs
- All computation is local
- No authentication, OAuth, or API key requirements

## Data Storage

**Databases:**
- None - Package operates on in-memory R objects and local files

**File Storage:**
- Local filesystem only
- Neuroimaging file formats via `neuroim2`:
  - NIfTI (.nii, .nii.gz) - Standard neuroimaging format
  - File paths passed to `neuroim2::read_vec()` for lazy loading

**Caching:**
- Internal Sinkhorn dual variable cache (`.dkge_sinkhorn_cache`)
  - LRU-style caching with max 64 entries
  - Cleared via `dkge_clear_sinkhorn_cache()`
- No external caching service (Redis, memcached, etc.)

## R Package Ecosystem Integrations

**Neuroimaging Stack (bbuchsbaum ecosystem):**

| Package | Integration Point | Purpose |
|---------|------------------|---------|
| `neuroim2` | `dkge_cluster_ts()`, `dkge_cluster_betas()`, `dkge_neuro_loader()` | NeuroVec/NeuroVol data handling, ClusteredNeuroVec support |
| `fmridesign` | `dkge_neuro_loader()` | Design matrix construction from fmridesign objects |
| `fmrireg` | `dkge_cluster_betas()`, `dkge_freedman_lane()` | GLM fitting via `fmri_ols_fit()` |
| `fmriAR` | Suggested dependency | Autoregressive noise modeling (optional) |

**Key Integration Functions:**
- `dkge_neuro_loader()` - Creates streaming loader from neuroim2 objects
  - Location: `R/dkge-neuroim2.R`
  - Accepts: `NeuroVec`, `ClusteredNeuroVec`, file paths
  - Returns: Loader with `n()`, `X(s)`, `B(s)`, `Omega(s)` methods

**Statistical/ML Packages:**

| Package | Integration Point | Purpose |
|---------|------------------|---------|
| `glmnet` | `cv.glmnet` imported | Regularized regression in classifiers |
| `FNN` | K-nearest neighbor mapper | Transport/mapping operations |
| `Matrix` | Sparse matrix support | Memory-efficient operations |
| `multivarious` | `multiblock_biprojector` class | Multiblock analysis interface compatibility |

**Visualization:**
- `ggplot2` - Full import for plotting functions
  - Custom theme: `theme_dkge()`
  - Plot functions: `dkge_plot_scree()`, `dkge_plot_effect_loadings()`, `dkge_plot_subject_contrib()`, etc.
- `patchwork` (suggested) - Multi-panel layouts in `dkge_plot_suite()`
- `ggrepel` (suggested) - Label positioning

## Authentication & Identity

**Auth Provider:**
- None required
- Package is purely computational with no user authentication

## Monitoring & Observability

**Error Tracking:**
- None (no Sentry, Bugsnag, etc.)

**Logs:**
- Standard R `message()` for progress reporting
- No structured logging framework

## CI/CD & Deployment

**Hosting:**
- GitHub repository
- GitHub Pages for documentation site

**CI Pipeline:**
- GitHub Actions workflow: `.github/workflows/pkgdown.yaml`
  - Triggers: push to main/master, pull requests, releases
  - Actions:
    1. Setup R with public RSPM
    2. Install dependencies
    3. Build pkgdown site
    4. Deploy to GitHub Pages (gh-pages branch)
  - Uses: `r-lib/actions/setup-r@v2`, `r-lib/actions/setup-r-dependencies@v2`

**No R CMD check workflow detected** - Consider adding standard R package CI

## Environment Configuration

**Required env vars:**
- None for core functionality
- `GITHUB_PAT` used in CI for GitHub API access (from `secrets.GITHUB_TOKEN`)

**Secrets location:**
- GitHub repository secrets (CI only)
- No local secrets required

## Webhooks & Callbacks

**Incoming:**
- None

**Outgoing:**
- None

## Data Format Integrations

**Input Formats:**
- Numeric matrices (R base)
- `neuroim2::NeuroVec` - 4D time-series neuroimaging data
- `neuroim2::NeuroVol` - 3D label/mask volumes
- `neuroim2::ClusteredNeuroVec` - Pre-clustered neuroimaging data
- `fmridesign` objects - Design matrix specifications

**Output Formats:**
- R objects (dkge, dkge_contrasts, dkge_inference, etc.)
- NIfTI files via `neuroim2` write functions
- Data frames via `as.data.frame()` methods

## Parallel Processing

**Backend:**
- `future.apply` - Provides `future_lapply()` for parallel subject processing
- OpenMP (optional) - C++ level parallelization in Sinkhorn solver
  - Enabled via `#pragma omp parallel for` when available
  - Compile-time detection via `#ifdef _OPENMP`

**Usage Pattern:**
```r
# User sets parallel backend
future::plan(future::multisession)
# Then runs DKGE functions that use future_lapply internally
```

## Service Dependencies Summary

| Service Type | Status | Notes |
|--------------|--------|-------|
| External APIs | None | Fully local computation |
| Databases | None | In-memory + filesystem |
| Authentication | None | No auth required |
| Caching | Internal only | Sinkhorn warm-start cache |
| CI/CD | GitHub Actions | pkgdown site deployment |
| Monitoring | None | Standard R messaging |
| File Storage | Local | NIfTI via neuroim2 |

---

*Integration audit: 2026-01-19*
