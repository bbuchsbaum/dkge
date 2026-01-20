# Phase 4: Numerical Edge Cases - Research

**Researched:** 2026-01-19
**Domain:** Numerical stability, rank-deficiency handling, NaN prevention, partial effect overlap
**Confidence:** HIGH

## Summary

This phase addresses numerical robustness of the dkge package for degenerate and edge-case inputs. The research examined the existing codebase to understand current handling of rank-deficiency, near-singular matrices, NaN/Inf propagation, and partial effect overlap across subjects.

**Key findings:**
1. The codebase already has defensive numerical patterns (pmax for eigenvalues, symmetrization, jitter) but lacks explicit rank-deficiency detection and graceful degradation messaging
2. The `kernel_roots()` function has eigenvalue clamping but limited condition number checking
3. Effect alignment in `.dkge_align_subjects_to_union()` fills missing effects with zeros, which creates a well-defined mathematical behavior but may need warnings for sparse subjects
4. NaN/Inf handling is post-hoc (replace with 0 in power iteration) rather than proactive with user feedback

**Primary recommendation:** Add early detection at data construction (`dkge_subject()`/`dkge_data()`) with informative warnings, implement graceful degradation with rank reduction, and track/report metadata about exclusions and degradations in the returned fit object.

## Standard Stack

This phase requires no new dependencies - it extends existing numerical patterns in base R.

### Core
| Function | Location | Purpose | Why Standard |
|----------|----------|---------|--------------|
| `eigen()` | base R | Eigendecomposition | Already used throughout; symmetric=TRUE ensures real eigenvalues |
| `svd()` | base R | Condition number via singular values | Standard for rank estimation |
| `qr()` | base R | QR factorization for rank detection | Standard for rank-deficient detection |
| `chol()` | base R | Cholesky factorization | Already used in `.dkge_compute_shared_ruler()` |
| `kappa()` | base R | Condition number estimation | R's built-in for matrix conditioning |

### Supporting Utilities Already in Codebase
| Function | Location | Purpose |
|----------|----------|---------|
| `pmax()` | base R | Eigenvalue floor at small positive values |
| `(M + t(M))/2` | Pattern | Symmetrization for numerical stability |
| `diag() + jitter` | Pattern | Ridge stabilization |

### Recommended New Internal Functions
| Function | Purpose | Location |
|----------|---------|----------|
| `.dkge_check_rank()` | Early rank detection with warnings | `R/dkge-data.R` or `R/dkge-utils.R` |
| `.dkge_voxel_exclusion_mask()` | Track NaN/Inf exclusions | `R/dkge-fit-core.R` |
| `.dkge_condition_number()` | Condition number check with threshold | `R/dkge-utils.R` |

**Installation:** No new packages required.

## Architecture Patterns

### Current Numerical Safety Patterns

The codebase already uses several numerical safety patterns:

**1. Eigenvalue flooring in `.dkge_kernel_roots()`:**
```r
# Source: R/dkge-fit.R:47-57
.dkge_kernel_roots <- function(K) {
  Ksym <- (K + t(K)) / 2
  eigK <- eigen(Ksym, symmetric = TRUE)
  vals <- pmax(eigK$values, 1e-10)  # Floor at 1e-10
  # ...
}
```

**2. Eigenvalue flooring in `.dkge_fit_solve()`:**
```r
# Source: R/dkge-fit-core.R:213-218
pos_idx <- eig_values > 1e-12
if (!all(pos_idx)) {
  eig_vectors <- eig_vectors[, pos_idx, drop = FALSE]
  eig_values <- eig_values[pos_idx]
  rank <- length(eig_values)
}
```

**3. NaN replacement in power iteration:**
```r
# Source: R/dkge-fit.R:69-70
if (!all(is.finite(X))) {
  X[!is.finite(X)] <- 0
}
```

**4. Symmetric square root with flooring:**
```r
# Source: R/design-kernel.R:241-253
vals <- pmax(vals, jitter)  # Clamp at jitter threshold
n_clamped <- sum(pinned)
if (n_clamped > 0 && n_clamped / max(n_total, 1) > 0.1) {
  warning(sprintf("%.0f%% of eigenvalues were clamped during kernel root computation.", ...))
}
```

### Recommended Project Structure

```
R/
├── dkge-data.R          # Add .dkge_check_rank() call in dkge_subject()
├── dkge-fit-core.R      # Add voxel exclusion tracking, degradation metadata
├── dkge-utils.R         # Add .dkge_condition_number(), .dkge_check_finite()
└── dkge-diagnostics.R   # (optional) Consolidated numerical diagnostics
```

### Pattern 1: Early Rank Detection

**What:** Check matrix rank at data construction time
**When to use:** In `dkge_subject()` and `dkge_data()` to detect issues early
**Example:**
```r
# Source: Recommended new pattern based on existing codebase analysis
.dkge_check_rank <- function(beta, design, subject_id = NULL) {
  # Check design matrix rank via QR
  qr_design <- qr(design)
  design_rank <- qr_design$rank
  expected_rank <- ncol(design)

  if (design_rank < expected_rank) {
    warning(sprintf(
      "Subject '%s': design matrix is rank-deficient (rank %d < %d columns). Effects may be aliased.",
      subject_id %||% "(unnamed)", design_rank, expected_rank
    ), call. = FALSE)
  }

  # Check beta matrix for rank (less critical but informative)
  beta_rank <- qr(beta)$rank
  if (beta_rank < nrow(beta)) {
    warning(sprintf(
      "Subject '%s': beta matrix has reduced rank (%d < %d effects).",
      subject_id %||% "(unnamed)", beta_rank, nrow(beta)
    ), call. = FALSE)
  }

  list(design_rank = design_rank, beta_rank = beta_rank)
}
```

### Pattern 2: Graceful Rank Reduction with Metadata

**What:** When eigenvalues are too small, reduce rank rather than error
**When to use:** In `.dkge_fit_solve()` when requested rank exceeds available rank
**Example:**
```r
# Source: Recommended enhancement to R/dkge-fit-core.R
# After eigendecomposition
effective_rank <- sum(eig_values > 1e-12)
if (rank > effective_rank) {
  warning(sprintf(
    "Requested rank %d exceeds effective rank %d. Reducing to %d components.",
    rank, effective_rank, effective_rank
  ), call. = FALSE)
  rank <- effective_rank
}
# Store in fit object
fit$effective_rank <- effective_rank
fit$rank_reduced <- (rank < rank_requested)
```

### Pattern 3: NaN Exclusion Tracking

**What:** Track which voxels are excluded due to NaN/Inf and report in metadata
**When to use:** When processing beta matrices containing NA/NaN/Inf
**Example:**
```r
# Source: Recommended pattern for R/dkge-fit-core.R
.dkge_voxel_exclusion_mask <- function(B_list, subject_ids = NULL) {
  S <- length(B_list)
  excluded_voxels <- vector("list", S)
  excluded_counts <- integer(S)

  for (s in seq_len(S)) {
    B <- B_list[[s]]
    bad_cols <- which(colSums(!is.finite(B)) > 0)
    excluded_voxels[[s]] <- bad_cols
    excluded_counts[s] <- length(bad_cols)

    if (length(bad_cols) > 0) {
      pct <- 100 * length(bad_cols) / ncol(B)
      warning(sprintf(
        "Subject '%s': %d voxels (%.1f%%) excluded due to NA/NaN/Inf values.",
        subject_ids[s] %||% s, length(bad_cols), pct
      ), call. = FALSE)
    }
  }

  list(
    excluded_voxels = excluded_voxels,
    excluded_counts = excluded_counts,
    total_excluded = sum(excluded_counts)
  )
}
```

### Pattern 4: Condition Number Checking

**What:** Check condition number against threshold (1e8 per CONTEXT.md)
**When to use:** After computing pooled Gram matrix, before Cholesky
**Example:**
```r
# Source: Recommended pattern for R/dkge-fit.R
.dkge_check_condition <- function(M, threshold = 1e8, name = "matrix") {
  cond <- kappa(M, exact = FALSE)
  if (cond > threshold) {
    warning(sprintf(
      "%s is ill-conditioned (condition number: %.2e > %.2e threshold). Results may be numerically unstable.",
      name, cond, threshold
    ), call. = FALSE)
  }
  cond
}
```

### Anti-Patterns to Avoid
- **Silent NaN propagation:** Don't let NaN silently corrupt results - either exclude and warn, or error
- **Hard failures on rank deficiency:** Don't error on rank-deficient matrices when graceful reduction is possible
- **Cryptic error messages:** Don't produce errors like "matrix not positive definite" without identifying which subject/effect caused it
- **Unchecked inversions:** Don't use `solve()` without checking condition number first

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Rank estimation | Manual SVD thresholding | `qr()$rank` or `Matrix::rankMatrix()` | Already handles tolerances correctly |
| Condition number | Manual eigenvalue ratio | `kappa()` or `rcond()` | Handles edge cases, uses LAPACK |
| Symmetric eigendecomposition | Custom iteration | `eigen(M, symmetric = TRUE)` | Guaranteed real eigenvalues, optimized |
| Cholesky with fallback | Manual positive-definiteness fix | Existing `.dkge_kernel_roots()` pattern with clamping | Already handles near-singular cases |
| NaN detection | Column-by-column checks | `colSums(!is.finite(M)) > 0` | Vectorized, efficient |

**Key insight:** R's base linear algebra functions (backed by LAPACK) already handle most numerical edge cases correctly. The main gap is _user-facing feedback_ about what went wrong and _metadata_ about degradations.

## Common Pitfalls

### Pitfall 1: Silent Rank Reduction
**What goes wrong:** Eigenvalue flooring reduces effective rank without user notification
**Why it happens:** The `pmax(vals, 1e-10)` pattern silently clamps small eigenvalues
**How to avoid:** Count clamped eigenvalues and warn if significant portion (>10%) are clamped
**Warning signs:** Number of returned components less than requested, unexpectedly low explained variance

### Pitfall 2: Design Matrix Singularity from Collinear Effects
**What goes wrong:** Cholesky of pooled Gram fails because design matrices are rank-deficient
**Why it happens:** User provides collinear or constant columns in design matrix
**How to avoid:** Check `qr(design)$rank` in `dkge_subject()` before proceeding
**Warning signs:** Error "matrix not positive definite" from `chol()`

### Pitfall 3: NaN from Zero Voxels
**What goes wrong:** Division by zero or log(0) creates NaN in downstream computations
**Why it happens:** Voxels with all-zero or constant betas have zero variance
**How to avoid:** Pre-filter or weight down zero-variance voxels; check `colSums(B^2) > 0`
**Warning signs:** NaN values in weights, projections, or contrast values

### Pitfall 4: Ill-Conditioning Amplification in Inverse Square Root
**What goes wrong:** `Kihalf` computation amplifies noise when K is ill-conditioned
**Why it happens:** `1/sqrt(tiny_eigenvalue)` produces huge values
**How to avoid:** Already handled by eigenvalue flooring, but should warn user
**Warning signs:** Very large values in Kihalf, condition number > 1e8

### Pitfall 5: Partial Effect Overlap Creates Zero Rows/Columns
**What goes wrong:** Subjects with few effects create sparse Gram matrices
**Why it happens:** `.dkge_align_subjects_to_union()` fills missing effects with zeros
**How to avoid:** Track and warn for subjects with >50% missing effects per CONTEXT.md
**Warning signs:** Low `pair_counts` diagonal values in provenance, many zero rows in aligned beta

### Pitfall 6: Power Iteration Non-Convergence
**What goes wrong:** Leading singular value estimation fails or gives wrong result
**Why it happens:** Zero matrix, all-NaN matrix, or numerical issues
**How to avoid:** Already handled with `if (!is.finite(w_norm) || w_norm == 0) break`
**Warning signs:** Very small or zero subject weights despite non-zero data

## Code Examples

### Existing Pattern: Eigenvalue Clamping in kernel_roots()
```r
# Source: R/design-kernel.R:234-278
kernel_roots <- function(K, jitter = 1e-10) {
  if (!isTRUE(isSymmetric(K, tol = 1e-8))) {
    warning("Kernel matrix not symmetric; applying symmetrization via (K + t(K))/2.")
  }
  Ks <- (K + t(K)) / 2
  ee <- eigen(Ks, symmetric = TRUE)
  vals <- ee$values
  V <- ee$vectors

  if (is.null(jitter)) {
    pinned <- vals <= 0
    tiny <- 1e-10
    vals[pinned] <- tiny
  } else {
    pinned <- vals < jitter
    vals <- pmax(vals, jitter)
  }

  n_clamped <- sum(pinned)
  n_total <- length(vals)
  if (n_clamped > 0 && n_clamped / max(n_total, 1) > 0.1) {
    warning(sprintf("%.0f%% of eigenvalues were clamped during kernel root computation.",
                    100 * n_clamped / n_total))
  }
  # ... rest of function
}
```

### Existing Pattern: Effect Union with Zero Fill
```r
# Source: R/dkge-align-data.R:31-47
embed_beta <- function(beta, current_ids) {
  stopifnot(nrow(beta) == length(current_ids))
  out <- matrix(0, n_union, ncol(beta),
                dimnames = list(union_ids, colnames(beta)))
  idx <- match(current_ids, union_ids)
  out[idx, ] <- beta
  out
}
```

### Existing Pattern: NaN Safety in Power Iteration
```r
# Source: R/dkge-fit.R:67-97
.dkge_leading_sv_squared <- function(X, tol = 1e-6, max_iter = 50) {
  X <- as.matrix(X)
  if (!all(is.finite(X))) {
    X[!is.finite(X)] <- 0  # Replace NaN/Inf with 0
  }
  # ... rest of power iteration
  if (!is.finite(v_norm) || v_norm == 0) {
    return(0)  # Graceful fallback
  }
  # ...
}
```

### Existing Pattern: Positive Eigenvalue Selection in Solve
```r
# Source: R/dkge-fit-core.R:213-218
pos_idx <- eig_values > 1e-12
if (!all(pos_idx)) {
  eig_vectors <- eig_vectors[, pos_idx, drop = FALSE]
  eig_values <- eig_values[pos_idx]
  rank <- length(eig_values)
}
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Silent eigenvalue clamping | Warning when >10% clamped | Already in `kernel_roots()` | Users know when kernel is ill-conditioned |
| Hard error on Cholesky failure | Not yet implemented | Target for Phase 4 | Need graceful degradation |
| No NaN tracking | NaN replaced with 0 silently | Current | Need user notification |
| No condition number check | Not yet implemented | Target for Phase 4 | Need 1e8 threshold warning |

**Current gaps:**
- No early rank detection in `dkge_subject()`/`dkge_data()`
- No condition number checking against 1e8 threshold
- No metadata in fit object tracking degradations/exclusions
- No sparse subject warning (>50% missing effects threshold)

## Open Questions

Things that couldn't be fully resolved:

1. **Minimum viable rank threshold heuristic**
   - What we know: Should error if rank drops below some minimum
   - What's unclear: Exact formula (e.g., `max(2, requested_rank * 0.5)`?)
   - Recommendation: Start with `min(2, requested_rank)` and refine based on testing

2. **Eigenvalue "effectively zero" threshold**
   - What we know: Currently uses 1e-10 (kernel roots) and 1e-12 (fit solve)
   - What's unclear: Whether these should be unified or remain context-specific
   - Recommendation: Keep context-specific but document rationale (kernel vs covariance)

3. **Per-subject vs aggregate rank checking**
   - What we know: Individual subjects may have rank issues
   - What's unclear: When to drop a subject vs reduce global rank
   - Recommendation: Per CONTEXT.md - drop subject with warning, don't fail entire fit

## Sources

### Primary (HIGH confidence)
- R/dkge-fit.R: Reviewed lines 1-312 for numerical patterns
- R/dkge-fit-core.R: Reviewed lines 1-497 for solve/assemble patterns
- R/design-kernel.R: Reviewed lines 230-302 for kernel_roots implementation
- R/dkge-align-data.R: Reviewed lines 1-127 for effect alignment patterns
- R/dkge-align-effects.R: Reviewed lines 1-391 for kernel alignment with partial overlap
- R/dkge-analytic.R: Reviewed lines 1-418 for perturbation handling

### Secondary (MEDIUM confidence)
- R/dkge-data.R: Reviewed lines 1-465 for data construction patterns
- tests/testthat/test-data-validation.R: Reviewed existing edge case test patterns
- tests/testthat/test-kernel-invariants.R: Reviewed numerical tolerance patterns

### Context Document (PRIMARY constraint)
- .planning/phases/04-numerical-edge-cases/04-CONTEXT.md: User decisions on error vs degradation policy, NaN handling, tolerance thresholds, partial overlap behavior

## Metadata

**Confidence breakdown:**
- Numerical patterns in codebase: HIGH - Direct code review of all relevant files
- Existing tolerance thresholds: HIGH - Documented in code and CONTEXT.md
- Recommended new patterns: MEDIUM - Based on codebase style but not yet implemented
- Minimum rank heuristics: LOW - Needs empirical validation during testing

**Research date:** 2026-01-19
**Valid until:** 30 days (stable numerical patterns, no external dependencies)
