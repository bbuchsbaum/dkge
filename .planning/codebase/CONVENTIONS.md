# Coding Conventions

**Analysis Date:** 2026-01-19

## Naming Patterns

**Files:**
- Source files: `dkge-{module}.R` (hyphen-separated, lowercase)
- Examples: `R/dkge-fit.R`, `R/dkge-contrast.R`, `R/dkge-weights.R`
- Test files: `test-{module}.R` (hyphen-separated)
- Helper files: `helper-{purpose}.R` for shared test utilities

**Functions:**
- Public functions: `dkge_{verb}_{noun}()` pattern with underscores
- Examples: `dkge_fit()`, `dkge_contrast()`, `dkge_predict_loadings()`
- Aliases: Some have shorter forms (e.g., `dkge()` as alias for `dkge_fit()`)
- Internal functions: `.dkge_{description}()` with leading dot and underscores
- Examples: `.dkge_fit_prepare()`, `.dkge_contrast_loso()`, `.dkge_apply()`

**Variables:**
- snake_case for local variables: `subject_ids`, `contrast_list`, `voxel_weights`
- Single uppercase letters for mathematical objects: `K`, `U`, `R`, `X`, `B`
- List suffixes: `_list` for subject-indexed lists: `B_list`, `X_list`, `Omega_list`
- Index suffixes: `_idx`, `_s` (subject index): `train_idx`, `Bs`, `fs`

**Types/Classes:**
- S3 classes: `dkge`, `dkge_data`, `dkge_subject`, `dkge_weights`, `dkge_contrasts`
- snake_case with `dkge_` prefix for all custom classes

**Constants/Parameters:**
- Mathematical parameters: single letters with subscripts via underscore
- Examples: `w_method`, `w_tau`, `cpca_T`, `cpca_ridge`

## Code Style

**Formatting:**
- No explicit .Rprofile or .editorconfig detected
- Implied 2-space indentation based on source files
- Line length appears flexible (~80-100 characters)

**Linting:**
- No `.lintr` configuration file present
- R CMD check is the primary quality gate

**Documentation:**
- roxygen2 with `@export`, `@param`, `@return`, `@examples`, `@seealso`
- Internal functions use `@keywords internal` and `@noRd`
- Block comments use `# ----` separators for sections within files

## Import Organization

**Order:**
1. Base R (`stats`, `utils`) - imported via `@importFrom` in `R/dkge-package.R`
2. External packages (Rcpp, Matrix, ggplot2, etc.)
3. Internal dependencies via `:::` for internal functions

**Path Aliases:**
- None; direct package qualification used where needed (e.g., `future.apply::future_lapply()`)

**Namespace Management:**
- `@importFrom stats aggregate contr.helmert contr.sum p.adjust pt setNames`
- `@importFrom utils head`
- Global variables declared in `R/globals.R` for ggplot2 NSE symbols

## Error Handling

**Patterns:**
- `stopifnot()` for argument validation:
```r
stopifnot(is.matrix(K), nrow(K) == q, ncol(K) == q)
stopifnot(inherits(fit, "dkge"))
```

- `stop()` with `call. = FALSE` for user-facing errors:
```r
stop("Row names of beta matrix must match design column names (effects).", call. = FALSE)
stop("Provide `B_list` via `newdata` or the `B_list` argument.", call. = FALSE)
```

- `warning()` with `call. = FALSE` for deprecation:
```r
warning("`transport` argument to dkge_contrast() is deprecated; use `dkge_transport_contrasts_to_medoid()`.",
        call. = FALSE)
```

- Custom error conditions for specific cases:
```r
cond <- structure(list(message = msg, call = sys.call()),
                 class = c("dkge_transport_needed", "error", "condition"))
stop(cond)
```

- `tryCatch()` for graceful degradation:
```r
K_effects <- tryCatch(
  .dkge_weights_resolve_k(weights, kernel_info),
  error = function(e) {
    warning(sprintf("k-energy weights unavailable (%s); using uniform weights", e$message))
    NULL
  }
)
```

## Logging

**Framework:** Base R `message()` and `warning()`

**Patterns:**
- Verbose progress via `message()`:
```r
if (verbose_flag) {
  message(sprintf("Computing %d contrast(s) via LOSO for %d subjects", n_contrasts, S))
}
```

- Global option toggle:
```r
.dkge_verbose <- function(verbose) {
  isTRUE(verbose) && isTRUE(getOption("dkge.verbose", TRUE))
}
```

## Comments

**When to Comment:**
- Section headers using `# ---- Title ----` delimiter pattern
- Mathematical explanations for non-obvious operations
- Brief descriptions of complex conditional logic

**JSDoc/TSDoc:**
- Not applicable (R package)

**Roxygen2 Documentation:**
- All exported functions fully documented
- `@details` section for algorithmic explanation
- `@examples` with `\dontrun{}` for examples requiring external dependencies

## Function Design

**Size:**
- Public functions decomposed into internal helpers (prepare/accumulate/solve/assemble)
- Internal functions kept focused on single responsibilities
- Example: `dkge_fit()` calls `.dkge_fit_prepare()`, `.dkge_fit_accumulate()`, `.dkge_fit_solve()`, `.dkge_fit_assemble()`

**Parameters:**
- Defaults provided for most parameters
- `match.arg()` for string options:
```r
w_method <- match.arg(w_method)
cpca_part <- match.arg(cpca_part)
```

- `%||%` null-coalescing operator for defaults:
```r
subject_labels <- fit$subject_ids %||% paste0("subject", seq_len(S))
```

**Return Values:**
- Named lists with `structure(..., class = "classname")` for S3 objects
- Rich metadata attached to return objects
- Consistent field names across similar object types

## Module Design

**Exports:**
- Public API functions exported via `@export` roxygen tag
- One primary entry point per module (e.g., `dkge_fit()`, `dkge_contrast()`)
- Convenience wrappers for common workflows

**Barrel Files:**
- Not used; R package NAMESPACE handles exports
- `R/dkge-package.R` provides package-level documentation

**Internal Organization:**
- File header comment describing module purpose
- Section separators: `# ---- Section Name ----`
- Related internal helpers grouped in same file

## S3 Method Conventions

**print methods:**
- Concise summary of object contents
- Return `invisible(x)` at end
- Use `cat()` for output formatting

**as.data.frame methods:**
- Support `stringsAsFactors` parameter
- Include metadata columns (method, correction, etc.)
- Handle empty cases gracefully

**as.matrix methods:**
- Check for compatible dimensions
- Raise informative errors when conversion not possible

## Mathematical Code Conventions

**Matrix Operations:**
- Explicit symmetrization: `(M + t(M)) / 2`
- Small ridge for numerical stability: `diag(G_pool) <- diag(G_pool) + jitter`
- Use `pmax()` to clamp small/negative eigenvalues
- Named results from eigen decomposition: `eigK$values`, `eigK$vectors`

**Precision Handling:**
- Default jitter: `1e-10` to `1e-8` range
- Small epsilon to avoid division by zero: `1e-12`

---

*Convention analysis: 2026-01-19*
