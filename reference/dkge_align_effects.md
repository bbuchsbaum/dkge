# Align effect kernels across subjects with partial overlap

Builds a fold-aware common effect geometry using masked averaging over
the training subjects, then completes each subject's kernel on the
shared effect index via Nyström, shrinkage, or intersection.

## Usage

``` r
dkge_align_effects(
  K_list,
  effects,
  subject_ids = NULL,
  folds = NULL,
  mode = c("nystrom", "shrinkage", "intersection"),
  weights = NULL,
  ridge = 1e-06,
  alpha = 0.25,
  ensure_psd = TRUE,
  psd_tol = 1e-10,
  min_train_coverage = 1L,
  intersection_scope = c("all_subjects", "train_only"),
  effect_prior = NULL,
  prior_weight = 0,
  verbose = FALSE
)
```

## Arguments

- K_list:

  List of per-subject symmetric kernels; \`K_list\[\[s\]\]\` has
  dimensions \`\|O_s\| x \|O_s\|\`.

- effects:

  List of character vectors. \`effects\[\[s\]\]\` gives the effect IDs
  associated with the rows and columns of \`K_list\[\[s\]\]\`.

- subject_ids:

  Optional character vector naming subjects; defaults to the list names
  or sequential labels.

- folds:

  Optional cross-validation structure defining held-out subjects per
  fold (see Details).

- mode:

  Completion mode. One of `"nystrom"`, `"shrinkage"`, or
  `"intersection"`.

- weights:

  Optional numeric vector of subject weights used when pooling the
  training kernels.

- ridge:

  Ridge factor used when inverting training blocks for Nyström.

- alpha:

  Shrinkage weight applied in `mode = "shrinkage"`.

- ensure_psd:

  Logical; when \`TRUE\` (default) project pooled and completed matrices
  to the PSD cone.

- psd_tol:

  Eigenvalue floor expressed as a fraction of the largest eigenvalue
  when projecting to PSD.

- min_train_coverage:

  Drop effects observed by fewer than this many training subjects when
  forming the union.

- intersection_scope:

  When `mode = "intersection"`, restrict the intersection to
  `"all_subjects"` (default) or `"train_only"`.

- effect_prior:

  Optional PSD matrix indexed by effect IDs used to seed zero-coverage
  entries of the group kernel.

- prior_weight:

  Blend factor in \[0, 1\] applied to \`effect_prior\` when available.

- verbose:

  Logical; emit messages when \`TRUE\`.

## Value

When \`folds = NULL\`, a list with fields - \`K_aligned\`: list of
aligned \`n x n\` kernels per subject - \`effect_ids\`: character vector
of shared effect IDs - \`G\`: pooled training kernel (when applicable) -
\`obs_mask\`: list of logical vectors indicating observed effects per
subject - \`pair_counts\`: integer matrix of training coverage per
effect pair - \`coverage\`: data frame summarising training coverage per
effect - \`mode\`: completion mode used.

When folds are supplied, returns \`list(folds = list(...))\` where each
fold entry includes the same fields along with \`train_idx\` and
\`test_idx\`.

## Details

The \`folds\` argument accepts: \`NULL\` (single context), a
\`dkge_folds\` object, a data frame with columns \`subject\` and
\`fold\`, or a list whose elements name the held-out subjects.
Fold-specific results are returned under \`result\$folds\[\[f\]\]\` with
training/test indices attached.

## Examples

``` r
K_list <- list(s1 = diag(5), s2 = diag(4), s3 = diag(5))
effects <- list(
  s1 = c("a", "b", "c", "d", "e"),
  s2 = c("a", "b", "c", "d"),
  s3 = c("b", "c", "d", "e", "f")
)
aligned <- dkge_align_effects(K_list, effects, mode = "intersection")
length(aligned$K_aligned)
#> [1] 3
```
