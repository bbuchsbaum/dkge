# DKGE workflow pain points (SDAM testbed)

Notes captured while porting the plsrri SDAM example to DKGE
(`inst/examples/sdam_firstlevel_task_dkge.R`) and rendering the Quarto
report on the full mask. Each item is something a user had to write
themselves but probably shouldn't have, and each one suggests a concrete API
tweak.

These are mostly **not** ergonomic nice-to-haves — they expose missing
*boundary objects* in the workflow: inputs from maps/NIfTIs, outputs back to
voxel space, and tidy second-level inference results.

## Priority ordering

The items are grouped by priority, not chronological discovery. Numbered IDs
within a section don't imply within-priority order.

### P1 — Workflow spine: inputs and outputs

These four turn the SDAM script from a bespoke pipeline into an idiomatic
DKGE workflow.

#### P1.1 `dkge_subject_from_maps(maps, cell_labels, id = NULL, omega = NULL)`

For post-MVPA / RSA-style data each subject contributes a fixed number of
maps, not a time-series GLM. Today the user has to:

* stack the maps into a `q x P` matrix in canonical cell order
* manufacture `X_s = I_q` for each subject as a stand-in design
* hand-build `dkge_data(B_list, designs = make_design_list(B_list), ...)`

A native constructor that accepts a named list (or `q x P` matrix) of cell
maps and returns a `dkge_subject` with the implicit identity design would
eliminate ~30 lines of `make_sdam_design_list()`-style boilerplate. See
`R/dkge-data.R:120`.

#### P1.2 `dkge_subject_from_nifti(files, cell_labels, mask, id = NULL, omega = NULL)`

`dkge_subject.NeuroVec` exists for *4D timeseries* and the
`ClusteredNeuroVec` overload assumes a fmrireg GLM. Neither covers the most
common practical case for downstream applications: a list of `NeuroVol`s
(one per cell), already in a common space, masked to a shared mask, that
need to become rows of a `q x P` matrix.

A helper that

1. reads each NIfTI via `neuroim2::read_vol()`
2. validates that all share `dim(mask)`
3. returns a `dkge_subject` with `beta = stacked rows`

combined with P1.1, shrinks the SDAM example from ~80 lines of plumbing to
roughly ten.

#### P1.3 `dkge_values_to_neurovol(values, mask)`

Round-tripping a voxel vector back into a `NeuroVol` with a mask is the same
five-line block every script writes (`sdam_values_to_neurovol()` here).
Export it; it's the natural pair to P1.2.

#### P1.4 "From NIfTIs" end-to-end vignette

`vignette("dkge-workflow")` jumps straight to `dkge_subject(matrix, ...)`
with simulated betas. A short companion vignette walking through
NIfTIs -> mask -> `B_list` -> kernel -> fit -> contrasts -> between-subject
test -> brain map closes the documentation gap.
`inst/examples/sdam_firstlevel_task_dkge.R` is the natural seed.

### P2 — Compact statistical summaries

These two ship together. Both call into existing DKGE machinery and just
hide the per-script glue users currently write by hand.

#### P2.1 `dkge_contrast_voxel_summary(contrasts_obj, ...)`

After `dkge_contrast(fit, C, method="loso")` the per-voxel parametric
reduction (mean / SE / one-sample z across subjects) is the same four lines
every time:

```r
M  <- do.call(rbind, contrasts_obj$values$task)   # subjects x voxels
mu <- colMeans(M)
se <- apply(M, 2, sd) / sqrt(nrow(M))
z  <- mu / se
```

Wrap it as a tidy parametric contrast summary:

```r
dkge_contrast_voxel_summary(contrasts_obj,
                            terms = NULL,
                            statistic = c("z", "mean", "se"),
                            group = NULL,    # optional factor for diff maps
                            data  = NULL)    # data with subject/group cols
```

Returns a long data frame (or named list of voxel vectors) ready for
masking into a `NeuroVol` via P1.3. This is the analytic-summary helper;
the bootstrap-ratio analogue is P2.2.

#### P2.2 `dkge_bsr(fit, contrast, B = 1000L, method = c("loso", "qspace"), ...)`

This is the actual PLS-`bsr()` analogue: a single call that returns a
per-voxel bootstrap z (mean / SD across resamples) for a given contrast.
Today it's `dkge_loso_contrast()` per subject, then
`dkge_bootstrap_projected()`, then `mean / sd` from the result. Wrap once.

#### P2.3 `$medoid` is the wrong name when no transport is involved

`dkge_bootstrap_projected()` puts its summary under `$medoid` even when the
caller is just averaging in-place subject vectors. Add `$projected` as a
backward-compatible alias (don't rename); call sites stay valid, new code
reads less misleadingly. See `R/dkge-bootstrap.R:96`.

### P3 — Tidy between-subject outputs

#### P3.1 `as.data.frame.dkge_between_permutation(x, what = c("summary", "features"))`

`between$perm$feature_tests[["group"]]` is a heterogeneous list:
`term_map`, `statistic`, `p`, `p_adjusted`, `null_max`, `feature_ids`.
Passing it to `kable()` produced a mangled multi-block table in the first
report draft; users have to know to assemble
`data.frame(feature, statistic, p, p_adjusted)` themselves.

The right surface is on the parent permutation object, not on each per-term
list:

```r
as.data.frame(perm, what = "summary")    # one row per term
as.data.frame(perm, what = "features")   # one row per term x feature
```

That gives a single tidy entry point so `kable(as.data.frame(perm, what =
"features"))` and `tibble::as_tibble(as.data.frame(perm, ...))` work
without per-script glue. `summary()`, `print()`, and a direct `kable(perm)`
would still need their own dispatched methods if we want that to work
unwrapped — worth doing as a follow-up but out of scope of the
data-frame coercion itself. See `R/dkge-between-infer.R:122` and
`R/dkge-between-infer.R:266`.

### P4 — Group needs to enter earlier

Yes, the current SDAM report treats group as a late-stage afterthought.
`dkge_make_target(type = "component_scores")` explicitly reduces each
subject to component means, so spatial group information is discarded
*before* RRR — today the report caches a 32 x 3 target at full mask,
throwing away ~218k voxels of spatial information for the group test.

The intended fix is a voxel-level `transported_maps` target plus a
between-subject RRR with `~ group`. The initial report needed workarounds
for two API gaps; both are now part of the core cleanup plan and the
example uses the direct path.

**API gap 1 — silent contrast collapse.** When `centroids = NULL` (the
common case for already-aligned MNI data), `dkge_make_target(type =
"transported_maps", contrast = C)` used to ignore all but the first
contrast and return only `contrast_obj$values[[1]]`. The fix is to stack
all contrast blocks as `subject_ids x (voxels x contrasts)` feature blocks
in the no-transport path, matching the transported path.

**API gap 2 — subject id binding.** `dkge_subject_model(~ group, data)`
defaulted to `data$subject_id` for row matching, but the participants table
in this dataset has `subject` (not `subject_id`). The fix is an explicit
`subject_id_col` argument.

With those fixes, the headline analysis is:

```r
target <- dkge_make_target(fit, type = "transported_maps",
                           contrast = sdam_contrast_matrix(),
                           crossfit = "loso")
design <- dkge_subject_model(~ group, participants,
                             subject_id_col = "subject")
fit_b  <- dkge_between_rrr(target, design, rank = 1)
perm   <- dkge_between_permute(fit_b, terms = "group",
                               B = 999, scope = "both",
                               feature_adjust = "maxT")
```

That is what should be the SDAM example's headline group test, not the
component-score reduction. See `R/dkge-target.R:67`.

#### P4b — Stratified bases for *latent-space* group differences

The `transported_maps` path above catches the most common case of
"shared cell-direction structure, group-specific spatial loadings." When
the *latent design space itself* may differ by group, the diagnostic is
stratified DKGE plus K-Procrustes:

```r
fit_ctl  <- fit_sdam_dkge(B_list[ctl_ids],  kernel)
fit_sdam <- fit_sdam_dkge(B_list[sdam_ids], kernel)
ang <- dkge_principal_angles_K(fit_ctl$U, fit_sdam$U, fit$K)
```

Interpret the angles against a stability baseline, not in isolation:
bootstrap or split-half each group and compare between-group angles against
within-group resampling angles. Small between-group angles relative to
within-group instability -> the pooled fit is defensible. Large
between-group angles that exceed within-group instability -> the latent
space may itself depend on group, and the group-aware encoding (P4c) becomes
the appropriate next step.

#### P4c — Group as an explicit kernel effect (q -> 2q)

If P4b shows large angles, the right answer is to make group a first-class
design factor in the kernel. The first-class implementation is not an
SDAM branch: it is a partial global effect space. `dkge_effect_grid()`
declares `group x task x measure`, marks `group` as between-subject and
`task`/`measure` as within-subject, and `dkge_data()` carries
`observed_rows`, `obs_mask`, and `pair_counts` after union alignment.

The conservative SDAM path uses `block_factors = "group"` so that
`K_group = I_group` replicates the within-subject task/measure kernel in
independent group blocks. Fit-time `missingness = "mask"` then keeps
zero-filled rendering rows from becoming load-bearing in the training
`Chat`; subject weights and per-subject contributions are also computed on
observed rows. The q = 8 contrasts are named 8-vectors:

```r
C <- sdam_group_contrast_matrix()
colnames(C)
# group, group:task, group:measure, group:task:measure
```

The report now includes an observed-vs-completed validation table for these
between/mixed contrasts. Generic supervised CP/Tucker tensor models remain
a research extension, not a core DKGE dependency.

## Lower-priority items

#### L1 `dkge_between_rrr` ergonomics: don't make alignment load-bearing

`dkge_between_rrr()` already reorders an existing `dkge_subject_model` by
`target$subject_ids`, so the original framing ("doesn't join by id") was
imprecise. The real footgun is constructing the model matrix without
reliable subject row names — users still have to round-trip
`participants$subject_id <- participants$subject` and ensure column order
themselves before the design is well-formed.

A cleaner API would push the join below `dkge_subject_model`:

```r
dkge_between_rrr(target, ~ group, data = participants, subject_id = "subject")
# or
dkge_subject_model(~ group, participants, align_to = target$subject_ids)
```

#### L2 `factorial_kernel()` shortcut

`design_kernel()` already defaults to all main effects plus the full
interaction when `terms = NULL`. A wrapper

```r
factorial_kernel(factors, rho = c(main = 1, interaction = 0.5),
                 basis = "cell")
```

mainly buys naming conventions, cell labels (see L3), and the
`rho = c(main=..., interaction=...)` shorthand. Real but lower-priority
than the input/output helpers.

#### L3 Cell-basis kernels need cell labels

`design_kernel(basis = "cell")` does not populate `info$blocks` or
`dimnames(K)`, so plot helpers fall back to `effect1`, `effect2`, ...
even when cells have natural names. Two fixes:

* Have `design_kernel()` set `dimnames(K)` from `info$cell_labels` (or
  a `cell_labels` argument).
* Allow `dkge_plot_effect_loadings(fit, labels = c(...))` overrides.

#### L4 Standalone PLS-style plots before a `pls_style = TRUE` flag

I had to write `plot_sdam_design_heatmap`, `plot_sdam_design_interaction`,
`plot_sdam_subject_scores`, and a multi-slice `plot_sdam_bsr_volume` —
~80 lines of plotting code that should be a few exported helpers, not one
mega-flag on `dkge_plot_suite()`. Each panel needs extra metadata the
helpers don't currently see (cell table, group labels, mask/volume
geometry, BSR values), so the right design is small focused helpers that
take those metadata arguments explicitly.

#### L5 Export the diverging colour scales

`.dkge_scale_fill_diverging()` / `.dkge_scale_colour_diverging()` are
internal, so user code that follows dkge's idiom redefines the same
gradient. Export them under non-dot names:

```r
dkge_scale_fill_diverging()
dkge_scale_colour_diverging()
```

`theme_dkge()` is already exported; pair these with it.
