# DKGE Architecture Guide

This guide describes the current DKGE implementation for contributors.
The fit lifecycle is already split into internal stages, and public
constructors already exist for pipeline specifications and services. The
important boundary is therefore not “old monolith versus future
refactor,” but which interfaces are public, which stages are internal,
and which service objects the pipeline actually consumes today.

## Where are the public entry points?

| Task | Current public API | Main implementation |
|----|----|----|
| Harmonize subject data | [`dkge_subject()`](https://bbuchsbaum.github.io/dkge/reference/dkge_subject.md), [`dkge_data()`](https://bbuchsbaum.github.io/dkge/reference/dkge_data.md) | `R/dkge-data.R`, `R/dkge-align-data.R` |
| Define effects and kernels | [`dkge_effect_grid()`](https://bbuchsbaum.github.io/dkge/reference/dkge_effect_grid.md), [`design_kernel()`](https://bbuchsbaum.github.io/dkge/reference/design_kernel.md) | `R/dkge-effect-grid.R`, `R/design-kernel.R` |
| Fit a group embedding | [`dkge()`](https://bbuchsbaum.github.io/dkge/reference/dkge.md), [`dkge_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit.md) | `R/dkge-data.R`, `R/dkge-fit.R`, `R/dkge-fit-core.R` |
| Compute held-out contrasts | [`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md) | `R/dkge-contrast.R`, `R/dkge-loso.R`, `R/dkge-kfold.R`, `R/dkge-analytic.R` |
| Align subject maps | [`dkge_prepare_transport()`](https://bbuchsbaum.github.io/dkge/reference/dkge_prepare_transport.md), [`dkge_transport_contrasts_to_medoid()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_contrasts_to_medoid.md) | `R/dkge-transport.R`, `R/dkge-mapper.R` |
| Test aligned maps | [`dkge_signflip_maxT()`](https://bbuchsbaum.github.io/dkge/reference/dkge_signflip_maxT.md), [`dkge_infer()`](https://bbuchsbaum.github.io/dkge/reference/dkge_infer.md) | `R/dkge-inference.R` |
| Run the common workflow | [`dkge_pipeline()`](https://bbuchsbaum.github.io/dkge/reference/dkge_pipeline.md) | `R/dkge-pipeline.R`, `R/dkge-services.R` |
| Score new subjects | [`dkge_predict_subjects()`](https://bbuchsbaum.github.io/dkge/reference/dkge_predict_subjects.md), [`dkge_predict()`](https://bbuchsbaum.github.io/dkge/reference/dkge_predict.md) | `R/dkge-predict.R` |
| Classify latent outcomes | [`dkge_classify()`](https://bbuchsbaum.github.io/dkge/reference/dkge_classify.md), [`dkge_classification_spec()`](https://bbuchsbaum.github.io/dkge/reference/dkge_classification_spec.md) | classification and latent-classifier modules |

These functions are exported in `NAMESPACE`. Functions beginning with
`.dkge_` are implementation details even when they provide useful seams
for unit tests. Vignettes and downstream packages should not call them
with `dkge:::`.

## How does `dkge_fit()` work now?

[`dkge_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit.md)
validates its user-facing options and orchestrates four internal stages
in order:

``` text
.dkge_fit_prepare()
        |
        v
.dkge_fit_accumulate()
        |
        v
.dkge_fit_solve()
        |
        v
.dkge_fit_assemble()
```

All four stages live in `R/dkge-fit-core.R`;
[`dkge_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit.md)
itself lives in `R/dkge-fit.R`.

### 1. Prepare

`.dkge_fit_prepare()` accepts a `dkge_data` bundle or raw beta/design
lists. It validates finite, common-q inputs; aligns kernel labels with
the global effect order; builds the pooled design ruler `R`; computes
row-standardized `Btil`; derives `Khalf` and `Kihalf`; and resolves
spatial and effect-weight specifications.

The output is an internal payload containing harmonized data and
everything needed by later stages. It is not a partially fitted public
object.

### 2. Accumulate

`.dkge_fit_accumulate()` resolves subject weights and effect precision,
builds per-subject raw effect moments, applies analytic or split-half
debiasing when requested, pools observed effect pairs under the selected
missingness policy, and transforms the pooled moment into `Chat`.

The order is a scientific contract: observation masks, precision
weighting, and debiasing operate in raw effect space before the pooled
ruler and kernel mix effect rows. The fitted object retains
`effect_moment`, `pair_counts`, `pair_weight`, `pair_ess`, and
`moment_diagnostics` so this step remains auditable.

The hidden check records the complete-coverage boundary: `none` pools a
sum, whereas `rescale` returns the observed-pair mean; with full
coverage and the shown settings, `shrink` reduces to that same mean.

### 3. Solve

`.dkge_fit_solve()` applies optional CPCA partitioning and ridge
regularization, then dispatches to the pooled eigensolver or
joint-diagonalization solver. It reduces a requested rank when the
transformed moment has fewer positive directions and maps retained
eigenvectors back through `Kihalf` to form `U`.

The JD branch deliberately fails closed when effect weighting,
missingness normalization, or debiasing is requested; those combinations
are not currently implemented by that solver.

### 4. Assemble

`.dkge_fit_assemble()` attaches diagnostics, provenance, kernel
metadata, and the stable fit fields. It also assigns one of two
representation contracts:

- `block_biprojector` for an exact unregularized block factorization.
  These fits can expose physical block loadings and, with
  `keep_X = TRUE`, the concatenated training matrix.
- `qspace_moment` when pair normalization, debiasing, ridge, CPCA, JD,
  or other moment transformations remove an exact subject-by-feature
  factorization. These fits set `v` and `X_concat` to `NULL` and fail if
  `keep_X = TRUE` is requested.

Contributors must preserve this distinction. A valid q-space eigensystem
is not evidence that physical block loadings exist.

## What does `dkge_pipeline()` integrate?

[`dkge_pipeline()`](https://bbuchsbaum.github.io/dkge/reference/dkge_pipeline.md)
currently performs five steps: obtain a fit, compute contrasts,
optionally transport them, optionally classify, and run sign-flip
inference. It returns all five results plus fit diagnostics.

The configuration surface is uneven but explicit:

| Concern | Objects that exist | Pipeline behavior today |
|----|----|----|
| Contrast | [`dkge_contrast_service()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast_service.md) | The pipeline constructs this service internally from `method` and `ridge`; callers cannot pass a contrast-service object as a pipeline argument. |
| Transport | [`dkge_transport_spec()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_spec.md), [`dkge_transport_service()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_service.md) | The pipeline accepts a list, spec, or transport service and runs the typed service. |
| Inference | [`dkge_inference_spec()`](https://bbuchsbaum.github.io/dkge/reference/dkge_inference_spec.md), [`dkge_inference_service()`](https://bbuchsbaum.github.io/dkge/reference/dkge_inference_service.md) | The pipeline accepts a list, spec, or inference service and runs the typed service. |
| Classification | [`dkge_classification_spec()`](https://bbuchsbaum.github.io/dkge/reference/dkge_classification_spec.md) | The pipeline accepts a classification result, spec/list, or direct target and calls [`dkge_classify()`](https://bbuchsbaum.github.io/dkge/reference/dkge_classify.md) directly. There is no classification-service class. |

[`dkge_run_service()`](https://bbuchsbaum.github.io/dkge/reference/dkge_run_service.md)
is the exported generic for manually executing the three service classes
that do exist: contrast, transport, and inference. This makes them
reusable outside
[`dkge_pipeline()`](https://bbuchsbaum.github.io/dkge/reference/dkge_pipeline.md),
but it does not imply that the pipeline has a general service graph or
accepts an arbitrary service sequence.

Likewise,
[`dkge_predict_subjects()`](https://bbuchsbaum.github.io/dkge/reference/dkge_predict_subjects.md)
is already public. It harmonizes matrices, `dkge_subject` objects,
lists, or a `dkge_data` bundle and forwards the named beta list to
[`dkge_predict()`](https://bbuchsbaum.github.io/dkge/reference/dkge_predict.md).
It is current API, not a proposed wrapper.

## What is genuinely unfinished?

Future work should be described narrowly from the current seams:

- allow callers to inject a contrast service instead of having the
  pipeline always construct one;
- decide whether classification needs a service abstraction, then
  implement and test it before documenting one; and
- if a general service graph is introduced, define typed inputs and
  outputs, ordering, cache ownership, and backward compatibility with
  the present pipeline arguments.

Streaming accumulators, alternative solvers, or new backends may be
reasonable experiments, but they are not implied public contracts. Add
them only behind tests that preserve the raw-moment ordering, K-metric
algebra, representation classification, and held-out inference
boundaries.

## How should contributors change the architecture?

1.  Put input reconciliation in the prepare stage, moment-estimand
    changes in accumulation, eigensystem changes in solve, and
    output-contract changes in assembly.
2.  Add focused stage tests and at least one public integration test.
    The fit algebra tests should continue to verify `Chat`,
    K-orthonormality, and the fail-closed `qspace_moment` contract.
3.  Keep services and specs honest about current pipeline support. Do
    not document an object merely because a constructor name has been
    proposed.
4.  Preserve subject as the inference unit and reuse transport operators
    across resamples only when the intended inferential procedure
    conditions on that mapping.
5.  Update user workflow vignettes only for exported behavior; keep
    internal stage details here and in source comments.

For the user-facing analysis path, see
[`vignette("dkge-workflow")`](https://bbuchsbaum.github.io/dkge/articles/dkge-workflow.md).
For the statistical meaning of the fitted moment, see
[`vignette("dkge-concepts")`](https://bbuchsbaum.github.io/dkge/articles/dkge-concepts.md)
and
[`vignette("dkge-partial-effect-spaces")`](https://bbuchsbaum.github.io/dkge/articles/dkge-partial-effect-spaces.md).
