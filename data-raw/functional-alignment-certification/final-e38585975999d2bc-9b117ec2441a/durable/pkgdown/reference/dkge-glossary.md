# DKGE glossary

Definitions for the terms the DKGE documentation uses most. Each one
states what the package does, not what the word means elsewhere.

## Effect space

The shared coordinate system in which DKGE works. Each subject's beta
matrix has one row per named design effect, so the rows mean the same
thing across subjects even when the columns (clusters) do not. That row
space is the effect space, and its dimension `q` is the number of design
effects. The design kernel `K` supplies the metric on it.

## Design kernel

A `q` by `q` positive semidefinite matrix stating which design effects
should be treated as similar. It expresses a scientific belief you can
name, such as adjacent ordinal levels being related. It is not a tuning
matrix, and `diag(q)` is the honest choice when the structure is
uncertain. See
[`design_kernel()`](https://bbuchsbaum.github.io/dkge/reference/design_kernel.md).

## Salience

The design-weighted component coordinates, \\K U\\. Rows are effects or
design cells and columns are latent components, so reading down a column
shows which effects a component combines. See
[`dkge_component_saliences()`](https://bbuchsbaum.github.io/dkge/reference/dkge_component_saliences.md).

## Cross-fitting

Evaluating a subject through a basis their own data did not shape. The
leave-one-subject-out path removes subject `s` from the pooled
compressed covariance, re-derives the basis by eigendecomposition of
what remains, and reads that subject's values through it. K-fold
cross-fitting does the same with folds instead of single subjects. It
limits basis-reuse optimism; it does not supply a population p-value.

## Reference support

The coordinates, labels, topology, and optional decoder on which aligned
values are represented. A subject parcellation or fixed anatomical/MNI
grid can supply support. Coordinates alone contain no functional
correspondence and no participant's activation values. See
[`dkge_reference_support()`](https://bbuchsbaum.github.io/dkge/reference/dkge_reference_support.md).

## Reference subject and medoid

A reference subject supplies a real subject parcellation as display
support. It is a *medoid* only when selected by a stated cohort
criterion, such as held-out symmetric functional reconstruction or
geometry-only loss. A caller-fixed subject is recorded as an explicit
reference, not relabelled a medoid. See
[`dkge_select_reference_subject()`](https://bbuchsbaum.github.io/dkge/reference/dkge_select_reference_subject.md).

## Functional correspondence

The fitted operator that maps a subject's parcels to one reference
support. Functional costs must come from typed response-signature
features; geometry may regularize or select support but does not become
functional data. Every subject, including a reference subject, is fitted
through the same mapper policy.

## Functional template

Target response signatures and masses attached to a reference support. A
selected subject can initialize a template without becoming privileged
functional truth.
[`dkge_fit_functional_template()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit_functional_template.md)
instead pools information across subjects with mass-aware updates and
fixed feature scale.

## Aligned subject rows

One contrast row per subject on the same identified support. These rows,
not parcels, are the observations for group inference. Fit-level MFA
weights are pooling weights for learning the DKGE basis and are not
silently reused as contrast inverse-variance weights. Obtain
operator-bound rows from
[`dkge_transport_contrasts_to_reference()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_contrasts_to_reference.md)
or
[`dkge_align_to_template()`](https://bbuchsbaum.github.io/dkge/reference/dkge_align_to_template.md);
the raw
[`dkge_aligned_maps()`](https://bbuchsbaum.github.io/dkge/reference/dkge_aligned_maps.md)
constructor is deliberately descriptive only.

## Group inference

Subject-level inference performed after correspondence has been fixed.
[`dkge_infer_aligned()`](https://bbuchsbaum.github.io/dkge/reference/dkge_infer_aligned.md)
never learns a mapper. It fails closed for ineligible alignments and
requires an explicit override for objects labelled approximate by the
calibration court.

## Rendering

Display or decoding of an already aligned statistic. A
[`dkge_renderer()`](https://bbuchsbaum.github.io/dkge/reference/dkge_renderer.md)
owns no subject correspondence and cannot turn a bare MNI grid into a
functional template.

## Transport

Applying fitted correspondence to move per-cluster values from each
subject's parcellation onto a common support. Sinkhorn, ridge, or OLS
mappers may implement the operator. Transport is separate from both the
shared effect-space metric `K` and downstream rendering.

## LOSO

Leave-one-subject-out. See Cross-fitting.
