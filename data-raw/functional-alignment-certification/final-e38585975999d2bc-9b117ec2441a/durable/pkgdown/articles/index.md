# Articles

### Start here

- [Getting Started with
  DKGE](https://bbuchsbaum.github.io/dkge/articles/dkge.md):

  Fit a first DKGE model from per-subject GLM betas, read its
  components, and find the page that answers your next question.

- [DKGE Workflow: From Subject Betas to a Shared
  Map](https://bbuchsbaum.github.io/dkge/articles/dkge-workflow.md):

  Take subject beta matrices to a shared map: fit, read components,
  obtain held-out contrast fields, and transport them to one
  parcellation.

- [DKGE Concepts: What Is Estimated and What Can Be
  Claimed?](https://bbuchsbaum.github.io/dkge/articles/dkge-concepts.md):

  Separate what DKGE estimates from what you may claim about it, and
  decide the estimand before fitting.

### Core analysis

- [Design Kernels and Model
  Tuning](https://bbuchsbaum.github.io/dkge/articles/dkge-design-kernels.md):

  Choose a design kernel that states a scientific belief, and compare
  candidates on held-out subjects rather than on fit.

- [Contrasts and
  Inference](https://bbuchsbaum.github.io/dkge/articles/dkge-contrasts-inference.md):

  Build a prespecified contrast, evaluate it out of sample with LOSO,
  and choose between the analytic and bootstrap routes to inference.

- [Components and
  Interpretability](https://bbuchsbaum.github.io/dkge/articles/dkge-components.md):

  Read a fitted component in effect space and on shared clusters, and
  rotate or project it without treating description as inference.

- [Classification with
  DKGE](https://bbuchsbaum.github.io/dkge/articles/dkge-classification.md):

  Cross-validate a latent classifier on DKGE components, choose between
  cell and cell_cross modes, and read the decoding result honestly.

### Study design extensions

- [Partial Effect
  Spaces](https://bbuchsbaum.github.io/dkge/articles/dkge-partial-effect-spaces.md):

  Declare a partial global effect grid, align subject-local rows to it,
  and choose a coverage policy before any fit.

- [Unbalanced Trialwise
  Designs](https://bbuchsbaum.github.io/dkge/articles/dkge-unbalanced-trialwise.md):

  Fit an unbalanced 3 x 5 x 4 trialwise design with count weighting,
  analytic debiasing, and a split-half alternative.

- [Between-Subject Multivariate
  Modeling](https://bbuchsbaum.github.io/dkge/articles/dkge-between-subjects.md):

  Ask whether group, trait, or covariate effects are distributed across
  a brain-derived feature vector per subject.

### Weighting

- [Weighting Strategies in
  DKGE](https://bbuchsbaum.github.io/dkge/articles/dkge-weighting.md):

  Separate location metrics, model-level spatial regularization, subject
  weights, and transport masses, and see which object each one changes.

- [Adaptive Voxel Weighting in
  DKGE](https://bbuchsbaum.github.io/dkge/articles/dkge-adaptive-weighting.md):

  Define, inspect, and structurally check a training-fold voxel
  weighting rule with
  [`dkge_weights()`](https://bbuchsbaum.github.io/dkge/reference/dkge_weights.md).

### Spatial mapping

- [Model-Level Spatial Regularization in
  DKGE](https://bbuchsbaum.github.io/dkge/articles/dkge-spatial-regularization.md):

  Build a graph Laplacian with adjoin, let smoothing change the fitted
  DKGE solution, tune its strength on raw held-out fields, and carry it
  into inference.

- [Transport and Dense
  Rendering](https://bbuchsbaum.github.io/dkge/articles/dkge-dense-rendering.md):

  Transport subject-specific cluster values onto a shared anchor space
  and render them as a dense field.

- [Functional alignment across
  subjects](https://bbuchsbaum.github.io/dkge/articles/dkge-functional-alignment.md):

- [Feature-Anchored
  DKGE](https://bbuchsbaum.github.io/dkge/articles/dkge-anchors.md):

  Align subjects whose items differ by using feature anchors instead of
  a shared discrete effect grid.

- [Mapper Customization and
  Performance](https://bbuchsbaum.github.io/dkge/articles/dkge-performance.md):

  Choose and cache a spatial mapper, or add your own, when transport
  dominates runtime.

### Extras

- [Plotting DKGE
  Fits](https://bbuchsbaum.github.io/dkge/articles/dkge-plotting.md):

  Read a fitted DKGE model through its scree, salience, contribution,
  and stability panels, and know what each one does not establish.

- [CPCA
  Filtering](https://bbuchsbaum.github.io/dkge/articles/dkge-cpca.md):

  Split the compressed covariance into design-aligned and residual parts
  before eigendecomposition.

- [DKGE versus Partial Least
  Squares](https://bbuchsbaum.github.io/dkge/articles/dkge-vs-pls.md):

  Compare DKGE with PLS on estimand, input contract, and inference route
  rather than on accuracy.
