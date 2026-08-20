# Build a DKGE subject from trialwise beta maps

Fits the second-stage model \`Y = X B + E\` without constructing any
voxel-by-voxel covariance. By default the estimator is IID OLS.
Supplying a trial covariance or precision uses GLS and stores \`(X' W
X)^-1\` as the unit-spatial-variance covariance of the effect estimates.
The returned subject retains q-space sufficient statistics and the
shared \`T x q\` design, but never the trialwise \`y\`. Use
\[dkge_trial_subject_chunks()\] when \`T x P\` itself is too large to
materialize.

## Usage

``` r
dkge_trial_subject(
  y,
  design,
  id = NULL,
  omega = NULL,
  effect_precision = NULL,
  trial_covariance = NULL,
  trial_precision = NULL,
  effect_noise_cov = NULL,
  residual_variance = NULL,
  noise_trace = NULL,
  split = c("none", "within_cell", "alternate", "run", "explicit"),
  split_labels = NULL,
  run_labels = NULL,
  split_independent = FALSE,
  tol = 1e-08
)
```

## Arguments

- y:

  Numeric \`T x P\` matrix of trialwise beta maps.

- design:

  Numeric \`T x q\` second-stage design matrix. This may be a one-hot
  cell design or a general full-rank basis design.

- id:

  Optional subject identifier.

- omega:

  Optional voxel/parcel weighting passed to \[dkge_subject()\].

- effect_precision:

  Optional direct q-vector of effect precisions, or \`"split_half"\` to
  derive bounded effect precision from the stored halves with
  \[dkge_split_effect_precision()\].

- trial_covariance:

  Optional symmetric positive-definite \`T x T\` relative covariance of
  trial errors. Mutually exclusive with \`trial_precision\`.

- trial_precision:

  Optional symmetric positive-definite \`T x T\` trial precision, or a
  function mapping a \`T x k\` matrix to its precision-weighted
  counterpart. A function is materialised only in trial space for
  validation; no \`P x P\` matrix is formed.

- effect_noise_cov:

  Optional directly supplied q-by-q covariance multiplier for the
  estimated effects. When supplied it overrides the covariance implied
  by the OLS or GLS information matrix.

- residual_variance:

  Optional externally estimated length-P spatial noise variances. The
  default is the (weighted for GLS) residual sum of squares divided by
  the residual degrees of freedom.

- noise_trace:

  Optional externally supplied unweighted spatial noise trace. By
  default this is \`sum(residual_variance)\` when estimable. A supplied
  value is used by analytic debiasing when no extra spatial weights are
  applied.

- split:

  Optional split-half sufficient statistic. \`"within_cell"\`
  deterministically balances trials within each one-hot cell;
  \`"alternate"\` alternates all trial rows; \`"run"\` assigns whole
  runs to halves; and \`"explicit"\` uses \`split_labels\`. Both
  half-designs must remain full rank.

- split_labels:

  Optional length-T vector containing exactly two explicit half labels.
  Supplying it selects \`split = "explicit"\`.

- run_labels:

  Optional length-T run labels. With \`split = "run"\`, whole runs are
  assigned deterministically to halves. With another split mode, these
  labels are used to audit whether the halves are run-disjoint.

- split_independent:

  Logical; explicitly declare a non-run-disjoint partition independent.
  Without run-disjoint labels or this declaration, split moments are
  recorded as exploratory and a warning is emitted.

- tol:

  Numerical tolerance used for rank and one-hot checks.

## Value

A \`dkge_subject\` carrying second-stage sufficient statistics.

## Details

The analytic noise correction assumes a common trial-dependence matrix
across spatial columns and a separable error covariance: the effect
covariance is \`(X' W X)^-1\`, scaled by each spatial column's residual
variance. Trialwise beta maps from overlapping or temporally correlated
first-level estimates generally violate IID errors; supply
\`trial_covariance\` or \`trial_precision\` in that case. Using the IID
default when trial errors are correlated mis-scales both the analytic
subtraction and any precision derived from it.

\`split = "within_cell"\` and \`split = "alternate"\` only create
deterministic partitions; they do not establish independence. A split
cross-moment is labelled independent only for run-disjoint halves or
when \`split_independent = TRUE\` is explicitly justified by the caller.

## Examples

``` r
set.seed(1)
X <- model.matrix(~ 0 + factor(rep(1:3, each = 4)))
Y <- X %*% matrix(rnorm(3 * 5), 3, 5) + matrix(rnorm(12 * 5), 12, 5)
subject <- dkge_trial_subject(Y, X, id = "s1")
```
