# Build a trialwise DKGE subject from feature chunks

Reduces successive \`T x P_block\` trialwise matrices without retaining
a full \`T x P\` response matrix. The shared trial design and error
model are applied to every block, and the resulting q-space coefficients
and per-feature residual statistics are concatenated into one
\[dkge_subject()\] record. The returned object retains one \`T x q\`
design and one \`q x P\` beta matrix.

## Usage

``` r
dkge_trial_subject_chunks(
  chunks,
  design,
  id = NULL,
  omega = NULL,
  effect_precision = NULL,
  trial_covariance = NULL,
  trial_precision = NULL,
  effect_noise_cov = NULL,
  split = c("none", "within_cell", "alternate", "run", "explicit"),
  split_labels = NULL,
  run_labels = NULL,
  split_independent = FALSE,
  tol = 1e-08,
  max_chunks = 100000L
)
```

## Arguments

- chunks:

  Either a list of numeric \`T x P_block\` matrices or a function called
  as \`chunks(i)\` that returns the next matrix and returns \`NULL\`
  after the last block. A function source is the memory-bounded path.

- design:

  Numeric \`T x q\` second-stage design matrix. This may be a one-hot
  cell design or a general full-rank basis design.

- id:

  Optional subject identifier.

- omega:

  Optional final voxel/parcel weighting of length total \`P\`, or a
  total-\`P\` square matrix. It is applied only after all chunks are
  reduced.

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

- max_chunks:

  Safety limit for function sources that fail to terminate.

## Value

A \`dkge_subject\` with the sufficient statistics needed by fitting. The
full trialwise response and the reconstructible q-by-P effect score
(\`effect_information

## Examples

``` r
set.seed(1)
X <- model.matrix(~ 0 + factor(rep(1:3, each = 4)))
blocks <- list(matrix(rnorm(12 * 2), 12, 2),
               matrix(rnorm(12 * 3), 12, 3))
subject <- dkge_trial_subject_chunks(blocks, X)
```
