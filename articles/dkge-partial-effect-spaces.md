# Unbalanced and Partial Effect Spaces

DKGE can represent one global effect space when subjects contribute
different cells or estimate the same cells with very different trial
counts. The output is still a common q-dimensional group basis, but
three distinct problems must be handled separately:

| Problem | What changes | DKGE mechanism |
|----|----|----|
| A cell is absent for a subject | Whether an effect pair is observed | `observed_rows` and `missingness` |
| A cell has 2 trials for one subject and 20 for another | Precision of an observed effect | [`dkge_effect_weights()`](https://bbuchsbaum.github.io/dkge/reference/dkge_effect_weights.md) |
| Cell estimates contain finite-trial estimation noise | Expected diagonal noise in the second moment | `debias` |

Coverage, precision weighting, and debiasing are complementary. A count
weight does not remove estimation-noise bias; analytic subtraction does
not decide which subject should contribute more; and zero-filling alone
does not record that a cell was unobserved.

The key contract is:

- [`dkge_data()`](https://bbuchsbaum.github.io/dkge/reference/dkge_data.md)
  aligns subject-local effect rows to the union of all effect labels.
- Missing rows are rendered as zero rows in the aligned matrices.
- The original observation pattern is retained as `observed_rows`,
  `provenance$obs_mask`, and `provenance$pair_counts`.
- `dkge_fit(missingness = ...)` controls how the partial coverage is
  handled when the q-space covariance is accumulated.

## How do you declare a partial global grid?

Use
[`dkge_effect_grid()`](https://bbuchsbaum.github.io/dkge/reference/dkge_effect_grid.md)
to declare the global effect cells and which factors are within- versus
between-subject.

``` r

grid <- dkge_effect_grid(
  factors = list(
    group = c("control", "patient"),
    task = c("A", "B"),
    measure = c("low", "high")
  ),
  scope = c(group = "between", task = "within", measure = "within"),
  block_factors = "group"
)

grid$cell_labels
#> [1] "control:A:low"  "control:A:high" "control:B:low"  "control:B:high"
#> [5] "patient:A:low"  "patient:A:high" "patient:B:low"  "patient:B:high"
```

`block_factors = "group"` requests independent group blocks in the
design kernel. Terms such as `task` and `measure` are replicated within
each group block rather than coupled across groups by an all-ones group
factor.

``` r

kernel <- design_kernel(
  grid,
  terms = list(
    "group", "task", "measure",
    c("group", "task"),
    c("group", "measure"),
    c("task", "measure"),
    c("group", "task", "measure")
  ),
  basis = "cell",
  normalize = "none"
)

kernel$info$term_scope
#>              group               task            measure         group:task 
#>          "between"           "within"           "within"            "mixed" 
#>      group:measure       task:measure group:task:measure 
#>            "mixed"           "within"            "mixed"
kernel$info$block_factors
#> [1] "group"
```

The kernel metadata classifies terms as `within`, `between`, or `mixed`.
Contrast helpers use this metadata to recommend the appropriate
inference route.

## How are subject-local rows aligned?

Each subject supplies only the four rows for their own group.
[`dkge_data()`](https://bbuchsbaum.github.io/dkge/reference/dkge_data.md)
embeds those rows into the eight-row union and records which global rows
were actually observed.

``` r

cell4 <- c("A:low", "A:high", "B:low", "B:high")
subject_info <- data.frame(
  subject_id = paste0("s", 1:4),
  group = c("control", "control", "patient", "patient")
)

make_subject_beta <- function(group) {
  rows <- paste(group, sub(":.*$", "", cell4), sub("^.*:", "", cell4), sep = ":")
  B <- matrix(rnorm(4 * 12), nrow = 4,
              dimnames = list(rows, paste0("feature", 1:12)))
  X <- diag(4)
  colnames(X) <- rows
  list(B = B, X = X)
}

subjects <- lapply(subject_info$group, make_subject_beta)
betas <- lapply(subjects, `[[`, "B")
designs <- lapply(subjects, `[[`, "X")

dat <- dkge_data(betas, designs, subject_ids = subject_info$subject_id)
dat$effects
#> [1] "control:A:low"  "control:A:high" "control:B:low"  "control:B:high"
#> [5] "patient:A:low"  "patient:A:high" "patient:B:low"  "patient:B:high"
dat$observed_rows
#> [[1]]
#> [1] 1 2 3 4
#> 
#> [[2]]
#> [1] 1 2 3 4
#> 
#> [[3]]
#> [1] 5 6 7 8
#> 
#> [[4]]
#> [1] 5 6 7 8
dat$provenance$pair_counts
#>                control:A:low control:A:high control:B:low control:B:high
#> control:A:low              2              2             2              2
#> control:A:high             2              2             2              2
#> control:B:low              2              2             2              2
#> control:B:high             2              2             2              2
#> patient:A:low              0              0             0              0
#> patient:A:high             0              0             0              0
#> patient:B:low              0              0             0              0
#> patient:B:high             0              0             0              0
#>                patient:A:low patient:A:high patient:B:low patient:B:high
#> control:A:low              0              0             0              0
#> control:A:high             0              0             0              0
#> control:B:low              0              0             0              0
#> control:B:high             0              0             0              0
#> patient:A:low              2              2             2              2
#> patient:A:high             2              2             2              2
#> patient:B:low              2              2             2              2
#> patient:B:high             2              2             2              2
```

Controls observe rows 1-4 and patients observe rows 5-8. Cross-group row
pairs have zero pair counts.

## How should partial coverage enter the fit?

The default `missingness = "none"` preserves the historical zero-filled
accumulation. That compatibility setting is appropriate when every
subject observes every effect. For a genuinely partial global space,
choose a policy that uses the recorded coverage:

- `"mask"` zeros entries whose applicable coverage measure is below
  `min_pairs`: observed pair mass without effect weights, and Kish pair
  ESS when effect-precision weights are active.
- `"rescale"` returns a per-pair mean: it divides by observed pair mass
  without effect weights and by total pair precision when those weights
  are active.
- `"shrink"` blends that per-pair mean toward its diagonal according to
  pair coverage (pair ESS on the precision-weighted branch).

``` r

fit <- dkge_fit(
  dat,
  K = kernel,
  rank = 2,
  w_method = "none",
  missingness = "mask",
  miss_args = list(min_pairs = 1)
)

fit$missingness
#> [1] "mask"
fit$pair_counts[1:4, 5:8]
#>                patient:A:low patient:A:high patient:B:low patient:B:high
#> control:A:low              0              0             0              0
#> control:A:high             0              0             0              0
#> control:B:low              0              0             0              0
#> control:B:high             0              0             0              0
```

Here too, `w_method = "none"` isolates the coverage policy. With the
default MFA subject scaling, zero filling still occurs first, but the
later transformed block energies also affect each subject’s scalar
contribution.

### What the zero-filled rows do, and do not, do

Coverage policies act in *raw* effect space, before the pooled ruler `R`
and the design kernel mix rows. A subject’s raw second moment `B_s B_s'`
therefore has exact zero rows and columns wherever that subject observed
nothing: placeholders contribute no energy of their own. Zero filling
also precedes subject-weight derivation, but a non-`"none"` `w_method`
scores the transformed block after the pooled ruler `R` and `Khalf` have
mixed effect coordinates. Subject weights therefore do not evaluate each
raw zero row in isolation.

They are not, however, insulated from the metric. The fit embeds each
moment as `K^{1/2} R' (B_s B_s') R K^{1/2}`, and if `K` couples an
observed cell to an unobserved one, `K^{1/2}` will place some of the
observed energy on the unobserved coordinate. That is the K-metric doing
its job — it is the same smoothing that makes an ordinal or circular
kernel useful — not a leak of values the subject never supplied.

If you want strict separation, say so in the kernel rather than in the
accumulator. `block_factors = "group"`, as used above, makes `K`
block-diagonal across groups, so `K^{1/2}` is block-diagonal too and no
control-group energy can reach a patient-group coordinate:

``` r

Khalf <- kernel_roots(kernel$K)$Khalf
# Cross-group blocks of K^{1/2} are numerically zero.
max(abs(Khalf[1:4, 5:8]))
#> [1] 0
```

## How do you fit an unbalanced 3 x 5 x 4 trialwise design?

Now consider a fully within-subject `condition x delay x response`
design with 3, 5, and 4 levels: 60 possible cells per subject. The
response is ordinal, and trial counts vary by subject and cell.
[`dkge_effect_grid()`](https://bbuchsbaum.github.io/dkge/reference/dkge_effect_grid.md)
pins the global row order while
[`design_kernel()`](https://bbuchsbaum.github.io/dkge/reference/design_kernel.md)
says that adjacent response levels are more similar than distant ones.

``` r

grid60 <- dkge_effect_grid(
  factors = list(
    condition = c("c1", "c2", "c3"),
    delay = paste0("d", 1:5),
    response = list(L = 4, type = "ordinal",
                    levels = as.character(1:4), l = 1)
  )
)

kernel60 <- design_kernel(
  grid60,
  terms = list(
    "condition", "delay", "response",
    c("condition", "delay"), c("condition", "response"),
    c("delay", "response"), c("condition", "delay", "response")
  ),
  basis = "cell",
  normalize = "unit_trace"
)

c(q = length(grid60$cell_labels), kernel_rows = nrow(kernel60$K))
#>           q kernel_rows 
#>          60          60
```

The trialwise constructor fits $`Y_s = X_s B_s + E_s`$ for each subject.
It retains $`B_s`$, $`(X_s^\top X_s)^{-1}`$, residual variances, and
one-hot cell counts, but not the full trial-by-feature response. Here
every cell has at least three trials so that the within-cell split
leaves both halves estimable.

|     | total_trials | min_cell | median_cell | max_cell |
|:----|-------------:|---------:|------------:|---------:|
| s1  |          389 |        3 |           6 |       10 |
| s2  |          406 |        3 |           7 |       10 |
| s3  |          411 |        3 |           7 |       10 |
| s4  |          392 |        3 |           7 |       10 |
| s5  |          364 |        3 |           5 |       10 |

Each subject contributes a different count profile, but all are aligned
to the same 60 labels. The group fit below makes three choices explicit:

1.  `effect_scaling = "none"` keeps cell means in their common beta
    units rather than applying the pooled design ruler.
2.  `effect_weights = dkge_effect_weights("count")` gives more influence
    to better-estimated subject-by-cell rows.
3.  `debias = "analytic"` subtracts the expected finite-trial noise
    moment before pooling.

`w_method = "none"` is deliberate in this diagnostic example: it
isolates cell-level precision weighting and debiasing from the package’s
default subject-level MFA scaling. It is not a general recommendation to
disable subject weighting.

``` r

fit_analytic <- dkge_fit(
  trial_data,
  K = kernel60,
  rank = 3,
  w_method = "none",
  effect_scaling = "none",
  effect_weights = dkge_effect_weights("count"),
  debias = "analytic",
  missingness = "none"
)

fit_analytic$rank
#> [1] 3
```

For cells $`c`$ and $`c'`$, count weighting uses pair reliability
$`\sqrt{n_{sc} n_{sc'}}`$. DKGE first computes the precision-weighted
mean for each pair, then restores the cohort scale. `fit$pair_ess` is
Kish’s effective number of contributing subjects; it falls when one or
two subjects dominate a cell pair even if every subject observed it.

``` r

pair_diagnostics <- data.frame(
  diagnostic = c("minimum pair ESS", "maximum pair ESS",
                 "negative raw-effect mass", "negative transformed mass"),
  value = c(min(fit_analytic$pair_ess), max(fit_analytic$pair_ess),
            fit_analytic$moment_diagnostics$effect$negative_mass,
            fit_analytic$moment_diagnostics$transformed$negative_mass)
)
knitr::kable(pair_diagnostics, digits = 3)
```

| diagnostic                |   value |
|:--------------------------|--------:|
| minimum pair ESS          |   3.699 |
| maximum pair ESS          |   4.997 |
| negative raw-effect mass  | 145.840 |
| negative transformed mass |   1.364 |

### What exactly does debiasing change?

Finite-trial noise can be addressed in either of two ways:

- `debias = "analytic"` subtracts `noise_trace * (X'X)^{-1}` per
  subject. The noise trace includes residual variance and any diagonal
  spatial weights.
- `debias = "split_half"` replaces the raw second moment with the
  symmetrized cross-product of two stored half estimates. Independent
  half-errors then have zero expected cross-product.

The chunked constructor uses the same weighted analytic correction as
the dense constructor. In particular, a non-unit spatial `omega` is
included when the noise trace is reconstructed from per-feature residual
variances; an unweighted cached trace is not reused as if it were
already weighted.

``` r

set.seed(19411)
X_chunk <- model.matrix(~ 0 + factor(rep(1:2, each = 6)))
colnames(X_chunk) <- c("e1", "e2")
omega_chunk <- c(0.2, 1, 3, 0.5)
make_chunk_y <- function() {
  truth <- matrix(c(1, -0.5, 0.25, 2, 0.4, -1, 0.7, 0.1), 2, 4)
  Y <- X_chunk %*% truth +
    matrix(rnorm(nrow(X_chunk) * 4, sd = 0.35), nrow(X_chunk), 4)
  colnames(Y) <- paste0("v", 1:4)
  Y
}
Y_chunk <- list(make_chunk_y(), make_chunk_y())
dense_subjects <- lapply(seq_along(Y_chunk), function(s) {
  dkge_trial_subject(Y_chunk[[s]], X_chunk, id = paste0("s", s),
                     omega = omega_chunk)
})
chunked_subjects <- lapply(seq_along(Y_chunk), function(s) {
  Y <- Y_chunk[[s]]
  dkge_trial_subject_chunks(
    list(Y[, 1:2, drop = FALSE], Y[, 3:4, drop = FALSE]),
    X_chunk,
    id = paste0("s", s),
    omega = omega_chunk
  )
})
K_chunk <- diag(2)
dimnames(K_chunk) <- list(colnames(X_chunk), colnames(X_chunk))
fit_dense <- dkge_fit(
  dkge_data(dense_subjects), K = K_chunk, rank = 1,
  w_method = "none", effect_scaling = "none", debias = "analytic"
)
fit_chunked <- dkge_fit(
  dkge_data(chunked_subjects), K = K_chunk, rank = 1,
  w_method = "none", effect_scaling = "none", debias = "analytic"
)
chunked_check <- data.frame(
  weighted_noise_trace = sum(
    omega_chunk * dense_subjects[[1]]$residual_variance
  ),
  max_abs_Chat_difference = max(abs(fit_dense$Chat - fit_chunked$Chat))
)
knitr::kable(chunked_check, digits = 12)
```

| weighted_noise_trace | max_abs_Chat_difference |
|---------------------:|------------------------:|
|            0.5501382 |                       0 |

Because the subjects above stored within-cell halves, the alternative
fit is runnable with the same data:

``` r

fit_split <- dkge_fit(
  trial_data,
  K = kernel60,
  rank = 3,
  w_method = "none",
  effect_scaling = "none",
  effect_weights = dkge_effect_weights("count"),
  debias = "split_half"
)

c(analytic = fit_analytic$moment_diagnostics$effect$negative_mass,
  split_half = fit_split$moment_diagnostics$effect$negative_mass)
#>   analytic split_half 
#>   145.8396   193.1882
```

Here the split-half estimate has more negative spectral mass (about 193
versus 146), but that ordering is not a performance score. Both numbers
diagnose the finite-sample indefiniteness of their respective moment
estimators; choosing between them depends on whether the split errors
are credibly independent and whether the analytic covariance model is
credible.

Analytic subtraction and pair normalization can produce an indefinite
q-by-q estimate. DKGE therefore uses a symmetric eigendecomposition,
retains the leading positive eigenpairs, and exposes negative spectral
mass through `fit$moment_diagnostics`; an SVD would incorrectly turn
negative directions into positive components.

The constructor’s `split = "within_cell"` alternates trials within each
cell; it does not prove that the two half-errors are independent. If
runs, sessions, or temporal dependence define independence in your
experiment, use an appropriately constructed split outside this
convenience path or prefer the analytic estimator with a justified
covariance model.

## Which contrasts are estimable within subject?

Contrasts named after kernel terms are tagged with the term scope.

``` r

task4 <- c(-0.5, -0.5, 0.5, 0.5)
contrasts <- list(
  task = rep(task4, 2),
  group = c(rep(-0.25, 4), rep(0.25, 4)),
  "group:task" = c(-task4, task4)
)

task_res <- dkge_contrast(fit, contrasts["task"], method = "loso", align = FALSE)
knitr::kable(task_res$metadata$contrast_estimability)
```

| contrast | estimability | recommended_inference     |
|:---------|:-------------|:--------------------------|
| task     | within       | LOSO/k-fold cross-fitting |

`task` is a within-subject contrast, so LOSO/k-fold contrast inference
is the right family. `group` and `group:task` are between or mixed
effects. DKGE will still return descriptive cross-fitted maps, but the
confirmatory route is the between-subject RRR layer: construct a subject
model with
[`dkge_subject_model()`](https://bbuchsbaum.github.io/dkge/reference/dkge_subject_model.md),
then call
[`dkge_between_rrr()`](https://bbuchsbaum.github.io/dkge/reference/dkge_between_rrr.md)
and assess it with
[`dkge_between_permute()`](https://bbuchsbaum.github.io/dkge/reference/dkge_between_permute.md).

``` r

group_res <- dkge_contrast(
  fit, contrasts[c("group", "group:task")],
  method = "loso", align = FALSE
)
#> Warning: Contrast(s) 'group', 'group:task' are between/mixed effects; loso
#> results are descriptive. Use subject-label permutation or dkge_between_*
#> inference for group-effect testing.
knitr::kable(group_res$metadata$contrast_estimability)
```

| contrast | estimability | recommended_inference |
|:---|:---|:---|
| group | between | subject-label permutation or dkge_between\_\* inference |
| group:task | mixed | subject-label permutation or dkge_between\_\* inference |

## What is deliberately outside this workflow?

A supervised tensor or CP/Tucker model could be useful for a research
extension, especially when the goal is to discover a latent
group-conditioned expression pattern. That is not part of the DKGE core
contract here. The core feature is partial global effect spaces with
explicit coverage, q-space missingness policies, and term-scope-aware
inference recommendations.

For group, trait, or mixed-effect tests on the resulting subject
representations, continue with
[`vignette("dkge-between-subjects")`](https://bbuchsbaum.github.io/dkge/articles/dkge-between-subjects.md).
For a conceptual map of component-, contrast-, and feature-level claims,
see
[`vignette("dkge-concepts")`](https://bbuchsbaum.github.io/dkge/articles/dkge-concepts.md).
