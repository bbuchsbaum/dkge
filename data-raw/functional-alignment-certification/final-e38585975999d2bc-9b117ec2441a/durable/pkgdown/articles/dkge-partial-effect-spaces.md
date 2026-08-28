# Partial Effect Spaces

DKGE can represent one global effect space when subjects contribute
different cells or estimate the same cells with very different trial
counts. The output is still a common q-dimensional group basis, but
three distinct problems must be handled separately:

This page is for the case where the usual rectangular data picture is
false. Perhaps controls and patients occupy different rows of a global
grid, or every subject attempted the same factorial design but some
cells have many more trials than others. If every subject has every
effect with comparable precision, you do not need these controls; use
[`vignette("dkge-workflow")`](https://bbuchsbaum.github.io/dkge/articles/dkge-workflow.md).

The guiding rule is simple: **absence, imprecision, and estimation noise
are not synonyms**. The rest of the vignette shows where each enters the
fit.

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
- The original observation pattern is retained on the returned bundle as
  `observed_rows`, plus `obs_mask` and `pair_counts` under `provenance`.
- `dkge_fit(missingness = ...)` controls how the partial coverage is
  handled when the q-space covariance is accumulated.

``` text
local observed rows -> global labelled grid -> coverage-aware raw moment -> precision/debiasing -> kernel transform
```

The first example isolates missing cells. The larger trialwise example
then adds unequal precision and debiasing. Keeping those examples
separate prevents one mechanism from appearing to solve all three
problems.

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

bundle_partial <- dkge_data(betas, designs, subject_ids = subject_info$subject_id)
bundle_partial$effects
#> [1] "control:A:low"  "control:A:high" "control:B:low"  "control:B:high"
#> [5] "patient:A:low"  "patient:A:high" "patient:B:low"  "patient:B:high"
bundle_partial$observed_rows
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
bundle_partial$provenance$pair_counts
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

![Binary heatmap showing that control subjects observe the first four
global effect cells and patient subjects observe the last
four.](dkge-partial-effect-spaces_files/figure-html/coverage-map-1.png)

The empty half of each row is structural absence, not a measured beta of
zero. The observation mask preserves that distinction after the aligned
matrices are expanded to eight rows.

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
fit_partial <- dkge_fit(
  bundle_partial,
  K = kernel,
  rank = 2,
  w_method = "none",
  missingness = "mask",
  miss_args = list(min_pairs = 1)
)

fit_partial$missingness
#> [1] "mask"
fit_partial$pair_counts[1:4, 5:8]
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

## Which contrasts are estimable within subject?

Contrasts named after kernel terms are tagged with the term scope.

``` r
task4 <- c(-0.5, -0.5, 0.5, 0.5)
contrasts <- list(
  task = rep(task4, 2),
  group = c(rep(-0.25, 4), rep(0.25, 4)),
  "group:task" = c(-task4, task4)
)

task_res <- dkge_contrast(fit_partial, contrasts["task"], method = "loso", align = FALSE)
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
  fit_partial, contrasts[c("group", "group:task")],
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

## Where to go next

- [`vignette("dkge-unbalanced-trialwise")`](https://bbuchsbaum.github.io/dkge/articles/dkge-unbalanced-trialwise.md)
  fits a 3 x 5 x 4 trialwise design on the contract this page defines.
- [`vignette("dkge-between-subjects")`](https://bbuchsbaum.github.io/dkge/articles/dkge-between-subjects.md)
  is the confirmatory route for the between and mixed contrasts named
  above.
- [`vignette("dkge-concepts")`](https://bbuchsbaum.github.io/dkge/articles/dkge-concepts.md)
  maps component-, contrast-, and feature-level claims onto what each
  analysis supports.
