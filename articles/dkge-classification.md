# Classification with DKGE

``` r
library(dkge)
set.seed(12)
```

The DKGE classification helpers allow you to reuse the cross-fitted
group basis to decode subject-level contrasts with statistical rigor.
This vignette demonstrates a complete workflow: data preparation → model
fit → target specification → cross-validated decoding → statistical
inference.

The approach is fully cross-validated: by default
[`dkge_classify()`](https://bbuchsbaum.github.io/dkge/reference/dkge_classify.md)
uses leave-one-subject-out (LOSO) cross-fitting, re-estimating the group
basis without the held-out subject at each fold. This ensures that the
reported accuracies are unbiased estimates of generalisation
performance.

## Setup

We define two experimental factors — condition (A vs B) and time (4 time
points) — and construct a design kernel that encodes their structure:

``` r
factors <- list(cond = list(L = 2), time = list(L = 4))
kern    <- design_kernel(factors, basis = "effect")
q       <- nrow(kern$K)   # number of effect-coded columns
q
#> [1] 7
```

We simulate 10 subjects, each with `q` design effects and 60 voxels. A
detectable signal distinguishing the two conditions is embedded in the
condition contrast:

``` r
n_subjects <- 10L
v <- 60L

make_subject <- function(id) {
  # design: identity matrix over effects (one trial per effect)
  design <- diag(q)
  colnames(design) <- rownames(kern$K)

  # signal: condition contrast (first effect column) drives voxels 1-10
  signal <- matrix(0, nrow = q, ncol = v)
  signal[1, 1:10] <- 0.5    # positive signal for condition contrast

  beta <- signal + matrix(rnorm(q * v, sd = 1.0), nrow = q)
  dkge_subject(beta, design = design, id = paste0("sub", id))
}

subjects <- lapply(seq_len(n_subjects), make_subject)
```

Fit the DKGE model using the design kernel.
[`design_kernel()`](https://bbuchsbaum.github.io/dkge/reference/design_kernel.md)
returns a list with `$K` (the q×q kernel matrix) and `$info` (factor
metadata used by
[`dkge_targets()`](https://bbuchsbaum.github.io/dkge/reference/dkge_targets.md)):

``` r
fit <- dkge(subjects, K = kern, rank = 2)
fit
#> Multiblock Bi-Projector object:
#>   Projection matrix dimensions:  600 x 2 
#>   Block indices:
#>     Block 1: 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60
#>     Block 2: 61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120
#>     Block 3: 121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180
#>     Block 4: 181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240
#>     Block 5: 241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300
#>     Block 6: 301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360
#>     Block 7: 361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408,409,410,411,412,413,414,415,416,417,418,419,420
#>     Block 8: 421,422,423,424,425,426,427,428,429,430,431,432,433,434,435,436,437,438,439,440,441,442,443,444,445,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,471,472,473,474,475,476,477,478,479,480
#>     Block 9: 481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,496,497,498,499,500,501,502,503,504,505,506,507,508,509,510,511,512,513,514,515,516,517,518,519,520,521,522,523,524,525,526,527,528,529,530,531,532,533,534,535,536,537,538,539,540
#>     Block 10: 541,542,543,544,545,546,547,548,549,550,551,552,553,554,555,556,557,558,559,560,561,562,563,564,565,566,567,568,569,570,571,572,573,574,575,576,577,578,579,580,581,582,583,584,585,586,587,588,589,590,591,592,593,594,595,596,597,598,599,600
```

## Classification targets

[`dkge_targets()`](https://bbuchsbaum.github.io/dkge/reference/dkge_targets.md)
maps design formula terms to classification contrasts using the factor
metadata stored in the kernel. Each term becomes a binary (or
multi-class) decoding problem over the group embedding.

``` r
targets <- dkge_targets(fit, ~ cond + time)
length(targets)                      # one target per formula term
#> [1] 2
targets[[1]]$name                    # "cond"
#> [1] "cond"
nrow(targets[[1]]$weight_matrix)     # 2 classes (A vs B)
#> [1] 2
```

The weight matrix transforms subject-level beta vectors into
class-specific pattern vectors for the downstream classifier.

## Cross-validated decoding

``` r
cls <- dkge_classify(
  fit,
  targets = targets,
  method  = "lda",     # "lda" (default) or "logit"
  n_perm  = 199,       # permutation test for significance
  seed    = 99
)
print(cls)
#> DKGE Classification
#> --------------------
#> Targets: 2
#> Classifier: lda
#> Metrics: accuracy, logloss
#> Permutations: 199
#>   cond: accuracy=1.000, logloss=0.000
#>     accuracy p=0.005, logloss p=0.005
#>   time: accuracy=0.500, logloss=0.678
#>     accuracy p=0.010, logloss p=0.005
```

The result is a `dkge_classification` object. Each entry of `$results`
corresponds to one target:

``` r
res <- cls$results[["cond"]]
res$mode       # decoding mode: "cell" for within-subject targets
#> [1] "cell"
res$metrics    # cross-validated accuracy and log-loss
#>     accuracy      logloss 
#> 1.000000e+00 9.999779e-13
```

Convert to a tidy data frame for plotting or downstream analysis:

``` r
df <- as.data.frame(cls)
head(df)
#>   target   metric        value p_value n_perm
#> 1   cond accuracy 1.000000e+00   0.005    199
#> 2   cond  logloss 9.999779e-13   0.005    199
#> 3   time accuracy 5.000000e-01   0.010    199
#> 4   time  logloss 6.781927e-01   0.005    199
```

### Permutation p-values

When `n_perm > 0`, each target gets an empirical p-value from the
sign-flip max-T permutation distribution:

``` r
res$p_values    # named numeric vector per metric (NULL when n_perm = 0)
#> accuracy  logloss 
#>    0.005    0.005
```

### Subject-level predictions

Individual predictions and fold assignments live in `$row_data`:

``` r
head(res$row_data[, c("subject_label", "class_label", "fold")])
#>   subject_label class_label fold
#> 1          sub1           1    1
#> 2          sub1           2    1
#> 3          sub2           1    2
#> 4          sub2           2    2
#> 5          sub3           1    3
#> 6          sub3           2    3
```

Full predicted probabilities are in `res$probabilities` (rows = subject
× class combinations, columns = class labels).

### Per-fold diagnostics

Confusion matrices and class counts per fold are in
`res$diagnostics$folds`:

``` r
diag_fold1 <- res$diagnostics$folds[[1]]
diag_fold1$confusion           # confusion matrix for fold 1
#>    
#>     1 2
#>   1 1 0
#>   2 0 1
diag_fold1$class_counts_test   # observed class counts in test set
#> 1 2 
#> 1 1
```

## Two modes: `cell` vs `cell_cross`

The `mode` argument controls how the group basis is applied at test
time:

| Mode           | Basis used                      | Use case                                                 |
|----------------|---------------------------------|----------------------------------------------------------|
| `"cell"`       | Global `fit$U`                  | Fast; mild basis-step leakage — `fit$U` saw all subjects |
| `"cell_cross"` | Fold-specific LOSO `U_fold`     | Fully cross-validated, more conservative                 |
| `"delta"`      | Global `fit$U` + subject labels | Subject-level binary test; requires `y` argument         |

`"auto"` (default) selects `"cell"` for within-subject targets and
`"delta"` for between-subject targets.

``` r
# Fully cross-validated basis — recommended for publication-quality results
cls_strict <- dkge_classify(fit, targets = targets, mode = "cell_cross", n_perm = 0)
```

## Multiple betas per condition

When you have multiple beta estimates per condition (e.g., separate
scanner runs), stack the run-specific effects as extra rows in the
design matrix and supply a weight matrix that combines them into a
single class pattern.

``` r
# 4 effects: A_run1, A_run2, B_run1, B_run2
make_subject_multi <- function(id) {
  design <- diag(4)
  colnames(design) <- c("A_run1", "A_run2", "B_run1", "B_run2")

  signal <- matrix(0, nrow = 4, ncol = v)
  signal[1:2, 1:10] <-  0.5   # both A runs share signal
  signal[3:4, 1:10] <- -0.5   # both B runs share the opposite signal

  beta <- signal + matrix(rnorm(4 * v, sd = 1.0), nrow = 4)
  dkge_subject(beta, design = design, id = paste0("sub", id))
}

subjects_multi <- lapply(seq_len(n_subjects), make_subject_multi)
fit_multi      <- dkge(subjects_multi, K = diag(4), rank = 2)

# Weight matrix: rows = classes, columns = effects.
# Each class averages its two runs.
W_runs <- matrix(0, nrow = 2, ncol = 4)
W_runs[1, 1:2] <- 0.5    # class A = mean(A_run1, A_run2)
W_runs[2, 3:4] <- 0.5    # class B = mean(B_run1, B_run2)
rownames(W_runs) <- c("A", "B")
```

A plain matrix is accepted directly as `targets` — DKGE wraps it
automatically:

``` r
cls_multi <- dkge_classify(fit_multi,
                           targets = W_runs,   # plain matrix dispatch
                           n_perm  = 99,
                           seed    = 101)

res_multi <- cls_multi$results[[1]]
res_multi$metrics
#>     accuracy      logloss 
#> 1.000000e+00 6.953406e-06
as.data.frame(cls_multi)
#>    target   metric        value p_value n_perm
#> 1 target1 accuracy 1.000000e+00    0.01     99
#> 2 target1  logloss 6.953406e-06    0.01     99
```

The cross-validation remains unbiased: run averaging happens inside the
weight matrix, so each LOSO fold observes exactly one pattern per
condition.

## Hyperdesign inputs and fold bridges

Labs using `multidesign::hyperdesign()` can supply design kernels and
fold assignments through coercion S3 generics. Any object implementing
[`as_dkge_kernel()`](https://bbuchsbaum.github.io/dkge/reference/as_dkge_kernel.md)
and
[`as_dkge_folds()`](https://bbuchsbaum.github.io/dkge/reference/as_dkge_folds.md)
flows through the existing pipeline without altering the core solvers.

``` r
# Pseudocode — requires multidesign package and a user-defined hyperdesign object
library(multidesign)

hd     <- make_demo_hyperdesign()                # user-supplied helper
Kobj   <- as_dkge_kernel(hd, basis = "effect")   # list(K = ..., info = ...)

fit_hd <- dkge(
  betas       = dkge_data_from_hd(hd),
  K           = Kobj,
  keep_inputs = TRUE
)

folds  <- as_dkge_folds(fold_over(hd, over = "subject", k = 5, seed = 1), fit_hd)

res_hd <- dkge_contrast(fit_hd,
                        contrasts = c(1, -1, 0, 0),
                        method    = "kfold",
                        folds     = folds)
```

Matrices, lists with a `$K` element, and existing `dkge_folds` objects
continue to work as before — only packages with richer design structures
need to implement these generics.

## Practical tips

**Target specification.** Define the full factorial structure with
formula notation (e.g., `~ cond + time + cond:time`) and reuse the same
target definitions across both classification and contrast analyses for
methodological consistency.

**Backend choice.** LDA (`method = "lda"`) is fast and usually
sufficient. Use `method = "logit"` with `class_weights = "balanced"`
when experimental conditions have unequal numbers of trials.

**Permutation testing.** Increase `n_perm` for precise p-values near
significance thresholds. Set `n_perm = 0` during exploratory analysis to
skip the permutation loop.

**Spatial alignment.** When subjects have different voxel grids, apply
transport to a common medoid parcellation via
[`dkge_transport_contrasts_to_medoid()`](https://bbuchsbaum.github.io/dkge/reference/dkge_transport_contrasts_to_medoid.md)
before classification.

**Input caching.** Keep `keep_inputs = TRUE` (default) to enable
[`dkge_update_weights()`](https://bbuchsbaum.github.io/dkge/reference/dkge_update_weights.md)
without re-running the full fit.

## Next steps

- **Adaptive weighting**
  ([`vignette("dkge-adaptive-weighting")`](https://bbuchsbaum.github.io/dkge/articles/dkge-adaptive-weighting.md)):
  Emphasize spatially reliable voxels before computing the group
  embedding.
- **[`dkge_pipeline()`](https://bbuchsbaum.github.io/dkge/reference/dkge_pipeline.md)**:
  Orchestrate fit → contrast → transport → inference → classification in
  a single call.
- **Diagnostics**: Inspect `fit$evals` (variance explained) and per-fold
  confusion matrices to identify underperforming subjects or poorly
  separated conditions.
