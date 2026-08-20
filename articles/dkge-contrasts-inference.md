# Contrasts, Inference, and Bootstraps

This vignette follows one estimand from construction to interpretation:
the subject-level effect-1 minus effect-2 contrast on a shared simulated
spatial support. It then distinguishes that contrast estimand from
component-level inference and shows how a subject bootstrap quantifies
sampling variation conditional on the fitted basis and transport
operators.

## Demo Dataset

The first two effect rows contain opposite halves of a smooth spatial
signal, so their raw contrast equals `planted`. The remaining effects
contain noise. Every subject uses the same cluster ordering with small
coordinate jitter, giving the transport step a known target while
retaining subject variation.

``` r

library(dkge)
S <- 8; q <- 4; P <- 25; T <- 60
effect_names <- paste0("effect", seq_len(q))
cluster_axis <- seq(-1, 1, length.out = P)
planted <- exp(-5 * (cluster_axis + 0.35)^2) -
  0.65 * exp(-7 * (cluster_axis - 0.45)^2)
truth <- matrix(0, q, P, dimnames = list(effect_names, NULL))
truth[1, ] <- planted / 2
truth[2, ] <- -planted / 2

betas <- lapply(seq_len(S), function(s) {
  truth + matrix(rnorm(q * P, sd = 0.30), q, P)
})
designs <- lapply(seq_len(S), function(s) {
  X <- matrix(rnorm(T * q), T, q)
  X <- qr.Q(qr(X))
  colnames(X) <- effect_names
  X
})
reference_centroids <- cbind(x = 10 * cluster_axis, y = 0, z = 0)
centroids <- lapply(seq_len(S), function(s) {
  reference_centroids + matrix(rnorm(P * 3, sd = 0.02), P, 3)
})
subjects <- lapply(seq_len(S), function(s) dkge_subject(betas[[s]], designs[[s]], id = paste0("sub", s)))
bundle <- dkge_data(subjects)
fit <- dkge(bundle, K = diag(q), rank = 3)
fit$centroids <- centroids  # attach for transport helpers; pass explicitly in production
```

## Leave-One-Subject-Out (LOSO) Contrasts

The LOSO approach refits the latent basis without each held-out subject
before computing that subject’s contrast field. This limits basis-reuse
optimism; it does not by itself guarantee unbiased population inference.
The public `dkge_contrast(..., method = "loso")` path returns those
subject-specific fields, after which one fixed set of transport
operators can place them on a common reference for the inferential
procedure defined below.

``` r

contrast_vec <- c(1, -1, 0, 0)  # effect1 minus effect2
loso <- dkge_contrast(fit, contrast_vec, method = "loso")

transported <- dkge_transport_contrasts_to_medoid(
  fit,
  loso,
  medoid = 1L,
  centroids = centroids,
  mapper = dkge_mapper_spec("sinkhorn", epsilon = 0.05, lambda_spa = 1,
                            max_iter = 5000, tol = 1e-4)
)

str(loso$values, max.level = 1)
#> List of 1
#>  $ contrast1:List of 8
```

The resulting data structure contains several important components that
capture different stages of the contrast computation process. The
`loso$values` component stores the raw subject-level contrast vectors
before any spatial transport has been applied. After transport to the
medoid space, `transported[[1]]$value` contains the contrast values that
have been aligned to the medoid parcellation and combined by the
coordinate-wise median across subjects. Additionally,
`transported[[1]]$subj_values` preserves the individual-subject
contributions after transport, which can be useful for detecting
outliers or understanding between-subject variability.

We summarize the region where the planted raw contrast is largest. DKGE
projection and transport change the numerical scale, so agreement in
direction is meaningful here; equality to the raw planted amplitude is
not expected.

``` r

active_region <- which(planted >= stats::quantile(planted, 0.80))
active_subject_values <- rowMeans(
  transported[[1]]$subj_values[, active_region, drop = FALSE]
)
contrast_summary <- data.frame(
  target = "mean transported contrast in planted-positive region",
  estimate = mean(active_subject_values),
  subject_sd = sd(active_subject_values),
  positive_subjects = sum(active_subject_values > 0),
  n_subjects = length(active_subject_values)
)
contrast_summary
#>                                                 target estimate subject_sd
#> 1 mean transported contrast in planted-positive region    0.808      0.185
#>   positive_subjects n_subjects
#> 1                 8          8
```

In this seeded example, the active-region estimate is 0.808, and 8 of 8
subjects have the planted positive direction. This is recovery evidence
for the simulation, not a population effect size or external validation
result.

We can visualize the LOSO medoid map to quickly inspect the spatial
pattern of contrast values.

``` r

loso_medoid <- transported[[1]]$value
plot(loso_medoid, type = "h", lwd = 2, main = "LOSO contrast on medoid clusters",
     xlab = "Medoid cluster", ylab = "Contrast value")
abline(h = 0, col = "grey60", lty = 2)
```

![Lollipop plot of LOSO contrast values on medoid
clusters.](dkge-contrasts-inference_files/figure-html/loso-plot-1.png)

## Analytic Inference on Components

[`dkge_component_stats()`](https://bbuchsbaum.github.io/dkge/reference/dkge_component_stats.md)
answers a different question from the contrast above: it transports
subject component loadings and tests their coordinate-wise group means.
The parametric option below uses a t reference distribution, then
applies the requested p-value adjustment across returned
component-coordinate tests.

``` r

comp_stats <- dkge_component_stats(
  fit,
  components = 1:2,
  mapper = dkge_mapper_spec("sinkhorn", epsilon = 0.05, lambda_spa = 0.5,
                            max_iter = 2000, tol = 1e-4),
  centroids = centroids,
  inference = list(type = "parametric"),
  medoid = 1L
)
head(comp_stats$summary)
#>   component cluster   stat        p   p_adj significant
#> 1         1       1  6.100 4.91e-04 0.00136        TRUE
#> 2         1       2  1.838 1.09e-01 0.13924       FALSE
#> 3         1       3  0.798 4.51e-01 0.50110       FALSE
#> 4         1       4  4.763 2.05e-03 0.00331        TRUE
#> 5         1       5  4.950 1.66e-03 0.00296        TRUE
#> 6         1       6 12.473 4.90e-06 0.00013        TRUE
```

The table reports component, reference-cluster index, t statistic, raw
p-value, FDR-adjusted p-value, and the declared significance flag. These
are component-loading tests under the mapper and parametric assumptions
above; they are not tests of the effect-1 minus effect-2 contrast.

## Bootstrap Inference

The following nonparametric subject bootstrap resamples the already
transported contrast vectors and recomputes their mean. It estimates
sampling variation conditional on the fitted basis, chosen medoid, and
fixed transport operators; it does not propagate uncertainty from those
earlier steps.

``` r

subject_medoid <- lapply(seq_len(nrow(transported[[1]]$subj_values)), function(idx) {
  transported[[1]]$subj_values[idx, ]
})

boot <- dkge_bootstrap_projected(
  subject_medoid,
  B = 200,
  seed = 2
)

head(boot$medoid$mean)
#> [1] 0.4953 0.1480 0.0644 0.5279 0.5701 0.8248
```

For the declared active-region summary, the bootstrap distribution is:

``` r

active_boot <- rowMeans(boot$medoid$boot[, active_region, drop = FALSE])
bootstrap_summary <- data.frame(
  bootstrap_mean = mean(active_boot),
  lower_95 = unname(stats::quantile(active_boot, 0.025)),
  upper_95 = unname(stats::quantile(active_boot, 0.975))
)
bootstrap_summary
#>   bootstrap_mean lower_95 upper_95
#> 1          0.807    0.714     0.94
```

The resulting conditional 95% percentile interval is \[0.714, 0.94\].
Because the earlier basis and transport operators are held fixed, this
interval should not be described as unconditional uncertainty for the
full DKGE pipeline.

[`dkge_bootstrap_projected()`](https://bbuchsbaum.github.io/dkge/reference/dkge_bootstrap_projected.md)
also returns coordinate-wise means, standard deviations, z ratios,
intervals, and (optionally) the draws used above. When the estimand
lives in q-space instead,
[`dkge_bootstrap_qspace()`](https://bbuchsbaum.github.io/dkge/reference/dkge_bootstrap_qspace.md)
performs multiplier resampling there and can reuse a declared transport
cache.

## Practical Recommendations

Choose the inferential target before choosing a procedure. The component
table tests transported component loadings; the projected bootstrap
above summarizes an explicit contrast after transport. Parametric
inference relies on its reference-distribution assumptions. Subject
resampling avoids that particular reference distribution but still
requires representative, independent subjects and enough observations
for the empirical distribution to be informative.

Before reporting results, state the contrast, analysis space,
aggregation rule, resampling unit, mapper, medoid, and multiplicity
correction. Inspect subject-level values as a diagnostic, but do not
remove unusual subjects without a prespecified rule and sensitivity
analysis.
