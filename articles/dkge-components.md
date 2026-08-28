# Components and Interpretability

A component number has no intrinsic scientific meaning. Once DKGE
returns a fitted basis you still have to say what each direction
represents, and that requires reading it in more than one place.

You need to connect three views of the same fitted direction: the
effects that define it, the spatial units on which subjects express it,
and any target used to rotate or label it.

This vignette asks how a fitted component is expressed in effect space
and on a shared toy cluster grid. It then shows rotation and projection
without treating those descriptive operations as inferential results.

``` text
effect saliences (what) -> subject cluster scores (where, within support) -> rotation or new-data projection (reuse)
```

## Example Fit

``` r

library(dkge)
S <- 4; q <- 5; P <- 18; T <- 70
effect_names <- paste0("effect", seq_len(q))
cluster_names <- paste0("cluster", seq_len(P))
cluster_axis <- seq(-1, 1, length.out = P)
shared_signal <- rbind(
  exp(-4 * (cluster_axis + 0.45)^2),
  -exp(-4 * (cluster_axis - 0.45)^2),
  cluster_axis,
  rep(0, P),
  rep(0, P)
)
betas <- lapply(seq_len(S), function(s) {
  B <- shared_signal + matrix(rnorm(q * P, sd = 0.35), q, P)
  dimnames(B) <- list(effect_names, cluster_names)
  B
})
designs <- replicate(S, {
  X <- matrix(rnorm(T * q), T, q)
  X <- qr.Q(qr(X))
  colnames(X) <- effect_names
  X
}, simplify = FALSE)
subjects <- lapply(seq_len(S), function(s) dkge_subject(betas[[s]], designs[[s]], id = paste0("sub", s)))
bundle <- dkge_data(subjects)
fit <- dkge(bundle, K = diag(q), rank = 3)
```

By default
[`dkge()`](https://bbuchsbaum.github.io/dkge/reference/dkge.md)
accumulates subject contributions in the shared effect space and
diagonalizes the pooled moment. To request joint diagonalization of the
subject matrices instead, set `solver = "jd"` (and optionally pass
`jd_control = dkge_jd_control(...)` to configure the optimizer). Treat
agreement or disagreement with the pooled solution as a model
comparison; the solver choice does not guarantee greater stability.

``` r

fit_jd <- dkge(bundle,
               K = diag(q),
               rank = 3,
               solver = "jd",
               jd_control = dkge_jd_control(maxit = 200, tol = 1e-8))
```

Both solvers return the same object structure (`fit$U`, `fit$Chat`,
projections), so downstream code does not change. subsequent workflow in
this vignette applies unchanged regardless of the solver choice.

## Which effects define each component?

``` r

c(effects = nrow(fit$U), components = ncol(fit$U), subjects = length(fit$Btil))
#>    effects components   subjects 
#>          5          3          4
```

The compact summary establishes the object shapes without printing
internal matrices. The interpretable effect-space quantity is the
design-weighted salience (K U):

``` r

dkge_plot_effect_loadings(fit, comps = 1:3)
```

![Heatmap showing how five design effects contribute to the first three
DKGE
components.](dkge-components_files/figure-html/loading-heatmap-1.png)

Read each column as a relative pattern. Effects with opposite colors
pull the component in opposite directions; a column dominated by one row
is closer to a single-effect axis. Component signs are arbitrary, so the
within-column pattern matters more than whether red or blue is positive.

The salience matrix $`K U`$ is a dual read-out basis; it is not the
contrast matrix that isolates fitted components. Passing those columns
to
[`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md)
would apply `K` again after the pooled-design ruler. Use
[`dkge_component_contrasts()`](https://bbuchsbaum.github.io/dkge/reference/dkge_component_contrasts.md)
instead. It returns $`C = R U`$, for which the full-data component
coefficients are exactly the identity:

``` r

component_contrasts <- dkge_component_contrasts(fit, comps = 1:3)
component_alpha <- crossprod(
  fit$U,
  fit$K %*% backsolve(fit$R, component_contrasts)
)
round(component_alpha, 10)
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
```

Those columns may be supplied directly to
[`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md).
Under LOSO or K-fold evaluation, genuine changes in a held-out basis can
still mix component coordinates relative to the full-data reference; the
helper prevents the additional, unintended mixing caused by using
$`K U`$ as an input contrast.

`fit$Btil` stores the row-standardized subject beta blocks. Multiplying
those blocks by the component directions produces spatial-unit scores,
which answer the next question.

## Projecting Subjects into Component Space

``` r

subject_scores <- dkge_project_btil(fit, fit$Btil)
str(subject_scores, max.level = 1)
#> List of 4
#>  $ sub1: num [1:18, 1:3] 2.339 2.324 0.844 1.03 1.314 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ sub2: num [1:18, 1:3] 1.78 1.62 1.73 0.98 1.94 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ sub3: num [1:18, 1:3] 1.715 2.515 0.924 1.741 2.142 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ sub4: num [1:18, 1:3] 2.818 2.517 1.179 1.312 0.685 ...
#>   ..- attr(*, "dimnames")=List of 2
```

Each subject score matrix has dimensions clusters × components. Cluster
rows are only comparable across subjects when they already refer to the
same spatial support. This simulation deliberately uses a shared
18-cluster indexing scheme, so a row-wise mean is meaningful. With
subject-specific parcellations, first transport the scores to an
identified reference, anchor, or voxel space; equal matrix dimensions
alone do not establish anatomical correspondence.

``` r

avg_scores <- Reduce("+", subject_scores) / length(subject_scores)
avg_df <- as.data.frame(avg_scores)
component_cols <- paste0("Component ", seq_len(ncol(avg_df)))
names(avg_df) <- component_cols
avg_df$Cluster <- seq_len(nrow(avg_df))
avg_long <- tidyr::pivot_longer(avg_df,
                                cols = seq_along(component_cols),
                                names_to = "Component",
                                values_to = "Score")

ggplot(avg_long, aes(x = Cluster, y = Score, colour = Component)) +
  geom_line(linewidth = 1.1) +
  labs(title = "Average component scores", y = "Score") +
  scale_x_continuous(breaks = seq_len(nrow(avg_scores))) +
  theme(legend.position = "top")
```

![Line plot showing average DKGE component scores across
clusters.](dkge-components_files/figure-html/score-plot-1.png)

The plot summarizes expression on this shared toy grid. It is
descriptive: a large value identifies a cluster-component pairing to
inspect, not a significant spatial effect.

## Rotating Components

[`dkge_procrustes_K()`](https://bbuchsbaum.github.io/dkge/reference/dkge_procrustes_K.md)
rotates components toward a target basis in the $`K`$ metric: canonical
contrasts, or loadings from an earlier study. The rotation stays inside
the design-kernel metric, so the solution’s properties are preserved.

``` r

# target basis: identity for first two effects
B_target <- dkge_k_orthonormalize(diag(1, q)[, 1:2], fit$K)
rot <- dkge_procrustes_K(B_target, fit$U[, 1:2], fit$K)
round(rot$U_aligned, 3)
#>        [,1]   [,2]
#> [1,]  0.796 -0.207
#> [2,] -0.207  0.757
#> [3,] -0.540 -0.615
#> [4,] -0.053 -0.029
#> [5,]  0.174 -0.071
```

The returned `rot$U_aligned` is the fitted basis rotated toward
`B_target`. `rot$d` reports the objective actually achieved by `rot$R`;
when reflections are forbidden, it can be smaller than
`rot$unconstrained_d`.

## Projecting New Data

``` r

new_beta <- matrix(rnorm(q * P), q, P)
projected <- dkge_project_clusters(fit, new_beta)
head(projected)
#>           [,1]    [,2]     [,3]
#> [1,]  0.109688  0.0776  0.12167
#> [2,]  0.000515 -0.3380  0.53417
#> [3,] -0.273249  0.2488  0.07049
#> [4,] -0.204818 -0.0151  0.15044
#> [5,]  0.018553  0.0610 -0.10772
#> [6,] -0.063557 -0.0554  0.00225
```

[`dkge_project_clusters()`](https://bbuchsbaum.github.io/dkge/reference/dkge_project_clusters.md)
returns component scores for each cluster in the new data, ready for
contrast testing or transport to a reference parcellation.

## Interpreting Components

Several outputs answer different interpretive questions. `fit$weights`
contains one fitting weight per subject, not a component-specific
contribution. For component-specific participation, inspect
`fit$contribs`,
[`dkge_plot_subject_contrib()`](https://bbuchsbaum.github.io/dkge/reference/dkge_plot_subject_contrib.md),
or
[`dkge_subject_component_projections()`](https://bbuchsbaum.github.io/dkge/reference/dkge_subject_component_projections.md).
`fit$v` concatenates subject-block cluster loadings; interpret its rows
spatially only when the blocks share a declared support or after
transport to a common space. The deprecated
[`dkge_component_stats()`](https://bbuchsbaum.github.io/dkge/reference/dkge_component_stats.md)
can provide a descriptive legacy summary, but it refuses inference
because its correspondence is learned from the full-fit loadings. For
spatial inference, construct typed aligned subject rows and call
[`dkge_infer_aligned()`](https://bbuchsbaum.github.io/dkge/reference/dkge_infer_aligned.md).

## Where to go next

- [`vignette("dkge-plotting")`](https://bbuchsbaum.github.io/dkge/articles/dkge-plotting.md)
  — how to display saliences, subject contributions and subspace
  stability without over-reading a polished figure.
- [`vignette("dkge-contrasts-inference")`](https://bbuchsbaum.github.io/dkge/articles/dkge-contrasts-inference.md)
  — what it takes to attach uncertainty to a component or a contrast,
  rather than describing one.
- [`vignette("dkge-concepts")`](https://bbuchsbaum.github.io/dkge/articles/dkge-concepts.md)
  — which claims a fitted component supports and which it does not.
