# Mapper Customisation and Performance

This vignette provides guidance on selecting between kNN and Sinkhorn
mappers based on your analysis needs, demonstrates how to implement
custom transport methods, and shows how to leverage warm starts to
improve computational performance.

## Mapper Factory Recap

``` r
library(dkge)
dkge_mapper("knn", k = 8, sigx = 3)
#> $method
#> [1] "knn"
#> 
#> $pars
#> $pars$k
#> [1] 8
#> 
#> $pars$sigx
#> [1] 3
#> 
#> 
#> attr(,"class")
#> [1] "dkge_mapper_knn" "dkge_mapper"
dkge_mapper("sinkhorn", epsilon = 0.05, lambda_xyz = 1, lambda_feat = 0)
#> $method
#> [1] "sinkhorn"
#> 
#> $pars
#> $pars$epsilon
#> [1] 0.05
#> 
#> $pars$lambda_xyz
#> [1] 1
#> 
#> $pars$lambda_feat
#> [1] 0
#> 
#> 
#> attr(,"class")
#> [1] "dkge_mapper_sinkhorn" "dkge_mapper"
```

The
[`dkge_mapper()`](https://bbuchsbaum.github.io/dkge/reference/dkge_mapper.md)
function creates a mapper configuration that stores the backend
parameters for your chosen transport method. When you call
[`fit_mapper()`](https://bbuchsbaum.github.io/dkge/reference/fit_mapper.md)
on this configuration, it produces subject-specific transport plans that
include cached computational state, enabling efficient reuse across
multiple mapping operations.

## Speed Comparison

``` r
S <- 3; q <- 3; P <- 200
centroids <- replicate(S, matrix(rnorm(P * 3), P, 3), simplify = FALSE)
anchors <- matrix(rnorm(1500 * 3), 1500, 3)
values <- lapply(seq_len(S), function(s) rnorm(P))

knn_mapper <- dkge_mapper("knn", k = 12, sigx = 5)
knn_fit <- lapply(seq_len(S), function(s) fit_mapper(knn_mapper, centroids[[s]], anchors))

sink_mapper <- dkge_mapper("sinkhorn", epsilon = 0.02, lambda_xyz = 1)
sink_fit <- lapply(seq_len(S), function(s) fit_mapper(sink_mapper, centroids[[s]], anchors))

microbenchmark::microbenchmark(
  knn = lapply(seq_len(S), function(s) apply_mapper(knn_fit[[s]], values[[s]])),
  sinkhorn = lapply(seq_len(S), function(s) apply_mapper(sink_fit[[s]], values[[s]])),
  times = 10L
)
#> Unit: microseconds
#>      expr  min   lq mean median   uq   max neval
#>       knn  331  336  569    379  433  2297    10
#>  sinkhorn 5173 5245 7928   8314 9366 10799    10
```

The performance characteristics of these two mapping approaches differ
significantly. kNN mapping operates through purely local neighborhoods
and scales linearly with the number of neighbors `k`, making it
computationally efficient for most applications. In contrast, Sinkhorn
mapping offers greater expressiveness by supporting feature-based costs
and soft matching between points, but this flexibility comes with
increased computational overhead. However, the optimized C++ solver
implementation with warm start capabilities makes repeated Sinkhorn
calls much more tractable for iterative analyses.

## Warm Starts and Dual Caching

The example below calls `dkge:::sinkhorn_plan_cpp()`, an **internal**
C++ function. It is shown here for illustration only — the DKGE pipeline
manages warm starts automatically through its renderer objects. Do not
rely on this function’s signature in production code; use
[`fit_mapper()`](https://bbuchsbaum.github.io/dkge/reference/fit_mapper.md)
/
[`apply_mapper()`](https://bbuchsbaum.github.io/dkge/reference/apply_mapper.md)
instead.

``` r
C <- as.matrix(dist(anchors[1:200, ]))
mu <- rep(1/nrow(C), nrow(C))
nu <- rep(1/ncol(C), ncol(C))
# Internal C++ entry point — exposed here for illustration only
res1 <- dkge:::sinkhorn_plan_cpp(C, mu, nu, epsilon = 0.05, max_iter = 400, tol = 1e-7)
res2 <- dkge:::sinkhorn_plan_cpp(C, mu, nu, epsilon = 0.05, max_iter = 400, tol = 1e-7,
                                 log_u_init = res1$log_u, log_v_init = res1$log_v)
res1$iterations
res2$iterations
```

The second optimization run converges significantly faster because it
leverages the saved dual variables from the previous computation. This
warm start mechanism is automatically handled by the DKGE rendering
pipeline, which stores these dual variables in the renderer object for
subsequent use. When working directly with the lower-level
`sinkhorn_plan_cpp()` function for experimentation, you can manually
provide these initialization values to achieve similar performance
gains.

## Adding Custom Mappers

To integrate a custom mapping approach into the DKGE framework, you need
to implement two key S3 methods that follow the established naming
convention recognized by
[`dkge_mapper()`](https://bbuchsbaum.github.io/dkge/reference/dkge_mapper.md).

``` r
fit_mapper.dkge_mapper_mytransport <- function(mapper, subj_points, anchor_points, ...) {
  # return object with class 'dkge_mapper_fit_mytransport'
}

apply_mapper.dkge_mapper_fit_mytransport <- function(fitted_mapper, values, ...) {
  # return length(anchor_points) vector
}
```

As long as you maintain consistent output dimensions that match the
expected anchor point structure, DKGE will seamlessly integrate your
custom mapper into the existing analysis pipeline without requiring
additional configuration.

## Practical Guidance

The choice between mapping methods should be guided by your specific
analysis requirements and data characteristics. For quick exploratory
analyses or situations where subject clusters already show good
anatomical alignment, kNN mapping provides an efficient and
straightforward solution.

However, when dealing with substantial anatomical shifts between
subjects or when latent features should guide the matching process,
Sinkhorn mapping becomes the preferred choice. In these cases, you can
tune the `epsilon` parameter upward to achieve faster convergence and
smoother transport plans, though this comes with a trade-off in matching
precision.

To maximize computational efficiency across multiple analyses, make sure
to reuse renderer objects when performing bootstraps or computing
contrasts. This practice allows you to exploit the cached dual variables
and avoid the computational overhead of recomputing transport plans from
scratch.
