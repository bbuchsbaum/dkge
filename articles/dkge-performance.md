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
function configures the dense spatial pipeline. A fitted Sinkhorn mapper
keeps the joint plan, the target-conditional application operator, and
convergence diagnostics. Reusing that fitted object makes later
[`apply_mapper()`](https://bbuchsbaum.github.io/dkge/reference/apply_mapper.md)
calls matrix multiplications rather than new OT solves.

## Speed Comparison

``` r

knitr::kable(performance, digits = 5)
```

| stage | mapper   | seconds_per_call |
|:------|:---------|-----------------:|
| fit   | kNN      |          0.00100 |
| fit   | Sinkhorn |          0.00400 |
| apply | kNN      |          0.00011 |
| apply | Sinkhorn |          0.00025 |

The table separates the solve from application. kNN fitting constructs
local neighbourhoods; Sinkhorn fitting constructs a dense cost matrix
and iterates to match both marginals, so its fit is usually the
expensive stage. Once fitted, both backends reuse their mappings.
Measure both stages at your own parcel and anchor counts rather than
extrapolating from apply-only timings.

## Warm Starts and Dual Caching

DKGE caches duals only after a solve satisfies its requested marginal
tolerance. The key includes every cost and mass entry, so warm starts
cannot cross between merely similar problems. The cache is
process-local; fitted renderer objects store plans and operators, not
dual vectors.

``` r

dkge_clear_sinkhorn_cache()
cached_mapper <- dkge_mapper("sinkhorn", epsilon = 0.05, lambda_xyz = 1)
cold_fit <- fit_mapper(cached_mapper, centroids, anchors)
warm_fit <- fit_mapper(cached_mapper, centroids, anchors)

rbind(
  cold = unlist(cold_fit$stats$diagnostics[c("iterations", "cache_hit")]),
  repeated = unlist(warm_fit$stats$diagnostics[c("iterations", "cache_hit")])
)
#>          iterations cache_hit
#> cold             10         0
#> repeated          0         1
```

For this exact repeat, the cached converged plan is returned with zero
new iterations. If you request a tighter tolerance, DKGE instead uses
the cached duals as an initialization and continues solving.

## Adding Custom Mappers

To integrate a custom mapping approach into the DKGE framework, you need
to implement two key S3 methods that follow the established naming
convention recognized by
[`dkge_mapper()`](https://bbuchsbaum.github.io/dkge/reference/dkge_mapper.md).

``` r

fit_mapper.dkge_mapper_mytransport <- function(mapper, subj_points,
                                                anchor_points, ...) {
  structure(list(Q = nrow(anchor_points)),
            class = "dkge_mapper_fit_mytransport")
}

apply_mapper.dkge_mapper_fit_mytransport <- function(fitted_mapper, values, ...) {
  rep(mean(values), fitted_mapper$Q)
}

custom <- dkge_mapper("mytransport")
custom_fit <- fit_mapper(custom, centroids, anchors)
head(apply_mapper(custom_fit, values))
#> [1] 0.172 0.172 0.172 0.172 0.172 0.172
```

The constructor accepts the custom identifier and S3 dispatch finds the
two methods. A renderer-compatible extension must also validate the
standard point inputs and return one value per anchor, as this example
does.

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

For repeated maps, reuse renderer or fitted-mapper objects so the
already-solved plans and application operators are retained. Warm-start
caching is most useful when an identical fit is requested again; check
`stats$diagnostics` instead of assuming that a cache was used.
