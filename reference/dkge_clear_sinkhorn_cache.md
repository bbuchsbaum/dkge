# Clear cached dual variables for Sinkhorn warm-starts

Releases memory held by the internal Sinkhorn cache. Call this after
large transport batches or before serialising objects.

## Usage

``` r
dkge_clear_sinkhorn_cache()
```

## Value

Logical \`TRUE\` invisibly.

## Examples

``` r
dkge_clear_sinkhorn_cache()
```
