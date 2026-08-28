# Construct contrasts that isolate fitted DKGE components

Returns the effect-space contrast matrix `C = R %*% U`, where `R` is the
pooled-design ruler and `U` is the fitted K-orthonormal basis. These are
the contrasts that isolate fitted component coordinates because
`crossprod(U, K %*% backsolve(R, C))` equals the identity.

## Usage

``` r
dkge_component_contrasts(fit, comps = NULL)
```

## Arguments

- fit:

  A fitted `dkge` object.

- comps:

  Components to include. Defaults to all fitted components and otherwise
  accepts numeric component indices.

## Value

A numeric effects-by-components matrix suitable for the `contrasts`
argument of
[`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md).

## Details

This is deliberately different from
[`dkge_component_saliences()`](https://bbuchsbaum.github.io/dkge/reference/dkge_component_saliences.md),
which returns the dual read-out basis `K %*% U`. Saliences describe how
effects load on components; passing `K %*% U` back to
[`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md)
applies the kernel a second time and generally does not isolate
components. With LOSO or K-fold contrasts, genuine fold-to-fold basis
variation can still mix coordinates relative to the full-data basis.

## See also

[`dkge_component_saliences()`](https://bbuchsbaum.github.io/dkge/reference/dkge_component_saliences.md),
[`dkge_contrast()`](https://bbuchsbaum.github.io/dkge/reference/dkge_contrast.md)

## Examples

``` r
toy <- dkge_sim_toy(
  factors = list(cond = list(L = 3)), active_terms = "cond",
  S = 3, P = 10, snr = 4
)
fit <- dkge(toy$B_list, toy$X_list, K = toy$K, rank = 2)
C <- dkge_component_contrasts(fit)
round(crossprod(fit$U, fit$K %*% backsolve(fit$R, C)), 10)
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
```
