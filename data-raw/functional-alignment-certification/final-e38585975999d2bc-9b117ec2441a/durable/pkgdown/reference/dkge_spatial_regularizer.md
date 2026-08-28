# Construct a model-level DKGE spatial regularizer

Builds a sparse graph-Laplacian regularizer for the spatial columns of
each subject beta block. Coordinate-based construction delegates to
[`adjoin::spatial_laplacian()`](https://rdrr.io/pkg/adjoin/man/spatial_laplacian.html);
DKGE validates the resulting operator, binds it to the beta-column
domain, and smooths the beta block with the Tikhonov resolvent
\$\$H\_\lambda = (I + \lambda L)^{-1},\$\$ replacing \\\widetilde B_s\\
by \\\widetilde B_s H\_\lambda\\ in the fitted moment and in every
downstream component or contrast map.

## Usage

``` r
dkge_spatial_regularizer(
  coords = NULL,
  laplacian = NULL,
  lambda = 1,
  dthresh = 1.42,
  nnk = 27L,
  weight_mode = c("binary", "heat"),
  sigma = dthresh/2,
  normalized = FALSE,
  handle_isolates = c("keep_zero", "self_loop"),
  domain = NULL
)
```

## Arguments

- coords:

  Numeric coordinate matrix (`P x d`) or a list of such matrices, one
  per subject. Mutually exclusive with `laplacian`.

- laplacian:

  Precomputed graph Laplacian (`P x P`) or a list of Laplacians. It must
  be a symmetric combinatorial graph Laplacian: finite, non-negative
  diagonal, non-positive off-diagonal, and zero row sums.

- lambda:

  Non-negative Tikhonov smoothing strength. `0` is an exact identity
  operation and reproduces an unsmoothed DKGE fit.

- dthresh, nnk, weight_mode, sigma, handle_isolates:

  Arguments forwarded to
  [`adjoin::spatial_laplacian()`](https://rdrr.io/pkg/adjoin/man/spatial_laplacian.html)
  when `coords` is supplied. `dthresh` is a distance cutoff in the units
  of `coords` and `nnk` caps the number of neighbours; see the "Choosing
  `dthresh`" section, as the defaults assume unit-spaced voxel indices
  rather than millimetres. DKGE does not allow
  `handle_isolates = "drop"`, because dropping rows would change the
  beta column domain.

- normalized:

  Must currently be `FALSE`. DKGE deliberately uses a combinatorial
  Laplacian whose null space contains the constant field;
  degree-normalized Laplacians generally do not preserve constants under
  the resolvent used here.

- domain:

  Optional spatial-unit labels. Supply one character vector for a shared
  matrix or a list aligned with subject-specific matrices. Coordinate
  row names or Laplacian dimnames are used when `domain` is omitted.

## Value

An object of class `dkge_spatial_regularizer`, suitable for the
`spatial` argument of
[`dkge()`](https://bbuchsbaum.github.io/dkge/reference/dkge.md) or
[`dkge_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit.md).

## Effect on the pooled moment

Because the resolvent is applied to the beta block and the moment is
then formed from the smoothed block, \\H\_\lambda\\ enters the \\q
\times q\\ moment **twice**: \$\$M_s(\lambda) = \widetilde B_s
H\_\lambda \Omega_s H\_\lambda \widetilde B_s^{\mathsf T},\$\$ which
reduces to \\\widetilde B_s H\_\lambda^2 \widetilde B_s^{\mathsf T}\\
under an identity spatial metric. The effective moment-level smoothing
is therefore the *squared* resolvent, not \\H\_\lambda\\; keep that in
mind when comparing `lambda` against a target smoothing kernel.

Two per-unit weightings sit on **opposite sides** of the smoother.
Writing \\W\\ for
[`dkge_weights()`](https://bbuchsbaum.github.io/dkge/reference/dkge_weights.md)
voxel weights and \\\Omega_s\\ for the spatial metric in `Omega_list`,
the moment DKGE actually accumulates is \$\$M_s(\lambda) = (\widetilde
B_s W^{1/2}) H\_\lambda \Omega_s H\_\lambda (\widetilde B_s
W^{1/2})^{\mathsf T}.\$\$ Voxel weights are applied to the beta columns
*before* smoothing, so a down-weighted unit is shrunk first and then
diffuses its reduced value to its neighbours – a zero-weight unit acts
as a hole in the field rather than being removed from the graph.
`Omega_list` is applied *after* smoothing, as a metric on the smoothed
field. To exclude a unit completely, remove the matching beta column,
graph node, domain label, and every aligned per-unit input; changing
only `coords`, `laplacian`, or a voxel weight is not a domain exclusion
operation.

## Choosing `dthresh`

`dthresh` is expressed in the units of `coords`. The defaults
(`dthresh = 1.42`, `nnk = 27`) assume **unit-spaced voxel indices**,
where 1.42 admits face and edge neighbours of a 3x3x3 neighbourhood.
Passing millimetre coordinates with the default threshold isolates every
unit, which makes `L` the zero matrix and the resolvent the identity:
the fit is then bit-identical to an unregularized one despite a positive
`lambda`. Construction warns when a positive-`lambda` graph ends up with
no edges. Spatial status is `inactive` for `lambda = 0`, `inert` when
all graphs are edgeless, `partial` when only some subject graphs are
edgeless, and `active` when every graph can smooth. Both
[`print()`](https://rdrr.io/r/base/print.html) and
`dkge_diagnostics(fit)$spatial$diagnostics` report the status and edge
count.

Supply either `coords` or `laplacian`. A single matrix defines one
shared spatial domain and therefore requires every subject to have the
same number and ordering of spatial units. A list allows
subject-specific domains and may be named with the fitted subject IDs.
When calling
[`dkge_fit()`](https://bbuchsbaum.github.io/dkge/reference/dkge_fit.md)
on raw lists, use `dkge_data(..., subject_ids = ...)` to declare those
IDs explicitly. When both beta columns and `domain` are named, DKGE
matches them by name and fails if they do not describe the same units.

## Examples

``` r
coords <- cbind(x = 0:11, y = 0, z = 0)
spatial <- dkge_spatial_regularizer(
  coords = coords,
  lambda = 0.5,
  dthresh = 1.01,
  nnk = 3,
  weight_mode = "binary"
)
spatial
#> <dkge_spatial_regularizer>
#>   source   : adjoin 
#>   domains  : shared 
#>   lambda   : 0.5 
#>   status   : active 
#>   units    : 12 
#>   edges    : 11 
```
