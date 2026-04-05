# Fit a Design-Kernel Group Embedding (DKGE) model

Fit a Design-Kernel Group Embedding (DKGE) model

## Usage

``` r
dkge_fit(
  data,
  designs = NULL,
  K = NULL,
  Omega_list = NULL,
  w_method = c("mfa_sigma1", "energy", "none"),
  w_tau = 0.3,
  ridge = 0,
  rank = NULL,
  keep_X = FALSE,
  cpca_blocks = NULL,
  cpca_T = NULL,
  cpca_part = c("none", "design", "resid", "both"),
  cpca_ridge = 0,
  weights = NULL,
  solver = c("pooled", "jd"),
  jd_control = dkge_jd_control(),
  jd_mask = NULL,
  jd_init = NULL
)
```

## Arguments

- data:

  \`dkge_data\` bundle (preferred) or a list of subject beta matrices.

- designs:

  Per-subject design matrices (ignored when \`data\` already carries
  them). Columns must align with the effects encoded by \`K\`.

- K:

  qxq design kernel in effect space. Supply either a raw matrix or the
  list returned by \[design_kernel()\].

- Omega_list:

  Optional per-subject spatial reliabilities/weights. Each element may
  be \`NULL\` (no extra weighting), a numeric vector of length \`P_s\`
  (diagonal weights for clusters/voxels), or a full \`P_s x P_s\` matrix
  specifying custom covariance. These weights are applied both when
  accumulating the compressed covariance and when computing MFA/energy
  block normalisation.

- w_method:

  Subject-level weighting scheme. \* \`"mfa_sigma1"\` (default): inverse
  squared leading singular value of \\K^{1/2} Btil_s \Omega_s^{1/2}\\
  (Multiple Factor Analysis scaling). \* \`"energy"\`: inverse Frobenius
  norm squared of the same block. \* \`"none"\`: disable block scaling
  (all weights = 1).

- w_tau:

  Shrinkage parameter (0..1) toward equal weights. 0 keeps the raw
  weights, 1 forces equal weighting.

- ridge:

  Ridge regularization parameter (default 0).

- rank:

  Desired rank for the decomposition. If NULL, uses full rank.

- keep_X:

  Logical; when \`TRUE\`, store the concatenated training matrix used to
  build the multiblock projection (can be large).

- cpca_blocks:

  Optional integer vector specifying which effect rows span the CPCA
  design subspace. Ignored when \`cpca_part = "none"\` or when
  \`cpca_T\` is provided.

- cpca_T:

  Optional qxq0 matrix giving the CPCA design basis explicitly. When
  supplied it overrides \`cpca_blocks\`.

- cpca_part:

  Which CPCA-filtered component to estimate: \`"none"\` (default)
  performs the standard DKGE fit; \`"design"\` uses only the CPCA design
  part; \`"resid"\` uses the residual part; \`"both"\` fits the design
  part but also stores the residual basis for inspection.

- cpca_ridge:

  Optional ridge added to the CPCA-filtered matrices before
  eigen-decomposition.

- weights:

  Either \`NULL\` (default) to apply the weighting scheme specified by
  \`w_method\`, or a \[\`dkge_weights()\`\] specification controlling
  additional voxel/anchor-level weighting. When supplied, it must
  inherit from \`dkge_weights\` and is resolved via \[dkge_weights()\].

- solver:

  Solver used for the q-space problem. \`"pooled"\` keeps the original
  eigen-decomposition; \`"jd"\` performs joint diagonalisation using
  \[dkge_jd_control()\] settings.

- jd_control:

  Control parameters for the JD solver.

- jd_mask:

  Optional mask (matrix or list of matrices) applied to the off-diagonal
  penalty during JD. When supplied as a single matrix it is recycled
  across subjects.

- jd_init:

  Optional orthogonal initialiser for the JD solver expressed in the
  whitened \\K^{1/2}\\ metric (qxq matrix).

## Value

A fitted \`dkge\` object. When the multivarious package is installed the
return value additionally inherits from \`multiblock_biprojector\`,
making it compatible with the multivarious multiblock interface.

## Examples

``` r
# Simulate toy data with known structure
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 4, P = 30, snr = 5
)

# Bundle into dkge_data and fit
data <- dkge_data(toy$B_list, toy$X_list)
fit <- dkge_fit(data, K = toy$K, rank = 2)
fit$rank
#> [1] 2
```
