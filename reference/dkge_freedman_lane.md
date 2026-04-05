# Freedman-Lane permutations for DKGE (scaffold)

This function orchestrates Freedman-Lane permutations at the
\*time-series\* level: for each subject, fit the reduced model (without
the effect of interest), permute residuals, reconstruct surrogate data,
refit the full GLM to get B\* betas, then re-run DKGE LOSO to obtain a
group statistic (e.g., max-\|t\| over medoid clusters). It requires the
caller to provide three adapter functions (or rely on
'fmrireg'/'neuroim2'): - fit_glm(Y_s, X_s, X0_s) -\> list(beta = qxP,
beta0 = q0xP, resid = TxP) - resample_resid(resid_s) -\> resid_s\* (TxP)
\[permute or phase-randomize per run\] - transport_and_stat(B_list,
X_list, K, c) -\> scalar (e.g., max-\|t\|)

## Usage

``` r
dkge_freedman_lane(
  Y_list,
  X_list,
  X0_list,
  K,
  c,
  B = 500,
  adapters,
  seed = 123L
)
```

## Arguments

- Y_list:

  list of neuroim2 BrainVectors (or TxP matrices) per subject

- X_list:

  list of Txq design matrices (full)

- X0_list:

  list of Txq0 reduced designs (null space of the contrast)

- K:

  design kernel (qxq)

- c:

  contrast vector (qx1)

- B:

  number of permutations

- adapters:

  list with functions: fit_glm, resample_resid, transport_and_stat

- seed:

  RNG seed for reproducibility

## Value

list with fields: stat_obs, stat_null (B-vector), p, details
