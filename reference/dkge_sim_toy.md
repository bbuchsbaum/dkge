# Simulate toy DKGE datasets with known factorial structure

Generates subject-level beta matrices \\B_s\\ (\`q x P_s\`) in effect
space with user-selected active design terms, controlled signal-to-noise
ratio (SNR), and per-subject noise levels. The routine also returns the
per-subject design 'rulers' (typically identity matrices) and the design
kernel \\K\\ so the synthetic dataset can be fed directly to
\[dkge_fit()\]. Components are planted inside chosen effect blocks, then
\\K\\-orthonormalised to provide a ground truth subspace.

## Usage

``` r
dkge_sim_toy(
  factors,
  terms = NULL,
  active_terms,
  r_per_term = NULL,
  S = 3,
  P = 10,
  snr = 8,
  noise_scales = NULL,
  seed = 1L,
  contrasts_type = c("helmert", "sum")
)
```

## Arguments

- factors:

  Named list describing the experimental factors (as for
  \[design_kernel()\]), e.g. \`list(A = list(L = 2), B = list(L = 3))\`.

- terms:

  Optional list of character vectors specifying which terms to include
  in the kernel. Defaults to the full factorial set used by
  \[design_kernel()\].

- active_terms:

  Character vector of term names to activate (e.g., \`c("A", "B",
  "A:B")\`). These must be present in the kernel.

- r_per_term:

  Named integer vector giving how many latent columns to draw within
  each active term. Defaults to one column per term.

- S:

  Number of subjects.

- P:

  Either a single integer (clusters/voxels per subject) or a
  length-\`S\` integer vector.

- snr:

  Target Frobenius SNR (signal / noise) per subject. Scalar or
  length-\`S\`.

- noise_scales:

  Optional length-\`S\` multipliers applied to the noise standard
  deviation for each subject.

- seed:

  RNG seed for reproducibility.

- contrasts_type:

  Factor-contrast system, one of \`"helmert"\` or \`"sum"\`.

## Value

A list with the following entries:

- B_list:

  List of \`q x P_s\` beta matrices.

- X_list:

  List of subject-level design rulers (identity matrices).

- K:

  \`q x q\` design kernel.

- info:

  Design metadata returned by \[design_kernel()\].

- U_true:

  Ground-truth \`q x r_true\` component basis (\\K\\-orthonormal).

- M_list:

  Per-subject spatial patterns (\`P_s x r_true\`).

- active_cols:

  Named list mapping each active term to the selected column indices in
  the effect basis.

- subject_ids:

  Character vector of subject identifiers.

## Details

Construction follows the model \$\$B_s = U\_{true} M_s^\top + E_s,\$\$
where \`U_true\` lies inside the requested design-term blocks, \`M_s\`
encodes the subject-specific spatial loadings, and \`E_s\` is Gaussian
noise adjusted to hit the requested Frobenius SNR.

## Examples

``` r
# Generate 3-subject toy data with 2x3 factorial design
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"),
  S = 3, P = 20, snr = 5
)
length(toy$B_list)   # 3 subjects
#> [1] 3
dim(toy$B_list[[1]]) # 5 effects x 20 clusters
#> [1]  5 20
dim(toy$K)           # 5x5 kernel
#> [1] 5 5
```
