# Project subjects onto DKGE components

Computes one \*\*signed\*\* scalar per subject and component: the
subject's mean cluster brain score along that component.

## Usage

``` r
dkge_subject_component_projections(
  fit,
  groups = NULL,
  mode = c("loso", "pooled"),
  comps = NULL,
  align = TRUE,
  ridge = 0
)
```

## Arguments

- fit:

  Fitted \`dkge\` object.

- groups:

  Optional vector of group labels, or a data frame with subject and
  group columns. Named vectors are matched by subject name and must
  cover every subject; unnamed vectors are matched positionally.

- mode:

  \`"loso"\` for held-out supplementary projections or \`"pooled"\` for
  descriptive projections on the pooled fit.

- comps:

  Components to include.

- align:

  Rotate each LOSO fold basis onto the pooled component axes
  (\`fit\$U\`) by K-Procrustes before scoring.

- ridge:

  Ridge passed to the LOSO fold basis builder.

## Value

Tidy data frame with subject, group, component, component_id,
projection, and mode.

## Details

For subject \\s\\ with whitened, voxel-weighted betas \\\tilde B_s\\
(\\q \times P_s\\) and group component \\u_j\\, the per-cluster brain
scores are the columns of \\a\_{s,j} = \tilde B_s^\top K u_j \in
R^{P_s}\\. The reported projection is their average over clusters,
\$\$\pi\_{s,j} = \frac{1}{P_s} \sum\_{p} a\_{s,j}\[p\] = \langle K u_j,
\bar b_s \rangle, \qquad \bar b_s = \frac{1}{P_s}\sum_p \tilde B_s\[,
p\].\$\$ Equivalently it is the K-inner product between the component
salience \\K u_j\\ and the subject's cluster-averaged effect-space
profile.

This quantity is \*\*linear\*\* in the subject's betas, so it is signed
and \\\pi\_{s,j}\\ flips sign when \\\tilde B_s\\ is negated. That is
the point of the function: it complements the sign-blind quadratic
participation measure \\\lVert \tilde B_s^\top K u_j \rVert^2\\ reported
by \[dkge_plot_subject_contrib()\], which cannot distinguish a subject
expressing a component from a subject expressing its mirror image.
Averaging over clusters is what makes the score subject-independent in
cluster space and therefore comparable across subjects with different
parcel counts; a subject whose cluster scores are large but evenly split
in sign will have a small projection and a large energy.

In \`"loso"\` mode the same definition is applied with the basis
\\U^{(-s)}\\ estimated without subject \\s\\, using the fold loaders'
voxel weighting. When \`align = TRUE\` each fold basis is rotated onto
the \*pooled\* basis \`fit\$U\` by K-Procrustes
(\[dkge_procrustes_K()\]), so component \\j\\ names the same pooled axis
for every subject. (The fold builder's own alignment uses fold 1 as its
reference, not \`fit\$U\`, and is deliberately not reused here.) With
\`align = FALSE\` component signs and order are only defined up to each
fold's own eigen-solve, so cross-subject comparison is not meaningful.

Voxel weights are normalised to mean 1 and applied as \\\sqrt{w}\\
column scales — the same convention as the training blocks — so a
uniform weight of any magnitude leaves the projection unchanged. This is
not the column-wise weighted mean \\\sum_p w_p b_p / \sum_p w_p\\; that
would change the scale relative to \[dkge_plot_subject_contrib()\]
energy, which uses the same \\\sqrt{w}\\ reweighting. Weights produced
by \[dkge_fit()\] already have mean 1.

## Examples

``` r
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 4, P = 15, snr = 5
)
fit <- dkge(toy$B_list, toy$X_list, K = toy$K, rank = 2)
dkge_subject_component_projections(fit, mode = "pooled", comps = 1:2)
#>   subject group component component_id   projection   mode
#> 1   sub01   all       LV1            1  0.365425776 pooled
#> 2   sub01   all       LV2            2 -0.003821853 pooled
#> 3   sub02   all       LV1            1 -0.027837086 pooled
#> 4   sub02   all       LV2            2  0.558071921 pooled
#> 5   sub03   all       LV1            1 -0.720802091 pooled
#> 6   sub03   all       LV2            2  0.081939297 pooled
#> 7   sub04   all       LV1            1 -0.383363113 pooled
#> 8   sub04   all       LV2            2  0.554849800 pooled
```
