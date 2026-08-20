# Plot subject projections onto DKGE components

Plot subject projections onto DKGE components

## Usage

``` r
dkge_plot_subject_component_projections(
  fit = NULL,
  groups = NULL,
  mode = c("loso", "pooled"),
  comps = NULL,
  align = TRUE,
  ridge = 0,
  projections = NULL
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

- projections:

  Optional data frame previously returned by
  \[dkge_subject_component_projections()\]. Supplying it skips
  recomputation, which matters for \`mode = "loso"\` because that path
  refits one basis per subject. It must contain \`subject\`, \`group\`,
  \`component\`, and a numeric \`projection\`. When supplied, \`fit\`,
  \`groups\`, \`comps\`, \`align\`, \`ridge\`, and \`mode\` are all
  ignored: the panel label is read from the frame's own \`mode\` column,
  and a frame without one (or mixing modes) is left unlabelled rather
  than being captioned with the \`mode\` default.

## Value

A ggplot object.

## Examples

``` r
toy <- dkge_sim_toy(
  factors = list(A = list(L = 2), B = list(L = 3)),
  active_terms = c("A", "B"), S = 4, P = 15, snr = 5
)
fit <- dkge(toy$B_list, toy$X_list, K = toy$K, rank = 2)
proj <- dkge_subject_component_projections(fit, mode = "pooled", comps = 1:2)
dkge_plot_subject_component_projections(fit, projections = proj)
```
