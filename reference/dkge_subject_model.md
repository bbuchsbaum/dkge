# Build a between-subject model matrix

Wraps formula handling for subject-level DKGE analyses. The result owns
the model matrix, term metadata, and subject identifiers used by
\[dkge_between_rrr()\].

## Usage

``` r
dkge_subject_model(
  formula,
  data,
  subject_ids = NULL,
  nuisance = NULL,
  contrasts.arg = NULL,
  na.action = stats::na.fail
)
```

## Arguments

- formula:

  A model formula such as \`~ group \* trait + age + sex\`.

- data:

  Data frame containing subject-level variables.

- subject_ids:

  Optional subject identifiers aligned with \`data\` rows. Defaults to
  the \`subject_id\` column when present, then \`rownames(data)\`, then
  sequential IDs. Identifiers are stored as \`rownames\` of \`data\`
  before \[stats::model.frame()\] so \`na.action = stats::na.omit\`
  keeps IDs in sync with the retained rows.

- nuisance:

  Optional character vector of term labels to treat as nuisance in
  downstream inference. When \`dkge_between_permute(terms = NULL)\`,
  these terms are excluded from the default test set.

- contrasts.arg:

  Optional contrasts passed to \[stats::model.matrix()\].

- na.action:

  NA handler for the model frame. Defaults to \[stats::na.fail()\].

## Value

Object of class \`dkge_subject_model\`.

## Examples

``` r
dat <- data.frame(
  subject_id = paste0("s", 1:8),
  group = factor(rep(c("A", "B"), each = 4)),
  trait = rnorm(8)
)
dkge_subject_model(~ group * trait, dat)
#> <dkge_subject_model>
#>   subjects : 8 
#>   columns  : 4 
#>   terms    : (Intercept), group, trait, group:trait 
```
